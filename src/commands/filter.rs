//! `filter`: allele-wise filtering, tag refilling and bcftools filtering
//! expressions, applied in parallel over genomic windows.
//!
//! Per record, in order:
//! 1. `-S/--samples-file` sample subsetting (as `bcftools view -S`, including
//!    recomputing INFO/AC and INFO/AN from the kept genotypes);
//! 2. `-I/--include-genotypes` or `-E/--exclude-genotypes`: set the GT of
//!    failing samples to missing, as `bcftools filter -S .` does with
//!    `-i`/`-e` (updating INFO/AC and INFO/AN in place when present);
//! 3. `--min-ac/--min-af` allele filtering (as `allelefilter`); a site with no
//!    passing ALT allele is dropped (even in soft-filter mode);
//! 4. `--fill-tags` recomputation (as `fill-tags` / `bcftools +fill-tags`);
//! 5. the `-i/--include` or `-e/--exclude` expression, evaluated by bcftools'
//!    own expression engine (vendored `filter.c`), with `bcftools filter`
//!    semantics for the FILTER column (`-s/--soft-filter`, `-m/--mode`).

use std::sync::Arc;

use anyhow::{Result, anyhow, bail};
use clap::Args;
use rust_htslib::bcf::Record;
use rust_htslib::htslib as raw;

use crate::args::IndexedInput;
use crate::commands::io;
use crate::filltags::{self, FillTags, TagSet};
use crate::hts::{self, HeaderView, SampleSubset};
use crate::processor::RecordFn;

#[derive(Args, Debug)]
pub struct FilterArgs {
    #[command(flatten)]
    pub indexed: IndexedInput,

    /// Output path. Indexed input writes bgzf VCF and requires this; stdin
    /// input writes uncompressed VCF to stdout when omitted or "-".
    #[arg(long = "output", short = 'o')]
    pub output: Option<String>,

    /// Window size in base pairs for parallel processing. Defaults to 16384 bp (one
    /// tabix/CSI index bin).
    #[arg(long)]
    pub window_size: Option<u64>,

    /// Keep only the samples listed in FILE (one per line), in that order;
    /// "^FILE" drops them instead. Applied before every other step, as
    /// `bcftools view -S`: INFO/AC and INFO/AN are recomputed from the kept
    /// genotypes, and removing every sample gives sites-only output.
    #[arg(long = "samples-file", short = 'S', value_name = "FILE")]
    pub samples_file: Option<String>,

    /// Set GT to missing for samples where the expression is NOT true
    /// (bcftools syntax, typically on FORMAT fields, e.g. 'FMT/DP>=3'). As
    /// `bcftools filter -S . -i`; a site-level expression masks all samples
    /// of failing sites. Applied after -S, before all other steps.
    #[arg(long = "include-genotypes", short = 'I', conflicts_with = "exclude_genotypes", value_name = "EXPR")]
    pub include_genotypes: Option<String>,

    /// Set GT to missing for samples where the expression is true, e.g.
    /// -E 'FMT/DP<3' (as `bcftools filter -S . -e`).
    #[arg(long = "exclude-genotypes", short = 'E', value_name = "EXPR")]
    pub exclude_genotypes: Option<String>,

    /// Exclude sites for which the expression is true (bcftools syntax).
    #[arg(long = "exclude", short = 'e', conflicts_with = "include")]
    pub exclude: Option<String>,

    /// Include only sites for which the expression is true (bcftools syntax).
    #[arg(long = "include", short = 'i')]
    pub include: Option<String>,

    /// Annotate the FILTER column of failing sites with NAME instead of
    /// dropping them. "+" picks a name automatically ("Filter1", ...).
    #[arg(long = "soft-filter", short = 's')]
    pub soft_filter: Option<String>,

    /// bcftools FILTER mode: "+" adds NAME to existing FILTER values instead
    /// of replacing them; "x" resets passing sites to PASS. May be combined.
    #[arg(long = "mode", short = 'm')]
    pub mode: Option<String>,

    /// Minimum allele count (inclusive), applied per ALT allele.
    #[arg(long)]
    pub min_ac: Option<u32>,

    /// Minimum allele frequency (inclusive), applied per ALT allele.
    #[arg(long)]
    pub min_af: Option<f64>,

    /// Recompute tags before evaluating the expression: comma-separated list
    /// (AN,AC,AF,MAF,NS,AC_Hom,AC_Het,AC_Hemi,HWE,ExcHet,F_MISSING,VAF,VAF1)
    /// or "all".
    #[arg(long = "fill-tags")]
    pub fill_tags: Option<String>,
}

/// Set failing samples' genotypes to missing, exactly as bcftools
/// vcffilter.c `set_genotypes` with `-S .`: every GT value (of any ploidy)
/// becomes an unphased missing allele, and INFO/AN and INFO/AC, if present
/// with the expected number of values, are decremented accordingly.
/// `masked(i)` says whether sample `i` is to be masked.
fn mask_genotypes(rec: &mut Record, masked: impl Fn(usize) -> bool) -> Result<()> {
    let n = hts::n_samples(rec.header());
    if n == 0 {
        return Ok(());
    }
    let Ok(mut gt) = hts::format_i32(rec, "GT") else { return Ok(()) };
    if gt.is_empty() {
        return Ok(());
    }
    let n_allele = rec.allele_count() as usize;
    let mut an = match hts::info_i32(rec, "AN") {
        Ok(v) if v.len() == 1 => Some(v[0]),
        _ => None,
    };
    let mut ac = match hts::info_i32(rec, "AC") {
        Ok(v) if v.len() == n_allele.saturating_sub(1) => Some(v),
        _ => None,
    };
    let width = gt.len() / n;
    let mut changed = false;
    for i in (0..n).filter(|&i| masked(i)) {
        for v in &mut gt[i * width..(i + 1) * width] {
            if *v == hts::INT32_VECTOR_END {
                break;
            }
            if !hts::gt_is_missing(*v) {
                let ial = hts::gt_allele(*v);
                if let Some(ac) = ac.as_mut()
                    && ial > 0
                    && (ial as usize) <= n_allele
                    && let Some(c) = ac.get_mut(ial as usize - 1)
                {
                    *c -= 1;
                }
                if let Some(an) = an.as_mut() {
                    *an -= 1;
                }
            }
            *v = 0; // bcf_gt_missing
            changed = true;
        }
    }
    if !changed {
        return Ok(());
    }
    hts::set_format_i32(rec, "GT", &gt)?;
    if let Some(an) = an {
        hts::set_info_i32(rec, "AN", &[an])?;
    }
    if let Some(ac) = ac {
        hts::set_info_i32(rec, "AC", &ac)?;
    }
    Ok(())
}

/// bcftools `filter` FILTER-column handling for one record (vcffilter.c).
fn apply_filter_column(rec: &mut Record, pass: bool, soft_id: Option<i32>, mode_add: bool, mode_reset: bool) -> bool {
    if pass {
        if mode_reset || hts::n_filters(rec) == 0 {
            hts::add_filter(rec, 0); // PASS
        }
        true
    } else if let Some(id) = soft_id {
        if mode_add {
            hts::add_filter(rec, id);
        } else {
            hts::set_filter(rec, id);
        }
        true
    } else {
        false
    }
}

pub fn run(args: FilterArgs) -> Result<()> {
    let mut mode_add = false;
    let mut mode_reset = false;
    for c in args.mode.as_deref().unwrap_or("").chars() {
        match c {
            '+' => mode_add = true,
            'x' => mode_reset = true,
            _ => bail!("the -m parameter \"{}\" is not recognised", args.mode.as_deref().unwrap_or("")),
        }
    }
    let tags = match &args.fill_tags {
        Some(spec) => filltags::parse_tags(spec)?,
        None => TagSet::default(),
    };
    let (expr, exclude) = match (&args.exclude, &args.include) {
        (Some(e), _) => (Some(e.clone()), true),
        (None, Some(i)) => (Some(i.clone()), false),
        (None, None) => (None, false),
    };
    let (gt_expr, gt_exclude) = match (&args.exclude_genotypes, &args.include_genotypes) {
        (Some(e), _) => (Some(e.clone()), true),
        (None, Some(i)) => (Some(i.clone()), false),
        (None, None) => (None, false),
    };
    let min_ac = args.min_ac;
    let min_af = args.min_af;

    // Soft filter name, resolved against the input header in `prepare`.
    let soft_name: Arc<std::sync::OnceLock<Option<String>>> = Arc::new(std::sync::OnceLock::new());

    let prepare = {
        let tags = tags.clone();
        let expr = expr.clone();
        let soft = args.soft_filter.clone();
        let soft_name = Arc::clone(&soft_name);
        move |h: &HeaderView| -> Result<Vec<String>> {
            let mut lines = tags.header_lines();
            let name = match (&soft, &expr) {
                (Some(s), Some(e)) => {
                    let name = if s == "+" {
                        (1..)
                            .map(|i| format!("Filter{i}"))
                            .find(|n| !hts::header_has(h, raw::BCF_HL_FLT, n))
                            .unwrap()
                    } else {
                        s.clone()
                    };
                    let desc = format!(
                        "Set if {}: {}",
                        if exclude { "true" } else { "not true" },
                        e.replace('"', "\\\"")
                    );
                    lines.push(format!("##FILTER=<ID={name},Description=\"{desc}\">"));
                    Some(name)
                }
                _ => None,
            };
            soft_name.set(name).map_err(|_| anyhow!("prepare called twice"))?;
            Ok(lines)
        }
    };

    let factory = move |h: &Arc<HeaderView>| -> Result<RecordFn> {
        let mut fill = if tags.is_empty() { None } else { Some(FillTags::new(&tags, h)?) };
        let mut filter = expr.as_deref().map(|e| hts::Filter::new(h, e)).transpose()?;
        let mut gt_filter = gt_expr.as_deref().map(|e| hts::Filter::new(h, e)).transpose()?;
        let soft_id = match soft_name.get().cloned().flatten() {
            Some(name) => Some(hts::filter_id(h, &name).ok_or_else(|| anyhow!("FILTER {name} missing from header"))?),
            None => None,
        };
        Ok(Box::new(move |rec: &mut Record| -> Result<bool> {
            if let Some(f) = gt_filter.as_mut() {
                let (site, per_sample) = f.test_samples(rec);
                match per_sample {
                    // sample-level expression: mask samples where it fails
                    Some(p) => mask_genotypes(rec, |i| p[i] == gt_exclude)?,
                    // site-level expression: mask every sample if the site fails
                    None if site == gt_exclude => mask_genotypes(rec, |_| true)?,
                    None => {}
                }
            }
            if (min_ac.is_some() || min_af.is_some()) && !crate::filter::filter_alleles_at_site(rec, min_af, min_ac)? {
                return Ok(false);
            }
            if let Some(fill) = fill.as_mut() {
                fill.fill(rec)?;
            }
            // FILTER-column handling (incl. setting PASS) is `bcftools filter`'s
            // and applies only when an expression is given.
            match filter.as_mut() {
                Some(f) => {
                    let pass = f.test(rec) != exclude;
                    Ok(apply_filter_column(rec, pass, soft_id, mode_add, mode_reset))
                }
                None => Ok(true),
            }
        }))
    };

    io::run_records(
        "filter",
        &args.indexed,
        &args.output,
        args.window_size,
        args.samples_file.map(|file| SampleSubset { file }),
        prepare,
        factory,
    )
}
