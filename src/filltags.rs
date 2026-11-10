//! Recompute per-site INFO/FORMAT tags from genotypes, mirroring
//! `bcftools +fill-tags` (bcftools 1.23 `plugins/fill-tags.c`).
//!
//! The logic is a line-by-line port operating on the same htslib data
//! (`bcf_get_*` / `bcf_update_*`), so values, update order and quirks match
//! bcftools:
//!
//! * per sample, each *distinct* allele in the genotype is counted once; a
//!   homozygous genotype contributes 2 to `AC_Hom` regardless of ploidy, a
//!   heterozygous genotype contributes 1 to `AC_Het` for each distinct allele,
//!   a haploid call (or a half-missing call such as `./1`) is hemizygous;
//! * `AN`/`AC`/`AF`/`MAF` are derived from those counts;
//! * `HWE`/`ExcHet` use the Wigginton et al. (2005) exact test per ALT allele;
//! * `F_MISSING` is evaluated with the bcftools expression engine, exactly as
//!   `+fill-tags` does;
//! * `VAF`/`VAF1` are FORMAT tags computed from `FORMAT/AD`.

use std::sync::Arc;

use anyhow::{Result, anyhow, bail};
use rust_htslib::bcf::Record;

use crate::hts::{self, HeaderView};

/// A tag that can be (re)computed.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum Tag {
    An,
    Ac,
    Ns,
    AcHom,
    AcHet,
    AcHemi,
    Af,
    Maf,
    Hwe,
    ExcHet,
    FMissing,
    /// FORMAT/VAF: per-sample ALT read fractions from FORMAT/AD.
    Vaf,
    /// FORMAT/VAF1: per-sample summed-ALT read fraction from FORMAT/AD.
    Vaf1,
}

impl Tag {
    /// Every supported tag, in bcftools' header-declaration order.
    pub const ALL: [Tag; 13] = [
        Tag::FMissing,
        Tag::An,
        Tag::Ac,
        Tag::Ns,
        Tag::AcHom,
        Tag::AcHet,
        Tag::AcHemi,
        Tag::Af,
        Tag::Maf,
        Tag::Hwe,
        Tag::ExcHet,
        Tag::Vaf,
        Tag::Vaf1,
    ];

    pub fn id(self) -> &'static str {
        match self {
            Tag::An => "AN",
            Tag::Ac => "AC",
            Tag::Ns => "NS",
            Tag::AcHom => "AC_Hom",
            Tag::AcHet => "AC_Het",
            Tag::AcHemi => "AC_Hemi",
            Tag::Af => "AF",
            Tag::Maf => "MAF",
            Tag::Hwe => "HWE",
            Tag::ExcHet => "ExcHet",
            Tag::FMissing => "F_MISSING",
            Tag::Vaf => "VAF",
            Tag::Vaf1 => "VAF1",
        }
    }

    /// The header line bcftools appends for this tag.
    fn header_line(self) -> &'static str {
        match self {
            Tag::An => r#"##INFO=<ID=AN,Number=1,Type=Integer,Description="Total number of alleles in called genotypes">"#,
            Tag::Ac => r#"##INFO=<ID=AC,Number=A,Type=Integer,Description="Allele count in genotypes">"#,
            Tag::Ns => r#"##INFO=<ID=NS,Number=1,Type=Integer,Description="Number of samples with data">"#,
            Tag::AcHom => r#"##INFO=<ID=AC_Hom,Number=A,Type=Integer,Description="Allele counts in homozygous genotypes">"#,
            Tag::AcHet => r#"##INFO=<ID=AC_Het,Number=A,Type=Integer,Description="Allele counts in heterozygous genotypes">"#,
            Tag::AcHemi => r#"##INFO=<ID=AC_Hemi,Number=A,Type=Integer,Description="Allele counts in hemizygous genotypes">"#,
            Tag::Af => r#"##INFO=<ID=AF,Number=A,Type=Float,Description="Allele frequency">"#,
            Tag::Maf => r#"##INFO=<ID=MAF,Number=1,Type=Float,Description="Frequency of the second most common allele">"#,
            Tag::Hwe => r#"##INFO=<ID=HWE,Number=A,Type=Float,Description="HWE test (PMID:15789306); 1=good, 0=bad">"#,
            Tag::ExcHet => r#"##INFO=<ID=ExcHet,Number=A,Type=Float,Description="Test excess heterozygosity; 1=good, 0=bad">"#,
            Tag::FMissing => r#"##INFO=<ID=F_MISSING,Number=1,Type=Float,Description="Added by +fill-tags expression F_MISSING:1=F_MISSING">"#,
            Tag::Vaf => r#"##FORMAT=<ID=VAF,Number=A,Type=Float,Description="The fraction of reads with alternate allele (nALT/nSumAll)">"#,
            Tag::Vaf1 => r#"##FORMAT=<ID=VAF1,Number=1,Type=Float,Description="The fraction of reads with alternate alleles (nSumALT/nSumAll)">"#,
        }
    }
}

/// Parsed set of tags to fill, plus options.
#[derive(Clone, Debug, Default)]
pub struct TagSet {
    tags: Vec<Tag>,
    /// bcftools `-d/--drop-missing`: do not count half-missing genotypes
    /// (e.g. `./1`) as hemizygous.
    pub drop_missing: bool,
}

impl TagSet {
    pub fn new(tags: Vec<Tag>) -> Self {
        Self { tags, drop_missing: false }
    }

    pub fn is_empty(&self) -> bool {
        self.tags.is_empty()
    }

    pub fn contains(&self, tag: Tag) -> bool {
        self.tags.contains(&tag)
    }

    pub fn tags(&self) -> &[Tag] {
        &self.tags
    }

    fn needs_counts(&self) -> bool {
        self.tags.iter().any(|t| !matches!(t, Tag::FMissing | Tag::Vaf | Tag::Vaf1))
    }

    /// Header lines to append (bcftools order; existing IDs are kept as they
    /// are by htslib).
    pub fn header_lines(&self) -> Vec<String> {
        self.tags.iter().map(|t| t.header_line().to_string()).collect()
    }
}

/// Parse a bcftools-style tag list, e.g. `AN,AC,AF` or `all`. Names are
/// case-insensitive and may carry an `INFO/` or `FORMAT/` prefix.
pub fn parse_tags(spec: &str) -> Result<TagSet> {
    let mut tags = Vec::new();
    let push = |t: Tag, tags: &mut Vec<Tag>| {
        if !tags.contains(&t) {
            tags.push(t);
        }
    };
    for raw in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
        if raw.eq_ignore_ascii_case("all") {
            for t in Tag::ALL {
                push(t, &mut tags);
            }
            continue;
        }
        let name = if raw.len() > 5 && raw[..5].eq_ignore_ascii_case("INFO/") {
            &raw[5..]
        } else if raw.len() > 7 && raw[..7].eq_ignore_ascii_case("FORMAT/") {
            &raw[7..]
        } else {
            raw
        };
        let tag = Tag::ALL
            .iter()
            .copied()
            .find(|t| t.id().eq_ignore_ascii_case(name))
            .ok_or_else(|| {
                anyhow!(
                    "unsupported tag \"{raw}\" (supported: {}, or \"all\")",
                    Tag::ALL.iter().map(|t| t.id()).collect::<Vec<_>>().join(",")
                )
            })?;
        push(tag, &mut tags);
    }
    // Keep bcftools' declaration order regardless of the order given.
    tags.sort_by_key(|t| Tag::ALL.iter().position(|x| x == t));
    Ok(TagSet::new(tags))
}

/// Per-allele genotype counts, as accumulated by bcftools `set_counts`.
#[derive(Clone, Copy, Default, Debug)]
struct Counts {
    nhet: i32,
    nhom: i32,
    nhemi: i32,
    nac: i32,
}

impl Counts {
    fn total(&self) -> i32 {
        self.nhet + self.nhom + self.nhemi + self.nac
    }
}

/// Per-worker tag filler (holds the compiled F_MISSING expression, which is
/// bound to the worker's header).
pub struct FillTags {
    tags: TagSet,
    n_samples: usize,
    f_missing: Option<hts::Filter>,
    warned_ad: bool,
}

impl FillTags {
    /// `header` must already contain [`TagSet::header_lines`].
    pub fn new(tags: &TagSet, header: &Arc<HeaderView>) -> Result<Self> {
        let f_missing = if tags.contains(Tag::FMissing) {
            Some(hts::Filter::new(header, "F_MISSING")?)
        } else {
            None
        };
        Ok(Self { tags: tags.clone(), n_samples: hts::n_samples(header), f_missing, warned_ad: false })
    }

    /// Fill the requested tags on one record (bcftools `process()`).
    pub fn fill(&mut self, rec: &mut Record) -> Result<()> {
        if let Some(f) = self.f_missing.as_mut() {
            // ftf_filter_expr with Number=1, Type=Float
            let v = f.values(rec);
            let x = match v.first() {
                Some(d) if !is_double_missing(*d) => *d as f32,
                _ => hts::float_missing(),
            };
            hts::set_info_f32(rec, "F_MISSING", &[x])?;
        }
        if self.tags.needs_counts() || self.tags.contains(Tag::Vaf) || self.tags.contains(Tag::Vaf1) {
            if self.tags.needs_counts() {
                self.process_fmt(rec)?;
                self.process_info_af(rec)?;
            }
            self.process_vaf_vaf1(rec)?;
        }
        Ok(())
    }

    fn process_fmt(&mut self, rec: &mut Record) -> Result<()> {
        let Ok(gt) = hts::format_i32(rec, "GT") else {
            return Ok(()); // no GT tag
        };
        let n_allele = rec.allele_count() as usize;
        let nsmpl = self.n_samples;
        if nsmpl == 0 {
            return Ok(());
        }
        let width = gt.len() / nsmpl;
        let tags = &self.tags;

        let mut counts = vec![Counts::default(); n_allele];
        let mut ns: i32 = 0;
        let mut seen = vec![false; n_allele];
        let mut distinct: Vec<usize> = Vec::with_capacity(4);
        for i in 0..nsmpl {
            let p = &gt[i * width..(i + 1) * width];
            let mut nals = 0usize;
            let mut ial = 0usize;
            distinct.clear();
            while ial < width {
                let v = p[ial];
                if v == hts::INT32_VECTOR_END {
                    break; // smaller ploidy
                }
                ial += 1;
                if hts::gt_is_missing(v) {
                    continue; // missing allele
                }
                let idx = hts::gt_allele(v) as usize;
                nals += 1;
                if idx >= n_allele {
                    bail!("incorrect allele (\"{idx}\") at {}", hts::rec_locus(rec));
                }
                if !seen[idx] {
                    seen[idx] = true;
                    distinct.push(idx);
                }
            }
            for &k in &distinct {
                seen[k] = false;
            }
            if nals == 0 {
                continue; // missing genotype
            }
            let is_hom = distinct.len() == 1;
            let (is_hemi, is_half) = if nals != ial {
                if tags.drop_missing { (false, true) } else { (true, false) }
            } else if nals == 1 {
                (true, false)
            } else {
                (false, false)
            };
            for &k in &distinct {
                let c = &mut counts[k];
                if is_half {
                    c.nac += 1;
                } else if !is_hom {
                    c.nhet += 1;
                } else if !is_hemi {
                    c.nhom += 2;
                } else {
                    c.nhemi += 1;
                }
            }
            ns += 1;
        }

        if tags.contains(Tag::Ns) {
            hts::set_info_i32(rec, "NS", &[ns])?;
        }
        if tags.contains(Tag::An) {
            let an: i32 = counts.iter().map(Counts::total).sum();
            hts::set_info_i32(rec, "AN", &[an])?;
        }
        if tags.contains(Tag::Af) || tags.contains(Tag::Maf) {
            let mut farr = vec![0f32; n_allele];
            let mut an: i32 = 0;
            if n_allele > 1 {
                for j in 0..n_allele {
                    farr[j] = counts[j].total() as f32;
                    an += farr[j] as i32;
                }
                if an != 0 {
                    for f in farr.iter_mut() {
                        *f /= an as f32;
                    }
                } else {
                    farr.iter_mut().for_each(|f| *f = hts::float_missing());
                }
            }
            if tags.contains(Tag::Af) {
                hts::set_info_f32(rec, "AF", &farr[1.min(farr.len())..n_allele.max(1)])?;
            }
            if n_allele > 1 && tags.contains(Tag::Maf) {
                if an != 0 {
                    farr.sort_by(|a, b| b.partial_cmp(a).unwrap_or(std::cmp::Ordering::Equal));
                }
                hts::set_info_f32(rec, "MAF", &farr[1..2])?;
            }
        }
        if tags.contains(Tag::Ac) {
            let v: Vec<i32> = (1..n_allele).map(|j| counts[j].total()).collect();
            hts::set_info_i32(rec, "AC", &v)?;
        }
        if tags.contains(Tag::AcHet) {
            let v: Vec<i32> = (1..n_allele).map(|j| counts[j].nhet).collect();
            hts::set_info_i32(rec, "AC_Het", &v)?;
        }
        if tags.contains(Tag::AcHom) {
            let v: Vec<i32> = (1..n_allele).map(|j| counts[j].nhom).collect();
            hts::set_info_i32(rec, "AC_Hom", &v)?;
        }
        if tags.contains(Tag::AcHemi) && n_allele > 1 {
            let v: Vec<i32> = (1..n_allele).map(|j| counts[j].nhemi).collect();
            hts::set_info_i32(rec, "AC_Hemi", &v)?;
        }
        if tags.contains(Tag::Hwe) || tags.contains(Tag::ExcHet) {
            let n_alt = n_allele.saturating_sub(1);
            let mut fhwe = vec![0f32; n_alt];
            let mut fexc = vec![0f32; n_alt];
            if n_allele > 1 {
                let mut nref_tot = counts[0].nhom;
                for c in &counts {
                    nref_tot += c.nhet; // NB: neglects multiallelic genotypes, as bcftools
                }
                for j in 1..n_allele {
                    let nref = nref_tot - counts[j].nhet;
                    let nalt = counts[j].nhet + counts[j].nhom;
                    let nhet = counts[j].nhet;
                    if nref > 0 && nalt > 0 {
                        let (h, e) = calc_hwe(nref, nalt, nhet).map_err(|e| anyhow!("{e} at {}", hts::rec_locus(rec)))?;
                        fhwe[j - 1] = h;
                        fexc[j - 1] = e;
                    } else {
                        fhwe[j - 1] = 1.0;
                        fexc[j - 1] = 1.0;
                    }
                }
            }
            if tags.contains(Tag::Hwe) {
                hts::set_info_f32(rec, "HWE", &fhwe)?;
            }
            if tags.contains(Tag::ExcHet) {
                hts::set_info_f32(rec, "ExcHet", &fexc)?;
            }
        }
        Ok(())
    }

    /// Sites-only VCFs: AF from INFO/AN and INFO/AC.
    fn process_info_af(&mut self, rec: &mut Record) -> Result<()> {
        if !self.tags.contains(Tag::Af) || self.n_samples > 0 {
            return Ok(());
        }
        let Ok(an) = hts::info_i32(rec, "AN") else { return Ok(()) };
        if an.len() != 1 || an[0] == 0 {
            return Ok(());
        }
        let an = an[0];
        let Ok(ac) = hts::info_i32(rec, "AC") else { return Ok(()) };
        if ac.len() != rec.allele_count() as usize - 1 {
            return Ok(());
        }
        let af: Vec<f32> = ac.iter().map(|&c| (c as f64 / an as f64) as f32).collect();
        hts::set_info_f32(rec, "AF", &af)
    }

    fn process_vaf_vaf1(&mut self, rec: &mut Record) -> Result<()> {
        let (want_vaf, want_vaf1) = (self.tags.contains(Tag::Vaf), self.tags.contains(Tag::Vaf1));
        if !want_vaf && !want_vaf1 {
            return Ok(());
        }
        let n_allele = rec.allele_count() as usize;
        if n_allele <= 1 {
            return Ok(());
        }
        let ad = match hts::format_i32(rec, "AD") {
            Ok(v) if !v.is_empty() => v,
            _ => {
                if !self.warned_ad {
                    eprintln!(
                        "Warning: cannot add the VAF/VAF1 annotations, the required FORMAT/AD tag is missing at {}.\n         (This warning is printed only once.)",
                        hts::rec_locus(rec)
                    );
                    self.warned_ad = true;
                }
                return Ok(());
            }
        };
        let nsmpl = self.n_samples;
        if nsmpl == 0 || ad.len() != nsmpl * n_allele {
            return Ok(()); // incorrect number of values (possibly all missing)
        }
        let nval = ad.len() / nsmpl;
        for (mode_vaf, wanted) in [(true, want_vaf), (false, want_vaf1)] {
            if !wanted {
                continue;
            }
            let nval1 = if mode_vaf { n_allele - 1 } else { 1 };
            let mut farr = vec![0f32; nval1 * nsmpl];
            for i in 0..nsmpl {
                let src = &ad[i * nval..(i + 1) * nval];
                let dst = &mut farr[i * nval1..(i + 1) * nval1];
                let mut sum = 0f32;
                let mut j = 0;
                while j < nval {
                    if src[j] == hts::INT32_MISSING || src[j] == hts::INT32_VECTOR_END {
                        break;
                    }
                    sum += src[j] as f32;
                    j += 1;
                }
                if j != nval {
                    dst[0] = hts::float_missing();
                    for d in dst.iter_mut().skip(1) {
                        *d = hts::float_vector_end();
                    }
                    continue;
                }
                if !mode_vaf {
                    dst[0] = if sum != 0.0 { (sum - src[0] as f32) / sum } else { 0.0 };
                    continue;
                }
                for (k, d) in dst.iter_mut().enumerate() {
                    *d = if sum != 0.0 { src[k + 1] as f32 / sum } else { 0.0 };
                }
            }
            hts::set_format_f32(rec, if mode_vaf { "VAF" } else { "VAF1" }, &farr)?;
        }
        Ok(())
    }
}

/// `bcf_double_is_missing`
fn is_double_missing(d: f64) -> bool {
    d.to_bits() == 0x7FF0_0000_0000_0001
}

/// Wigginton et al. 2005 (PMID 15789306), as implemented in bcftools
/// `calc_hwe`. Returns `(p_hwe, p_exc_het)`.
pub fn calc_hwe(nref: i32, nalt: i32, nhet: i32) -> Result<(f32, f32)> {
    let ngt = (nref + nalt) / 2;
    let nrare = nref.min(nalt);
    if (nrare & 1) ^ (nhet & 1) != 0 {
        bail!("nrare/nhet should be both odd or even: nrare={nrare} nref={nref} nalt={nalt} nhet={nhet}");
    }
    if nrare < nhet {
        bail!("fewer rare alleles than hets? nrare={nrare} nref={nref} nalt={nalt} nhet={nhet}");
    }
    if (nref + nalt) & 1 != 0 {
        bail!("expected diploid genotypes: nref={nref} nalt={nalt}");
    }
    let nrare_u = nrare as usize;
    let mut probs = vec![0f64; nrare_u + 1];

    let mut mid = (nrare as f64 * (nref + nalt - nrare) as f64 / (nref + nalt) as f64) as i32;
    if (nrare & 1) ^ (mid & 1) != 0 {
        mid += 1;
    }

    let mut hom_r = (nrare - mid) / 2;
    let mut hom_c = ngt - mid - hom_r;
    probs[mid as usize] = 1.0;
    let mut sum = 1.0f64;

    let mut het = mid;
    while het > 1 {
        let h = het as usize;
        probs[h - 2] = probs[h] * het as f64 * (het as f64 - 1.0)
            / (4.0 * (hom_r as f64 + 1.0) * (hom_c as f64 + 1.0));
        sum += probs[h - 2];
        hom_r += 1;
        hom_c += 1;
        het -= 2;
    }

    hom_r = (nrare - mid) / 2;
    hom_c = ngt - mid - hom_r;
    het = mid;
    while het <= nrare - 2 {
        let h = het as usize;
        probs[h + 2] = probs[h] * 4.0 * hom_r as f64 * hom_c as f64
            / ((het as f64 + 2.0) * (het as f64 + 1.0));
        sum += probs[h + 2];
        hom_r -= 1;
        hom_c -= 1;
        het += 2;
    }

    for p in probs.iter_mut() {
        *p /= sum;
    }

    let nhet_u = nhet as usize;
    let mut prob = probs[nhet_u];
    for p in &probs[nhet_u + 1..] {
        prob += p;
    }
    let p_exc_het = prob;

    let mut prob = 0.0;
    for &p in &probs {
        if p > probs[nhet_u] {
            continue;
        }
        prob += p;
    }
    if prob > 1.0 {
        prob = 1.0;
    }
    Ok((prob as f32, p_exc_het as f32))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn parse_all_and_subset() {
        let t = parse_tags("all").unwrap();
        assert_eq!(t.tags().len(), Tag::ALL.len());
        let t = parse_tags("INFO/af,AN").unwrap();
        assert_eq!(t.tags(), &[Tag::An, Tag::Af]);
        assert!(parse_tags("END").is_err());
        assert_eq!(parse_tags("FORMAT/VAF").unwrap().tags(), &[Tag::Vaf]);
    }

    #[test]
    fn hwe_perfect_and_excess_het() {
        // 1 AA, 2 AB, 1 BB: nref=4, nalt=4, nhet=2
        let (h, e) = calc_hwe(4, 4, 2).unwrap();
        assert!(h > 0.99 && h <= 1.0, "{h}");
        assert!(e > 0.0 && e < 1.0, "{e}");
        // All hets: 4 AB: nref=4 nalt=4 nhet=4
        let (h2, e2) = calc_hwe(4, 4, 4).unwrap();
        assert!(h2 < h);
        assert!(e2 < e);
    }
}
