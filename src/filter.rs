//! Per-site allele filtering by minimum allele count and/or frequency.
//!
//! Counts alleles directly from per-sample genotype calls (the INFO/AC and
//! INFO/AF values in the input are *not* trusted), drops alleles that fail
//! either threshold, and rebuilds the site with htslib's
//! `bcf_remove_allele_set` (which re-slices INFO/FORMAT fields declared
//! Number=A/R/G and renumbers genotypes). Samples whose genotype touched a
//! removed allele are first set to fully missing (all FORMAT fields). The
//! recomputable summary fields (AC, AF, AN, F_MISSING) are then updated if
//! the header declares them.

use std::os::raw::{c_int, c_ulong};

use anyhow::{Result, bail};
use rust_htslib::bcf::Record;
use rust_htslib::htslib as raw;

use crate::hts::{self, HeaderView};

/// Apply minimum AC and AF thresholds independently — an allele is kept only
/// when it satisfies *both* thresholds (when both are supplied). When neither
/// is supplied, every allele is kept (the count fields are still recomputed).
///
/// Returns `Ok(false)` when no ALT allele passes — the entire site should be
/// dropped from the output.
pub fn filter_alleles_at_site(rec: &mut Record, min_af: Option<f64>, min_ac: Option<u32>) -> Result<bool> {
    let n_total = rec.allele_count() as usize;
    let n_alt = n_total.saturating_sub(1);

    // Sites with no ALT alleles can't be filtered, but we still refresh the
    // per-site count fields so the output is internally consistent.
    if n_alt == 0 {
        let counts = count_alleles(rec, n_total.max(1));
        update_count_info_fields(rec, &counts, 0)?;
        return Ok(true);
    }

    let counts = count_alleles(rec, n_total);
    let keep: Vec<bool> = (0..n_alt)
        .map(|i| {
            let ac = counts.allele_counts[i + 1];
            let af = if counts.n_called_alleles > 0 {
                ac as f64 / counts.n_called_alleles as f64
            } else {
                0.0
            };
            min_ac.is_none_or(|m| ac >= m) && min_af.is_none_or(|m| af >= m)
        })
        .collect();
    let n_kept = keep.iter().filter(|&&k| k).count();
    if n_kept == 0 {
        return Ok(false);
    }
    if n_kept == n_alt {
        update_count_info_fields(rec, &counts, n_alt)?;
        return Ok(true);
    }

    // removed[allele_index]
    let mut removed = vec![false; n_total];
    for (i, &k) in keep.iter().enumerate() {
        removed[i + 1] = !k;
    }
    blank_samples_with_removed_alleles(rec, &removed)?;
    remove_alleles(rec, &removed)?;

    let n_total_new = n_kept + 1;
    let counts = count_alleles(rec, n_total_new);
    update_count_info_fields(rec, &counts, n_kept)?;
    Ok(true)
}

#[derive(Debug)]
struct AlleleCounts {
    /// Called allele instances per allele index (0 = REF).
    allele_counts: Vec<u32>,
    /// Total number of called allele instances across all samples.
    n_called_alleles: u32,
    n_samples: u32,
    /// Samples with at least one missing allele in their genotype.
    n_samples_with_missing: u32,
}

fn count_alleles(rec: &Record, n_total: usize) -> AlleleCounts {
    let n_samples = hts::n_samples(rec.header());
    let mut allele_counts = vec![0u32; n_total];
    let mut n_called_alleles = 0u32;
    let mut n_samples_with_missing = 0u32;
    match hts::format_i32(rec, "GT") {
        Ok(gt) if n_samples > 0 && !gt.is_empty() => {
            let width = gt.len() / n_samples;
            for s in 0..n_samples {
                let mut missing = false;
                for &v in &gt[s * width..(s + 1) * width] {
                    if v == hts::INT32_VECTOR_END {
                        break;
                    }
                    if hts::gt_is_missing(v) {
                        missing = true;
                        continue;
                    }
                    let p = hts::gt_allele(v) as usize;
                    if p < n_total {
                        allele_counts[p] += 1;
                        n_called_alleles += 1;
                    } else {
                        missing = true;
                    }
                }
                if missing {
                    n_samples_with_missing += 1;
                }
            }
        }
        _ => n_samples_with_missing = n_samples as u32,
    }
    AlleleCounts { allele_counts, n_called_alleles, n_samples: n_samples as u32, n_samples_with_missing }
}

/// Set every FORMAT value of samples whose genotype includes a removed
/// allele to missing; their GT keeps its ploidy with all alleles missing.
fn blank_samples_with_removed_alleles(rec: &mut Record, removed: &[bool]) -> Result<()> {
    let n_samples = hts::n_samples(rec.header());
    let Ok(gt) = hts::format_i32(rec, "GT") else { return Ok(()) };
    if n_samples == 0 || gt.is_empty() {
        return Ok(());
    }
    let width = gt.len() / n_samples;
    let lost: Vec<usize> = (0..n_samples)
        .filter(|&s| {
            gt[s * width..(s + 1) * width].iter().take_while(|&&v| v != hts::INT32_VECTOR_END).any(|&v| {
                !hts::gt_is_missing(v) && removed.get(hts::gt_allele(v) as usize).copied().unwrap_or(false)
            })
        })
        .collect();
    if lost.is_empty() {
        return Ok(());
    }
    let gt_id = gt_header_id(rec.header());
    let line = hts::rec_ptr(rec);
    unsafe {
        let n_fmt = (*line).n_fmt() as usize;
        for f in 0..n_fmt {
            let fmt = &mut *(*line).d.fmt.add(f);
            if fmt.p.is_null() || fmt.n <= 0 {
                continue;
            }
            let is_gt = Some(fmt.id) == gt_id;
            for &s in &lost {
                let base = fmt.p.add(s * fmt.size as usize);
                blank_sample_field(base, fmt.type_, fmt.n as usize, is_gt)?;
            }
        }
        (*line).d.indiv_dirty = 1;
    }
    Ok(())
}

fn gt_header_id(h: &HeaderView) -> Option<c_int> {
    let c = std::ffi::CString::new("GT").ok()?;
    let id = unsafe { raw::bcf_hdr_id2int(h.as_ptr(), raw::BCF_DT_ID as c_int, c.as_ptr()) };
    (id >= 0).then_some(id)
}

/// Write missing values into one sample's slot of a FORMAT field.
unsafe fn blank_sample_field(base: *mut u8, ty: c_int, n: usize, is_gt: bool) -> Result<()> {
    // BCF_BT_* type codes
    const INT8: c_int = 1;
    const INT16: c_int = 2;
    const INT32: c_int = 3;
    const FLOAT: c_int = 5;
    const CHAR: c_int = 7;
    unsafe {
        match ty {
            INT8 | INT16 | INT32 => {
                let (missing, vend, width): (i64, i64, usize) = match ty {
                    INT8 => (i8::MIN as i64, i8::MIN as i64 + 1, 1),
                    INT16 => (i16::MIN as i64, i16::MIN as i64 + 1, 2),
                    _ => (i32::MIN as i64, i32::MIN as i64 + 1, 4),
                };
                let read = |j: usize| -> i64 {
                    let p = base.add(j * width);
                    match width {
                        1 => (p as *const i8).read_unaligned() as i64,
                        2 => (p as *const i16).read_unaligned() as i64,
                        _ => (p as *const i32).read_unaligned() as i64,
                    }
                };
                let write = |j: usize, v: i64| {
                    let p = base.add(j * width);
                    match width {
                        1 => (p as *mut i8).write_unaligned(v as i8),
                        2 => (p as *mut i16).write_unaligned(v as i16),
                        _ => (p as *mut i32).write_unaligned(v as i32),
                    }
                };
                if is_gt {
                    // keep ploidy; every allele becomes unphased missing (0)
                    for j in 0..n {
                        if read(j) == vend {
                            break;
                        }
                        write(j, 0);
                    }
                } else {
                    write(0, missing);
                    for j in 1..n {
                        write(j, vend);
                    }
                }
            }
            FLOAT => {
                let p = base as *mut u32;
                p.write_unaligned(0x7F80_0001);
                for j in 1..n {
                    p.add(j).write_unaligned(0x7F80_0002);
                }
            }
            CHAR => {
                *base = b'.';
                for j in 1..n {
                    *base.add(j) = 0;
                }
            }
            _ => bail!("unexpected FORMAT value type {ty}"),
        }
    }
    Ok(())
}

/// `bcf_remove_allele_set` with a bitset of the alleles to remove.
fn remove_alleles(rec: &mut Record, removed: &[bool]) -> Result<()> {
    // kbitset_t { size_t n, n_max; unsigned long b[n + 1]; } — only
    // kbs_exists() is used on it.
    const BITS: usize = c_ulong::BITS as usize;
    let nwords = removed.len().div_ceil(BITS);
    let mut buf: Vec<c_ulong> = vec![0; 2 + nwords + 1];
    buf[0] = nwords as c_ulong;
    buf[1] = nwords as c_ulong;
    for (i, &r) in removed.iter().enumerate() {
        if r {
            buf[2 + i / BITS] |= 1 << (i % BITS);
        }
    }
    buf[2 + nwords] = c_ulong::MAX; // sentinel, as kbs_init sets b[n] non-zero
    let set = buf.as_ptr() as *const raw::kbitset_t;
    let h = unsafe { rec.header().as_ptr() };
    let r = unsafe { raw::bcf_remove_allele_set(h, hts::rec_ptr(rec), set) };
    if r != 0 {
        bail!("failed to remove alleles at {}", hts::rec_locus(rec));
    }
    Ok(())
}

fn update_count_info_fields(rec: &mut Record, counts: &AlleleCounts, n_alt: usize) -> Result<()> {
    let has = |rec: &Record, id: &str| hts::header_has(rec.header(), raw::BCF_HL_INFO, id);
    if has(rec, "AN") {
        hts::set_info_i32(rec, "AN", &[counts.n_called_alleles as i32])?;
    }
    if has(rec, "AC") && n_alt > 0 {
        let ac: Vec<i32> = (1..=n_alt).map(|i| counts.allele_counts[i] as i32).collect();
        hts::set_info_i32(rec, "AC", &ac)?;
    }
    if has(rec, "AF") && n_alt > 0 {
        let total = counts.n_called_alleles as f64;
        let af: Vec<f32> = (1..=n_alt)
            .map(|i| if total > 0.0 { (counts.allele_counts[i] as f64 / total) as f32 } else { 0.0 })
            .collect();
        hts::set_info_f32(rec, "AF", &af)?;
    }
    if has(rec, "F_MISSING") {
        let f = if counts.n_samples > 0 {
            counts.n_samples_with_missing as f32 / counts.n_samples as f32
        } else {
            0.0
        };
        hts::set_info_f32(rec, "F_MISSING", &[f])?;
    }
    Ok(())
}
