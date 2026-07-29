//! Per-site allele filtering by minimum allele count and/or frequency.
//!
//! Counts alleles directly from per-sample genotype calls (the INFO/AC and
//! INFO/AF values in the input are *not* trusted), drops alleles that fail
//! either threshold, rebuilds the site, sets every genotype that touched a
//! removed allele to fully missing, and updates per-record INFO/FORMAT fields
//! whose definitions in the header are dependent on the allele list (Number=A,
//! R, G) plus the recomputable summary fields (AC, AF, AN, F_MISSING).

use std::collections::HashMap;

use noodles_vcf as vcf;
use vcf::Header;
use vcf::header::record::value::map::format::Number as FormatNumber;
use vcf::header::record::value::map::info::Number as InfoNumber;
use vcf::variant::RecordBuf;
use vcf::variant::record_buf::info::field::Value as InfoValue;
use vcf::variant::record_buf::info::field::value::Array as InfoArray;
use vcf::variant::record::Samples as _;
use vcf::variant::record_buf::samples::Samples;
use vcf::variant::record_buf::samples::sample::Value as SampleValue;
use vcf::variant::record_buf::samples::sample::value::Array as SampleArray;
use vcf::variant::record_buf::samples::sample::value::Genotype;
use vcf::variant::record_buf::samples::sample::value::genotype::Allele;


/// Apply minimum AC and AF thresholds independently — an allele is kept only
/// when it satisfies *both* thresholds (when both are supplied). When neither
/// is supplied, every allele is kept (the function still recomputes count
/// fields).
///
/// Returns `None` when no ALT allele passes — the entire site should be
/// dropped from the output.
pub fn filter_alleles_at_site(
    mut buf: RecordBuf,
    header: &Header,
    min_af: Option<f64>,
    min_ac: Option<u32>,
) -> Option<RecordBuf> {
    let n_alt = buf.alternate_bases().as_ref().len();

    // Sites with no ALT alleles can't be filtered, but we still refresh the
    // per-site count fields so the output is internally consistent.
    if n_alt == 0 {
        let counts = count_alleles(&buf, 1);
        update_count_info_fields(&mut buf, header, &counts, n_alt);
        return Some(buf);
    }

    let n_total_old = n_alt + 1;
    let counts_pre = count_alleles(&buf, n_total_old);

    let kept_alt_idx: Vec<usize> = (0..n_alt)
        .filter(|&i| {
            let ac = counts_pre.allele_counts[i + 1];
            let af = if counts_pre.n_called_alleles > 0 {
                ac as f64 / counts_pre.n_called_alleles as f64
            } else {
                0.0
            };
            let pass_ac = min_ac.is_none_or(|m| ac >= m);
            let pass_af = min_af.is_none_or(|m| af >= m);
            pass_ac && pass_af
        })
        .collect();

    if kept_alt_idx.is_empty() {
        return None;
    }

    if kept_alt_idx.len() < n_alt {
        // At least one allele was removed — rebuild the site.
        let n_kept_alt = kept_alt_idx.len();
        let n_total_new = n_kept_alt + 1;

        // remap[orig_idx] = Some(new_idx) when kept, None when removed.
        let mut remap: Vec<Option<usize>> = vec![None; n_total_old];
        remap[0] = Some(0);
        for (new_pos, &orig_alt) in kept_alt_idx.iter().enumerate() {
            remap[orig_alt + 1] = Some(new_pos + 1);
        }

        // inv_remap[new_idx] = orig_idx — needed for slicing G fields.
        let mut inv_remap: Vec<usize> = vec![0; n_total_new];
        for (orig, slot) in remap.iter().enumerate() {
            if let Some(new) = *slot {
                inv_remap[new] = orig;
            }
        }

        // Slice ALT.
        let new_alts: Vec<String> = kept_alt_idx
            .iter()
            .map(|&i| buf.alternate_bases().as_ref()[i].clone())
            .collect();
        *buf.alternate_bases_mut().as_mut() = new_alts;

        let format_numbers: Vec<FormatNumber> = buf.samples().keys().as_ref().iter()
            .map(|key| {
                header.formats().get(key)
                    .map(|m| m.number())
                    .unwrap_or(FormatNumber::Count(1))
            })
            .collect();

        rebuild_samples(
            &mut buf,
            &remap,
            &inv_remap,
            n_total_old,
            n_total_new,
            &format_numbers,
        );
        rebuild_info(
            &mut buf,
            header,
            &kept_alt_idx,
            &inv_remap,
            n_total_old,
            n_total_new,
        );

        let counts_post = count_alleles(&buf, n_total_new);
        update_count_info_fields(&mut buf, header, &counts_post, n_kept_alt);
    } else {
        update_count_info_fields(&mut buf, header, &counts_pre, n_alt);
    }

    Some(buf)
}

#[derive(Debug)]
struct AlleleCounts {
    /// `allele_counts[i]` = number of called allele *instances* for original
    /// allele index `i` (0 = REF). Sum of all entries equals
    /// `n_called_alleles`.
    allele_counts: Vec<u32>,
    /// Total number of called allele instances across all samples (i.e.
    /// genotype slots whose value was a known allele index).
    n_called_alleles: u32,
    /// Number of samples in the record.
    n_samples: u32,
    /// Number of samples with at least one missing allele in their genotype.
    n_samples_with_missing: u32,
}

fn count_alleles(buf: &RecordBuf, n_total: usize) -> AlleleCounts {
    let samples = buf.samples();
    let gt_idx = samples.keys().as_ref().get_index_of("GT");
    let mut allele_counts = vec![0u32; n_total];
    let mut n_called_alleles = 0u32;
    let mut n_samples_with_missing = 0u32;
    let n_samples = samples.len() as u32;

    for sample in samples.values() {
        let mut sample_has_missing = false;

        if let Some(idx) = gt_idx {
            match sample.values().get(idx).and_then(|v| v.as_ref()) {
                Some(SampleValue::Genotype(g)) => {
                    for allele in g.as_ref() {
                        match allele.position() {
                            Some(p) if p < n_total => {
                                allele_counts[p] += 1;
                                n_called_alleles += 1;
                            }
                            _ => sample_has_missing = true,
                        }
                    }
                }
                _ => sample_has_missing = true,
            }
        } else {
            sample_has_missing = true;
        }

        if sample_has_missing {
            n_samples_with_missing += 1;
        }
    }

    AlleleCounts {
        allele_counts,
        n_called_alleles,
        n_samples,
        n_samples_with_missing,
    }
}

fn rebuild_samples(
    buf: &mut RecordBuf,
    remap: &[Option<usize>],
    inv_remap: &[usize],
    n_total_old: usize,
    n_total_new: usize,
    format_numbers: &[FormatNumber],
) {
    let gt_idx = buf.samples().keys().as_ref().get_index_of("GT");
    let mut new_values: Vec<Vec<Option<SampleValue>>> = Vec::with_capacity(buf.samples().len());

    for row in buf.samples().values() {
        let row_values = row.values();
        let (gt_lost, ploidy) = inspect_genotype(row_values, gt_idx, remap);

        if gt_lost {
            let mut new_row = Vec::with_capacity(row_values.len());
            for (i, _) in row_values.iter().enumerate() {
                if Some(i) == gt_idx {
                    let missing = (0..ploidy)
                        .map(|_| Allele::new(None, vcf::variant::record::samples::series::value::genotype::Phasing::Unphased))
                        .collect::<Vec<_>>();
                    let mut g = Genotype::default();
                    *g.as_mut() = missing;
                    new_row.push(Some(SampleValue::Genotype(g)));
                } else {
                    new_row.push(None);
                }
            }
            new_values.push(new_row);
            continue;
        }

        let mut new_row: Vec<Option<SampleValue>> = Vec::with_capacity(row_values.len());
        for (i, cell) in row_values.iter().enumerate() {
            if Some(i) == gt_idx {
                new_row.push(remap_genotype_cell(cell, remap));
                continue;
            }
            new_row.push(reslice_sample_cell(
                cell,
                format_numbers[i],
                inv_remap,
                n_total_old,
                n_total_new,
                ploidy,
            ));
        }
        new_values.push(new_row);
    }

    let keys = buf.samples().keys().clone();
    *buf.samples_mut() = Samples::new(keys, new_values);
}

fn inspect_genotype(
    row: &[Option<SampleValue>],
    gt_idx: Option<usize>,
    remap: &[Option<usize>],
) -> (bool, usize) {
    let Some(idx) = gt_idx else {
        return (false, 0);
    };
    match row.get(idx).and_then(|v| v.as_ref()) {
        Some(SampleValue::Genotype(g)) => {
            let alleles = g.as_ref();
            let ploidy = alleles.len();
            let lost = alleles.iter().any(|a| match a.position() {
                Some(p) => p < remap.len() && remap[p].is_none(),
                None => false, // already missing; not "lost"
            });
            (lost, ploidy)
        }
        _ => (false, 0),
    }
}

fn remap_genotype_cell(
    cell: &Option<SampleValue>,
    remap: &[Option<usize>],
) -> Option<SampleValue> {
    let SampleValue::Genotype(g) = cell.as_ref()? else {
        return cell.clone();
    };
    let new_alleles: Vec<Allele> = g
        .as_ref()
        .iter()
        .map(|a| {
            let new_pos = a.position().and_then(|p| remap.get(p).copied().flatten());
            Allele::new(new_pos, a.phasing())
        })
        .collect();
    let mut new_g = Genotype::default();
    *new_g.as_mut() = new_alleles;
    Some(SampleValue::Genotype(new_g))
}

fn reslice_sample_cell(
    cell: &Option<SampleValue>,
    number: FormatNumber,
    inv_remap: &[usize],
    n_total_old: usize,
    n_total_new: usize,
    ploidy: usize,
) -> Option<SampleValue> {
    let v = cell.as_ref()?;
    let SampleValue::Array(arr) = v else {
        return Some(v.clone());
    };
    match number {
        FormatNumber::AlternateBases => Some(SampleValue::Array(slice_sample_array(
            arr.clone(),
            &a_indices(inv_remap),
        ))),
        FormatNumber::ReferenceAlternateBases => Some(SampleValue::Array(slice_sample_array(
            arr.clone(),
            &r_indices(inv_remap),
        ))),
        FormatNumber::Samples => slice_sample_array_g(arr.clone(), n_total_old, n_total_new, ploidy, inv_remap)
            .map(SampleValue::Array),
        _ => Some(SampleValue::Array(arr.clone())),
    }
}

fn rebuild_info(
    buf: &mut RecordBuf,
    header: &Header,
    kept_alt_idx: &[usize],
    inv_remap: &[usize],
    n_total_old: usize,
    n_total_new: usize,
) {
    let info = buf.info_mut();
    let keys: Vec<String> = info.keys().cloned().collect();
    for key in keys {
        // AC/AF/AN/F_MISSING are recomputed wholesale below; skip them here.
        if matches!(key.as_str(), "AC" | "AF" | "AN" | "F_MISSING") {
            continue;
        }
        let Some(map) = header.infos().get(&key) else {
            continue;
        };
        let number = map.number();
        let cell = match info.get_mut(&key) {
            Some(c) => c,
            None => continue,
        };
        match number {
            InfoNumber::AlternateBases => {
                if let Some(InfoValue::Array(arr)) = cell.take() {
                    *cell = Some(InfoValue::Array(slice_info_array(arr, kept_alt_idx)));
                }
            }
            InfoNumber::ReferenceAlternateBases => {
                if let Some(InfoValue::Array(arr)) = cell.take() {
                    let r_idx: Vec<usize> = std::iter::once(0)
                        .chain(kept_alt_idx.iter().map(|&i| i + 1))
                        .collect();
                    *cell = Some(InfoValue::Array(slice_info_array(arr, &r_idx)));
                }
            }
            InfoNumber::Samples => {
                // Per-genotype INFO fields are rare; ploidy is undefined at the
                // INFO level. Slice using diploid as best effort, drop on length
                // mismatch.
                if let Some(InfoValue::Array(arr)) = cell.take() {
                    *cell = slice_info_array_g(arr, n_total_old, n_total_new, 2, inv_remap)
                        .map(InfoValue::Array);
                }
            }
            _ => {}
        }
    }
}

fn update_count_info_fields(
    buf: &mut RecordBuf,
    header: &Header,
    counts: &AlleleCounts,
    n_alt: usize,
) {
    let info = buf.info_mut();

    if header.infos().contains_key("AN") {
        info.insert(
            "AN".into(),
            Some(InfoValue::Integer(counts.n_called_alleles as i32)),
        );
    }

    if header.infos().contains_key("AC") && n_alt > 0 {
        let ac: Vec<Option<i32>> = (1..=n_alt)
            .map(|i| Some(counts.allele_counts[i] as i32))
            .collect();
        info.insert("AC".into(), Some(InfoValue::Array(InfoArray::Integer(ac))));
    }

    if header.infos().contains_key("AF") && n_alt > 0 {
        let total = counts.n_called_alleles as f64;
        let af: Vec<Option<f32>> = (1..=n_alt)
            .map(|i| {
                if total > 0.0 {
                    Some((counts.allele_counts[i] as f64 / total) as f32)
                } else {
                    Some(0.0)
                }
            })
            .collect();
        info.insert("AF".into(), Some(InfoValue::Array(InfoArray::Float(af))));
    }

    if header.infos().contains_key("F_MISSING") {
        let f_missing = if counts.n_samples > 0 {
            counts.n_samples_with_missing as f32 / counts.n_samples as f32
        } else {
            0.0
        };
        info.insert("F_MISSING".into(), Some(InfoValue::Float(f_missing)));
    }
}

// ---------- helpers: index sets and array slicing ----------

fn a_indices(inv_remap: &[usize]) -> Vec<usize> {
    // ALT array is indexed 0..n_alt (i.e. original allele index minus one).
    inv_remap
        .iter()
        .skip(1) // skip REF
        .map(|&orig| orig - 1)
        .collect()
}

fn r_indices(inv_remap: &[usize]) -> Vec<usize> {
    // R array is indexed 0..n_total (allele index).
    inv_remap.to_vec()
}

fn slice_info_array(arr: InfoArray, indices: &[usize]) -> InfoArray {
    match arr {
        InfoArray::Integer(v) => InfoArray::Integer(slice_indexed(&v, indices)),
        InfoArray::Float(v) => InfoArray::Float(slice_indexed(&v, indices)),
        InfoArray::Character(v) => InfoArray::Character(slice_indexed(&v, indices)),
        InfoArray::String(v) => InfoArray::String(slice_indexed(&v, indices)),
    }
}

fn slice_sample_array(arr: SampleArray, indices: &[usize]) -> SampleArray {
    match arr {
        SampleArray::Integer(v) => SampleArray::Integer(slice_indexed(&v, indices)),
        SampleArray::Float(v) => SampleArray::Float(slice_indexed(&v, indices)),
        SampleArray::Character(v) => SampleArray::Character(slice_indexed(&v, indices)),
        SampleArray::String(v) => SampleArray::String(slice_indexed(&v, indices)),
    }
}

fn slice_indexed<T: Clone>(src: &[Option<T>], indices: &[usize]) -> Vec<Option<T>> {
    indices
        .iter()
        .map(|&i| src.get(i).cloned().unwrap_or(None))
        .collect()
}

fn slice_info_array_g(
    arr: InfoArray,
    n_total_old: usize,
    n_total_new: usize,
    ploidy: usize,
    inv_remap: &[usize],
) -> Option<InfoArray> {
    let map = build_g_index_map(arr_len(&arr), n_total_old, n_total_new, ploidy, inv_remap)?;
    Some(match arr {
        InfoArray::Integer(v) => InfoArray::Integer(remap_g(&v, &map)),
        InfoArray::Float(v) => InfoArray::Float(remap_g(&v, &map)),
        InfoArray::Character(v) => InfoArray::Character(remap_g(&v, &map)),
        InfoArray::String(v) => InfoArray::String(remap_g(&v, &map)),
    })
}

fn slice_sample_array_g(
    arr: SampleArray,
    n_total_old: usize,
    n_total_new: usize,
    ploidy: usize,
    inv_remap: &[usize],
) -> Option<SampleArray> {
    let map = build_g_index_map(sample_arr_len(&arr), n_total_old, n_total_new, ploidy, inv_remap)?;
    Some(match arr {
        SampleArray::Integer(v) => SampleArray::Integer(remap_g(&v, &map)),
        SampleArray::Float(v) => SampleArray::Float(remap_g(&v, &map)),
        SampleArray::Character(v) => SampleArray::Character(remap_g(&v, &map)),
        SampleArray::String(v) => SampleArray::String(remap_g(&v, &map)),
    })
}

fn arr_len(arr: &InfoArray) -> usize {
    match arr {
        InfoArray::Integer(v) => v.len(),
        InfoArray::Float(v) => v.len(),
        InfoArray::Character(v) => v.len(),
        InfoArray::String(v) => v.len(),
    }
}

fn sample_arr_len(arr: &SampleArray) -> usize {
    match arr {
        SampleArray::Integer(v) => v.len(),
        SampleArray::Float(v) => v.len(),
        SampleArray::Character(v) => v.len(),
        SampleArray::String(v) => v.len(),
    }
}

/// Returns `Some(map)` where `map[new_g_idx] = old_g_idx`, or `None` when the
/// observed array length doesn't match the expected `C(n_total_old + ploidy - 1, ploidy)`
/// (in which case the field can't be safely re-sliced and the caller should drop it).
fn build_g_index_map(
    observed_len: usize,
    n_total_old: usize,
    n_total_new: usize,
    ploidy: usize,
    inv_remap: &[usize],
) -> Option<Vec<usize>> {
    if ploidy == 0 {
        return None;
    }
    let old_layout = enumerate_genotypes(n_total_old, ploidy);
    if old_layout.len() != observed_len {
        return None;
    }
    let new_layout = enumerate_genotypes(n_total_new, ploidy);
    let old_idx: HashMap<Vec<usize>, usize> = old_layout
        .into_iter()
        .enumerate()
        .map(|(i, t)| (t, i))
        .collect();
    let mut map = Vec::with_capacity(new_layout.len());
    for new_g in &new_layout {
        let mut old_g: Vec<usize> = new_g.iter().map(|&a| inv_remap[a]).collect();
        old_g.sort_unstable();
        let i = *old_idx.get(&old_g).expect("inv_remap is in-range");
        map.push(i);
    }
    Some(map)
}

fn remap_g<T: Clone>(src: &[Option<T>], map: &[usize]) -> Vec<Option<T>> {
    map.iter().map(|&i| src[i].clone()).collect()
}

/// Enumerate all unordered genotypes for the given allele count and ploidy in
/// VCF Number=G order: outermost loop on the largest allele, inner loops on
/// the next smaller, etc.
fn enumerate_genotypes(n_alleles: usize, ploidy: usize) -> Vec<Vec<usize>> {
    let mut out = Vec::new();
    if ploidy == 0 || n_alleles == 0 {
        return out;
    }
    let mut current = vec![0usize; ploidy];
    enumerate_inner(&mut out, &mut current, ploidy, n_alleles - 1);
    out
}

fn enumerate_inner(
    out: &mut Vec<Vec<usize>>,
    current: &mut Vec<usize>,
    pos: usize,
    max_val: usize,
) {
    if pos == 0 {
        out.push(current.clone());
        return;
    }
    for v in 0..=max_val {
        current[pos - 1] = v;
        enumerate_inner(out, current, pos - 1, v);
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn enumerate_genotypes_diploid_three_alleles() {
        let g = enumerate_genotypes(3, 2);
        assert_eq!(
            g,
            vec![
                vec![0, 0],
                vec![0, 1],
                vec![1, 1],
                vec![0, 2],
                vec![1, 2],
                vec![2, 2],
            ]
        );
    }

    #[test]
    fn enumerate_genotypes_triploid_two_alleles() {
        let g = enumerate_genotypes(2, 3);
        assert_eq!(
            g,
            vec![vec![0, 0, 0], vec![0, 0, 1], vec![0, 1, 1], vec![1, 1, 1]]
        );
    }

    #[test]
    fn build_g_map_diploid_drop_middle_alt() {
        // Original alleles: 0=REF, 1=A, 2=C, 3=G. We drop allele 2 (C).
        // Old layout (n=4, P=2): 10 genotypes
        //   indices: 00, 01, 11, 02, 12, 22, 03, 13, 23, 33
        // New layout (n=3, P=2): 0=REF, 1=A, 2=G
        //   indices: 00, 01, 11, 02, 12, 22
        // map[new] = old:
        //   (0,0) old (0,0) -> 0
        //   (0,1) old (0,1) -> 1
        //   (1,1) old (1,1) -> 2
        //   (0,2) old (0,3) -> 6
        //   (1,2) old (1,3) -> 7
        //   (2,2) old (3,3) -> 9
        let inv_remap = vec![0, 1, 3]; // new 0->old 0, new 1->old 1, new 2->old 3
        let map = build_g_index_map(10, 4, 3, 2, &inv_remap).unwrap();
        assert_eq!(map, vec![0, 1, 2, 6, 7, 9]);
    }
}
