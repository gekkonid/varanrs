use std::io::Write;

use anyhow::Result;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use histr::StreamHist;
use indexmap::IndexMap;
use rust_htslib::bcf::Record;

use crate::hts::{self, HeaderView};

struct InfoHistogram {
    name: String,
    /// Value index used for multi-valued tags (Number=R: first ALT).
    target_idx: usize,
    is_float: bool,
    hist: StreamHist,
}

impl InfoHistogram {
    fn insert(&mut self, rec: &Record) {
        let v = if self.is_float {
            hts::info_f32(rec, &self.name).ok().and_then(|v| {
                let x = *v.get(self.target_idx)?;
                (!hts::is_float_missing(x) && !hts::is_float_vector_end(x)).then_some(x as f64)
            })
        } else {
            hts::info_i32(rec, &self.name).ok().and_then(|v| {
                let x = *v.get(self.target_idx)?;
                (x != hts::INT32_MISSING && x != hts::INT32_VECTOR_END).then_some(x as f64)
            })
        };
        if let Some(v) = v
            && v.is_finite()
        {
            self.hist.insert(v);
        }
    }
}

pub struct SketchAccumulator {
    pub sample_ids: Vec<String>,
    pub n_sites: usize,
    positions: Vec<(u32, u32)>,
    contig_rank: IndexMap<String, usize>,
    calls: Vec<BitVec<u64, Lsb0>>,
    any_alt: Vec<BitVec<u64, Lsb0>>,
    hom_alt: Vec<BitVec<u64, Lsb0>>,
    dp_sums: Vec<u64>,
    dp_counts: Vec<u64>,
    info_histograms: Vec<InfoHistogram>,
}

impl SketchAccumulator {
    pub fn new(sample_ids: Vec<String>, contig_rank: IndexMap<String, usize>) -> Self {
        let n = sample_ids.len();
        Self {
            sample_ids,
            n_sites: 0,
            positions: Vec::new(),
            contig_rank,
            calls: vec![BitVec::new(); n],
            any_alt: vec![BitVec::new(); n],
            hom_alt: vec![BitVec::new(); n],
            dp_sums: vec![0u64; n],
            dp_counts: vec![0u64; n],
            info_histograms: Vec::new(),
        }
    }

    /// Append another accumulator's sites (from the same samples, contig
    /// ranks and histogram setup). Pairwise and per-sample statistics don't
    /// depend on site order; genotype output is sorted by position.
    pub fn merge(&mut self, other: SketchAccumulator) {
        assert_eq!(self.calls.len(), other.calls.len(), "sample count mismatch");
        self.n_sites += other.n_sites;
        self.positions.extend(other.positions);
        for (a, b) in self.calls.iter_mut().zip(other.calls) {
            a.extend_from_bitslice(&b);
        }
        for (a, b) in self.any_alt.iter_mut().zip(other.any_alt) {
            a.extend_from_bitslice(&b);
        }
        for (a, b) in self.hom_alt.iter_mut().zip(other.hom_alt) {
            a.extend_from_bitslice(&b);
        }
        for (a, b) in self.dp_sums.iter_mut().zip(other.dp_sums) {
            *a += b;
        }
        for (a, b) in self.dp_counts.iter_mut().zip(other.dp_counts) {
            *a += b;
        }
        for (a, b) in self.info_histograms.iter_mut().zip(other.info_histograms) {
            debug_assert_eq!(a.name, b.name);
            a.hist.merge(b.hist);
        }
    }

    pub fn n_samples(&self) -> usize {
        self.calls.len()
    }

    pub fn sample_called(&self, s_idx: usize, site: usize) -> bool {
        self.calls[s_idx][site]
    }

    pub fn sample_any_alt(&self, s_idx: usize, site: usize) -> bool {
        self.any_alt[s_idx][site]
    }

    pub fn sample_hom_alt(&self, s_idx: usize, site: usize) -> bool {
        self.hom_alt[s_idx][site]
    }

    pub fn sample_call_count(&self, s_idx: usize) -> u64 {
        self.calls[s_idx].count_ones() as u64
    }

    pub fn process_record(&mut self, rec: &Record) {
        let contig_rank = self
            .contig_rank
            .get(&hts::rec_contig(rec))
            .copied()
            .unwrap_or(u32::MAX as usize);
        let pos = (rec.pos() + 1) as u32;
        let n_samples = self.calls.len();

        let gt = hts::format_i32(rec, "GT").ok().filter(|v| !v.is_empty());
        let dp = hts::format_i32(rec, "DP").ok().filter(|v| !v.is_empty());
        let gt_width = gt.as_ref().map_or(0, |v| v.len() / n_samples.max(1));
        let dp_width = dp.as_ref().map_or(0, |v| v.len() / n_samples.max(1));

        for s_idx in 0..n_samples {
            let mut has_missing = true;
            let mut n_alt: u8 = 0;
            let mut ploidy: u8 = 0;
            if let Some(gt) = gt.as_ref() {
                let alleles = gt[s_idx * gt_width..(s_idx + 1) * gt_width]
                    .iter()
                    .take_while(|&&v| v != hts::INT32_VECTOR_END);
                has_missing = false;
                for &v in alleles {
                    ploidy += 1;
                    if hts::gt_is_missing(v) {
                        has_missing = true;
                    } else if hts::gt_allele(v) != 0 {
                        n_alt += 1;
                    }
                }
            }

            if has_missing || ploidy == 0 {
                self.calls[s_idx].push(false);
                self.any_alt[s_idx].push(false);
                self.hom_alt[s_idx].push(false);
            } else {
                self.calls[s_idx].push(true);
                self.any_alt[s_idx].push(n_alt > 0);
                self.hom_alt[s_idx].push(n_alt == ploidy);
            }

            // Number=1 DP only (one value per sample)
            if let Some(dp) = dp.as_ref()
                && dp_width == 1
            {
                let v = dp[s_idx];
                if v != hts::INT32_MISSING && v != hts::INT32_VECTOR_END {
                    self.dp_sums[s_idx] += v as u64;
                    self.dp_counts[s_idx] += 1;
                }
            }
        }

        self.positions.push((contig_rank as u32, pos));
        self.n_sites += 1;

        for entry in &mut self.info_histograms {
            entry.insert(rec);
        }
    }

    pub fn init_info_histograms(&mut self, header: &HeaderView, bin_count: usize) {
        for (name, number, ty) in hts::header_infos(header) {
            if number == "G" {
                continue;
            }
            let is_float = match ty.as_str() {
                "Integer" => false,
                "Float" => true,
                _ => continue,
            };
            self.info_histograms.push(InfoHistogram {
                name,
                target_idx: if number == "R" { 1 } else { 0 },
                is_float,
                hist: StreamHist::with_capacity(bin_count),
            });
        }
    }

    pub fn compute_pair_stats(&self, i: usize, j: usize) -> (u64, u64, f64) {
        let ci = self.calls[i].as_raw_slice();
        let cj = self.calls[j].as_raw_slice();
        let ai = self.any_alt[i].as_raw_slice();
        let aj = self.any_alt[j].as_raw_slice();
        let hi = self.hom_alt[i].as_raw_slice();
        let hj = self.hom_alt[j].as_raw_slice();

        let n_words = ci.len();
        let mut n_common: u64 = 0;
        let mut diff1: u64 = 0;
        let mut diff2: u64 = 0;

        for w in 0..n_words {
            let common = ci[w] & cj[w];
            let a_xor = ai[w] ^ aj[w];
            let h_xor = hi[w] ^ hj[w];
            n_common += common.count_ones() as u64;
            diff1 += (common & (a_xor ^ h_xor)).count_ones() as u64;
            diff2 += (common & a_xor & h_xor).count_ones() as u64;
        }

        let n_diff = diff1 + 2 * diff2;
        let distance = if n_common > 0 {
            n_diff as f64 / (2.0 * n_common as f64)
        } else {
            f64::NAN
        };

        (n_diff, n_common, distance)
    }

    pub fn write_pairs_csv<W: Write + Send>(&self, writer: W) -> Result<()> {
        use rayon::prelude::*;
        use std::sync::Mutex;

        let n = self.n_samples();
        let mut w = csv::Writer::from_writer(writer);

        w.write_record([
            "sample_i",
            "sample_j",
            "sample_i_id",
            "sample_j_id",
            "n_diff",
            "n_common",
            "distance",
        ])?;

        if n < 2 {
            return Ok(());
        }

        let w = Mutex::new(w);

        (0..n).into_par_iter().try_for_each(|i| -> Result<()> {
            let id_i = self.sample_ids[i].clone();
            let i_s = i.to_string();
            let ci = self.calls[i].as_raw_slice();
            let ai = self.any_alt[i].as_raw_slice();
            let hi = self.hom_alt[i].as_raw_slice();
            let n_words = ci.len();

            struct Row {
                j_s: String,
                id_j: String,
                n_diff_s: String,
                n_common_s: String,
                dist_s: String,
            }

            let batch_size = n - i - 1;
            let mut batch: Vec<Row> = Vec::with_capacity(batch_size);

            for j in (i + 1)..n {
                let cj = self.calls[j].as_raw_slice();
                let aj = self.any_alt[j].as_raw_slice();
                let hj = self.hom_alt[j].as_raw_slice();

                let mut n_common: u64 = 0;
                let mut diff1: u64 = 0;
                let mut diff2: u64 = 0;

                for w_idx in 0..n_words {
                    let common = ci[w_idx] & cj[w_idx];
                    let a_xor = ai[w_idx] ^ aj[w_idx];
                    let h_xor = hi[w_idx] ^ hj[w_idx];
                    n_common += common.count_ones() as u64;
                    diff1 += (common & (a_xor ^ h_xor)).count_ones() as u64;
                    diff2 += (common & a_xor & h_xor).count_ones() as u64;
                }

                let n_diff = diff1 + 2 * diff2;
                let dist = if n_common > 0 {
                    n_diff as f64 / (2.0 * n_common as f64)
                } else {
                    f64::NAN
                };

                batch.push(Row {
                    j_s: j.to_string(),
                    id_j: self.sample_ids[j].clone(),
                    n_diff_s: n_diff.to_string(),
                    n_common_s: n_common.to_string(),
                    dist_s: format!("{:.6}", dist),
                });
            }

            let mut gw = w.lock().unwrap();
            for row in &batch {
                gw.write_record([
                    i_s.as_str(),
                    row.j_s.as_str(),
                    id_i.as_str(),
                    row.id_j.as_str(),
                    row.n_diff_s.as_str(),
                    row.n_common_s.as_str(),
                    row.dist_s.as_str(),
                ])?;
            }

            Ok(())
        })?;

        w.into_inner().unwrap().flush()?;
        Ok(())
    }

    pub fn write_sample_stats_csv<W: Write>(&self, writer: W) -> Result<()> {
        let n = self.n_samples();
        let mut w = csv::Writer::from_writer(writer);

        w.write_record([
            "sample_id", "n_missing", "n_total", "miss_rate",
            "n_het", "het_rate", "avg_dp",
        ])?;

        let n_total = self.n_sites as u64;
        for s_idx in 0..n {
            let n_called = self.calls[s_idx].count_ones() as u64;
            let n_missing = n_total.saturating_sub(n_called);
            let miss_rate = if n_total > 0 {
                n_missing as f64 / n_total as f64
            } else {
                f64::NAN
            };

            let n_alt = self.any_alt[s_idx].count_ones() as u64;
            let n_hom_alt = self.hom_alt[s_idx].count_ones() as u64;
            let n_het = n_alt.saturating_sub(n_hom_alt);
            let het_rate = if n_called > 0 {
                n_het as f64 / n_called as f64
            } else {
                f64::NAN
            };

            let avg_dp = if self.dp_counts[s_idx] > 0 {
                self.dp_sums[s_idx] as f64 / self.dp_counts[s_idx] as f64
            } else {
                f64::NAN
            };

            let avg_dp_str;
            let avg_dp_col = if avg_dp.is_nan() {
                avg_dp_str = String::new();
                avg_dp_str.as_str()
            } else {
                avg_dp_str = format!("{:.1}", avg_dp);
                avg_dp_str.as_str()
            };

            w.write_record([
                self.sample_ids[s_idx].as_str(),
                &n_missing.to_string(),
                &n_total.to_string(),
                &format!("{:.6}", miss_rate),
                &n_het.to_string(),
                &format!("{:.6}", het_rate),
                avg_dp_col,
            ])?;
        }

        w.flush()?;
        Ok(())
    }

    pub fn write_genotypes_csv<W: Write>(&self, writer: W) -> Result<()> {
        let n = self.n_samples();
        let n_sites = self.n_sites;
        let mut w = csv::Writer::from_writer(writer);

        let mut header: Vec<String> = Vec::with_capacity(1 + n);
        header.push("site".to_string());
        for id in &self.sample_ids {
            header.push(id.clone());
        }
        w.write_record(header)?;

        let mut order: Vec<usize> = (0..n_sites).collect();
        order.sort_by_key(|&site| self.positions[site]);

        for (out_idx, &site) in order.iter().enumerate() {
            let mut row: Vec<String> = Vec::with_capacity(1 + n);
            row.push(out_idx.to_string());
            for s_idx in 0..n {
                let gt = if !self.calls[s_idx][site] {
                    ".".to_string()
                } else if self.hom_alt[s_idx][site] {
                    "2".to_string()
                } else if self.any_alt[s_idx][site] {
                    "1".to_string()
                } else {
                    "0".to_string()
                };
                row.push(gt);
            }
            w.write_record(row)?;
        }

        w.flush()?;
        Ok(())
    }

    pub fn write_info_stats_csv<W: Write>(&self, writer: W) -> Result<()> {
        let mut w = csv::Writer::from_writer(writer);
        w.write_record(["field", "bin_mean", "count"])?;
        for entry in &self.info_histograms {
            for bin in &entry.hist.bins {
                let (mean, count): (f64, u64) = bin.into();
                w.write_record([
                    entry.name.as_str(),
                    &format!("{:.6}", mean),
                    &count.to_string(),
                ])?;
            }
        }
        w.flush()?;
        Ok(())
    }
}
