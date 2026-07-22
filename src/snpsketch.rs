use std::io::Write;

use anyhow::Result;
use bitvec::order::Lsb0;
use bitvec::vec::BitVec;
use indexmap::IndexMap;
use noodles_vcf as vcf;
use vcf::variant::record_buf::samples::sample::Value as SampleValue;
use vcf::variant::RecordBuf;

pub struct SketchAccumulator {
    pub sample_ids: Vec<String>,
    pub n_sites: usize,
    positions: Vec<(u32, u32)>,
    contig_rank: IndexMap<String, usize>,
    calls: Vec<BitVec<u64, Lsb0>>,
    any_alt: Vec<BitVec<u64, Lsb0>>,
    hom_alt: Vec<BitVec<u64, Lsb0>>,
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

    pub fn process_record(&mut self, buf: &RecordBuf) {
        let contig_rank = self
            .contig_rank
            .get(buf.reference_sequence_name())
            .copied()
            .unwrap_or(u32::MAX as usize);
        let pos = buf
            .variant_start()
            .map(|p| p.get() as u32)
            .unwrap_or(0);

        let samples = buf.samples();
        let keys = samples.keys().as_ref();
        let gt_idx = keys.get_index_of("GT");
        let n_samples = self.calls.len();

        let mut processed = 0usize;

        for (s_idx, sample) in samples.values().enumerate() {
            if s_idx >= n_samples {
                break;
            }
            processed = s_idx + 1;

            let mut has_missing = true;
            let mut n_alt: u8 = 0;
            let mut ploidy: u8 = 0;

            if let Some(idx) = gt_idx
                && let Some(SampleValue::Genotype(g)) =
                    sample.values().get(idx).and_then(|v| v.as_ref())
            {
                let alleles = g.as_ref();
                ploidy = alleles.len() as u8;
                if ploidy > 0 {
                    has_missing = false;
                    for allele in alleles {
                        match allele.position() {
                            Some(0) => {}
                            Some(_) => n_alt += 1,
                            None => has_missing = true,
                        }
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
        }

        for s_idx in processed..n_samples {
            self.calls[s_idx].push(false);
            self.any_alt[s_idx].push(false);
            self.hom_alt[s_idx].push(false);
        }

        self.positions.push((contig_rank as u32, pos));
        self.n_sites += 1;
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

    pub fn write_missingness_csv<W: Write>(&self, writer: W) -> Result<()> {
        let n = self.n_samples();
        let mut w = csv::Writer::from_writer(writer);

        w.write_record(["sample_id", "n_missing", "n_total", "miss_rate"])?;

        let n_total = self.n_sites as u64;
        for s_idx in 0..n {
            let n_called = self.calls[s_idx].count_ones() as u64;
            let n_missing = n_total.saturating_sub(n_called);
            let rate = if n_total > 0 {
                n_missing as f64 / n_total as f64
            } else {
                f64::NAN
            };
            w.write_record([
                self.sample_ids[s_idx].as_str(),
                &n_missing.to_string(),
                &n_total.to_string(),
                &format!("{:.6}", rate),
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
}
