//! Per-window averages of numeric INFO fields.
//!
//! Each record is assigned to window `floor(POS0 / window_size)` of its
//! contig (by start position). For every requested tag, the mean over the
//! window's sites with a non-missing value is reported. Multi-valued tags
//! (Number=A/R/.) use the value for the first ALT allele (as `snpsketch
//! --info-stats`).

use std::collections::HashMap;
use std::io::Write;

use anyhow::{Result, anyhow, bail};
use rust_htslib::bcf::Record;

use crate::hts::{self, HeaderView};

/// A numeric INFO field to average.
#[derive(Clone, Debug)]
pub struct TagSpec {
    pub name: String,
    is_float: bool,
    /// Value index used for multi-valued tags (Number=R: first ALT).
    idx: usize,
}

impl TagSpec {
    /// Resolve comma-separated tag names against the header's INFO
    /// definitions. Each must be a declared Integer or Float field (not
    /// Number=G).
    pub fn parse(spec: &str, header: &HeaderView) -> Result<Vec<TagSpec>> {
        let infos = hts::header_infos(header);
        let mut out = Vec::new();
        for name in spec.split(',').map(str::trim).filter(|s| !s.is_empty()) {
            let name = name.strip_prefix("INFO/").unwrap_or(name);
            let (_, number, ty) = infos
                .iter()
                .find(|(id, _, _)| id == name)
                .ok_or_else(|| anyhow!("INFO/{name} is not defined in the header (use --fill-tags to compute it?)"))?;
            let is_float = match ty.as_str() {
                "Integer" => false,
                "Float" => true,
                _ => bail!("INFO/{name} has Type={ty}; only Integer and Float fields can be averaged"),
            };
            if number == "G" {
                bail!("INFO/{name} has Number=G, which cannot be averaged per site");
            }
            out.push(TagSpec { name: name.to_string(), is_float, idx: if number == "R" { 1 } else { 0 } });
        }
        if out.is_empty() {
            bail!("no tags given");
        }
        Ok(out)
    }

    fn value(&self, rec: &Record) -> Option<f64> {
        let v = if self.is_float {
            hts::info_f32(rec, &self.name).ok().and_then(|v| {
                let x = *v.get(self.idx)?;
                (!hts::is_float_missing(x) && !hts::is_float_vector_end(x)).then_some(x as f64)
            })
        } else {
            hts::info_i32(rec, &self.name).ok().and_then(|v| {
                let x = *v.get(self.idx)?;
                (x != hts::INT32_MISSING && x != hts::INT32_VECTOR_END).then_some(x as f64)
            })
        };
        v.filter(|x| x.is_finite())
    }
}

#[derive(Clone, Debug)]
struct Bin {
    n_sites: u64,
    sum: Vec<f64>,
    n: Vec<u64>,
}

impl Bin {
    fn new(n_tags: usize) -> Self {
        Self { n_sites: 0, sum: vec![0.0; n_tags], n: vec![0; n_tags] }
    }
}

/// Accumulates per-window sums; mergeable across workers.
#[derive(Clone, Debug)]
pub struct WindowStats {
    window_size: u64,
    tags: Vec<TagSpec>,
    bins: HashMap<(String, u64), Bin>,
}

impl WindowStats {
    pub fn new(window_size: u64, tags: Vec<TagSpec>) -> Self {
        Self { window_size: window_size.max(1), tags, bins: HashMap::new() }
    }

    pub fn process_record(&mut self, rec: &Record) {
        let key = (hts::rec_contig(rec), rec.pos().max(0) as u64 / self.window_size);
        let n_tags = self.tags.len();
        let bin = self.bins.entry(key).or_insert_with(|| Bin::new(n_tags));
        bin.n_sites += 1;
        for (i, t) in self.tags.iter().enumerate() {
            if let Some(v) = t.value(rec) {
                bin.sum[i] += v;
                bin.n[i] += 1;
            }
        }
    }

    pub fn merge(&mut self, other: WindowStats) {
        let n_tags = self.tags.len();
        for (key, b) in other.bins {
            let a = self.bins.entry(key).or_insert_with(|| Bin::new(n_tags));
            a.n_sites += b.n_sites;
            for i in 0..n_tags {
                a.sum[i] += b.sum[i];
                a.n[i] += b.n[i];
            }
        }
    }

    /// Write `contig,start,end,n_sites,<tag>...` rows (1-based, inclusive
    /// coordinates). Contigs are written in the order given, then any
    /// contigs with data that were not listed (sorted by name). Every window
    /// up to a contig's length is written (empty ones with n_sites=0 and
    /// empty means); without a length, up to its last window with data.
    pub fn write_csv<W: Write>(&self, writer: W, contigs: &[(String, Option<u64>)]) -> Result<()> {
        let mut w = csv::Writer::from_writer(writer);
        let mut header = vec!["contig".to_string(), "start".into(), "end".into(), "n_sites".into()];
        header.extend(self.tags.iter().map(|t| t.name.clone()));
        w.write_record(&header)?;

        // Last window index with data, per contig.
        let mut last: HashMap<&str, u64> = HashMap::new();
        for (c, k) in self.bins.keys() {
            let e = last.entry(c.as_str()).or_insert(*k);
            *e = (*e).max(*k);
        }
        let mut order: Vec<(&str, Option<u64>)> = contigs.iter().map(|(c, l)| (c.as_str(), *l)).collect();
        let mut extra: Vec<&str> = last.keys().copied().filter(|c| !contigs.iter().any(|(n, _)| n == c)).collect();
        extra.sort_unstable();
        order.extend(extra.into_iter().map(|c| (c, None)));

        let ws = self.window_size;
        let empty = Bin::new(self.tags.len());
        let mut row: Vec<String> = Vec::with_capacity(header.len());
        for (contig, len) in order {
            let by_len = len.filter(|&l| l > 0).map(|l| l.div_ceil(ws));
            let by_data = last.get(contig).map(|&k| k + 1);
            let n_windows = match (by_len, by_data) {
                (Some(a), Some(b)) => a.max(b),
                (a, b) => a.or(b).unwrap_or(0),
            };
            for k in 0..n_windows {
                let bin = self.bins.get(&(contig.to_string(), k)).unwrap_or(&empty);
                let start = k * ws + 1;
                let mut end = (k + 1) * ws;
                if let Some(l) = len.filter(|&l| l >= start) {
                    end = end.min(l);
                }
                row.clear();
                row.push(contig.to_string());
                row.push(start.to_string());
                row.push(end.to_string());
                row.push(bin.n_sites.to_string());
                for i in 0..self.tags.len() {
                    row.push(if bin.n[i] > 0 { format!("{:.6}", bin.sum[i] / bin.n[i] as f64) } else { String::new() });
                }
                w.write_record(&row)?;
            }
        }
        w.flush()?;
        Ok(())
    }
}
