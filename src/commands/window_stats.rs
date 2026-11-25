//! `window-stats`: average numeric INFO fields over fixed-size genomic
//! windows, optionally recomputing tags first (as `filter --fill-tags`).

use std::io::Write as _;
use std::sync::{Arc, Mutex};

use anyhow::{Context, Result};
use clap::Args;
use rust_htslib::bcf::Record;

use crate::args::{self, IndexedInput};
use crate::commands::io;
use crate::filltags::{self, FillTags, TagSet};
use crate::hts::{self, HeaderView, VcfReader};
use crate::processor::{ParallelVariantWindowProcessor, RecordFn};
use crate::windowstats::{TagSpec, WindowStats};

#[derive(Args, Debug)]
pub struct WindowStatsArgs {
    #[command(flatten)]
    pub indexed: IndexedInput,

    /// Comma-separated numeric INFO fields to average per window, e.g.
    /// "AC_Het,AF". Multi-valued fields use the first ALT allele's value.
    #[arg(long = "tags", short = 't', required = true)]
    pub tags: String,

    /// Window size in base pairs.
    #[arg(long = "window-size", short = 'w', default_value = "100000")]
    pub window_size: u64,

    /// Recompute tags before averaging: comma-separated list
    /// (AN,AC,AF,MAF,NS,AC_Hom,AC_Het,AC_Hemi,HWE,ExcHet,F_MISSING,VAF,VAF1)
    /// or "all".
    #[arg(long = "fill-tags")]
    pub fill_tags: Option<String>,

    /// Output CSV path (stdout when omitted or "-").
    #[arg(long = "output", short = 'o')]
    pub output: Option<String>,
}

pub fn run(args: WindowStatsArgs) -> Result<()> {
    anyhow::ensure!(args.window_size > 0, "--window-size must be positive");
    let fill_set = match &args.fill_tags {
        Some(spec) => filltags::parse_tags(spec)?,
        None => TagSet::default(),
    };
    let lines = fill_set.header_lines();
    let input = args.indexed.input.as_ref().filter(|p| p.as_os_str() != "-");

    let mut reader = match input {
        Some(p) => VcfReader::open_indexed(p)?,
        None => VcfReader::open("-")?,
    };
    reader.append_header_lines(&lines)?;
    let tags = TagSpec::parse(&args.tags, reader.header())?;

    // Output contig order/lengths: header, then .fai.
    let mut contigs = hts::header_contigs(reader.header());
    if let Some(fai_path) = args.indexed.fai.as_ref() {
        let fai = args::parse_fai(fai_path)?;
        for (name, len) in &fai.ordered {
            match contigs.iter_mut().find(|(n, _)| n == name) {
                Some((_, l)) => *l = l.or(Some(*len)),
                None => contigs.push((name.clone(), Some(*len))),
            }
        }
    }
    if !args.indexed.contig.is_empty() {
        contigs.retain(|(n, _)| args.indexed.contig.contains(n));
    }

    let new_fill = {
        let fill_set = fill_set.clone();
        move |h: &Arc<HeaderView>| -> Result<Option<FillTags>> {
            if fill_set.is_empty() { Ok(None) } else { Ok(Some(FillTags::new(&fill_set, h)?)) }
        }
    };

    let stats = match input {
        None => {
            let mut fill = new_fill(reader.header())?;
            let mut acc = WindowStats::new(args.window_size, tags);
            let filter = &args.indexed.contig;
            // Alternate two record buffers: bcftools' filter.c caches
            // per-record data keyed by the record's address.
            let mut recs = [reader.empty_record(), reader.empty_record()];
            let mut cur = 0usize;
            loop {
                cur ^= 1;
                let rec = &mut recs[cur];
                if !reader.read(rec)? {
                    break;
                }
                if !filter.is_empty() && !filter.contains(&hts::rec_contig(rec)) {
                    continue;
                }
                if let Some(f) = fill.as_mut() {
                    f.fill(rec).with_context(|| format!("at {}", hts::rec_locus(rec)))?;
                }
                acc.process_record(rec);
            }
            acc
        }
        Some(path) => {
            // Fail on bad --fill-tags before spawning workers.
            drop(new_fill(reader.header())?);
            drop(reader);
            let workers: Arc<Mutex<Vec<Arc<Mutex<WindowStats>>>>> = Arc::new(Mutex::new(Vec::new()));
            let factory = {
                let workers = Arc::clone(&workers);
                let tags = tags.clone();
                let ws = args.window_size;
                move |h: &Arc<HeaderView>| -> Result<RecordFn> {
                    let mut fill = new_fill(h)?;
                    let acc = Arc::new(Mutex::new(WindowStats::new(ws, tags.clone())));
                    workers.lock().unwrap().push(Arc::clone(&acc));
                    Ok(Box::new(move |rec: &mut Record| -> Result<bool> {
                        if let Some(f) = fill.as_mut() {
                            f.fill(rec)?;
                        }
                        acc.lock().unwrap().process_record(rec);
                        Ok(false)
                    }))
                }
            };
            let mut builder = ParallelVariantWindowProcessor::builder()
                .input(path)
                .worker_threads(args.indexed.resolve_threads())
                .window_size(args.window_size)
                .header_lines(lines)
                .record_factory(factory)
                .progress_callback(|bp, _, _| {
                    eprint!("\r  {} Mbp processed", bp / 1_000_000);
                    std::io::stderr().flush().ok();
                });
            if !args.indexed.contig.is_empty() {
                builder = builder.contigs(args.indexed.contig.clone());
            }
            if let Some(fai_path) = args.indexed.fai.as_ref() {
                builder = builder.contig_lengths(args::parse_fai(fai_path)?.lengths);
            }
            let result = builder.run();
            eprintln!();
            result?;
            let mut merged = WindowStats::new(args.window_size, tags);
            for w in workers.lock().unwrap().drain(..) {
                let a = std::mem::replace(&mut *w.lock().unwrap(), WindowStats::new(1, Vec::new()));
                merged.merge(a);
            }
            merged
        }
    };

    stats.write_csv(io::open_writer(&args.output)?, &contigs)
}
