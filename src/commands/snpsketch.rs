use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use anyhow::{Context, Result, anyhow};
use clap::Args;
use indexmap::IndexMap;

use crate::args::{self, IndexedInput};
use crate::hts::{self, VcfReader};
use crate::processor::{ParallelVariantWindowProcessor, RecordFn};
use crate::snpsketch::SketchAccumulator;
use crate::util::human_bp;

#[derive(Args, Debug)]
pub struct SketchArgs {
    #[command(flatten)]
    pub indexed: IndexedInput,

    #[arg(long, default_value = "100")]
    pub stride: usize,

    #[arg(long, default_value = "16384")]
    pub chunk: u64,

    #[arg(long, default_value = "pairs.csv")]
    pub pairs: PathBuf,

    #[arg(long, default_value = "sample_stats.csv")]
    pub sample_stats: PathBuf,

    #[arg(long)]
    pub genotypes: Option<PathBuf>,

    #[arg(long)]
    pub info_stats: Option<PathBuf>,

    #[arg(long, default_value = "100")]
    pub info_hist_bins: usize,
}

pub fn run(args: SketchArgs) -> Result<()> {
    let input = args.indexed.input.as_ref()
        .ok_or_else(|| anyhow!("snpsketch requires an indexed input file"))?;

    let reader = VcfReader::open_indexed(input)?;
    let header = reader.header();

    let sample_ids: Vec<String> = header
        .samples()
        .iter()
        .map(|s| String::from_utf8_lossy(s).into_owned())
        .collect();

    let (contig_rank, contig_lengths) = if let Some(ref fai_path) = args.indexed.fai {
        let fai = args::parse_fai(fai_path)?;
        let rank: IndexMap<String, usize> = fai
            .ordered
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.clone(), i))
            .collect();
        (rank, Some(fai.lengths))
    } else {
        let rank: IndexMap<String, usize> = hts::header_contigs(header)
            .into_iter()
            .enumerate()
            .map(|(i, (name, _))| (name, i))
            .collect();
        (rank, None)
    };

    if sample_ids.is_empty() {
        return Err(anyhow!("no samples found in VCF header"));
    }

    let n_samples = sample_ids.len();

    if n_samples > 5000 {
        eprintln!(
            "varanrs snpsketch: {} samples, output pairs.csv will have {} rows (~{:.1} GB).",
            n_samples,
            n_samples * (n_samples - 1) / 2,
            (n_samples * (n_samples - 1) / 2) as f64 * 45.0 / 1_000_000_000.0
        );
    }

    eprintln!(
        "varanrs snpsketch: {} samples, {} contig(s), stride={}, chunk={}bp",
        n_samples,
        if args.indexed.contig.is_empty() { "all" } else { "filtered" },
        args.stride,
        args.chunk,
    );


    let threads = args.indexed.resolve_threads();

    // One accumulator per worker (no lock contention); merged at the end.
    let info_bins = args.info_stats.is_some().then_some(args.info_hist_bins);
    let new_accumulator = {
        let sample_ids = sample_ids.clone();
        move |h: &rust_htslib::bcf::header::HeaderView| {
            let mut a = SketchAccumulator::new(sample_ids.clone(), contig_rank.clone());
            if let Some(bins) = info_bins {
                a.init_info_histograms(h, bins);
            }
            a
        }
    };
    let mut merged = new_accumulator(header);
    let workers: Arc<Mutex<Vec<Arc<Mutex<SketchAccumulator>>>>> = Arc::new(Mutex::new(Vec::new()));
    let factory = {
        let workers = Arc::clone(&workers);
        move |h: &Arc<rust_htslib::bcf::header::HeaderView>| -> Result<RecordFn> {
            let acc = Arc::new(Mutex::new(new_accumulator(h)));
            workers.lock().unwrap().push(Arc::clone(&acc));
            Ok(Box::new(move |rec: &mut rust_htslib::bcf::Record| {
                acc.lock().unwrap().process_record(rec);
                Ok(false)
            }) as RecordFn)
        }
    };

    let mut builder = ParallelVariantWindowProcessor::builder()
        .input(input)
        .worker_threads(threads)
        .window_size(args.chunk)
        .stride(args.stride)
        .record_factory(factory)
        .progress_callback(|cum_bp, done, total| {
            let pct = if total > 0 {
                done as f64 / total as f64 * 100.0
            } else {
                0.0
            };
            let (val, unit) = human_bp(cum_bp);
            eprint!(
                "\r  {:.1}% ({}/{})  {:.1}{}  ",
                pct, done, total, val, unit
            );
            use std::io::Write as _;
            std::io::stderr().flush().ok();
        });

    if !args.indexed.contig.is_empty() {
        builder = builder.contigs(args.indexed.contig);
    }

    if let Some(lengths) = contig_lengths {
        builder = builder.contig_lengths(lengths);
    }

    builder.run()?;
    eprintln!();

    for w in workers.lock().unwrap().drain(..) {
        let a = std::mem::replace(&mut *w.lock().unwrap(), SketchAccumulator::new(Vec::new(), IndexMap::new()));
        merged.merge(a);
    }
    let acc = &merged;

    eprintln!(
        "  {} sites processed across {} samples",
        acc.n_sites, n_samples
    );

    let pairs_file = BufWriter::new(
        File::create(&args.pairs)
            .with_context(|| format!("failed to create {}", args.pairs.display()))?,
    );
    eprintln!("  writing {}", args.pairs.display());
    acc.write_pairs_csv(pairs_file)
        .with_context(|| format!("failed to write {}", args.pairs.display()))?;

    let sample_stats_file = BufWriter::new(
        File::create(&args.sample_stats)
            .with_context(|| format!("failed to create {}", args.sample_stats.display()))?,
    );
    eprintln!("  writing {}", args.sample_stats.display());
    acc.write_sample_stats_csv(sample_stats_file)
        .with_context(|| format!("failed to write {}", args.sample_stats.display()))?;

    if let Some(ref genotypes_path) = args.genotypes {
        let genotypes_file = BufWriter::new(
            File::create(genotypes_path)
                .with_context(|| format!("failed to create {}", genotypes_path.display()))?,
        );
        eprintln!("  writing {}", genotypes_path.display());
        acc.write_genotypes_csv(genotypes_file)
            .with_context(|| format!("failed to write {}", genotypes_path.display()))?;
    }

    if let Some(ref info_stats_path) = args.info_stats {
        let info_stats_file = BufWriter::new(
            File::create(info_stats_path)
                .with_context(|| format!("failed to create {}", info_stats_path.display()))?,
        );
        eprintln!("  writing {}", info_stats_path.display());
        acc.write_info_stats_csv(info_stats_file)
            .with_context(|| format!("failed to write {}", info_stats_path.display()))?;
    }

    eprintln!("  done.");
    Ok(())
}
