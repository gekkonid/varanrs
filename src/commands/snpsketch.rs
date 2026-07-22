use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use anyhow::{Context, Result, anyhow};
use clap::Args;
use noodles_util::variant;

use crate::processor::ParallelVariantWindowProcessor;
use crate::snpsketch::SketchAccumulator;

#[derive(Args, Debug)]
pub struct SketchArgs {
    #[arg(long)]
    pub input: PathBuf,

    #[arg(long, default_value = "100")]
    pub stride: usize,

    #[arg(long, default_value = "16384")]
    pub chunk: u64,

    #[arg(long, default_value = "pairs.csv")]
    pub pairs: PathBuf,

    #[arg(long, default_value = "missing.csv")]
    pub missing: PathBuf,

    #[arg(long)]
    pub genotypes: Option<PathBuf>,

    #[arg(long)]
    pub threads: Option<usize>,

    #[arg(long)]
    pub contig: Vec<String>,
}

pub fn run(args: SketchArgs) -> Result<()> {
    let mut reader = variant::io::indexed_reader::Builder::default()
        .build_from_path(&args.input)
        .with_context(|| format!("failed to open indexed reader: {}", args.input.display()))?;

    let header = reader
        .read_header()
        .context("failed to read header")?;

    let sample_ids: Vec<String> = header.sample_names().iter().cloned().collect();

    let contig_rank: HashMap<String, usize> = header
        .contigs()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect();

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
        if args.contig.is_empty() { "all" } else { "filtered" },
        args.stride,
        args.chunk,
    );

    drop(reader);

    let threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
    });

    let accumulator = Arc::new(Mutex::new(SketchAccumulator::new(
        sample_ids.clone(),
        contig_rank,
    )));

    let acc = Arc::clone(&accumulator);
    let record_callback = move |buf| {
        acc.lock().unwrap().process_record(&buf);
        None
    };

    let mut builder = ParallelVariantWindowProcessor::builder()
        .input(args.input)
        .worker_threads(threads)
        .window_size(args.chunk)
        .stride(args.stride)
        .record_callback(record_callback)
        .progress_callback(|bp| {
            eprint!("\r  {} Mbp sampled", bp / 1_000_000);
        });

    if !args.contig.is_empty() {
        builder = builder.contigs(args.contig);
    }

    builder.run()?;
    eprintln!();

    let acc = accumulator.lock().unwrap();

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

    let missing_file = BufWriter::new(
        File::create(&args.missing)
            .with_context(|| format!("failed to create {}", args.missing.display()))?,
    );
    eprintln!("  writing {}", args.missing.display());
    acc.write_missingness_csv(missing_file)
        .with_context(|| format!("failed to write {}", args.missing.display()))?;

    if let Some(ref genotypes_path) = args.genotypes {
        let genotypes_file = BufWriter::new(
            File::create(genotypes_path)
                .with_context(|| format!("failed to create {}", genotypes_path.display()))?,
        );
        eprintln!("  writing {}", genotypes_path.display());
        acc.write_genotypes_csv(genotypes_file)
            .with_context(|| format!("failed to write {}", genotypes_path.display()))?;
    }

    eprintln!("  done.");
    Ok(())
}
