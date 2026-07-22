use std::fs::File;
use std::io::BufWriter;
use std::path::PathBuf;
use std::sync::{Arc, Mutex};

use anyhow::{Context, Result, anyhow};
use clap::Args;
use indexmap::IndexMap;
use noodles_util::variant;

use crate::args::{self, IndexedInput};
use crate::processor::ParallelVariantWindowProcessor;
use crate::snpsketch::SketchAccumulator;

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

    #[arg(long, default_value = "missing.csv")]
    pub missing: PathBuf,

    #[arg(long)]
    pub genotypes: Option<PathBuf>,
}

pub fn run(args: SketchArgs) -> Result<()> {
    let input = &args.indexed.input;

    let mut reader = variant::io::indexed_reader::Builder::default()
        .build_from_path(input)
        .with_context(|| format!("failed to open indexed reader: {}", input.display()))?;

    let header = reader
        .read_header()
        .context("failed to read header")?;

    let sample_ids: Vec<String> = header.sample_names().iter().cloned().collect();

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
        let rank: IndexMap<String, usize> = header
            .contigs()
            .iter()
            .enumerate()
            .map(|(i, (name, _))| (name.to_string(), i))
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

    drop(reader);

    let threads = args.indexed.resolve_threads();

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
        .input(input)
        .worker_threads(threads)
        .window_size(args.chunk)
        .stride(args.stride)
        .record_callback(record_callback)
        .progress_callback(|bp| {
            eprint!("\r  {} Mbp sampled", bp / 1_000_000);
        });

    if !args.indexed.contig.is_empty() {
        builder = builder.contigs(args.indexed.contig);
    }

    if let Some(lengths) = contig_lengths {
        builder = builder.contig_lengths(lengths);
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
