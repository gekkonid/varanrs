//! `upper`: force REF/ALT alleles to uppercase. GLnexus workaround.

use std::path::PathBuf;

use anyhow::Result;
use clap::Args;

use crate::args::{self, IndexedInput};
use crate::processor::ParallelVariantWindowProcessor;
use crate::util::uppercase_alleles;

#[derive(Args, Debug)]
pub struct UpperArgs {
    #[command(flatten)]
    pub indexed: IndexedInput,

    /// Output vcf.gz path.
    #[arg(long)]
    pub output: PathBuf,

    /// Window size in base pairs. Defaults to 1 Mbp.
    #[arg(long)]
    pub window_size: Option<u64>,
}

pub fn run(args: UpperArgs) -> Result<()> {
    let threads = args.indexed.resolve_threads();

    let mut builder = ParallelVariantWindowProcessor::builder()
        .input(args.indexed.input)
        .with_output_file(args.output)
        .worker_threads(threads)
        .record_callback(uppercase_alleles)
        .progress_callback(|bp| {
            eprint!("\r  {} Mbp processed", bp / 1_000_000);
        });

    if let Some(ws) = args.window_size {
        builder = builder.window_size(ws);
    }

    if !args.indexed.contig.is_empty() {
        builder = builder.contigs(args.indexed.contig);
    }

    if let Some(ref fai_path) = args.indexed.fai {
        let fai = args::parse_fai(fai_path)?;
        builder = builder.contig_lengths(fai.lengths);
    }

    let result = builder.run();
    eprintln!();
    result
}
