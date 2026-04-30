//! `upper`: force REF/ALT alleles to uppercase. GLnexus workaround.

use std::path::PathBuf;

use anyhow::Result;
use clap::Args;

use crate::processor::ParallelVariantWindowProcessor;
use crate::util::uppercase_alleles;

#[derive(Args, Debug)]
pub struct UpperArgs {
    /// Input variant file (must be indexed: vcf.gz+tbi or bcf+csi).
    #[arg(long)]
    pub input: PathBuf,

    /// Output vcf.gz path.
    #[arg(long)]
    pub output: PathBuf,

    /// Worker thread count. Defaults to the system parallelism.
    #[arg(long)]
    pub threads: Option<usize>,

    /// Window size in base pairs. Defaults to 1 Mbp.
    #[arg(long)]
    pub window_size: Option<u64>,
}

pub fn run(args: UpperArgs) -> Result<()> {
    let threads = args.threads.unwrap_or_else(|| {
        std::thread::available_parallelism()
            .map(usize::from)
            .unwrap_or(1)
    });

    let mut builder = ParallelVariantWindowProcessor::builder()
        .input(args.input)
        .with_output_file(args.output)
        .worker_threads(threads)
        .record_callback(uppercase_alleles)
        .progress_callback(|bp| {
            eprint!("\r  {} Mbp processed", bp / 1_000_000);
        });

    if let Some(ws) = args.window_size {
        builder = builder.window_size(ws);
    }

    let result = builder.run();
    eprintln!(); // newline after the in-place progress line.
    result
}
