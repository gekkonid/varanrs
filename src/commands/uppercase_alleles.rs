//! `uppercase-alleles`: force REF/ALT alleles to uppercase. GLnexus workaround.

use std::io::Write;

use anyhow::{Context, Result};
use clap::Args;
use noodles_util::variant;
use noodles_vcf as vcf;
use vcf::variant::RecordBuf;
use vcf::variant::io::Write as _;

use crate::args::{self, IndexedInput};
use crate::processor::ParallelVariantWindowProcessor;
use crate::util::uppercase_alleles;

#[derive(Args, Debug)]
pub struct UppercaseAllelesArgs {
    #[command(flatten)]
    pub indexed: IndexedInput,

    /// Output vcf.gz path. Writes to stdout when omitted or set to "-".
    #[arg(long = "output", short = 'o')]
    pub output: Option<String>,

    /// Window size in base pairs. Defaults to 1 Mbp.
    #[arg(long)]
    pub window_size: Option<u64>,
}

fn is_stdout(output: &Option<String>) -> bool {
    output.as_deref().is_none_or(|s| s == "-")
}

fn open_writer(output: &Option<String>) -> Result<Box<dyn Write>> {
    if is_stdout(output) {
        Ok(Box::new(std::io::stdout().lock()))
    } else {
        let path = output.as_ref().unwrap();
        Ok(Box::new(
            std::fs::File::create(path)
                .with_context(|| format!("failed to create output: {path}"))?,
        ))
    }
}

pub fn run(args: UppercaseAllelesArgs) -> Result<()> {
    let Some(input) = args.indexed.input.as_ref() else {
        // Stdin path: auto-detect VCF or uBCF, apply uppercase, write stdout
        let boxed: Box<dyn std::io::BufRead> =
            Box::new(std::io::BufReader::new(std::io::stdin()));
        let mut reader = variant::io::reader::Builder::default()
            .build_from_reader(boxed)
            .context("failed to build reader from stdin")?;
        let mut header = reader.read_header().context("failed to read header")?;
        crate::util::normalize_header_for_noodles(&mut header);

        let mut writer = vcf::io::Writer::new(open_writer(&args.output)?);
        writer.write_header(&header).context("failed to write header")?;

        for result in reader.records(&header) {
            let record = result.context("failed to read record")?;
            let buf = RecordBuf::try_from_variant_record(&header, record.as_ref())
                .context("failed to convert record")?;
            if let Some(out) = uppercase_alleles(buf) {
                writer.write_variant_record(&header, &out)?;
            }
        }
        return Ok(());
    };

    // Indexed path: parallel processor
    let threads = args.indexed.resolve_threads();

    let mut builder = ParallelVariantWindowProcessor::builder()
        .input(input)
        .worker_threads(threads)
        .record_callback(uppercase_alleles)
        .progress_callback(|bp, _, _| {
            eprint!("\r  {} Mbp processed", bp / 1_000_000);
            use std::io::Write as _;
            std::io::stderr().flush().ok();
        });

    if is_stdout(&args.output) {
        return Err(anyhow::anyhow!(
            "uppercase-alleles with indexed input requires --output (cannot write bgzf to stdout)"
        ));
    }
    builder = builder.with_output_file(args.output.as_ref().unwrap());

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
