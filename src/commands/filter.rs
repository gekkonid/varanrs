//! `filter`: per-site allele filtering by minimum AC and/or AF.

use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::Args;
use noodles_vcf as vcf;
use vcf::variant::io::Write as _;

use crate::filter::filter_alleles_at_site;

#[derive(Args, Debug)]
pub struct FilterArgs {
    /// Input VCF/BCF path.
    #[arg(long)]
    pub input: PathBuf,

    /// Output VCF path.
    #[arg(long)]
    pub output: PathBuf,

    /// Minimum allele count (inclusive).
    #[arg(long)]
    pub min_ac: Option<u32>,

    /// Minimum allele frequency (inclusive).
    #[arg(long)]
    pub min_af: Option<f64>,
}

pub fn run(args: FilterArgs) -> Result<()> {
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_path(&args.input)
        .with_context(|| format!("failed to open input: {}", args.input.display()))?;

    let header = reader
        .read_header()
        .with_context(|| "failed to read VCF header")?;

    let mut writer = vcf::io::Writer::new(
        std::fs::File::create(&args.output)
            .with_context(|| format!("failed to create output: {}", args.output.display()))?,
    );

    writer
        .write_header(&header)
        .with_context(|| "failed to write VCF header")?;

    for result in reader.record_bufs(&header) {
        let buf = result.with_context(|| "failed to read record")?;
        if let Some(out) = filter_alleles_at_site(buf, &header, args.min_af, args.min_ac) {
            writer
                .write_variant_record(&header, &out)
                .with_context(|| "failed to write record")?;
        }
    }

    Ok(())
}
