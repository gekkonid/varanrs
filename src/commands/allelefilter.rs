//! `allelefilter`: per-site allele filtering by minimum AC and/or AF.

use std::io::{BufWriter, Write};

use anyhow::{Context, Result};
use clap::Args;
use noodles_util::variant;
use noodles_vcf as vcf;
use vcf::variant::RecordBuf;
use vcf::variant::io::Write as _;

use crate::filter::filter_alleles_at_site;

#[derive(Args, Debug)]
pub struct AllelefilterArgs {
    /// Input VCF/BCF path. Reads stdin when omitted or set to "-".
    #[arg()]
    pub input: Option<String>,

    /// Output VCF path. Writes stdout when omitted or set to "-".
    #[arg(long = "output", short = 'o')]
    pub output: Option<String>,

    /// Minimum allele count (inclusive).
    #[arg(long)]
    pub min_ac: Option<u32>,

    /// Minimum allele frequency (inclusive).
    #[arg(long)]
    pub min_af: Option<f64>,
}

fn is_std(path: &Option<String>) -> bool {
    path.as_deref().is_none_or(|s| s == "-")
}

fn open_reader(
    input: &Option<String>,
) -> Result<variant::io::Reader<Box<dyn std::io::BufRead>>> {
    let boxed: Box<dyn std::io::BufRead> = if is_std(input) {
        Box::new(std::io::BufReader::new(std::io::stdin()))
    } else {
        let path = input.as_ref().unwrap();
        let file = std::fs::File::open(path)
            .with_context(|| format!("failed to open: {path}"))?;
        Box::new(std::io::BufReader::new(file))
    };
    variant::io::reader::Builder::default()
        .build_from_reader(boxed)
        .context("failed to build reader")
}

fn open_writer(output: &Option<String>) -> Result<Box<dyn Write>> {
    if is_std(output) {
        Ok(Box::new(BufWriter::new(std::io::stdout().lock())))
    } else {
        let path = output.as_ref().unwrap();
        Ok(Box::new(
            BufWriter::new(
                std::fs::File::create(path)
                    .with_context(|| format!("failed to create output: {path}"))?,
            ),
        ))
    }
}

pub fn run(args: AllelefilterArgs) -> Result<()> {
    let mut reader = open_reader(&args.input)?;

    let mut header = reader
        .read_header()
        .with_context(|| "failed to read VCF header")?;
    crate::util::normalize_header_for_noodles(&mut header);

    let mut writer = vcf::io::Writer::new(open_writer(&args.output)?);
    writer
        .write_header(&header)
        .with_context(|| "failed to write VCF header")?;

    for result in reader.records(&header) {
        let record = result.with_context(|| "failed to read record")?;
        let buf = RecordBuf::try_from_variant_record(&header, record.as_ref())
            .with_context(|| "failed to convert record")?;
        if let Some(out) = filter_alleles_at_site(buf, &header, args.min_af, args.min_ac) {
            writer
                .write_variant_record(&header, &out)
                .with_context(|| "failed to write record")?;
        }
    }

    Ok(())
}
