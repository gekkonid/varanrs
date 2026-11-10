//! `allelefilter`: per-site allele filtering by minimum AC and/or AF.

use anyhow::Result;
use clap::Args;
use rust_htslib::bcf::Record;

use crate::commands::io;
use crate::filter::filter_alleles_at_site;
use crate::processor::RecordFn;

#[derive(Args, Debug)]
pub struct AllelefilterArgs {
    /// Input VCF/BCF path (any compression). Reads stdin when omitted or "-".
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

pub fn run(args: AllelefilterArgs) -> Result<()> {
    let input = if io::is_std(&args.input) { "-".to_string() } else { args.input.clone().unwrap() };
    let (min_ac, min_af) = (args.min_ac, args.min_af);
    io::stream_records(
        &input,
        None,
        io::open_writer(&args.output)?,
        |_| Ok(Vec::new()),
        move |_| -> Result<RecordFn> {
            Ok(Box::new(move |rec: &mut Record| filter_alleles_at_site(rec, min_af, min_ac)))
        },
    )
}
