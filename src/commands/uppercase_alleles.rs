//! `uppercase-alleles`: force REF/ALT alleles to uppercase. GLnexus workaround.

use anyhow::Result;
use clap::Args;
use rust_htslib::bcf::Record;

use crate::args::IndexedInput;
use crate::commands::io;
use crate::processor::RecordFn;
use crate::util::uppercase_alleles;

#[derive(Args, Debug)]
pub struct UppercaseAllelesArgs {
    #[command(flatten)]
    pub indexed: IndexedInput,

    /// Output path. Indexed input writes bgzf VCF and requires this; stdin
    /// input writes uncompressed VCF to stdout when omitted or "-".
    #[arg(long = "output", short = 'o')]
    pub output: Option<String>,

    /// Window size in base pairs. Defaults to 16384 bp (one
    /// tabix/CSI index bin).
    #[arg(long)]
    pub window_size: Option<u64>,
}

pub fn run(args: UppercaseAllelesArgs) -> Result<()> {
    io::run_records(
        "uppercase-alleles",
        &args.indexed,
        &args.output,
        args.window_size,
        None,
        |_| Ok(Vec::new()),
        |_| -> Result<RecordFn> { Ok(Box::new(|rec: &mut Record| Ok(uppercase_alleles(rec)))) },
    )
}
