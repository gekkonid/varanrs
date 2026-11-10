//! `fill-tags`: recompute INFO/FORMAT tags from genotypes, compatible with
//! `bcftools +fill-tags`.

use std::sync::Arc;

use anyhow::Result;
use clap::Args;
use rust_htslib::bcf::Record;

use crate::args::IndexedInput;
use crate::commands::io;
use crate::filltags::{self, FillTags};
use crate::hts::HeaderView;
use crate::processor::RecordFn;

#[derive(Args, Debug)]
pub struct FillTagsArgs {
    #[command(flatten)]
    pub indexed: IndexedInput,

    /// Output path. Indexed input writes bgzf VCF and requires this; stdin
    /// input writes uncompressed VCF to stdout when omitted or "-".
    #[arg(long = "output", short = 'o')]
    pub output: Option<String>,

    /// Window size in base pairs for parallel processing. Defaults to 16384 bp (one
    /// tabix/CSI index bin).
    #[arg(long)]
    pub window_size: Option<u64>,

    /// Comma-separated tags to fill (AN,AC,AF,MAF,NS,AC_Hom,AC_Het,AC_Hemi,
    /// HWE,ExcHet,F_MISSING,VAF,VAF1) or "all".
    #[arg(long = "tags", short = 't', default_value = "all")]
    pub tags: String,

    /// Do not count half-missing genotypes (e.g. "./1") as hemizygous.
    #[arg(long = "drop-missing", short = 'd')]
    pub drop_missing: bool,
}

pub fn run(args: FillTagsArgs) -> Result<()> {
    let mut tags = filltags::parse_tags(&args.tags)?;
    tags.drop_missing = args.drop_missing;
    let lines = tags.header_lines();
    let factory = move |h: &Arc<HeaderView>| -> Result<RecordFn> {
        let mut fill = FillTags::new(&tags, h)?;
        Ok(Box::new(move |rec: &mut Record| -> Result<bool> {
            fill.fill(rec)?;
            Ok(true)
        }))
    };
    io::run_records("fill-tags", &args.indexed, &args.output, args.window_size, None, move |_| Ok(lines), factory)
}
