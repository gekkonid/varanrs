//! varanrs CLI entrypoint. Git-style subcommands.

use anyhow::Result;
use clap::{Parser, Subcommand};
use varanrs::commands;

#[derive(Parser, Debug)]
#[command(name = "varanrs", about = "Bioinformatics toolkit", version)]
struct Cli {
    #[command(subcommand)]
    command: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Per-site allele filtering by minimum AC and/or AF.
    Allelefilter(commands::allelefilter::AllelefilterArgs),
    /// Subsample an indexed VCF/BCF, estimate pairwise distances and per-sample missingness.
    Snpsketch(commands::snpsketch::SketchArgs),
    /// Force REF and ALT alleles to uppercase (GLnexus workaround).
    UppercaseAlleles(commands::uppercase_alleles::UppercaseAllelesArgs),
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Allelefilter(args) => commands::allelefilter::run(args),
        Command::Snpsketch(args) => commands::snpsketch::run(args),
        Command::UppercaseAlleles(args) => commands::uppercase_alleles::run(args),
    }
}
