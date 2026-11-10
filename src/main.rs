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
    /// Filter sites with bcftools-style expressions, allele filters and tag refilling (parallel).
    Filter(commands::filter::FilterArgs),
    /// Recompute INFO tags (AC, AN, AF, MAF, NS, HWE, ...) from genotypes, like bcftools +fill-tags.
    FillTags(commands::fill_tags::FillTagsArgs),
    /// Subsample an indexed VCF/BCF, estimate pairwise distances and per-sample missingness.
    Snpsketch(commands::snpsketch::SketchArgs),
    /// Force REF and ALT alleles to uppercase (GLnexus workaround).
    UppercaseAlleles(commands::uppercase_alleles::UppercaseAllelesArgs),
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Allelefilter(args) => commands::allelefilter::run(args),
        Command::Filter(args) => commands::filter::run(args),
        Command::FillTags(args) => commands::fill_tags::run(args),
        Command::Snpsketch(args) => commands::snpsketch::run(args),
        Command::UppercaseAlleles(args) => commands::uppercase_alleles::run(args),
    }
}
