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
    Filter(commands::filter::FilterArgs),
    /// Force REF and ALT alleles to uppercase (GLnexus workaround).
    Upper(commands::upper::UpperArgs),
}

fn main() -> Result<()> {
    match Cli::parse().command {
        Command::Filter(args) => commands::filter::run(args),
        Command::Upper(args) => commands::upper::run(args),
    }
}
