use std::path::PathBuf;

use anyhow::{Context, Result, anyhow};
use clap::Args;
use indexmap::IndexMap;

/// Shared CLI arguments for subcommands that operate on indexed VCF/BCF files.
#[derive(Args, Debug)]
pub struct IndexedInput {
    /// Variant file (must be indexed: vcf.gz+tbi or bcf+csi).
    /// Reads stdin when omitted or set to "-".
    #[arg()]
    pub input: Option<PathBuf>,

    /// Worker thread count. Defaults to system parallelism.
    #[arg(long)]
    pub threads: Option<usize>,

    /// Restrict processing to the given contig(s). Repeatable.
    #[arg(long)]
    pub contig: Vec<String>,

    /// FASTA index for contig names and lengths, used when the VCF header
    /// lacks contig length information or when custom contig ordering is
    /// desired.
    #[arg(long)]
    pub fai: Option<PathBuf>,
}

impl IndexedInput {
    pub fn resolve_threads(&self) -> usize {
        self.threads.unwrap_or_else(|| {
            std::thread::available_parallelism()
                .map(usize::from)
                .unwrap_or(1)
        })
    }
}

/// Parsed .fai data preserving line order.
pub struct FaiData {
    pub ordered: Vec<(String, u64)>,
    pub lengths: IndexMap<String, u64>,
}

pub fn parse_fai(path: &PathBuf) -> Result<FaiData> {
    let content = std::fs::read_to_string(path)
        .with_context(|| format!("failed to read fai: {}", path.display()))?;
    let mut ordered = Vec::new();
    let mut map = IndexMap::new();
    for line in content.lines() {
        let mut fields = line.split('\t');
        if let (Some(name), Some(len_str)) = (fields.next(), fields.next()) {
            let len = len_str.parse::<u64>().with_context(|| {
                format!("invalid length in fai line: {line}")
            })?;
            ordered.push((name.to_string(), len));
            map.insert(name.to_string(), len);
        }
    }
    if map.is_empty() {
        return Err(anyhow!("empty or invalid fai: {}", path.display()));
    }
    Ok(FaiData { ordered, lengths: map })
}

