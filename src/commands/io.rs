//! Shared input/output plumbing for record-transforming subcommands.
//!
//! A subcommand supplies:
//! * `prepare`: given the input header, the `##` lines to add to it;
//! * `factory`: given a (worker's) header that already contains those lines,
//!   a per-record function (`Ok(true)` keeps the record).
//!
//! The same pair drives either the streaming path (stdin, any VCF/BCF,
//! single-threaded, uncompressed VCF out) or the region-parallel path
//! (indexed input, [`ParallelVariantWindowProcessor`], bgzf VCF out). Records
//! are parsed and written by htslib, so output is byte-identical to bcftools.

use std::io::{BufWriter, Write};
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};

use crate::args::{self, IndexedInput};
use crate::hts::{self, HeaderView, KString, SampleSubset, VcfReader};
use crate::processor::{ParallelVariantWindowProcessor, RecordFn};

pub fn is_std(path: &Option<String>) -> bool {
    path.as_deref().is_none_or(|s| s == "-")
}

pub fn open_writer(output: &Option<String>) -> Result<Box<dyn Write>> {
    if is_std(output) {
        Ok(Box::new(BufWriter::new(std::io::stdout().lock())))
    } else {
        let path = output.as_ref().unwrap();
        Ok(Box::new(BufWriter::new(
            std::fs::File::create(path).with_context(|| format!("failed to create output: {path}"))?,
        )))
    }
}

/// Stream every record of `input` ("-" for stdin; VCF or BCF, any
/// compression) through the record function into `writer` as VCF text.
pub fn stream_records<W, P, F>(
    input: &str,
    samples: Option<&SampleSubset>,
    mut writer: W,
    prepare: P,
    factory: F,
) -> Result<()>
where
    W: Write,
    P: FnOnce(&HeaderView) -> Result<Vec<String>>,
    F: Fn(&Arc<HeaderView>) -> Result<RecordFn>,
{
    let mut reader = VcfReader::open(input)?;
    subset_samples(&mut reader, samples)?;
    let lines = prepare(reader.header())?;
    reader.append_header_lines(&lines)?;
    writer.write_all(&hts::format_header(reader.header())?)?;
    let mut f = factory(reader.header())?;
    // Alternate two record buffers: bcftools' filter.c caches per-record data
    // keyed by the record's address.
    let mut recs = [reader.empty_record(), reader.empty_record()];
    let mut cur = 0usize;
    let mut ks = KString::new();
    let mut out = Vec::with_capacity(1 << 16);
    loop {
        cur ^= 1;
        let rec = &mut recs[cur];
        if !reader.read(rec)? {
            break;
        }
        if f(rec).with_context(|| format!("at {}", hts::rec_locus(rec)))? {
            out.clear();
            hts::format_record(rec, &mut ks, &mut out)?;
            writer.write_all(&out)?;
        }
    }
    writer.flush()?;
    Ok(())
}

/// Apply sample subsetting (if any) to a reader, reporting warnings.
fn subset_samples(reader: &mut VcfReader, samples: Option<&SampleSubset>) -> Result<()> {
    if let Some(s) = samples {
        for w in reader.subset_samples(s)? {
            eprintln!("Warning: {w}");
        }
    }
    Ok(())
}

/// Run over stdin (streaming) or an indexed input (parallel).
pub fn run_records<P, F>(
    name: &str,
    indexed: &IndexedInput,
    output: &Option<String>,
    window_size: Option<u64>,
    samples: Option<SampleSubset>,
    prepare: P,
    factory: F,
) -> Result<()>
where
    P: FnOnce(&HeaderView) -> Result<Vec<String>>,
    F: Fn(&Arc<HeaderView>) -> Result<RecordFn> + Send + Sync + 'static,
{
    let Some(input) = indexed.input.as_ref().filter(|p| p.as_os_str() != "-") else {
        return stream_records("-", samples.as_ref(), open_writer(output)?, prepare, factory);
    };

    if is_std(output) {
        return Err(anyhow!(
            "{name} with indexed input requires --output (cannot write bgzf to stdout)"
        ));
    }

    let lines = {
        let mut reader = VcfReader::open_indexed(input)?;
        subset_samples(&mut reader, samples.as_ref())?;
        let lines = prepare(reader.header())?;
        // Build the record function once up front so invalid arguments (e.g.
        // an expression naming an undeclared tag) fail even when the input
        // has no records for any worker to process.
        reader.append_header_lines(&lines)?;
        drop(factory(reader.header())?);
        lines
    };

    let mut builder = ParallelVariantWindowProcessor::builder()
        .input(input)
        .worker_threads(indexed.resolve_threads())
        .header_lines(lines)
        .samples(samples)
        .record_factory(factory)
        .with_output_file(output.as_ref().unwrap())
        .progress_callback(|bp, _, _| {
            eprint!("\r  {} Mbp processed", bp / 1_000_000);
            std::io::stderr().flush().ok();
        });
    if let Some(ws) = window_size {
        builder = builder.window_size(ws);
    }
    if !indexed.contig.is_empty() {
        builder = builder.contigs(indexed.contig.clone());
    }
    if let Some(ref fai_path) = indexed.fai {
        let fai = args::parse_fai(fai_path)?;
        builder = builder.contig_lengths(fai.lengths);
    }
    let result = builder.run();
    eprintln!();
    result
}
