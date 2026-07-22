//! Parallel, regioned variant processor.
//!
//! Reads an indexed VCF/BCF, splits each contig into fixed-size windows, and
//! processes the windows in parallel: each worker opens its own indexed reader,
//! queries its assigned region, applies a user-supplied callback to every
//! record, and serializes the results to an in-memory uncompressed VCF buffer.
//! The main thread re-orders the per-window buffers and feeds them to a
//! `bgzf::MultithreadedWriter` so the output is correctly ordered and
//! bgzf-compressed in parallel.

use std::collections::{BTreeMap, HashSet};
use std::fs::File;
use std::io::Write as _;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};
use crossbeam_channel::bounded;
use indexmap::IndexMap;
use noodles_bgzf as bgzf;
use noodles_core::{Position, Region};
use noodles_util::variant;
use noodles_vcf as vcf;
use vcf::variant::RecordBuf;
use vcf::variant::io::Write as _;

#[allow(dead_code)]
fn dump_record_fields(record: &dyn vcf::variant::Record, header: &vcf::Header) -> String {
    use std::fmt::Write;

    let chrom = record
        .reference_sequence_name(header)
        .unwrap_or("<err>");
    let pos = match record.variant_start() {
        Some(Ok(p)) => p.get().to_string(),
        _ => String::from("."),
    };

    let mut ref_bases = String::new();
    for b in record.reference_bases().iter() {
        match b {
            Ok(c) => ref_bases.push(c as char),
            Err(_) => ref_bases.push_str("<err>"),
        }
    }
    if ref_bases.is_empty() {
        ref_bases.push('.');
    }

    let mut alts = String::new();
    for (i, a) in record.alternate_bases().iter().enumerate() {
        if i > 0 {
            alts.push(',');
        }
        match a {
            Ok(s) => alts.push_str(s),
            Err(_) => alts.push_str("<err>"),
        }
    }
    if alts.is_empty() {
        alts.push('.');
    }

    let qual = match record.quality_score() {
        Some(Ok(q)) => format!("{q}"),
        _ => String::from("."),
    };

    let mut buf = String::new();
    write!(buf, "record dump: {chrom}\t{pos}\t.\t{ref_bases}\t{alts}\t{qual}\t.\t").ok();

    let info = record.info();
    let iter = info.iter(header);
    {
        let mut first = true;
        for result in iter {
            match result {
                Ok((key, value)) => {
                    if !first { buf.push(';'); }
                    first = false;
                    match value {
                        Some(v) => write!(buf, "{key}={v:?}").ok(),
                        None => write!(buf, "{key}").ok(),
                    };
                }
                Err(e) => {
                    if !first { buf.push(';'); }
                    first = false;
                    write!(buf, "<info_field_err:{e}>").ok();
                }
            }
        }
        if first { buf.push('.'); }
    }

    // FORMAT + sample summary
    match record.samples() {
        Ok(samples) => {
            let keys = samples.column_names(header);
            write!(buf, "\t").ok();
            let mut first = true;
            for k in keys {
                if !first { buf.push(':'); }
                first = false;
                match k {
                    Ok(name) => { buf.push_str(name); }
                    Err(_) => { buf.push_str("<err>"); }
                }
            }
            write!(buf, "\t{} samples", samples.len()).ok();
        }
        Err(e) => {
            write!(buf, "\t<samples_err:{e}>").ok();
        }
    }

    buf
}

/// Default per-window size, in base pairs.
pub const DEFAULT_WINDOW_SIZE: u64 = 1_000_000;

type RecordCallback = dyn Fn(RecordBuf) -> Option<RecordBuf> + Send + Sync + 'static;
type ProgressCallback = dyn Fn(u64, usize, usize) + Send + Sync + 'static;
type HeaderCallback = Box<dyn FnOnce(&mut vcf::Header) + Send + 'static>;

/// One contiguous slice of a single contig.
#[derive(Clone, Debug)]
struct Window {
    idx: usize,
    contig: Box<str>,
    start: u64,
    end: u64,
}

impl Window {
    fn bp(&self) -> u64 {
        self.end - self.start + 1
    }
}

/// Builder for [`ParallelVariantWindowProcessor`].
#[derive(Default)]
pub struct ParallelVariantWindowProcessorBuilder {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    worker_threads: Option<usize>,
    window_size: Option<u64>,
    stride: Option<usize>,
    contig_filter: Option<HashSet<String>>,
    contig_lengths: Option<IndexMap<String, u64>>,
    record_callback: Option<Arc<RecordCallback>>,
    header_callback: Option<HeaderCallback>,
    progress_callback: Option<Arc<ProgressCallback>>,
}

impl ParallelVariantWindowProcessorBuilder {
    pub fn input(mut self, path: impl Into<PathBuf>) -> Self {
        self.input = Some(path.into());
        self
    }

    /// Set the output `.vcf.gz` path. Optional: when omitted, the processor
    /// runs in side-effect-only mode — workers apply `record_callback` to
    /// every variant but no VCF is written and records aren't serialized.
    pub fn with_output_file(mut self, path: impl Into<PathBuf>) -> Self {
        self.output = Some(path.into());
        self
    }

    pub fn worker_threads(mut self, n: usize) -> Self {
        self.worker_threads = Some(n);
        self
    }

    pub fn window_size(mut self, bp: u64) -> Self {
        self.window_size = Some(bp);
        self
    }

    /// Process every Nth window. Default 1 (process all windows).
    pub fn stride(mut self, n: usize) -> Self {
        self.stride = Some(n);
        self
    }

    /// Limit processing to the given contig(s). When omitted, all contigs are processed.
    pub fn contigs(mut self, contigs: Vec<String>) -> Self {
        self.contig_filter = Some(contigs.into_iter().collect());
        self
    }

    /// Supply contig name -> length pairs from an external source (e.g. .fai).
    /// These are used as fallback lengths for contigs that lack length info in
    /// the VCF header, and as a source of additional contigs not listed in the
    /// header.
    pub fn contig_lengths(mut self, lengths: IndexMap<String, u64>) -> Self {
        self.contig_lengths = Some(lengths);
        self
    }

    /// Required. Applied to each `RecordBuf`; `Some(buf)` is written, `None` drops the record.
    pub fn record_callback<F>(mut self, f: F) -> Self
    where
        F: Fn(RecordBuf) -> Option<RecordBuf> + Send + Sync + 'static,
    {
        self.record_callback = Some(Arc::new(f));
        self
    }

    /// Optional. Mutates the input header before it's written to the output.
    pub fn header_callback<F>(mut self, f: F) -> Self
    where
        F: FnOnce(&mut vcf::Header) + Send + 'static,
    {
        self.header_callback = Some(Box::new(f));
        self
    }

    /// Optional. Invoked on the main writer thread after each window's bytes
    /// have been written to the bgzf stream. Arguments: (cumulative bp,
    /// completed windows, total windows).
    pub fn progress_callback<F>(mut self, f: F) -> Self
    where
        F: Fn(u64, usize, usize) + Send + Sync + 'static,
    {
        self.progress_callback = Some(Arc::new(f));
        self
    }

    pub fn run(self) -> Result<()> {
        let input = self.input.ok_or_else(|| anyhow!("missing input path"))?;
        let worker_threads = self
            .worker_threads
            .ok_or_else(|| anyhow!("missing worker_threads"))?
            .max(1);
        let window_size = self.window_size.unwrap_or(DEFAULT_WINDOW_SIZE).max(1);
        let stride = self.stride.unwrap_or(1).max(1);
        let record_callback = self
            .record_callback
            .ok_or_else(|| anyhow!("missing record_callback"))?;

        ParallelVariantWindowProcessor {
            input,
            output: self.output,
            worker_threads,
            window_size,
            stride,
            contig_filter: self.contig_filter,
            contig_lengths: self.contig_lengths,
            record_callback,
            header_callback: self.header_callback,
            progress_callback: self.progress_callback,
        }
        .execute()
    }
}

pub struct ParallelVariantWindowProcessor {
    input: PathBuf,
    output: Option<PathBuf>,
    worker_threads: usize,
    window_size: u64,
    stride: usize,
    contig_filter: Option<HashSet<String>>,
    contig_lengths: Option<IndexMap<String, u64>>,
    record_callback: Arc<RecordCallback>,
    header_callback: Option<HeaderCallback>,
    progress_callback: Option<Arc<ProgressCallback>>,
}

impl ParallelVariantWindowProcessor {
    pub fn builder() -> ParallelVariantWindowProcessorBuilder {
        ParallelVariantWindowProcessorBuilder::default()
    }

    fn execute(self) -> Result<()> {
        let Self {
            input,
            output,
            worker_threads,
            window_size,
            stride,
            contig_filter,
            contig_lengths,
            record_callback,
            header_callback,
            progress_callback,
        } = self;

        // 1. Open input + read header. Capture the set of contigs the index
        //    actually knows about so we can skip header-declared contigs that
        //    aren't represented in the file (the query would otherwise error).
        //    Tabix exposes this via `index.header().reference_sequence_names()`;
        //    CSI (BCF) returns `None` and we fall back to no pre-filtering.
        let mut reader = variant::io::indexed_reader::Builder::default()
            .build_from_path(&input)
            .with_context(|| format!("failed to open indexed reader: {}", input.display()))?;
        let mut input_header = reader.read_header().context("failed to read header")?;
        crate::util::normalize_header_for_noodles(&mut input_header);
        let indexed_contigs: Option<std::collections::HashSet<Vec<u8>>> = reader
            .index()
            .header()
            .map(|h| h.reference_sequence_names().iter().map(|b| b.to_vec()).collect());
        drop(reader); // workers each open their own.

        // 2. Build output header.
        let mut output_header = input_header.clone();
        if let Some(cb) = header_callback {
            cb(&mut output_header);
        }
        let input_header = Arc::new(input_header);
        let output_header = Arc::new(output_header);

        // 3. Enumerate windows from the input header's contigs.
        let windows = enumerate_windows(&input_header, window_size, stride, indexed_contigs.as_ref(), contig_filter.as_ref(), contig_lengths.as_ref());

        // 4. Open output and write the header through a multithreaded bgzf
        //    writer. If no output path was provided, run in side-effect-only
        //    mode: workers still apply `record_callback` (so the user can
        //    accumulate state, count, push to a channel, etc.), but no bytes
        //    are serialized or written.
        let mut mt_writer = match output.as_ref() {
            Some(path) => {
                let file = File::create(path)
                    .with_context(|| format!("failed to create output: {}", path.display()))?;
                let mut w = bgzf::io::multithreaded_writer::Builder::default()
                    .set_worker_count(NonZeroUsize::new(worker_threads).unwrap())
                    .build_from_writer(file);
                {
                    let mut writer = vcf::io::Writer::new(&mut w);
                    writer
                        .write_header(&output_header)
                        .context("failed to write VCF header")?;
                }
                Some(w)
            }
            None => None,
        };
        let serialize = mt_writer.is_some();

        // 5. Set up channels.
        let (work_tx, work_rx) = bounded::<Window>(worker_threads * 4);
        type Out = (usize, u64, Result<Vec<u8>>);
        let (out_tx, out_rx) = bounded::<Out>(worker_threads * 4);

        // 6. Run the producer + workers + writer in a scope.
        std::thread::scope(|scope| -> Result<()> {
            // Move out_rx into the scope so it's dropped when the writer
            // finishes (or errors), unblocking workers' sends.
            let out_rx = out_rx;

            // Producer: enumerate windows.
            let work_tx_producer = work_tx.clone();
            let windows_for_producer = windows.clone();
            scope.spawn(move || {
                for w in windows_for_producer {
                    if work_tx_producer.send(w).is_err() {
                        break;
                    }
                }
            });
            drop(work_tx);

            // Workers.
            for _ in 0..worker_threads {
                let work_rx = work_rx.clone();
                let out_tx = out_tx.clone();
                let input = input.clone();
                let input_header = Arc::clone(&input_header);
                let output_header = Arc::clone(&output_header);
                let record_callback = Arc::clone(&record_callback);
        scope.spawn(move || {
            let mut reader = match variant::io::indexed_reader::Builder::default()
                .build_from_path(&input)
            {
            Ok(r) => r,
            Err(e) => {
                    // Send the error tagged to whatever window we were going to handle.
                    // Drain in the meantime so the producer doesn't block.
                    while let Ok(w) = work_rx.recv() {
                        let _ = out_tx.send((
                            w.idx,
                            w.bp(),
                            Err(anyhow!(
                                "worker failed to open indexed reader for {}: {e}",
                                input.display()
                            )),
                        ));
                    }
                    return;
                }
            };

            while let Ok(window) = work_rx.recv() {
                let bp = window.bp();
                let result = process_window(
                    &mut reader,
                    &input_header,
                    if serialize { Some(output_header.as_ref()) } else { None },
                    &window,
                    record_callback.as_ref(),
                );
                if out_tx.send((window.idx, bp, result)).is_err() {
                    break;
                }
            }
        });
            }
            drop(out_tx);

            // Writer / re-orderer / progress (main scope thread).
            let mut buffer: BTreeMap<usize, (u64, Vec<u8>)> = BTreeMap::new();
            let mut next_idx = 0usize;
            let mut cum_bp: u64 = 0;
            let mut window_errors: Vec<(usize, anyhow::Error)> = Vec::new();
            let total_windows = windows.len();
            while let Ok((idx, bp, result)) = out_rx.recv() {
                match result {
                    Ok(bytes) => {
                        buffer.insert(idx, (bp, bytes));
                    }
                    Err(e) => {
                        window_errors.push((idx, e));
                        buffer.insert(idx, (bp, Vec::new()));
                    }
                }
                while let Some((bp, bytes)) = buffer.remove(&next_idx) {
                    if let Some(w) = mt_writer.as_mut() {
                        w.write_all(&bytes)
                            .context("failed to write window bytes")?;
                    }
                    cum_bp += bp;
                    if let Some(cb) = progress_callback.as_ref() {
                        cb(cum_bp, next_idx + 1, total_windows);
                    }
                    next_idx += 1;
                }
            }

            if !window_errors.is_empty() {
                for (idx, e) in &window_errors {
                    eprintln!("varanrs: window {idx} failed: {e}");
                }
                return Err(anyhow!(
                    "{} of {} windows failed (first: window {}: {})",
                    window_errors.len(),
                    windows.len(),
                    window_errors[0].0,
                    window_errors[0].1,
                ));
            }

            if next_idx != windows.len() {
                return Err(anyhow!(
                    "writer terminated with gaps: wrote {next_idx} of {} windows",
                    windows.len()
                ));
            }
            Ok(())
        })?;

        if let Some(mut w) = mt_writer {
            w.finish().context("failed to finalize bgzf writer")?;
        }
        Ok(())
    }
}

fn enumerate_windows(
    header: &vcf::Header,
    window_size: u64,
    stride: usize,
    indexed_contigs: Option<&HashSet<Vec<u8>>>,
    contig_filter: Option<&HashSet<String>>,
    contig_lengths: Option<&IndexMap<String, u64>>,
) -> Vec<Window> {
    let mut out = Vec::new();
    let mut abs_idx = 0usize;

    let header_contig_names: HashSet<&str> =
        header.contigs().iter().map(|(n, _)| n.as_str()).collect();

    // Build ordered list of (name, length). Header contigs come first in header
    // order; then contigs from contig_lengths not in the header, in their
    // insertion order.
    let mut ordered: Vec<(String, u64)> = Vec::new();
    for (contig, map) in header.contigs().iter() {
        let length = map
            .length()
            .or_else(|| {
                contig_lengths
                    .and_then(|cl| cl.get(contig.as_str()).copied())
                    .and_then(|l| usize::try_from(l).ok())
            });
        if let Some(len) = length {
            ordered.push((contig.to_string(), len as u64));
        } else {
            eprintln!("varanrs: contig {contig} has no length in header or fai; skipping");
        }
    }
    if let Some(cl) = contig_lengths {
        for (name, len) in cl {
            if !header_contig_names.contains(name.as_str()) {
                ordered.push((name.clone(), *len));
            }
        }
    }

    for (contig, length) in &ordered {
        if let Some(filter) = contig_filter
            && !filter.contains(contig.as_str())
        {
            continue;
        }
        if let Some(set) = indexed_contigs
            && !set.contains(contig.as_bytes())
        {
            continue;
        }
        let length = *length;
        let mut start = 1u64;
        while start <= length {
            let end = (start + window_size - 1).min(length);
            if abs_idx.is_multiple_of(stride) {
                out.push(Window {
                    idx: out.len(),
                    contig: Box::from(contig.as_str()),
                    start,
                    end,
                });
            }
            abs_idx += 1;
            start = end + 1;
        }
    }
    out
}

fn process_window<R>(
    reader: &mut variant::io::IndexedReader<R>,
    input_header: &vcf::Header,
    output_header: Option<&vcf::Header>,
    window: &Window,
    record_callback: &RecordCallback,
) -> Result<Vec<u8>>
where
    R: bgzf::io::BufRead + bgzf::io::Seek,
{
    let start = Position::try_from(window.start as usize).context("invalid window start")?;
    let end = Position::try_from(window.end as usize).context("invalid window end")?;
    let region = Region::new(window.contig.as_ref(), start..=end);

    let query = reader
        .query(input_header, &region)
        .with_context(|| format!("query failed for {}:{}-{}", window.contig, window.start, window.end))?;

    let mut bytes = Vec::new();
    let mut skipped = 0usize;
    match output_header {
        Some(out_header) => {
            let mut writer = vcf::io::Writer::new(&mut bytes);
            for result in query {
                let record = match result {
                    Ok(r) => r,
                    Err(e) => {
                        skipped += 1;
                        eprintln!(
                            "varanrs: window {} {}:{}-{}: error reading record: {e}",
                            window.idx, window.contig, window.start, window.end
                        );
                        continue;
                    }
                };
                let _pos = match take_pos_if_in_window(record.as_ref(), window) {
                    Ok(Some(p)) => p,
                    Ok(None) => continue,
                    Err(e) => {
                        skipped += 1;
                        eprintln!(
                            "varanrs: window {} {}:{}-{}: skipping record ({e})",
                            window.idx, window.contig, window.start, window.end
                        );
                        continue;
                    }
                };
                let buf = match RecordBuf::try_from_variant_record(input_header, record.as_ref()) {
                    Ok(b) => b,
                    Err(e) => {
                        skipped += 1;
                        if skipped <= 5 {
                            eprintln!(
                                "varanrs: window {} {}:{}-{}: skipping record (RecordBuf conversion failed: {e})",
                                window.idx, window.contig, window.start, window.end
                            );
                        }
                        continue;
                    }
                };
                if let Some(out) = record_callback(buf) {
                    writer
                        .write_variant_record(out_header, &out)
                        .context("failed to serialize record")?;
                }
            }
        }
        None => {
            for result in query {
                let record = match result {
                    Ok(r) => r,
                    Err(e) => {
                        skipped += 1;
                        eprintln!(
                            "varanrs: window {} {}:{}-{}: error reading record: {e}",
                            window.idx, window.contig, window.start, window.end
                        );
                        continue;
                    }
                };
                let _pos = match take_pos_if_in_window(record.as_ref(), window) {
                    Ok(Some(p)) => p,
                    Ok(None) => continue,
                    Err(e) => {
                        skipped += 1;
                        eprintln!(
                            "varanrs: window {} {}:{}-{}: skipping record ({e})",
                            window.idx, window.contig, window.start, window.end
                        );
                        continue;
                    }
                };
                let buf = match RecordBuf::try_from_variant_record(input_header, record.as_ref()) {
                    Ok(b) => b,
                    Err(e) => {
                        skipped += 1;
                        if skipped <= 5 {
                            eprintln!(
                                "varanrs: window {} {}:{}-{}: skipping record (RecordBuf conversion failed: {e})",
                                window.idx, window.contig, window.start, window.end
                            );
                        }
                        continue;
                    }
                };
                let _ = record_callback(buf);
            }
        }
    }
    if skipped > 5 {
        eprintln!(
            "varanrs: window {} {}:{}-{}: {} further records skipped",
            window.idx, window.contig, window.start, window.end,
            skipped - 5
        );
    }
    Ok(bytes)
}

/// Boundary dedup: tabix queries return any record whose span overlaps the
/// region, so a record starting *before* this window will also be returned by
/// the previous window. Keep it only in the window where its POS lands.
fn take_pos_if_in_window(
    record: &dyn vcf::variant::Record,
    window: &Window,
) -> Result<Option<Position>> {
    let pos = match record.variant_start() {
        Some(r) => r.context("missing POS")?,
        None => return Ok(None), // unplaced record; skip.
    };
    if (pos.get() as u64) < window.start {
        return Ok(None);
    }
    Ok(Some(pos))
}
