//! Parallel, regioned variant processor.
//!
//! Reads an indexed VCF/BCF, splits each contig into fixed-size windows, and
//! processes the windows in parallel: each worker opens its own indexed reader,
//! queries its assigned region, applies a user-supplied callback to every
//! record, and serializes the results to an in-memory uncompressed VCF buffer.
//! The main thread re-orders the per-window buffers and feeds them to a
//! `bgzf::MultithreadedWriter` so the output is correctly ordered and
//! bgzf-compressed in parallel.

use std::collections::BTreeMap;
use std::fs::File;
use std::io::Write as _;
use std::num::NonZeroUsize;
use std::path::PathBuf;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow};
use crossbeam_channel::bounded;
use noodles_bgzf as bgzf;
use noodles_core::{Position, Region};
use noodles_util::variant;
use noodles_vcf as vcf;
use vcf::variant::RecordBuf;
use vcf::variant::io::Write as _;

/// Default per-window size, in base pairs.
pub const DEFAULT_WINDOW_SIZE: u64 = 1_000_000;

type RecordCallback = dyn Fn(RecordBuf) -> Option<RecordBuf> + Send + Sync + 'static;
type ProgressCallback = dyn Fn(usize) + Send + Sync + 'static;
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
    /// have been written to the bgzf stream. Argument is the cumulative number
    /// of bp covered by all written-and-completed windows so far.
    pub fn progress_callback<F>(mut self, f: F) -> Self
    where
        F: Fn(usize) + Send + Sync + 'static,
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
        let record_callback = self
            .record_callback
            .ok_or_else(|| anyhow!("missing record_callback"))?;

        ParallelVariantWindowProcessor {
            input,
            output: self.output,
            worker_threads,
            window_size,
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
        let input_header = reader.read_header().context("failed to read header")?;
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
        let windows = enumerate_windows(&input_header, window_size, indexed_contigs.as_ref());

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
            while let Ok((idx, bp, result)) = out_rx.recv() {
                let bytes = result.with_context(|| format!("worker failed on window {idx}"))?;
                buffer.insert(idx, (bp, bytes));
                while let Some((bp, bytes)) = buffer.remove(&next_idx) {
                    if let Some(w) = mt_writer.as_mut() {
                        w.write_all(&bytes)
                            .context("failed to write window bytes")?;
                    }
                    cum_bp += bp;
                    if let Some(cb) = progress_callback.as_ref() {
                        cb(cum_bp as usize);
                    }
                    next_idx += 1;
                }
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
    indexed_contigs: Option<&std::collections::HashSet<Vec<u8>>>,
) -> Vec<Window> {
    let mut out = Vec::new();
    for (contig, map) in header.contigs().iter() {
        let Some(length) = map.length() else {
            eprintln!("pygopus: contig {contig} has no length in header; skipping");
            continue;
        };
        if let Some(set) = indexed_contigs {
            if !set.contains(contig.as_bytes()) {
                continue;
            }
        }
        let length = length as u64;
        let mut start = 1u64;
        while start <= length {
            let end = (start + window_size - 1).min(length);
            out.push(Window {
                idx: out.len(),
                contig: Box::from(contig.as_str()),
                start,
                end,
            });
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
    match output_header {
        Some(out_header) => {
            let mut writer = vcf::io::Writer::new(&mut bytes);
            for result in query {
                let record = result.context("error reading record")?;
                let Some(pos) = take_pos_if_in_window(record.as_ref(), window)? else {
                    continue;
                };
                let _ = pos;
                let buf = RecordBuf::try_from_variant_record(input_header, record.as_ref())
                    .context("failed to materialize RecordBuf")?;
                if let Some(out) = record_callback(buf) {
                    writer
                        .write_variant_record(out_header, &out)
                        .context("failed to serialize record")?;
                }
            }
        }
        None => {
            // Side-effect-only mode: never allocate writer, never serialize.
            for result in query {
                let record = result.context("error reading record")?;
                let Some(pos) = take_pos_if_in_window(record.as_ref(), window)? else {
                    continue;
                };
                let _ = pos;
                let buf = RecordBuf::try_from_variant_record(input_header, record.as_ref())
                    .context("failed to materialize RecordBuf")?;
                let _ = record_callback(buf);
            }
        }
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
