//! Parallel, regioned variant processor.
//!
//! Reads an indexed VCF/BCF (htslib), splits each contig into fixed-size
//! windows, and processes the windows in parallel: each worker opens its own
//! indexed reader, queries its assigned region, applies a per-worker record
//! function to every record, and serialises the kept records to an in-memory
//! VCF text buffer (htslib `vcf_format`, so byte-identical to bcftools). The
//! main thread re-orders the per-window buffers and feeds them to a
//! multithreaded BGZF writer, so the output is correctly ordered and
//! compressed in parallel.
//!
//! Contigs are processed in file order (the order of their data in the
//! index), matching `bcftools view`.
//!
//! Memory is bounded: a window is only dispatched once it holds one of
//! `max_in_flight` permits (default 2 × worker threads), and the permit is
//! returned only after the window has been written in order. So at most
//! `max_in_flight` windows exist at once, whether queued, being processed, or
//! waiting in the re-order buffer, however slow any single window is.

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
use rust_htslib::bcf::Record;

use crate::hts::{self, HeaderView, KString, SampleSubset, VcfReader};

/// Default per-window size, in base pairs: 2^14, one tabix/CSI linear-index
/// bin (min_shift = 14).
pub const DEFAULT_WINDOW_SIZE: u64 = 16_384;

/// Contig positions are 1-based `i32` in VCF/BCF; windows beyond a contig's
/// declared length are searched (via the index) up to this position.
const MAX_POS: u64 = i32::MAX as u64;

/// A per-worker record function: `Ok(true)` keeps (writes) the record,
/// `Ok(false)` drops it, `Err` aborts the run.
pub type RecordFn = Box<dyn FnMut(&mut Record) -> Result<bool>>;
type RecordFactory = dyn Fn(&Arc<HeaderView>) -> Result<RecordFn> + Send + Sync;
type ProgressCallback = dyn Fn(u64, usize, usize) + Send + Sync + 'static;

/// One unit of work, in output order: either a window to query, or a span of
/// windows the index proves to be empty (nothing is read for it).
#[derive(Clone, Debug)]
struct Window {
    /// Sequential item index, used to re-order output.
    idx: usize,
    contig: Arc<str>,
    start: u64,
    end: u64,
    /// Number of (stride-selected) windows this item stands for.
    n_windows: usize,
    /// Base pairs of those windows, for progress reporting.
    bp: u64,
    /// True when the index shows no record can overlap this span.
    skip: bool,
}

/// A contig to process, with its window numbering.
#[derive(Clone, Debug)]
struct ContigPlan {
    name: Arc<str>,
    /// Declared length (header or .fai); 0 when unknown. Records beyond it
    /// are still processed, but do not count towards progress totals.
    length: u64,
    /// Absolute index (across all contigs) of this contig's first window,
    /// used for stride selection.
    first_abs: usize,
}

impl ContigPlan {
    fn n_windows(&self, window_size: u64) -> usize {
        self.length.div_ceil(window_size) as usize
    }
}

/// Builder for [`ParallelVariantWindowProcessor`].
#[derive(Default)]
pub struct ParallelVariantWindowProcessorBuilder {
    input: Option<PathBuf>,
    output: Option<PathBuf>,
    worker_threads: Option<usize>,
    max_in_flight: Option<usize>,
    window_size: Option<u64>,
    stride: Option<usize>,
    contig_filter: Option<HashSet<String>>,
    contig_lengths: Option<IndexMap<String, u64>>,
    record_factory: Option<Arc<RecordFactory>>,
    header_lines: Vec<String>,
    samples: Option<SampleSubset>,
    progress_callback: Option<Arc<ProgressCallback>>,
}

impl ParallelVariantWindowProcessorBuilder {
    pub fn input(mut self, path: impl Into<PathBuf>) -> Self {
        self.input = Some(path.into());
        self
    }

    /// Set the output `.vcf.gz` path. Optional: when omitted, the processor
    /// runs in side-effect-only mode — workers apply the record function to
    /// every variant but no VCF is written and records aren't serialised.
    pub fn with_output_file(mut self, path: impl Into<PathBuf>) -> Self {
        self.output = Some(path.into());
        self
    }

    pub fn worker_threads(mut self, n: usize) -> Self {
        self.worker_threads = Some(n);
        self
    }

    /// Maximum number of windows held in memory at once (queued, being
    /// processed, or awaiting in-order output). Default `2 * worker_threads`;
    /// values below `worker_threads` are raised to it.
    pub fn max_in_flight(mut self, n: usize) -> Self {
        self.max_in_flight = Some(n);
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

    /// Per-worker record function factory. Called once in each worker thread
    /// with that worker's header (which already contains `header_lines`).
    pub fn record_factory<F>(mut self, f: F) -> Self
    where
        F: Fn(&Arc<HeaderView>) -> Result<RecordFn> + Send + Sync + 'static,
    {
        self.record_factory = Some(Arc::new(f));
        self
    }

    /// Convenience for stateless functions: `true` keeps the record.
    pub fn record_callback<F>(self, f: F) -> Self
    where
        F: Fn(&mut Record) -> bool + Send + Sync + 'static,
    {
        let f = Arc::new(f);
        self.record_factory(move |_h| {
            let f = Arc::clone(&f);
            Ok(Box::new(move |r: &mut Record| Ok(f(r))) as RecordFn)
        })
    }

    /// Optional. `##` header lines appended to the output header (and to the
    /// header records are parsed with, so new INFO/FORMAT/FILTER IDs exist).
    pub fn header_lines(mut self, lines: Vec<String>) -> Self {
        self.header_lines = lines;
        self
    }

    /// Optional. Subset samples as each record is read, before the record
    /// function sees it (the header lines are then appended to the subset
    /// header).
    pub fn samples(mut self, samples: Option<SampleSubset>) -> Self {
        self.samples = samples;
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
        let max_in_flight = self
            .max_in_flight
            .unwrap_or(2 * worker_threads)
            .max(worker_threads);
        let window_size = self.window_size.unwrap_or(DEFAULT_WINDOW_SIZE).max(1);
        let stride = self.stride.unwrap_or(1).max(1);
        let record_factory = self
            .record_factory
            .ok_or_else(|| anyhow!("missing record_callback/record_factory"))?;

        ParallelVariantWindowProcessor {
            input,
            output: self.output,
            worker_threads,
            max_in_flight,
            window_size,
            stride,
            contig_filter: self.contig_filter,
            contig_lengths: self.contig_lengths,
            record_factory,
            header_lines: self.header_lines,
            samples: self.samples,
            progress_callback: self.progress_callback,
        }
        .execute()
    }
}

pub struct ParallelVariantWindowProcessor {
    input: PathBuf,
    output: Option<PathBuf>,
    worker_threads: usize,
    max_in_flight: usize,
    window_size: u64,
    stride: usize,
    contig_filter: Option<HashSet<String>>,
    contig_lengths: Option<IndexMap<String, u64>>,
    record_factory: Arc<RecordFactory>,
    header_lines: Vec<String>,
    samples: Option<SampleSubset>,
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
            max_in_flight,
            window_size,
            stride,
            contig_filter,
            contig_lengths,
            record_factory,
            header_lines,
            samples,
            progress_callback,
        } = self;

        // 1. Open the input, build the output header and plan contigs.
        let mut reader = VcfReader::open_indexed(&input)?;
        if let Some(s) = samples.as_ref() {
            reader.subset_samples(s)?;
        }
        reader.append_header_lines(&header_lines)?;
        let header_text = hts::format_header(reader.header())?;
        let contigs = plan_contigs(&reader, window_size, contig_filter.as_ref(), contig_lengths.as_ref());
        drop(reader); // workers and the producer each open their own.
        let total_abs: usize = contigs.last().map_or(0, |c| c.first_abs + c.n_windows(window_size));
        let total_windows = total_abs.div_ceil(stride);

        // 2. Open output and write the header through a multithreaded bgzf
        //    writer. Without an output path, run in side-effect-only mode.
        let mut mt_writer = match output.as_ref() {
            Some(path) => {
                let file = File::create(path)
                    .with_context(|| format!("failed to create output: {}", path.display()))?;
                let mut w = bgzf::io::multithreaded_writer::Builder::default()
                    .set_worker_count(NonZeroUsize::new(worker_threads).unwrap())
                    .build_from_writer(file);
                w.write_all(&header_text).context("failed to write VCF header")?;
                Some(w)
            }
            None => None,
        };
        let serialize = mt_writer.is_some();

        // 3. Channels. Every window in flight holds a permit, so the channels
        //    never need more capacity than `max_in_flight`.
        let (work_tx, work_rx) = bounded::<Window>(max_in_flight);
        type Out = (usize, usize, u64, Result<Vec<u8>>);
        let (out_tx, out_rx) = bounded::<Out>(max_in_flight);
        let (permit_tx, permit_rx) = bounded::<()>(max_in_flight);
        for _ in 0..max_in_flight {
            permit_tx.send(()).expect("permit channel has capacity");
        }

        // 4. Producer + workers + writer.
        std::thread::scope(|scope| -> Result<()> {
            // Move out_rx and permit_tx into the scope so they're dropped when
            // the writer finishes (or errors), unblocking workers' sends and
            // the producer's permit wait.
            let out_rx = out_rx;
            let permit_tx = permit_tx;

            // Producer: generate windows lazily, in order, one permit each.
            // Because dispatch is in order, the next window the writer needs
            // is always already in flight, so this cannot deadlock.
            let producer_input = input.clone();
            scope.spawn(move || {
                // Own reader for index lookups (only used to skip empty spans;
                // if it cannot be opened, every window is simply queried).
                let reader = VcfReader::open_indexed(&producer_input).ok();
                let mut next_idx = 0usize;
                let mut emit = |mut w: Window| -> bool {
                    if permit_rx.recv().is_err() {
                        return false;
                    }
                    w.idx = next_idx;
                    next_idx += 1;
                    work_tx.send(w).is_ok()
                };
                for contig in &contigs {
                    let wgen = WindowGen {
                        contig,
                        window_size,
                        stride,
                        tid: reader.as_ref().map(|r| r.tid(&contig.name)),
                        reader: reader.as_ref(),
                    };
                    if !wgen.run(&mut emit) {
                        break;
                    }
                }
            });

            // Workers.
            for _ in 0..worker_threads {
                let work_rx = work_rx.clone();
                let out_tx = out_tx.clone();
                let input = input.clone();
                let header_lines = header_lines.clone();
                let samples = samples.clone();
                let record_factory = Arc::clone(&record_factory);
                scope.spawn(move || {
                    let setup = || -> Result<(VcfReader, RecordFn)> {
                        let mut reader = VcfReader::open_indexed(&input)?;
                        if let Some(s) = samples.as_ref() {
                            reader.subset_samples(s)?;
                        }
                        reader.append_header_lines(&header_lines)?;
                        let f = record_factory(reader.header())?;
                        Ok((reader, f))
                    };
                    let (mut reader, mut f) = match setup() {
                        Ok(x) => x,
                        Err(e) => {
                            // Report the error for every window we receive so
                            // the producer is not blocked.
                            let msg = format!("{e:#}");
                            while let Ok(w) = work_rx.recv() {
                                let _ = out_tx.send((w.idx, w.n_windows, w.bp, Err(anyhow!("{msg}"))));
                            }
                            return;
                        }
                    };
                    // Two record buffers used alternately: bcftools' filter.c
                    // caches per-record data keyed by the record's address,
                    // which must therefore differ between consecutive records.
                    let mut recs = [reader.empty_record(), reader.empty_record()];
                    let mut ks = KString::new();
                    while let Ok(window) = work_rx.recv() {
                        let result = if window.skip {
                            Ok(Vec::new())
                        } else {
                            process_window(&mut reader, &mut recs, &mut ks, serialize, &window, &mut f)
                        };
                        if out_tx.send((window.idx, window.n_windows, window.bp, result)).is_err() {
                            break;
                        }
                    }
                });
            }
            drop(out_tx);

            // Writer / re-orderer / progress (main scope thread).
            let mut buffer: BTreeMap<usize, (usize, u64, Vec<u8>)> = BTreeMap::new();
            let mut next_idx = 0usize;
            let mut done_windows = 0usize;
            let mut cum_bp: u64 = 0;
            let mut first_error: Option<anyhow::Error> = None;
            while let Ok((idx, n, bp, result)) = out_rx.recv() {
                match result {
                    Ok(bytes) => {
                        buffer.insert(idx, (n, bp, bytes));
                    }
                    Err(e) => {
                        // Fail the run: never silently drop data.
                        first_error = Some(e);
                        break;
                    }
                }
                while let Some((n, bp, bytes)) = buffer.remove(&next_idx) {
                    if let Some(w) = mt_writer.as_mut() {
                        w.write_all(&bytes).context("failed to write window bytes")?;
                    }
                    cum_bp += bp;
                    done_windows += n;
                    // Items beyond the declared contig length (n == 0) carry
                    // no progress.
                    if n > 0
                        && let Some(cb) = progress_callback.as_ref()
                    {
                        cb(cum_bp, done_windows, total_windows);
                    }
                    next_idx += 1;
                    // Window written: release its permit.
                    let _ = permit_tx.send(());
                }
            }

            if let Some(e) = first_error {
                return Err(e);
            }
            if done_windows != total_windows {
                return Err(anyhow!(
                    "writer terminated with gaps: wrote {done_windows} of {total_windows} windows"
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

/// Decide which contigs to process, in file order, with declared lengths.
fn plan_contigs(
    reader: &VcfReader,
    window_size: u64,
    contig_filter: Option<&HashSet<String>>,
    contig_lengths: Option<&IndexMap<String, u64>>,
) -> Vec<ContigPlan> {
    // Candidates: header contigs, then .fai contigs, then (tabix) contigs that
    // only the index knows about — so records on contigs missing from the
    // header are not lost.
    let mut ordered: Vec<(String, u64)> = Vec::new();
    let mut seen: HashSet<String> = HashSet::new();
    let mut push = |name: String, len: Option<u64>, ordered: &mut Vec<(String, u64)>| {
        if seen.insert(name.clone()) {
            let len = len.or_else(|| contig_lengths.and_then(|cl| cl.get(&name).copied()));
            ordered.push((name, len.unwrap_or(0)));
        }
    };
    for (name, len) in hts::header_contigs(reader.header()) {
        push(name, len, &mut ordered);
    }
    if let Some(cl) = contig_lengths {
        for (name, len) in cl {
            push(name.clone(), Some(*len), &mut ordered);
        }
    }
    if let Some(names) = reader.tabix_seqnames() {
        for name in names {
            push(name, None, &mut ordered);
        }
    }

    // Keep requested contigs that the index has data for (others have no
    // records), ordered by the file offset of their first data.
    let mut keyed: Vec<(u64, usize, String, u64)> = Vec::new();
    for (i, (name, len)) in ordered.into_iter().enumerate() {
        if let Some(filter) = contig_filter
            && !filter.contains(name.as_str())
        {
            continue;
        }
        let Some(tid) = reader.tid(&name) else { continue };
        let Some(itr) = reader.query(tid, 0, MAX_POS as i64) else { continue };
        let Some(off) = itr.min_offset() else { continue };
        keyed.push((off, i, name, len));
    }
    keyed.sort_by_key(|k| (k.0, k.1));

    let mut out = Vec::new();
    let mut abs = 0usize;
    for (_, _, name, length) in keyed {
        let plan = ContigPlan { name: Arc::from(name.as_str()), length, first_abs: abs };
        abs += plan.n_windows(window_size);
        out.push(plan);
    }
    out
}

/// Lazily generates the windows of one contig, in order. Windows beyond the
/// declared contig length (a "tail" up to the largest VCF position) are also
/// covered, so records past a wrong or missing header length are never
/// silently lost. Spans the index shows to be empty are emitted as single
/// cheap skip items instead of one item per window, so a contig declared as
/// 2^31 bp with data in the first kb costs a handful of index lookups rather
/// than millions of queries.
struct WindowGen<'a> {
    contig: &'a ContigPlan,
    window_size: u64,
    stride: usize,
    reader: Option<&'a VcfReader>,
    /// `None`: no index reader; `Some(None)`: the index has no such contig.
    tid: Option<Option<i32>>,
}

impl WindowGen<'_> {
    const FANOUT: usize = 8;

    /// Windows within the declared length.
    fn n_declared(&self) -> usize {
        self.contig.n_windows(self.window_size)
    }

    /// All windows to visit: declared ones plus, when every window is
    /// processed (stride 1), the tail up to the largest VCF position.
    fn n_total(&self) -> usize {
        let d = self.n_declared();
        if self.stride == 1 && self.reader.is_some() {
            d.max(MAX_POS.div_ceil(self.window_size) as usize)
        } else {
            d
        }
    }

    /// Number of stride-selected declared windows among windows `[k0, k1)`.
    fn n_selected(&self, k0: usize, k1: usize) -> usize {
        let d = self.n_declared();
        let (k0, k1) = (k0.min(d), k1.min(d));
        let a0 = self.contig.first_abs + k0;
        let a1 = self.contig.first_abs + k1;
        a1.div_ceil(self.stride) - a0.div_ceil(self.stride)
    }

    fn window_bp(&self, k: usize) -> u64 {
        let start = k as u64 * self.window_size + 1;
        let end = ((k as u64 + 1) * self.window_size).min(self.contig.length);
        end - start + 1
    }

    /// Base pairs of the stride-selected declared windows in `[k0, k1)`.
    fn selected_bp(&self, k0: usize, k1: usize) -> u64 {
        let n = self.n_selected(k0, k1) as u64;
        if n == 0 {
            return 0;
        }
        let mut bp = n * self.window_size;
        let last = self.n_declared() - 1;
        if k0 <= last && last < k1 && (self.contig.first_abs + last).is_multiple_of(self.stride) {
            bp -= self.window_size - self.window_bp(last);
        }
        bp
    }

    /// Whether the index proves that no record overlaps `[start, end]`
    /// (1-based, inclusive).
    fn provably_empty(&self, start: u64, end: u64) -> bool {
        let Some(reader) = self.reader else { return false };
        match self.tid {
            Some(None) => true, // contig absent from the index: no records
            Some(Some(tid)) => match reader.query(tid, start as i64 - 1, end as i64) {
                Some(itr) => itr.is_empty(),
                None => false, // be conservative: query the windows
            },
            None => false,
        }
    }

    /// Emit all items for the contig; returns false if the consumer hung up.
    fn run(&self, emit: &mut impl FnMut(Window) -> bool) -> bool {
        self.span(0, self.n_total(), emit)
    }

    fn span(&self, k0: usize, k1: usize, emit: &mut impl FnMut(Window) -> bool) -> bool {
        let n_sel = self.n_selected(k0, k1);
        let n_tail = k1.saturating_sub(k0.max(self.n_declared()));
        if n_sel == 0 && n_tail == 0 {
            return true;
        }
        let start = k0 as u64 * self.window_size + 1;
        // Query the full nominal span even past the declared length (that only
        // affects progress accounting): a wrong header length must not hide
        // records.
        let end = (k1 as u64 * self.window_size).min(MAX_POS);
        let item = |skip: bool| Window {
            idx: 0,
            contig: Arc::clone(&self.contig.name),
            start,
            end,
            n_windows: n_sel,
            bp: self.selected_bp(k0, k1),
            skip,
        };
        if start > MAX_POS || self.provably_empty(start, end) {
            return emit(item(true));
        }
        if k1 - k0 == 1 {
            return emit(item(false));
        }
        let step = (k1 - k0).div_ceil(Self::FANOUT);
        let mut k = k0;
        while k < k1 {
            let next = (k + step).min(k1);
            if !self.span(k, next, emit) {
                return false;
            }
            k = next;
        }
        true
    }
}

fn process_window(
    reader: &mut VcfReader,
    recs: &mut [Record; 2],
    ks: &mut KString,
    serialize: bool,
    window: &Window,
    f: &mut RecordFn,
) -> Result<Vec<u8>> {
    let mut bytes = Vec::new();
    let Some(tid) = reader.tid(&window.contig) else {
        return Ok(bytes);
    };
    let mut itr = reader
        .query(tid, window.start as i64 - 1, window.end as i64)
        .ok_or_else(|| anyhow!("index query failed for {}:{}-{}", window.contig, window.start, window.end))?;
    let mut cur = 0usize;
    loop {
        cur ^= 1;
        let rec = &mut recs[cur];
        if !reader
            .itr_next(&mut itr, rec)
            .with_context(|| format!("in {}:{}-{}", window.contig, window.start, window.end))?
        {
            break;
        }
        // Boundary dedup: index queries return every record overlapping the
        // region, so a record starting before this window was already
        // handled by the window holding its POS.
        if (rec.pos() as u64) + 1 < window.start {
            continue;
        }
        if f(rec).with_context(|| format!("at {}", hts::rec_locus(rec)))? && serialize {
            hts::format_record(rec, ks, &mut bytes)?;
        }
    }
    Ok(bytes)
}
