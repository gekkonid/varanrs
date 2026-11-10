//! Integration tests for `ParallelVariantWindowProcessor`.

use std::fs::File;
use std::io::Write as _;
use std::path::Path;
use std::sync::Mutex;

use anyhow::Result;
use varanrs::processor::ParallelVariantWindowProcessor;

mod common;

const INDEX_TOOLS: &[&str] = &["bgzip", "tabix"];
const BCF_TOOLS: &[&str] = &["bgzip", "tabix", "bcftools"];

/// Create a real `.vcf.gz` (bgzipped + tabix-indexed) at `path`.
///
/// Uses the system `bgzip` and `tabix` binaries. The caller passes the *uncompressed* VCF text;
/// this writes it to `<path>.tmp.vcf`, bgzips to `<path>`, and indexes.
fn write_indexed_vcf(path: &Path, text: &str) {
    let raw = path.with_extension("tmp.vcf");
    File::create(&raw).unwrap().write_all(text.as_bytes()).unwrap();

    let status = std::process::Command::new("bgzip")
        .arg("-c")
        .arg(&raw)
        .stdout(File::create(path).unwrap())
        .status()
        .unwrap();
    assert!(status.success(), "bgzip failed");

    let status = std::process::Command::new("tabix")
        .arg("-p")
        .arg("vcf")
        .arg(path)
        .status()
        .unwrap();
    assert!(status.success(), "tabix failed");
}

fn read_records(path: &Path) -> Vec<(String, usize, String, Vec<String>)> {
    let (_h, records) = varanrs::hts::read_all(path.to_str().unwrap())
        .unwrap_or_else(|e| panic!("{path:?}: read error: {e}"));
    records
        .iter()
        .map(|r| {
            let alleles: Vec<String> = r.alleles().iter().map(|a| String::from_utf8_lossy(a).into_owned()).collect();
            (varanrs::hts::rec_contig(r), r.pos() as usize + 1, alleles[0].clone(), alleles[1..].to_vec())
        })
        .collect()
}

const HEADER: &str = "##fileformat=VCFv4.2\n\
##contig=<ID=chr1,length=5000>\n\
##contig=<ID=chr2,length=3000>\n\
##INFO=<ID=DUMMY,Number=1,Type=Integer,Description=\"Placeholder\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n";

fn build_simple_input(path: &Path) {
    let mut text = String::from(HEADER);
    // 30 records on chr1 spread across all 5 windows of size 1000.
    for i in 0..30 {
        let pos = 50 + i * 150; // 50..4400
        let alt = if i % 3 == 0 { "C" } else { "T" };
        text.push_str(&format!(
            "chr1\t{pos}\t.\tA\t{alt}\t.\t.\t.\tGT\t0/1\n"
        ));
    }
    // 10 records on chr2.
    for i in 0..10 {
        let pos = 100 + i * 250;
        text.push_str(&format!(
            "chr2\t{pos}\t.\tG\tA\t.\t.\t.\tGT\t0/0\n"
        ));
    }
    write_indexed_vcf(path, &text);
}

#[test]
fn identity_round_trip() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(|_| true)
        .run()?;

    let before = read_records(&input);
    let after = read_records(&output);
    assert_eq!(before, after);
    Ok(())
}

#[test]
fn filtering_drops_records() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(|rec| {
            // Drop records whose first ALT is "C".
            rec.alleles().get(1).is_none_or(|a| *a != b"C")
        })
        .run()?;

    let after = read_records(&output);
    assert!(after.iter().all(|(_, _, _, alts)| alts.first().map(|s| s != "C").unwrap_or(true)));
    let before = read_records(&input);
    let dropped = before.iter().filter(|(_, _, _, alts)| alts.first().map(|s| s == "C").unwrap_or(false)).count();
    assert_eq!(after.len(), before.len() - dropped);
    Ok(())
}

#[test]
fn header_callback_is_applied_to_output_header() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(2)
        .window_size(1000)
        .record_callback(|_| true)
        .header_lines(vec![
            "##INFO=<ID=PYGO_TEST,Number=1,Type=Integer,Description=\"added by test\">".to_string(),
        ])
        .run()?;

    let (header, _) = varanrs::hts::read_all(output.to_str().unwrap())?;
    assert!(varanrs::hts::header_has(&header, rust_htslib::htslib::BCF_HL_INFO, "PYGO_TEST"));
    Ok(())
}

#[test]
fn output_is_globally_ordered_under_parallelism() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(8)
        .window_size(137) // small windows, lots of reordering
        .record_callback(|_| true)
        .run()?;

    let recs = read_records(&output);
    // chr1 records first (sorted by POS), then chr2 records (sorted by POS),
    // mirroring header contig order.
    let chr1: Vec<_> = recs.iter().filter(|r| r.0 == "chr1").collect();
    let chr2: Vec<_> = recs.iter().filter(|r| r.0 == "chr2").collect();
    assert!(chr1.windows(2).all(|w| w[0].1 <= w[1].1), "chr1 not sorted");
    assert!(chr2.windows(2).all(|w| w[0].1 <= w[1].1), "chr2 not sorted");
    let chr1_end = recs.iter().rposition(|r| r.0 == "chr1").unwrap();
    let chr2_start = recs.iter().position(|r| r.0 == "chr2").unwrap();
    assert!(chr1_end < chr2_start, "contigs interleaved");
    Ok(())
}

#[test]
fn boundary_record_is_not_duplicated() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");

    // Three records: at POS 1000 (window 1 end), 1001 (window 2 start),
    // and a deletion spanning 999-1002 (REF=AAAA at POS=999) which tabix
    // returns for both window queries.
    let mut text = String::from(HEADER);
    text.push_str("chr1\t999\t.\tAAAA\tA\t.\t.\t.\tGT\t0/1\n");
    text.push_str("chr1\t1000\t.\tA\tT\t.\t.\t.\tGT\t0/1\n");
    text.push_str("chr1\t1001\t.\tA\tT\t.\t.\t.\tGT\t0/1\n");
    write_indexed_vcf(&input, &text);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(4)
        .window_size(1000) // boundary at POS 1000/1001
        .record_callback(|_| true)
        .run()?;

    let recs = read_records(&output);
    assert_eq!(recs.len(), 3, "expected exactly 3 records, got {recs:?}");
    let positions: Vec<usize> = recs.iter().map(|r| r.1).collect();
    assert_eq!(positions, vec![999, 1000, 1001]);
    Ok(())
}

#[test]
fn progress_callback_reports_cumulative_bp() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    let calls: std::sync::Arc<Mutex<Vec<u64>>> = Default::default();
    let calls_for_cb = std::sync::Arc::clone(&calls);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(|_| true)
        .progress_callback(move |bp, _, _| calls_for_cb.lock().unwrap().push(bp))
        .run()?;

    let calls = calls.lock().unwrap().clone();
    // chr1 (length 5000) → 5 windows of 1000 bp
    // chr2 (length 3000) → 3 windows of 1000 bp
    assert_eq!(calls.len(), 8);
    // Monotonic non-decreasing.
    assert!(calls.windows(2).all(|w| w[0] <= w[1]));
    // Final value = total bp covered by all windows.
    assert_eq!(*calls.last().unwrap(), 5000 + 3000);
    Ok(())
}

#[test]
fn iterates_without_writing_when_no_output_set() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    build_simple_input(&input);

    // Side-effect-only callback: count records, never write a VCF.
    let count = Arc::new(AtomicUsize::new(0));
    let count_for_cb = Arc::clone(&count);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(move |_rec| {
            count_for_cb.fetch_add(1, Ordering::Relaxed);
            // Keeping or dropping doesn't matter — nothing is written.
            true
        })
        .run()?;

    // Fixture: 30 chr1 + 10 chr2 = 40 records.
    assert_eq!(count.load(Ordering::Relaxed), 40);
    Ok(())
}

#[test]
fn no_output_still_fires_progress_in_window_order() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    build_simple_input(&input);

    let calls: std::sync::Arc<Mutex<Vec<u64>>> = Default::default();
    let calls_for_cb = std::sync::Arc::clone(&calls);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(|_| true)
        .progress_callback(move |bp, _, _| calls_for_cb.lock().unwrap().push(bp))
        .run()?;

    let calls = calls.lock().unwrap().clone();
    assert_eq!(calls.len(), 8);
    assert!(calls.windows(2).all(|w| w[0] <= w[1]));
    assert_eq!(*calls.last().unwrap(), 5000 + 3000);
    Ok(())
}

fn write_indexed_bcf(path: &Path, vcf_gz_path: &Path) {
    let status = std::process::Command::new("bcftools")
        .args(["convert", "-O", "b", "-o"])
        .arg(path)
        .arg(vcf_gz_path)
        .status()
        .unwrap();
    assert!(status.success(), "bcftools convert failed");

    let status = std::process::Command::new("bcftools")
        .args(["index"])
        .arg(path)
        .status()
        .unwrap();
    assert!(status.success(), "bcftools index failed");
}

#[test]
fn identity_round_trip_with_symbolic_alleles() -> Result<()> {
    if !common::check_tools(BCF_TOOLS) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let input_gz = dir.path().join("in.vcf.gz");
    let output_gz = dir.path().join("out.vcf.gz");
    let input_bcf = dir.path().join("in.bcf");
    let output_bcf = dir.path().join("out.bcf.vcf.gz");

    let text = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##ALT=<ID=DEL,Description=\"Deletion\">
##ALT=<ID=DUP,Description=\"Duplication\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\t.\tA\tG\t.\t.\t.\tGT\t0/1\t1/1
chr1\t200\t.\tC\tT\t.\t.\t.\tGT\t0/0\t0/1
chr1\t300\t.\tT\t<DEL>\t.\t.\t.\tGT\t0/1\t1/1
chr1\t400\t.\tG\tA\t.\t.\t.\tGT\t1/1\t0/0
chr1\t500\t.\tA\t<DUP>\t.\t.\t.\tGT\t0/1\t0/1
chr1\t600\t.\tC\tG\t.\t.\t.\tGT\t0/0\t0/0
chr1\t700\t.\tT\tC,A\t.\t.\t.\tGT\t0/1\t1/2
";

    // VCF.gz round-trip
    write_indexed_vcf(&input_gz, text);
    ParallelVariantWindowProcessor::builder()
        .input(&input_gz)
        .with_output_file(&output_gz)
        .worker_threads(2)
        .window_size(5000)
        .record_callback(|_| true)
        .run()?;

    let before = read_records(&input_gz);
    let after = read_records(&output_gz);
    assert_eq!(before, after, "VCF.gz round-trip with symbolic alleles");
    assert_eq!(before.len(), 7);

    // BCF round-trip
    write_indexed_bcf(&input_bcf, &input_gz);
    ParallelVariantWindowProcessor::builder()
        .input(&input_bcf)
        .with_output_file(&output_bcf)
        .worker_threads(2)
        .window_size(5000)
        .record_callback(|_| true)
        .run()?;

    let before_bcf = read_records(&output_bcf);
    assert_eq!(before_bcf.len(), before.len(),
        "BCF round-trip with symbolic alleles: record count");

    for (i, ((orig_chr, orig_pos, orig_ref, orig_alts), (bcf_chr, bcf_pos, bcf_ref, bcf_alts))) in
        before.iter().zip(before_bcf.iter()).enumerate()
    {
        assert_eq!(orig_chr, bcf_chr, "record {i}: chrom mismatch");
        assert_eq!(orig_pos, bcf_pos, "record {i}: pos mismatch");
        assert_eq!(orig_ref, bcf_ref, "record {i}: ref mismatch");
        assert_eq!(orig_alts, bcf_alts, "record {i}: alt mismatch");
    }

    Ok(())
}

#[test]
fn glnexus_rnc_normalize_header_fix() -> Result<()> {
    if !common::check_tools(BCF_TOOLS) {
        return Ok(());
    }

    // GLnexus declares RNC as Type=Character,Number=2 but encodes "II" as a
    // single string. htslib (like bcftools) reads and writes this as-is, so
    // no header workaround is needed: records and RNC values round-trip.
    let dir = tempfile::tempdir()?;
    let input_gz = dir.path().join("in.vcf.gz");
    let output_gz = dir.path().join("out.vcf.gz");

    let text = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=RNC,Number=2,Type=Character,Description=\"Reason for No Call\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\t.\tA\tG\t.\t.\t.\tGT:RNC\t0/1:II\t1/1:XX
chr1\t200\t.\tC\tT\t.\t.\t.\tGT:RNC\t0/0:..\t1/1:YY
";

    write_indexed_vcf(&input_gz, text);

    assert_eq!(read_records(&input_gz).len(), 2, "htslib reads the raw input");

    ParallelVariantWindowProcessor::builder()
        .input(&input_gz)
        .with_output_file(&output_gz)
        .worker_threads(2)
        .window_size(5000)
        .record_callback(|_| true)
        .run()?;

    let out = read_records(&output_gz);
    assert_eq!(out.len(), 2, "output should have 2 records");
    let text = String::from_utf8(
        std::process::Command::new("bgzip").arg("-dc").arg(&output_gz).output()?.stdout,
    )?;
    assert!(text.contains("GT:RNC\t0/1:II\t1/1:XX"), "RNC values preserved: {text}");

    // BCF round-trip
    let input_bcf = dir.path().join("in.bcf");
    let output_bcf = dir.path().join("out.bcf.vcf.gz");
    write_indexed_bcf(&input_bcf, &input_gz);
    ParallelVariantWindowProcessor::builder()
        .input(&input_bcf)
        .with_output_file(&output_bcf)
        .worker_threads(2)
        .window_size(5000)
        .record_callback(|_| true)
        .run()?;

    let out_bcf = read_records(&output_bcf);
    assert_eq!(out_bcf.len(), 2, "BCF output should have 2 records after fix");
    assert_eq!(out, out_bcf,
        "VCF.gz and BCF round-trip should produce identical records after header fix"
    );

    Ok(())
}

/// At most `max_in_flight` windows may be dispatched before the first
/// (deliberately slow) window is written, however many windows remain.
fn windows_started_before_first_write(max_in_flight: Option<usize>) -> Result<usize> {
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    let mut text = String::from(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=4000>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n",
    );
    for w in 0..40 {
        text.push_str(&format!("chr1\t{}\t.\tA\tC\t.\t.\t.\n", w * 100 + 50));
    }
    write_indexed_vcf(&input, &text);

    let started = Arc::new(Mutex::new(Vec::<usize>::new()));
    let snapshot = Arc::new(AtomicUsize::new(usize::MAX));
    let (s1, s2, snap) = (started.clone(), started.clone(), snapshot.clone());

    let mut b = ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(2)
        .window_size(100)
        .record_callback(move |rec| {
            let w = rec.pos() as usize / 100;
            s1.lock().unwrap().push(w);
            if w == 0 {
                std::thread::sleep(std::time::Duration::from_millis(300));
            }
            true
        })
        .progress_callback(move |_, done, _| {
            if done == 1 {
                snap.store(s2.lock().unwrap().len(), Ordering::SeqCst);
            }
        });
    if let Some(n) = max_in_flight {
        b = b.max_in_flight(n);
    }
    b.run()?;
    assert_eq!(started.lock().unwrap().len(), 40);
    Ok(snapshot.load(Ordering::SeqCst))
}

#[test]
fn in_flight_windows_bounded() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }
    // default: 2 x 2 threads
    let n = windows_started_before_first_write(None)?;
    assert!(n <= 4, "{n} windows started before window 0 was written (limit 4)");
    let n = windows_started_before_first_write(Some(2))?;
    assert!(n <= 2, "{n} windows started before window 0 was written (limit 2)");
    Ok(())
}

/// Contigs declared with huge lengths (here 2^31-1, as some tools emit) must
/// not make the processor enumerate or query millions of empty windows:
/// empty spans are skipped via the index, and memory stays flat.
#[test]
fn huge_declared_contig_lengths_are_cheap() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }
    use std::sync::Arc;
    use std::sync::atomic::{AtomicUsize, Ordering};

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    let text = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=2147483647>\n\
##contig=<ID=2,length=2147483647>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t105\t.\tA\tC\t.\t.\t.\n\
1\t250\t.\tA\tG\t.\t.\t.\n\
1\t300000000\t.\tC\tT\t.\t.\t.\n\
2\t1\t.\tG\tA\t.\t.\t.\n";
    write_indexed_vcf(&input, text);

    let last = Arc::new((AtomicUsize::new(0), AtomicUsize::new(0)));
    let calls = Arc::new(AtomicUsize::new(0));
    let (l2, c2) = (last.clone(), calls.clone());
    let t0 = std::time::Instant::now();
    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(2)
        .window_size(100)
        .record_callback(|_| true)
        .progress_callback(move |_, done, total| {
            l2.0.store(done, Ordering::SeqCst);
            l2.1.store(total, Ordering::SeqCst);
            c2.fetch_add(1, Ordering::SeqCst);
        })
        .run()?;
    let elapsed = t0.elapsed();

    let recs = read_records(&output);
    let pos: Vec<(String, usize)> = recs.iter().map(|r| (r.0.clone(), r.1)).collect();
    assert_eq!(
        pos,
        vec![("1".into(), 105), ("1".into(), 250), ("1".into(), 300000000), ("2".into(), 1)]
    );
    // every nominal window is accounted for in progress
    let total = 2 * 2147483647usize.div_ceil(100);
    assert_eq!(last.1.load(Ordering::SeqCst), total);
    assert_eq!(last.0.load(Ordering::SeqCst), total);
    // ...but via a small number of work items, quickly
    assert!(calls.load(Ordering::SeqCst) < 2000, "{} items", calls.load(Ordering::SeqCst));
    assert!(elapsed.as_secs() < 20, "took {elapsed:?}");
    Ok(())
}

/// Records beyond a contig's declared length, or on contigs missing from the
/// header, must not be silently dropped (htslib/bcftools keep them).
#[test]
fn records_beyond_declared_length_are_kept() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    let text = "##fileformat=VCFv4.2\n\
##contig=<ID=1,length=1>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
1\t22\t.\tA\tG\t.\t.\t.\n\
1\t70000\t.\tC\tT\t.\t.\t.\n\
2\t5\t.\tG\tA\t.\t.\t.\n";
    write_indexed_vcf(&input, text);
    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .with_output_file(&output)
        .worker_threads(2)
        .window_size(1000)
        .record_callback(|_| true)
        .run()?;
    let pos: Vec<(String, usize)> = read_records(&output).into_iter().map(|r| (r.0, r.1)).collect();
    assert_eq!(pos, vec![("1".into(), 22), ("1".into(), 70000), ("2".into(), 5)]);
    Ok(())
}

fn run_identity(input: &Path, output: &Path) -> Result<()> {
    ParallelVariantWindowProcessor::builder()
        .input(input)
        .with_output_file(output)
        .worker_threads(2)
        .window_size(100)
        .record_callback(|_| true)
        .run()
}

fn body_lines(path: &Path) -> Vec<String> {
    let out = std::process::Command::new("bgzip").arg("-dc").arg(path).output().unwrap();
    String::from_utf8(out.stdout).unwrap().lines().filter(|l| !l.starts_with('#')).map(String::from).collect()
}

/// Records that are unusual but accepted by htslib/bcftools (IUPAC REF,
/// REF '*', INFO keys starting with a digit, undeclared tags) round-trip
/// byte-for-byte.
#[test]
fn unusual_records_round_trip_unchanged() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    let records = [
        "chr1\t10\t.\tB\tBGC\t.\t.\t.",
        "chr1\t20\t.\tR\tA\t.\t.\t1000g2012apr_eur=0.5",
        "chr1\t30\t.\t*\tA\t.\t.\tUNDECLARED=x",
    ];
    let text = format!(
        "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##INFO=<ID=1000g2012apr_eur,Number=1,Type=Float,Description=\"x\">\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n{}\n",
        records.join("\n")
    );
    write_indexed_vcf(&input, &text);
    run_identity(&input, &output)?;
    assert_eq!(body_lines(&output), records);
    Ok(())
}

/// Output follows the file's contig order, not the header's.
#[test]
fn output_follows_file_contig_order() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    let text = "##fileformat=VCFv4.2\n##contig=<ID=chr2,length=1000>\n##contig=<ID=chr1,length=1000>\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n\
chr1\t5\t.\tA\tC\t.\t.\t.\nchr1\t500\t.\tA\tC\t.\t.\t.\nchr2\t7\t.\tG\tT\t.\t.\t.\n";
    write_indexed_vcf(&input, text);
    run_identity(&input, &output)?;
    let order: Vec<(String, usize)> = read_records(&output).into_iter().map(|r| (r.0, r.1)).collect();
    assert_eq!(order, vec![("chr1".into(), 5), ("chr1".into(), 500), ("chr2".into(), 7)]);
    Ok(())
}

/// A record htslib cannot parse fails the run instead of being dropped.
#[test]
fn unparseable_record_fails_the_run() -> Result<()> {
    if !common::check_tools(INDEX_TOOLS) {
        return Ok(());
    }
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    // second record has an invalid genotype (bcftools: "Couldn't read GT data")
    let text = "##fileformat=VCFv4.2\n##contig=<ID=chr1,length=1000>\n##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
chr1\t5\t.\tA\tC\t.\t.\t.\tGT\t0/1\nchr1\t50\t.\tA\tC\t.\t.\t.\tGT\t0/x\n";
    write_indexed_vcf(&input, text);
    assert!(run_identity(&input, &output).is_err());
    Ok(())
}
