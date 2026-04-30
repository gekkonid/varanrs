//! Integration tests for `ParallelVariantWindowProcessor`.

use std::fs::File;
use std::io::Write as _;
use std::path::Path;
use std::sync::Mutex;

use anyhow::Result;
use noodles_util::variant;
use noodles_vcf as vcf;
use pygopus::processor::ParallelVariantWindowProcessor;

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
    let mut reader = variant::io::reader::Builder::default()
        .build_from_path(path)
        .unwrap();
    let header = reader.read_header().unwrap();
    let mut out = Vec::new();
    for r in reader.records(&header) {
        let r = r.unwrap();
        let buf = vcf::variant::RecordBuf::try_from_variant_record(&header, r.as_ref()).unwrap();
        let chrom = buf.reference_sequence_name().to_string();
        let pos = buf.variant_start().unwrap().get();
        let reference = buf.reference_bases().to_string();
        let alts: Vec<String> = buf.alternate_bases().as_ref().to_vec();
        out.push((chrom, pos, reference, alts));
    }
    out
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
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .output(&output)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(Some)
        .run()?;

    let before = read_records(&input);
    let after = read_records(&output);
    assert_eq!(before, after);
    Ok(())
}

#[test]
fn filtering_drops_records() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .output(&output)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(|buf| {
            // Drop records whose first ALT is "C".
            let drop = buf
                .alternate_bases()
                .as_ref()
                .first()
                .map(|s| s == "C")
                .unwrap_or(false);
            if drop { None } else { Some(buf) }
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
    use vcf::header::record::value::{Map, map::Info};
    use vcf::header::record::value::map::info::{Number, Type};

    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .output(&output)
        .worker_threads(2)
        .window_size(1000)
        .record_callback(Some)
        .header_callback(|h| {
            h.infos_mut().insert(
                "PYGO_TEST".into(),
                Map::<Info>::new(Number::Count(1), Type::Integer, "added by test"),
            );
        })
        .run()?;

    let mut reader = variant::io::reader::Builder::default().build_from_path(&output)?;
    let header = reader.read_header()?;
    assert!(header.infos().contains_key("PYGO_TEST"));
    Ok(())
}

#[test]
fn output_is_globally_ordered_under_parallelism() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .output(&output)
        .worker_threads(8)
        .window_size(137) // small windows, lots of reordering
        .record_callback(Some)
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
        .output(&output)
        .worker_threads(4)
        .window_size(1000) // boundary at POS 1000/1001
        .record_callback(Some)
        .run()?;

    let recs = read_records(&output);
    assert_eq!(recs.len(), 3, "expected exactly 3 records, got {recs:?}");
    let positions: Vec<usize> = recs.iter().map(|r| r.1).collect();
    assert_eq!(positions, vec![999, 1000, 1001]);
    Ok(())
}

#[test]
fn progress_callback_reports_cumulative_bp() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");
    build_simple_input(&input);

    let calls: std::sync::Arc<Mutex<Vec<usize>>> = Default::default();
    let calls_for_cb = std::sync::Arc::clone(&calls);

    ParallelVariantWindowProcessor::builder()
        .input(&input)
        .output(&output)
        .worker_threads(4)
        .window_size(1000)
        .record_callback(Some)
        .progress_callback(move |bp| calls_for_cb.lock().unwrap().push(bp))
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
