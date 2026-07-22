//! End-to-end validation tests for the VCF processor.
//!
//! These tests compare the processor's output against `bcftools view` for
//! every `.vcf.gz` file found in `tests/testdata/`.  The comparison operates
//! at the VCF text level: `bcftools view input.vcf.gz` vs
//! `bgzip -d -c our_output.vcf.gz`, compared line-by-line.
//!
//! For large inputs the comparison streams lines rather than loading
//! everything into RAM at once.
//!
//! Tests automatically discover all `*.vcf.gz` files in the testdata
//! directory, so you can drop in new real-world VCFs (from different callers,
//! sample sizes, etc.) without modifying this code.

use std::fs::File;
use std::io::{BufRead, BufReader, Read, Write};
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};

use noodles_util::variant;
use noodles_vcf::variant::RecordBuf;
use varanrs::processor::ParallelVariantWindowProcessor;

/// Path to the testdata directory (relative to crate root).
const TESTDATA_DIR: &str = "tests/testdata";

/// (sorted header lines, data line iterator)
type ParsedVcf = (Vec<String>, Box<dyn Iterator<Item = String>>);

/// Parse a VCF file on disk into a pair: (sorted header lines, data line iterator).
///
/// The header (all `##` and `#CHROM` lines) is read into a Vec and sorted:
/// - The first line is always the `##fileformat=...` line (preserved at index 0)
/// - The `#CHROM` header line is preserved at the position it appears in the unsorted list
/// - All other `##` lines are sorted lexicographically
///
/// The data lines (non-header, non-empty) are returned as a `Box<dyn Iterator<Item = String>>`
/// so they can be streamed without loading into RAM.
///
/// This function reads only the header into RAM; data lines remain streamed.
fn parse_vcf_file(path: &Path) -> Result<ParsedVcf, String> {
    let mut headers: Vec<String> = Vec::new();

    let mut line = String::new();
    let file = File::open(path)
        .map_err(|e| format!("failed to open {}: {}", path.display(), e))?;
    let mut reader = BufReader::new(file); // file is consumed here

    loop {
        line.clear();
        let n = reader.read_line(&mut line).map_err(|e| e.to_string())?;
        if n == 0 {
            break;
        }

        // Strip trailing newline for storage
        let line = line.trim_end_matches('\n').trim_end_matches('\r').to_string();

        if line.starts_with("#CHROM") {
            headers.push(line);
            break;
        } else if line.starts_with('#') {
            headers.push(line);
        } else {
            break;
            
        }
    }

    // Sort headers: fileformat first, then ## lines alphabetically, then #CHROM last
    headers = sort_vcf_headers(headers);

    Ok((
        headers,
        Box::new(
            reader
                .lines()
                .map_while(Result::ok)
                .filter(|l| !l.is_empty()),
        ),
    ))
}

/// Sort VCF header lines in place.
/// - First line: `##fileformat=...` (always first)
/// - Last line: `#CHROM` (always last)
/// - All other `##` lines: sorted lexicographically
fn sort_vcf_headers(headers: Vec<String>) -> Vec<String> {
    let mut fileformat = None;
    let mut chrom_line = None;
    let mut other_headers: Vec<String> = Vec::new();

    for (i, h) in headers.into_iter().enumerate() {
        if i == 0 && h.starts_with("##fileformat=") {
            fileformat = Some(h);
        } else if h.starts_with("#CHROM") {
            chrom_line = Some(h);
        } else {
            other_headers.push(h);
        }
    }

    other_headers.sort();

    if let Some(ff) = fileformat {
        other_headers.insert(0, ff);
    }
    if let Some(ch) = chrom_line {
        other_headers.push(ch);
    }
    other_headers
}

/// Discover all `.vcf.gz` files in `tests/testdata/`.
fn discover_test_files() -> Vec<PathBuf> {
    let base = PathBuf::from(TESTDATA_DIR);
    if !base.exists() {
        return vec![];
    }
    let mut files: Vec<PathBuf> = Vec::new();
    if let Ok(entries) = std::fs::read_dir(&base) {
        for entry in entries.flatten() {
            let p = entry.path();
            if p.extension().is_some_and(|e| e == "gz")
                && let Some(stem) = p.file_stem()
                    && let Some(stem_str) = stem.to_str()
                        && stem_str.ends_with(".vcf") {
                            files.push(p);
                        }
        }
    }
    files.sort();
    files
}

/// Check whether `bcftools` is available on PATH.
fn bcftools_available() -> bool {
    Command::new("bcftools")
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Run `bcftools view <input>` and stream the uncompressed VCF lines to a
/// temporary file on disk, returning its path.
///
/// This avoids loading the entire output into memory.
fn run_bcftools_view_to_file(input: &Path, tmpdir: &Path) -> PathBuf {
    let out_path = tmpdir.join("bcftools_output.vcf");

    let mut child = Command::new("bcftools")
        .args(["view", "-o", "-"]) // stdout
        .arg(input)
        .stdout(Stdio::piped())
        .spawn()
        .expect("bcftools command failed to start");

    let stdout = child.stdout.take().expect("bcftools has no stdout");
    let mut reader = BufReader::new(stdout);
    let mut file = File::create(&out_path).expect("failed to create bcftools output file");
    let mut line = String::new();
    loop {
        line.clear();
        let n = reader.read_line(&mut line).expect("bcftools stdout read");
        if n == 0 {
            break;
        }
        file.write_all(line.as_bytes()).expect("bcftools output write");
    }
    let status = child.wait().expect("bcftools wait");
    assert!(status.success(), "bcftools view failed on {}", input.display());

    out_path
}

/// Run `bgzip -d -c <path>` and stream the uncompressed VCF lines to a
/// temporary file on disk, returning its path.
fn bgzip_decompress_to_file(path: &Path, tmpdir: &Path) -> PathBuf {
    let suffix = ".decompressed.vcf";
    let out_path = tmpdir.join(format!(
        "{}_{}",
        path.file_stem().unwrap().to_str().unwrap(),
        suffix
    ));

    let mut child = Command::new("bgzip")
        .args(["-d", "-c"])
        .arg(path)
        .stdout(Stdio::piped())
        .spawn()
        .expect("bgzip command failed to start");

    let mut stdout = child.stdout.take().expect("bgzip has no stdout");
    let mut file = File::create(&out_path).expect("failed to create bgzip output file");
    let mut buf = [0u8; 65536];
    loop {
        let n = stdout.read(&mut buf).expect("bgzip stdout read");
        if n == 0 {
            break;
        }
        file.write_all(&buf[..n]).expect("bgzip output write");
    }
    let status = child.wait().expect("bgzip wait");
    assert!(status.success(), "bgzip -d -c failed on {}", path.display());

    out_path
}

/// Run the varanrs processor with identity callback on `input`, writing output to `output`.
fn run_varanrs_identity(input: &Path, output: &Path, worker_threads: usize) {
    ParallelVariantWindowProcessor::builder()
        .input(input)
        .with_output_file(output)
        .worker_threads(worker_threads)
        .window_size(varanrs::processor::DEFAULT_WINDOW_SIZE)
        .record_callback(Some)
        .run()
        .unwrap_or_else(|e| {
            panic!("varanrs processor failed on {}: {}", input.display(), e)
        });
}

/// Run the varanrs processor with identity callback and custom window size.
fn run_varanrs_identity_with_window(
    input: &Path,
    output: &Path,
    worker_threads: usize,
    window_size: u64,
) {
    ParallelVariantWindowProcessor::builder()
        .input(input)
        .with_output_file(output)
        .worker_threads(worker_threads)
        .window_size(window_size)
        .record_callback(Some)
        .run()
        .unwrap_or_else(|e| {
            panic!(
                "varanrs processor failed on {} (window={}): {}",
                input.display(),
                window_size,
                e
            )
        });
}

/// Compare two VCF text files line-by-line, streaming rather than loading
/// everything at once.  Returns Ok(()) if identical, or an error describing
/// the first mismatch.
///
/// This is the core validation: it mirrors `cmp <(bcftools view $INPUT) \
/// <(zcat our_output.vcf.gz)` but within Rust.
fn compare_vcf_text_line_by_line(
    bcftools_path: &Path,
    our_path: &Path,
) -> Result<(), String> {
    let (_bcftools_headers, mut bcftools_iter) =
        parse_vcf_file(bcftools_path).map_err(|e| format!("parse bcftools file: {}", e))?;
    let (_our_headers, mut our_iter) =
        parse_vcf_file(our_path).map_err(|e| format!("parse varanrs file: {}", e))?;

    let mut i = 0usize;
    loop {
        match (bcftools_iter.next(), our_iter.next()) {
            (Some(exp), Some(act)) => {
                i += 1;
                if exp != act {
                    return Err(format!(
                        "Line {} differs:\n  bcftools: {}\n  varanrs:    {}",
                        i, exp, act
                    ));
                }
            }
            (None, None) => break, // both exhausted — identical
            (Some(_), None) => {
                return Err(format!(
                    "varanrs file is shorter: bcftools has more lines after line {}",
                    i
                ));
            }
            (None, Some(_)) => {
                return Err(format!(
                    "bcftools file is shorter: varanrs has more lines after line {}",
                    i
                ));
            }
        }
    }

    Ok(())
}

/// Compare VCF headers for structural equivalence.
/// Header lines are already sorted by `parse_vcf_file`, so we can do a
/// direct line-by-line comparison.  bcftools-injected lines (##FILTER= on
/// PASS-only files, ##bcftools_*, ##source=) are ignored.
fn compare_vcf_headers(
    bcftools_path: &Path,
    our_path: &Path,
) -> Result<(), String> {
    let (bcftools_headers, _) =
        parse_vcf_file(bcftools_path).map_err(|e| format!("parse bcftools file: {}", e))?;
    let (our_headers, _) =
        parse_vcf_file(our_path).map_err(|e| format!("parse varanrs file: {}", e))?;

    let is_bcf_artifact = |h: &&String| -> bool {
        h.starts_with("##bcftools_")
            || h.starts_with("##source=")
            || (h.starts_with("##FILTER=") && h.contains("PASS"))
    };

    let filtered_bcftools: Vec<&String> = bcftools_headers
        .iter()
        .filter(|h| !is_bcf_artifact(h))
        .collect();
    let filtered_ours: Vec<&String> = our_headers
        .iter()
        .filter(|h| !is_bcf_artifact(h))
        .collect();

    if filtered_bcftools != filtered_ours {
        let max_len = filtered_bcftools.len().max(filtered_ours.len());
        let mut diffs = Vec::new();
        for idx in 0..max_len {
            let b = filtered_bcftools.get(idx).map(|s| s.as_str()).unwrap_or("N/A");
            let o = filtered_ours.get(idx).map(|s| s.as_str()).unwrap_or("N/A");
            if b != o {
                diffs.push(format!("  line {}: bcftools: {}\n              varanrs:    {}", idx + 1, b, o));
            }
        }
        return Err(format!(
            "Headers differ at {} positions:\n{}",
            diffs.len(),
            diffs.join("\n")
        ));
    }

    Ok(())
}

// ============================================================================
// Tests
// ============================================================================

/// Primary end-to-end validation: run the processor with identity callback and
/// compare the uncompressed VCF text against `bcftools view`.
///
/// This is the main correctness test — if this passes, every record is
/// preserved with identical content and ordering.
#[test]
fn e2e_identity_matches_bcftools_view() {
    if !bcftools_available() {
        eprintln!("SKIP: bcftools not available");
        return;
    }

    let test_files = discover_test_files();
    if test_files.is_empty() {
        panic!("No .vcf.gz files found in {}", TESTDATA_DIR);
    }

    let tmpdir = tempfile::tempdir().unwrap();

    for test_file in &test_files {
        let basename = test_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        eprintln!("  e2e identity test: {} ...", basename);

        let varanrs_output = tmpdir.path().join(format!("{}_varanrs.vcf.gz", basename));

        // Run varanrs with identity callback (passthrough)
        run_varanrs_identity(test_file, &varanrs_output, 4);

        // Run bcftools view to a temp file on disk (streaming, no full RAM load)
        let bcftools_out = run_bcftools_view_to_file(test_file, tmpdir.path());

        // Decompress our output to a temp file on disk
        let varanrs_decompressed = bgzip_decompress_to_file(&varanrs_output, tmpdir.path());

        // Compare line-by-line (streaming, no full RAM load)
        compare_vcf_text_line_by_line(&bcftools_out, &varanrs_decompressed)
            .unwrap_or_else(|e| panic!("{}: {}", basename, e));
    }
}

/// Test that output records are globally sorted: contigs in header order,
/// then by position within each contig.
#[test]
fn e2e_output_is_globally_ordered() {
    let test_files = discover_test_files();
    if test_files.is_empty() {
        panic!("No .vcf.gz files found in {}", TESTDATA_DIR);
    }

    let tmpdir = tempfile::tempdir().unwrap();

    for test_file in &test_files {
        let basename = test_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        eprintln!("  e2e ordering test: {} ...", basename);

        let output = tmpdir.path().join(format!("{}_ordered.vcf.gz", basename));

        // Run with small window size to maximize reordering stress
        run_varanrs_identity_with_window(test_file, &output, 4, 1000);

        // Get contig order from input header
        let mut input_reader = if test_file.extension().is_some_and(|e| e == "gz") {
            let tmp_vcf = tmpdir.path().join("read_input.vcf");
            let status = Command::new("bgzip")
                .args(["-d", "-c"])
                .arg(test_file)
                .stdout(File::create(&tmp_vcf).unwrap())
                .status()
                .unwrap();
            assert!(status.success());
            variant::io::reader::Builder::default()
                .build_from_path(&tmp_vcf)
                .unwrap()
        } else {
            variant::io::reader::Builder::default()
                .build_from_path(test_file)
                .unwrap()
        };
        let header = input_reader.read_header().unwrap();
        let contig_order: Vec<String> = header
            .contigs()
            .iter()
            .map(|(name, _)| name.to_string())
            .collect();

        // Decompress output for reading
        let output_decompressed = bgzip_decompress_to_file(&output, tmpdir.path());
        let mut output_reader = variant::io::reader::Builder::default()
            .build_from_path(&output_decompressed)
            .unwrap();
        let out_header = output_reader.read_header().unwrap();

        let mut prev_contig_idx: Option<usize> = None;
        let mut prev_pos: usize = 0;
        let mut record_count = 0;

        for result in output_reader.records(&out_header) {
            let record = result.unwrap();
            let buf: RecordBuf =
                RecordBuf::try_from_variant_record(&out_header, record.as_ref()).unwrap();

            let chrom = buf.reference_sequence_name().to_string();
            let pos = buf.variant_start().map(|p| p.get()).unwrap_or(0);

            let contig_idx = contig_order
                .iter()
                .position(|c| c == &chrom)
                .unwrap_or(usize::MAX);

            if let Some(prev_ci) = prev_contig_idx {
                if contig_idx < prev_ci {
                    panic!(
                        "{}: record {} (contig '{}' idx {}) appeared after '{}' idx {} — not in header order",
                        basename, record_count, chrom, contig_idx,
                        contig_order.get(prev_ci).map(|s| s.as_str()).unwrap_or(""),
                        prev_ci
                    );
                }
                if contig_idx == prev_ci && pos < prev_pos {
                    panic!(
                        "{}: record at POS {} appeared after POS {} on same contig '{}'",
                        basename, pos, prev_pos, chrom
                    );
                }
            }

            prev_contig_idx = Some(contig_idx);
            prev_pos = pos;
            record_count += 1;
        }

        assert!(record_count > 0, "{}: expected at least 1 record", basename);
    }
}

/// Test that running with different thread counts produces byte-identical output.
#[test]
fn e2e_parallel_thread_count_equivalence() {
    let test_files = discover_test_files();
    if test_files.is_empty() {
        panic!("No .vcf.gz files found in {}", TESTDATA_DIR);
    }

    let thread_counts = [1, 2, 4, 8];
    let tmpdir = tempfile::tempdir().unwrap();

    for test_file in &test_files {
        let basename = test_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        eprintln!("  e2e thread equivalence test: {} ...", basename);

        // Run with 1 thread (baseline)
        let baseline = tmpdir.path().join(format!("{}_t1.vcf.gz", basename));
        run_varanrs_identity(test_file, &baseline, 1);
        let baseline_bytes = std::fs::read(&baseline).unwrap();

        // Run with each other thread count and compare byte-for-byte
        for &threads in &thread_counts {
            if threads == 1 {
                continue;
            }

            let output = tmpdir.path().join(format!("{}_t{}.vcf.gz", basename, threads));
            run_varanrs_identity(test_file, &output, threads);
            let output_bytes = std::fs::read(&output).unwrap();

            assert_eq!(
                baseline_bytes, output_bytes,
                "{}: output with {} threads differs from 1-thread baseline",
                basename, threads
            );
        }
    }
}

/// Test boundary handling: records at window boundaries should not be
/// duplicated across windows.  Uses a small window size to maximize boundary
/// crossings.
#[test]
fn e2e_no_duplicate_records_at_boundaries() {
    let test_files = discover_test_files();
    if test_files.is_empty() {
        panic!("No .vcf.gz files found in {}", TESTDATA_DIR);
    }

    let tmpdir = tempfile::tempdir().unwrap();

    for test_file in &test_files {
        let basename = test_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        eprintln!("  e2e boundary dedup test: {} ...", basename);

        // Parse input to count data lines (streaming)
        let input_decompressed = bgzip_decompress_to_file(test_file, tmpdir.path());
        let (_, input_iter) =
            parse_vcf_file(&input_decompressed).unwrap();
        let input_count = input_iter.count();

        // Run with very small window to maximize boundary crossings
        let output = tmpdir.path().join(format!("{}_boundary.vcf.gz", basename));
        run_varanrs_identity_with_window(test_file, &output, 4, 100);

        // Count records in output (streaming)
        let output_decompressed = bgzip_decompress_to_file(&output, tmpdir.path());
        let (_, output_iter) =
            parse_vcf_file(&output_decompressed).unwrap();
        let output_count = output_iter.count();

        assert_eq!(
            input_count, output_count,
            "{}: boundary dedup failed — input {} records, output {} records",
            basename, input_count, output_count
        );
    }
}

/// Test that the record count is preserved (no records lost or duplicated).
/// Uses default window size — a simpler sanity check.
#[test]
fn e2e_record_count_preserved() {
    let test_files = discover_test_files();
    if test_files.is_empty() {
        panic!("No .vcf.gz files found in {}", TESTDATA_DIR);
    }

    let tmpdir = tempfile::tempdir().unwrap();

    for test_file in &test_files {
        let basename = test_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        eprintln!("  e2e record count test: {} ...", basename);

        // Count input records (streaming)
        let input_decompressed = bgzip_decompress_to_file(test_file, tmpdir.path());
        let (_, input_iter) =
            parse_vcf_file(&input_decompressed).unwrap();
        let expected_count = input_iter.count();

        // Run varanrs
        let output = tmpdir.path().join(format!("{}_count.vcf.gz", basename));
        run_varanrs_identity(test_file, &output, 4);

        // Count output records (streaming)
        let output_decompressed = bgzip_decompress_to_file(&output, tmpdir.path());
        let (_, output_iter) =
            parse_vcf_file(&output_decompressed).unwrap();
        let actual_count = output_iter.count();

        assert_eq!(
            expected_count, actual_count,
            "{}: expected {} records, got {}",
            basename, expected_count, actual_count
        );
    }
}

/// Test that the header metadata is preserved through the processor.
/// Compares header lines (fileformat, contigs) between bcftools output and
/// varanrs output.
#[test]
fn e2e_header_metadata_preserved() {
    if !bcftools_available() {
        eprintln!("SKIP: bcftools not available");
        return;
    }

    let test_files = discover_test_files();
    if test_files.is_empty() {
        panic!("No .vcf.gz files found in {}", TESTDATA_DIR);
    }

    let tmpdir = tempfile::tempdir().unwrap();

    for test_file in &test_files {
        let basename = test_file
            .file_stem()
            .and_then(|s| s.to_str())
            .unwrap_or("unknown");

        eprintln!("  e2e header test: {} ...", basename);

        let varanrs_output = tmpdir.path().join(format!("{}_header.vcf.gz", basename));
        run_varanrs_identity(test_file, &varanrs_output, 4);

        // Get bcftools header for comparison
        let bcftools_out = run_bcftools_view_to_file(test_file, tmpdir.path());
        let varanrs_decompressed = bgzip_decompress_to_file(&varanrs_output, tmpdir.path());

        // Compare headers (sorted, so order is normalized)
        compare_vcf_headers(&bcftools_out, &varanrs_decompressed)
            .unwrap_or_else(|e| panic!("{}: {}", basename, e));
    }
}