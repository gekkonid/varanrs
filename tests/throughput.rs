//! Throughput regression tests for varanrs commands.
//!
//! All tests are `#[ignore]` — run explicitly with:
//!
//! ```bash
//! cargo test --release -- --ignored --test-threads=1
//! ```
//!
//! Single-threaded execution avoids I/O contention on the shared generated
//! VCF.  Thresholds are calibrated for `--release`; debug builds may need
//! `--skip throughput`.

use std::fmt::Write as _;
use std::fs::File;
use std::io::Write as _;
use std::path::{Path, PathBuf};
use std::process::{Command, Stdio};
use std::sync::LazyLock;
use std::time::{Duration, Instant};

use anyhow::{Context, Result};

mod common;

// ---------------------------------------------------------------------------
// Test parameters
// ---------------------------------------------------------------------------

const N_VARIANTS: usize = 67_000;
const N_SAMPLES: usize = 100;

// ---------------------------------------------------------------------------
// VCF generator (pure Rust, deterministic, no rand dependency)
// ---------------------------------------------------------------------------

fn generate_large_vcf(n_variants: usize, n_samples: usize) -> String {
    let n_chr1 = n_variants * 3 / 5;
    let n_chr2 = n_variants - n_chr1;

    let sample_strings_biallelic: [&str; 5] = [
        "0/0:20,0:20:99",
        "0/1:10,10:20:99",
        "1/1:0,20:20:99",
        "0/1:10,10:20:99",
        "./.:0,0:0:0",
    ];
    let sample_strings_triallelic: [&str; 5] = [
        "0/0:20,0,0:20:99",
        "0/1:10,10,0:20:99",
        "1/1:0,20,0:20:99",
        "0/2:10,0,10:20:99",
        "./.:0,0,0:0:0",
    ];

    let alts = ["C", "T", "G"];

    let mut buf = String::with_capacity(n_variants * (n_samples * 20 + 60));

    // Header
    buf.push_str("##fileformat=VCFv4.2\n");
    buf.push_str(&format!(
        "##contig=<ID=chr1,length={}>\n",
        n_chr1 as u64 * 1000 + 1
    ));
    buf.push_str(&format!(
        "##contig=<ID=chr2,length={}>\n",
        n_chr2 as u64 * 1000 + 1
    ));
    buf.push_str("##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n");
    buf.push_str("##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n");
    buf.push_str("##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n");
    buf.push_str(
        "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n",
    );
    buf.push_str(
        "##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth\">\n",
    );
    buf.push_str(
        "##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Read depth\">\n",
    );
    buf.push_str(
        "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">\n",
    );

    buf.push_str("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for i in 0..n_samples {
        write!(buf, "\tS{i:04}").unwrap();
    }
    buf.push('\n');

    write_contig(&mut buf, "chr1", n_chr1, n_samples, &alts, &sample_strings_biallelic, &sample_strings_triallelic, 0);
    write_contig(&mut buf, "chr2", n_chr2, n_samples, &alts, &sample_strings_biallelic, &sample_strings_triallelic, n_chr1 as u64 * 1000);

    buf
}

fn write_contig(
    buf: &mut String,
    chrom: &str,
    n_records: usize,
    n_samples: usize,
    alts: &[&str; 3],
    bin: &[&str; 5],
    tn: &[&str; 5],
    start_pos: u64,
) {
    let spacing: u64 = 1000;
    for i in 0..n_records {
        let pos = start_pos + 1 + i as u64 * spacing;
        let is_multiallelic = (pos % 7) == 0;
        let alt = alts[(i / 1000) % alts.len()];

        let sample_patterns = if is_multiallelic { tn } else { bin };

        write!(buf, "{chrom}\t{pos}\t.\tA\t").unwrap();
        if is_multiallelic {
            let alt2 = alts[((i / 1000) + 1) % alts.len()];
            write!(buf, "{alt},{alt2}").unwrap();
        } else {
            buf.push_str(alt);
        }
        buf.push_str("\t.\t.\t.\tGT:AD:DP:GQ");

        for s in 0..n_samples {
            let gt_code = (pos as usize + s) % 5;
            buf.push('\t');
            buf.push_str(sample_patterns[gt_code]);
        }
        buf.push('\n');
    }
}

// ---------------------------------------------------------------------------
// Shared test data (generated once, shared across all #[ignore] tests)
// ---------------------------------------------------------------------------

fn test_dir() -> &'static Path {
    static DIR: LazyLock<PathBuf> = LazyLock::new(|| {
        let dir = tempfile::tempdir().expect("tempdir");
        dir.keep()
    });
    &DIR
}

fn vcf_path() -> &'static Path {
    static PATH: LazyLock<PathBuf> = LazyLock::new(|| {
        let path = test_dir().join("large.vcf");
        let text = generate_large_vcf(N_VARIANTS, N_SAMPLES);
        std::fs::write(&path, &text).expect("write large.vcf");
        path
    });
    &PATH
}

fn vcf_text() -> &'static str {
    static TEXT: LazyLock<String> =
        LazyLock::new(|| generate_large_vcf(N_VARIANTS, N_SAMPLES));
    &TEXT
}

fn indexed_vcf_path() -> &'static Path {
    static PATH: LazyLock<PathBuf> = LazyLock::new(|| {
        let gz = test_dir().join("large.vcf.gz");
        let status = Command::new("bgzip")
            .arg("-c")
            .arg(vcf_path())
            .stdout(File::create(&gz).expect("create large.vcf.gz"))
            .status()
            .expect("bgzip");
        assert!(status.success(), "bgzip failed");
        let status = Command::new("tabix")
            .arg("-p")
            .arg("vcf")
            .arg(&gz)
            .status()
            .expect("tabix");
        assert!(status.success(), "tabix failed");
        gz
    });
    &PATH
}

// ---------------------------------------------------------------------------
// Helpers
// ---------------------------------------------------------------------------

fn varanrs_binary() -> String {
    std::env::var("CARGO_BIN_EXE_varanrs").unwrap_or_else(|_| "target/release/varanrs".into())
}

fn count_data_lines(text: &str) -> usize {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.is_empty())
        .count()
}

/// Assert elapsed time is below threshold, printing timings.
fn assert_throughput(label: &str, elapsed: Duration, threshold: Duration) {
    let secs = elapsed.as_secs_f64();
    let thresh = threshold.as_secs_f64();
    eprintln!(
        "throughput: {label}  {secs:.3}s  (limit {thresh:.3}s)"
    );
    assert!(
        elapsed < threshold,
        "throughput: {label} took {secs:.3}s, exceeding limit of {thresh:.3}s"
    );
}

/// Spawn varanrs as a subprocess, write `stdin_text` to its stdin,
/// return stdout as a String and the elapsed wall time.
fn subprocess_pipe(
    args: &[&str],
    stdin_text: &str,
) -> Result<(String, Duration)> {
    let mut child = Command::new(varanrs_binary())
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .context("spawn varanrs")?;

    let t0 = Instant::now();

    let mut stdin = child.stdin.take().context("take stdin")?;
    stdin
        .write_all(stdin_text.as_bytes())
        .context("write stdin")?;
    drop(stdin);

    let output = child.wait_with_output().context("wait")?;
    let elapsed = t0.elapsed();

    if !output.status.success() {
        let stderr = String::from_utf8_lossy(&output.stderr);
        return Err(anyhow::anyhow!(
            "varanrs exited with {}: {}",
            output.status,
            stderr
        ));
    }
    let stdout = String::from_utf8(output.stdout).context("utf8 stdout")?;
    Ok((stdout, elapsed))
}

// =========================================================================
// allelefilter — streaming (file path)
// =========================================================================

#[test]
#[ignore]
fn throughput_allelefilter_file() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let out = dir.path().join("out.vcf");

    let args = varanrs::commands::allelefilter::AllelefilterArgs {
        input: Some(vcf_path().display().to_string()),
        output: Some(out.display().to_string()),
        min_ac: None,
        min_af: None,
    };

    let t0 = Instant::now();
    varanrs::commands::allelefilter::run(args)?;
    let elapsed = t0.elapsed();

    let out_text = std::fs::read_to_string(&out)?;
    let n_out = count_data_lines(&out_text);
    assert_eq!(n_out, N_VARIANTS, "identity: all records should be kept");

    assert_throughput(
        "allelefilter_file",
        elapsed,
        Duration::from_secs(90),
    );
    Ok(())
}

// =========================================================================
// allelefilter — streaming (stdin pipe)
// =========================================================================

#[test]
#[ignore]
fn throughput_allelefilter_stdin() -> Result<()> {
    let text = vcf_text();
    let expected_lines = count_data_lines(text);

    let (stdout, elapsed) = subprocess_pipe(&["allelefilter"], text)?;

    let n_out = count_data_lines(&stdout);
    assert_eq!(
        n_out, expected_lines,
        "stdin identity: all records should be kept"
    );

    assert_throughput(
        "allelefilter_stdin",
        elapsed,
        Duration::from_secs(90),
    );
    Ok(())
}

// =========================================================================
// allelefilter — with --min-ac to exercise rebuild path
// =========================================================================

#[test]
#[ignore]
fn throughput_allelefilter_filter() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let out = dir.path().join("out.vcf");

    let args = varanrs::commands::allelefilter::AllelefilterArgs {
        input: Some(vcf_path().display().to_string()),
        output: Some(out.display().to_string()),
        min_ac: Some(50),
        min_af: None,
    };

    let t0 = Instant::now();
    varanrs::commands::allelefilter::run(args)?;
    let elapsed = t0.elapsed();

    let out_text = std::fs::read_to_string(&out)?;
    let n_out = count_data_lines(&out_text);
    // With min-ac 50 most biallelic sites pass (AC ~= 80) but some
    // triallelic sites have the second ALT dropped.  No site should be
    // completely filtered so we expect all records to survive.
    assert_eq!(
        n_out, N_VARIANTS,
        "with --min-ac 50 all sites should survive (biallelic AC=80 > 50)"
    );

    assert_throughput(
        "allelefilter_filter",
        elapsed,
        Duration::from_secs(120),
    );
    Ok(())
}

// =========================================================================
// uppercase-alleles — streaming (stdin pipe)
// =========================================================================

#[test]
#[ignore]
fn throughput_uppercase_streaming() -> Result<()> {
    let text = vcf_text();
    let expected_lines = count_data_lines(text);

    let (stdout, elapsed) = subprocess_pipe(&["uppercase-alleles"], text)?;

    let n_out = count_data_lines(&stdout);
    assert_eq!(n_out, expected_lines, "all records should be kept");

    assert_throughput(
        "uppercase_streaming",
        elapsed,
        Duration::from_secs(120),
    );
    Ok(())
}

// =========================================================================
// uppercase-alleles — parallel (indexed VCF.gz)
// =========================================================================

#[test]
#[ignore]
fn throughput_uppercase_parallel() -> Result<()> {
    if !common::check_tools(&["bgzip", "tabix"]) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let out = dir.path().join("out.vcf.gz");

    let args = varanrs::commands::uppercase_alleles::UppercaseAllelesArgs {
        indexed: varanrs::args::IndexedInput {
            input: Some(indexed_vcf_path().to_path_buf()),
            threads: Some(4),
            contig: vec![],
            fai: None,
        },
        output: Some(out.display().to_string()),
        window_size: None,
    };

    let t0 = Instant::now();
    varanrs::commands::uppercase_alleles::run(args)?;
    let elapsed = t0.elapsed();

    assert_throughput(
        "uppercase_parallel",
        elapsed,
        Duration::from_secs(120),
    );
    Ok(())
}

// =========================================================================
// snpsketch — parallel (indexed VCF.gz)
// =========================================================================

#[test]
#[ignore]
fn throughput_snpsketch() -> Result<()> {
    if !common::check_tools(&["bgzip", "tabix"]) {
        return Ok(());
    }

    let dir = tempfile::tempdir()?;

    let args = varanrs::commands::snpsketch::SketchArgs {
        indexed: varanrs::args::IndexedInput {
            input: Some(indexed_vcf_path().to_path_buf()),
            threads: Some(4),
            contig: vec!["chr1".into()],
            fai: None,
        },
        stride: 2,
        chunk: 16384,
        pairs: dir.path().join("pairs.csv"),
        sample_stats: dir.path().join("sample_stats.csv"),
        genotypes: None,
        info_stats: None,
        info_hist_bins: 100,
    };

    let t0 = Instant::now();
    varanrs::commands::snpsketch::run(args)?;
    let elapsed = t0.elapsed();

    assert_throughput(
        "snpsketch",
        elapsed,
        Duration::from_secs(120),
    );
    Ok(())
}
