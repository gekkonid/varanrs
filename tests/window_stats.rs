//! Tests for `window-stats`.

use std::fs::File;
use std::io::Write as _;
use std::path::Path;
use std::process::{Command, Stdio};

use varanrs::hts;
use varanrs::windowstats::{TagSpec, WindowStats};

mod common;

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=2500>
##contig=<ID=chr2,length=1000>
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
##INFO=<ID=AF,Number=A,Type=Float,Description=\"AF\">
##INFO=<ID=AD,Number=R,Type=Integer,Description=\"AD\">
##INFO=<ID=STR,Number=1,Type=String,Description=\"String\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t1\t.\tA\tC\t.\t.\tDP=10;AF=0.5;AD=1,2\tGT\t0/1\t0/1\t0/0
chr1\t1000\t.\tA\tC,G\t.\t.\tDP=20;AF=0.25,0.75;AD=3,4,5\tGT\t0/1\t1/2\t0/0
chr1\t1001\t.\tA\tC\t.\t.\tAF=.\tGT\t1/1\t0/1\t./.
chr1\t2400\t.\tA\tC\t.\t.\tDP=7\tGT\t0/1\t0/1\t0/1
chr3\t5\t.\tA\tC\t.\t.\tDP=1\tGT\t0/0\t0/0\t0/0
";

const EXPECTED: &str = "\
contig,start,end,n_sites,DP,AF,AD
chr1,1,1000,2,15.000000,0.375000,3.000000
chr1,1001,2000,1,,,
chr1,2001,2500,1,7.000000,,
chr2,1,1000,0,,,
chr3,1,1000,1,1.000000,,
";

fn varanrs() -> String {
    std::env::var("CARGO_BIN_EXE_varanrs").unwrap_or_else(|_| "target/release/varanrs".into())
}

fn run(args: &[&str], stdin: Option<&str>) -> std::process::Output {
    let mut child = Command::new(varanrs())
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let input = child.stdin.take().unwrap();
    if let Some(text) = stdin {
        let mut input = input;
        input.write_all(text.as_bytes()).unwrap();
    }
    child.wait_with_output().unwrap()
}

fn run_ok(args: &[&str], stdin: Option<&str>) -> String {
    let out = run(args, stdin);
    assert!(out.status.success(), "varanrs {args:?} failed: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

fn write_indexed_vcf(path: &Path, text: &str) {
    let raw = path.with_extension("tmp.vcf");
    File::create(&raw).unwrap().write_all(text.as_bytes()).unwrap();
    let status = Command::new("bgzip").arg("-c").arg(&raw).stdout(File::create(path).unwrap()).status().unwrap();
    assert!(status.success());
    let status = Command::new("tabix").args(["-p", "vcf"]).arg(path).status().unwrap();
    assert!(status.success());
}

#[test]
fn accumulator_means_per_window() {
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.vcf");
    std::fs::write(&path, VCF).unwrap();
    let (header, records) = hts::read_all(path.to_str().unwrap()).unwrap();
    let tags = TagSpec::parse("DP,INFO/AF,AD", &header).unwrap();
    // Split across two accumulators to exercise merge.
    let mut a = WindowStats::new(1000, tags.clone());
    let mut b = WindowStats::new(1000, tags);
    for (i, rec) in records.iter().enumerate() {
        if i % 2 == 0 { a.process_record(rec) } else { b.process_record(rec) }
    }
    a.merge(b);
    let mut out = Vec::new();
    a.write_csv(&mut out, &hts::header_contigs(&header)).unwrap();
    assert_eq!(String::from_utf8(out).unwrap(), EXPECTED);
}

#[test]
fn rejects_unknown_and_non_numeric_tags() {
    for tag in ["NOPE", "STR"] {
        let out = run(&["window-stats", "-t", tag, "-w", "1000"], Some(VCF));
        assert!(!out.status.success(), "-t {tag} should fail");
    }
}

#[test]
fn stdin_and_fill_tags() {
    assert_eq!(run_ok(&["window-stats", "-t", "DP,AF,AD", "-w", "1000"], Some(VCF)), EXPECTED);
    // AF is recomputed from genotypes; AC_Het is new.
    let out = run_ok(&["window-stats", "--fill-tags", "AF,AC_Het", "-t", "AC_Het,AF", "-w", "1000", "--contig", "chr1"], Some(VCF));
    assert_eq!(
        out,
        "\
contig,start,end,n_sites,AC_Het,AF
chr1,1,1000,2,2.000000,0.333333
chr1,1001,2000,1,1.000000,0.750000
chr1,2001,2500,1,3.000000,0.500000
"
    );
}

#[test]
fn parallel_matches_streaming() {
    if !common::check_tools(&["bgzip", "tabix"]) {
        return;
    }
    let dir = tempfile::tempdir().unwrap();
    let path = dir.path().join("t.vcf.gz");
    write_indexed_vcf(&path, VCF);
    let p = path.to_str().unwrap();
    let args = ["window-stats", "--fill-tags", "AC_Het", "-t", "DP,AF,AD,AC_Het", "-w", "1000"];
    let streamed = run_ok(&args, Some(VCF));
    for threads in ["1", "4"] {
        let mut a = args.to_vec();
        a.extend(["--threads", threads, p]);
        assert_eq!(run_ok(&a, None), streamed, "threads={threads}");
    }
}
