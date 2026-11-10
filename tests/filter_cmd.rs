//! End-to-end tests for the `filter` and `fill-tags` subcommands.

use std::fs::File;
use std::io::Write as _;
use std::path::Path;
use std::process::{Command, Stdio};

use anyhow::Result;

mod common;

const VCF: &str = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##contig=<ID=chr2,length=10000>
##FILTER=<ID=q10,Description=\"Quality below 10\">
##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
##INFO=<ID=DB,Number=0,Type=Flag,Description=\"dbSNP\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype quality\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t100\trs1\tA\tC\t50\tPASS\tDP=30;DB\tGT:GQ\t0/1:40\t0/0:50\t1/1:30\t0/0:20
chr1\t2500\t.\tG\tT,GA\t5\tq10\tDP=10\tGT:GQ\t0/2:10\t0/0:10\t0/0:10\t./.:.
chr1\t5200\t.\tT\tG\t.\t.\tDP=40\tGT:GQ\t0/0:99\t0/0:99\t0/0:99\t0/0:99
chr2\t300\t.\tAC\tA\t70\t.\t.\tGT:GQ\t0/1:5\t0/1:5\t0/1:5\t0/1:5
chr2\t7100\t.\tC\tT\t20\tPASS\tDP=5\tGT:GQ\t1:60\t0/1:60\t.:.\t1/1:60
";

fn varanrs() -> String {
    std::env::var("CARGO_BIN_EXE_varanrs").unwrap_or_else(|_| "target/release/varanrs".into())
}

/// Run varanrs on VCF text via stdin; return stdout.
fn run_stdin(args: &[&str], input: &str) -> String {
    let mut child = Command::new(varanrs())
        .args(args)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    child.stdin.take().unwrap().write_all(input.as_bytes()).unwrap();
    let out = child.wait_with_output().unwrap();
    assert!(out.status.success(), "varanrs {args:?} failed: {}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

fn body(text: &str) -> Vec<String> {
    text.lines().filter(|l| !l.starts_with('#')).map(str::to_string).collect()
}

fn positions(text: &str) -> Vec<String> {
    body(text).iter().map(|l| l.split('\t').nth(1).unwrap().to_string()).collect()
}

#[test]
fn exclude_expression_drops_matching_sites() {
    let out = run_stdin(&["filter", "-e", "QUAL<30 || DP>35"], VCF);
    // QUAL missing at 5200 is not < 30, but DP=40 > 35 excludes it.
    assert_eq!(positions(&out), vec!["100", "300"]);
    // bcftools filter sets PASS on passing sites with an empty FILTER
    assert!(body(&out)[1].contains("\tPASS\t"));
}

#[test]
fn include_expression_with_type_filter_and_on_the_fly_ac() {
    let out = run_stdin(&["filter", "-i", "TYPE=\"snp\" && FILTER=\"PASS\" && AC>=2"], VCF);
    assert_eq!(positions(&out), vec!["100", "7100"]);
}

#[test]
fn format_and_genotype_expressions_are_supported() {
    // site selected if any sample is het with GQ < 20
    let out = run_stdin(&["filter", "-i", "GT=\"het\" & FMT/GQ<20"], VCF);
    assert_eq!(positions(&out), vec!["2500", "300"]);
}

/// GT-indexed subscripts use per-record genotype data inside bcftools'
/// filter.c; consecutive records must each be evaluated on their own GTs.
#[test]
fn gt_indexed_subscripts_follow_each_record() {
    let vcf = "##fileformat=VCFv4.2\n##contig=<ID=1>\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"g\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"a\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n\
1\t1\t.\tA\tC\t.\t.\t.\tGT:AD\t0/0:10,0\n\
1\t2\t.\tA\tC\t.\t.\t.\tGT:AD\t1/1:10,0\n\
1\t3\t.\tA\tC\t.\t.\t.\tGT:AD\t1/1:0,10\n";
    let out = run_stdin(&["filter", "-i", "FMT/AD[GT]==10"], vcf);
    assert_eq!(positions(&out), vec!["1", "3"]);
}

#[test]
fn soft_filter_annotates_instead_of_dropping() {
    let out = run_stdin(&["filter", "-e", "QUAL<30", "-s", "LowQual"], VCF);
    assert!(out.contains("##FILTER=<ID=LowQual,Description=\"Set if true: QUAL<30\">"));
    let filters: Vec<String> = body(&out).iter().map(|l| l.split('\t').nth(6).unwrap().to_string()).collect();
    assert_eq!(filters, vec!["PASS", "LowQual", "PASS", "PASS", "LowQual"]);
}

#[test]
/// Expected values verified against `bcftools +fill-tags` 1.23.
fn fill_tags_matches_bcftools_semantics() {
    let out = run_stdin(&["fill-tags", "-t", "AN,AC,AF,NS,AC_Hemi,F_MISSING"], VCF);
    let lines = body(&out);
    let info: Vec<&str> = lines.iter().map(|l| l.split('\t').nth(7).unwrap()).collect();
    assert_eq!(info[0], "DP=30;DB;F_MISSING=0;NS=4;AN=8;AF=0.375;AC=3;AC_Hemi=0");
    // multiallelic with a missing sample
    assert_eq!(info[1], "DP=10;F_MISSING=0.25;NS=3;AN=6;AF=0,0.166667;AC=0,1;AC_Hemi=0,0");
    // haploid call is hemizygous
    assert_eq!(info[4], "DP=5;F_MISSING=0.25;NS=3;AN=5;AF=0.8;AC=4;AC_Hemi=1");
}

fn write_indexed_vcf(path: &Path, text: &str) {
    let raw = path.with_extension("tmp.vcf");
    File::create(&raw).unwrap().write_all(text.as_bytes()).unwrap();
    let status = Command::new("bgzip").arg("-c").arg(&raw).stdout(File::create(path).unwrap()).status().unwrap();
    assert!(status.success());
    let status = Command::new("tabix").args(["-p", "vcf"]).arg(path).status().unwrap();
    assert!(status.success());
}

fn bcftools_body(args: &[&str]) -> String {
    let out = Command::new("bcftools").args(args).output().unwrap();
    assert!(out.status.success(), "{}", String::from_utf8_lossy(&out.stderr));
    String::from_utf8(out.stdout).unwrap()
}

/// Region-parallel filter + soft-filter + fill-tags, compared byte-for-byte
/// with the equivalent bcftools pipeline.
#[test]
fn parallel_filter_matches_bcftools() -> Result<()> {
    if !common::check_tools(&["bgzip", "tabix", "bcftools"]) {
        return Ok(());
    }
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    write_indexed_vcf(&input, VCF);

    let ours = dir.path().join("ours.vcf.gz");
    varanrs::commands::filter::run(varanrs::commands::filter::FilterArgs {
        indexed: varanrs::args::IndexedInput { input: Some(input.clone()), threads: Some(3), contig: vec![], fai: None },
        output: Some(ours.display().to_string()),
        window_size: Some(1000),
        exclude: Some("MAF<0.2 || QUAL<10".into()),
        include: None,
        soft_filter: Some("Bad".into()),
        mode: None,
        min_ac: None,
        min_af: None,
        fill_tags: Some("all".into()),
        samples_file: None,
        include_genotypes: None,
        exclude_genotypes: None,
    })?;

    let theirs = dir.path().join("theirs.vcf");
    let fill = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "bcftools +fill-tags {} -Ou -- -t all | bcftools filter -e 'MAF<0.2 || QUAL<10' -s Bad -Ov -o {}",
            input.display(),
            theirs.display()
        ))
        .status()?;
    if !fill.success() {
        eprintln!("SKIP: bcftools +fill-tags plugin unavailable");
        return Ok(());
    }
    let a = bcftools_body(&["view", "-H", ours.to_str().unwrap()]);
    let b = bcftools_body(&["view", "-H", theirs.to_str().unwrap()]);
    assert_eq!(a, b);
    // and the raw text is identical too (htslib formatting)
    let raw = Command::new("bgzip").arg("-dc").arg(&ours).output()?.stdout;
    let raw_body = body(&String::from_utf8(raw)?).join("\n") + "\n";
    assert_eq!(raw_body, b);
    Ok(())
}

fn sample_list(dir: &Path, names: &[&str]) -> String {
    let p = dir.join("samples.txt");
    std::fs::write(&p, names.iter().map(|n| format!("{n}\n")).collect::<String>()).unwrap();
    p.display().to_string()
}

fn sample_columns(text: &str) -> Vec<String> {
    let chrom = text.lines().find(|l| l.starts_with("#CHROM")).unwrap();
    chrom.split('\t').skip(9).map(str::to_string).collect()
}

/// `-S` subsets (and reorders) samples first, recomputes AC/AN, and the
/// allele filter then counts only the kept samples.
#[test]
fn samples_file_subsets_before_allele_filter() {
    let dir = tempfile::tempdir().unwrap();
    let list = sample_list(dir.path(), &["S3", "S2"]);
    let out = run_stdin(&["filter", "-S", &list], VCF);
    assert_eq!(sample_columns(&out), vec!["S3", "S2"]);
    let first = &body(&out)[0];
    let cols: Vec<&str> = first.split('\t').collect();
    assert_eq!(cols[7], "DP=30;DB;AC=2;AN=4");
    assert_eq!(&cols[9..], &["1/1:30", "0/0:50"]);

    // 2500's ALT alleles are only carried by S1: dropped once S1 is gone.
    let out = run_stdin(&["filter", "-S", &list, "--min-ac", "1"], VCF);
    assert_eq!(positions(&out), vec!["100", "300", "7100"]);
}

#[test]
fn samples_file_exclusion_keeps_header_order() {
    let dir = tempfile::tempdir().unwrap();
    let list = sample_list(dir.path(), &["S3", "S1"]);
    let out = run_stdin(&["filter", "-S", &format!("^{list}")], VCF);
    assert_eq!(sample_columns(&out), vec!["S2", "S4"]);
}

#[test]
fn samples_file_rejects_unknown_samples() {
    let dir = tempfile::tempdir().unwrap();
    let list = sample_list(dir.path(), &["S1", "nope"]);
    let mut child = Command::new(varanrs())
        .args(["filter", "-S", &list])
        .stdin(Stdio::piped())
        .stdout(Stdio::null())
        .stderr(Stdio::piped())
        .spawn()
        .unwrap();
    let _ = child.stdin.take().unwrap().write_all(VCF.as_bytes());
    let out = child.wait_with_output().unwrap();
    assert!(!out.status.success());
    assert!(String::from_utf8_lossy(&out.stderr).contains("\"nope\""));
}

/// Region-parallel `-S` + filter, byte-identical to `bcftools view -S` piped
/// into `bcftools filter` (header included).
#[test]
fn parallel_samples_file_matches_bcftools() -> Result<()> {
    if !common::check_tools(&["bgzip", "tabix", "bcftools"]) {
        return Ok(());
    }
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    write_indexed_vcf(&input, VCF);
    let list = sample_list(dir.path(), &["S4", "S1"]);

    let ours = dir.path().join("ours.vcf.gz");
    varanrs::commands::filter::run(varanrs::commands::filter::FilterArgs {
        indexed: varanrs::args::IndexedInput { input: Some(input.clone()), threads: Some(3), contig: vec![], fai: None },
        output: Some(ours.display().to_string()),
        window_size: Some(1000),
        exclude: Some("AC=0".into()),
        include: None,
        soft_filter: None,
        mode: None,
        min_ac: None,
        min_af: None,
        fill_tags: None,
        samples_file: Some(list.clone()),
        include_genotypes: None,
        exclude_genotypes: None,
    })?;
    let theirs = Command::new("sh")
        .arg("-c")
        .arg(format!(
            "bcftools view --no-version -S {list} {} -Ou | bcftools filter --no-version -e 'AC=0' -Ov",
            input.display()
        ))
        .output()?;
    assert!(theirs.status.success(), "{}", String::from_utf8_lossy(&theirs.stderr));
    let raw = Command::new("bgzip").arg("-dc").arg(&ours).output()?.stdout;
    assert_eq!(String::from_utf8(raw)?, String::from_utf8(theirs.stdout)?);
    Ok(())
}

/// `-E` sets GT to missing (any ploidy) for samples where the expression is
/// true, and decrements INFO/AC and INFO/AN when present.
#[test]
fn exclude_genotypes_masks_failing_samples() {
    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"AC\">
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"AN\">
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4
chr1\t100\t.\tA\tC\t50\t.\tAC=4;AN=7\tGT:DP\t0|1:2\t1/1:10\t1:0\t0/0:3
chr1\t200\t.\tA\tC\t50\t.\t.\tGT:DP\t0/1:9\t1/1:.\t0/1:1\t./.:2
";
    let out = run_stdin(&["filter", "-E", "FMT/DP<3"], vcf);
    let b = body(&out);
    let cols: Vec<&str> = b[0].split('\t').collect();
    assert_eq!(cols[6], ".", "FILTER is left alone");
    assert_eq!(cols[7], "AC=2;AN=4");
    assert_eq!(&cols[9..], &["./.:2", "1/1:10", ".:0", "0/0:3"]);
    let cols: Vec<&str> = b[1].split('\t').collect();
    assert_eq!(&cols[9..], &["0/1:9", "1/1:.", "./.:1", "./.:2"]);

    // -I masks where the expression is not true, so a missing DP (for
    // which neither DP<3 nor DP>=3 is true) is masked by -I but not by -E.
    let inc = run_stdin(&["filter", "-I", "FMT/DP>=3"], vcf);
    assert_eq!(body(&inc)[0], b[0]);
    assert!(body(&inc)[1].ends_with("\t0/1:9\t./.:.\t./.:1\t./.:2"));
    // site-level expressions mask whole sites
    let site = run_stdin(&["filter", "-I", "POS>150"], vcf);
    assert!(body(&site)[0].ends_with("\t./.:2\t./.:10\t.:0\t./.:3"));
    assert!(body(&site)[0].contains("\tAC=0;AN=0\t"));
}
