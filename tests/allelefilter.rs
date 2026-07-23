//! End-to-end tests for the `filter` subcommand using text VCF I/O.

use std::fs;
use std::path::Path;

use anyhow::Result;

fn run_filter(input: &Path, output: &Path, min_ac: Option<u32>, min_af: Option<f64>) -> Result<()> {
    let args = varanrs::commands::allelefilter::AllelefilterArgs {
        input: Some(input.display().to_string()),
        output: Some(output.display().to_string()),
        min_ac,
        min_af,
    };
    varanrs::commands::allelefilter::run(args)
}

fn read_text(path: &Path) -> String {
    fs::read_to_string(path).unwrap()
}

/// Test that, given no thresholds, we recompute the INFO tags correctly
#[test]
fn no_threshold_recomputes_counts() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\t.\tGT\t0/1\t0/2\t1/1\t0/0\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, None, None)?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\tAN=8;AC=3,1;AF=0.375,0.125;F_MISSING=0\tGT\t0/1\t0/2\t1/1\t0/0\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}


/// Test that, given an AC threshold that one allele fails, we turn it back into a biallelic site
/// with missing genotypes that had the 2nd filtered-out ALT, and that this updates the INFO fields
#[test]
fn min_ac_drops_allele() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\t.\tGT\t0/1\t0/2\t1/1\t0/0\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, Some(3), None)?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC\t.\t.\tAN=6;AC=3;AF=0.5;F_MISSING=0.25\tGT\t0/1\t./.\t1/1\t0/0\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}


#[test]
fn min_af_drops_allele() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\t.\tGT\t0/1\t0/2\t1/1\t0/0\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, None, Some(0.25))?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC\t.\t.\tAN=6;AC=3;AF=0.5;F_MISSING=0.25\tGT\t0/1\t./.\t1/1\t0/0\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}

#[test]
fn both_thresholds_combined() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\t.\tGT\t0/1\t0/2\t1/1\t0/0\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, Some(3), Some(0.25))?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC\t.\t.\tAN=6;AC=3;AF=0.5;F_MISSING=0.25\tGT\t0/1\t./.\t1/1\t0/0\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}

#[test]
fn all_alt_dropped_returns_none() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\t.\tGT\t0/1\t0/2\t1/1\t0/0\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, Some(5), None)?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}

#[test]
fn ar_format_resliced() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\t.\tGT:AD\t0/1:10,5,2\t0/2:8,3,1\t1/1:6,4,0\t0/0:12,0,0\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, Some(3), None)?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
##FORMAT=<ID=AD,Number=R,Type=Integer,Description=\"Allelic depth\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC\t.\t.\tAN=6;AC=3;AF=0.5;F_MISSING=0.25\tGT:AD\t0/1:10,5\t./.:.	1/1:6,4\t0/0:12,0\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}

#[test]
fn info_ar_resliced() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
##INFO=<ID=DUMMY_R,Number=R,Type=Integer,Description=\"dummy R\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\tDUMMY_R=100,50,20\tGT\t0/1\t0/2\t1/1\t0/0\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, Some(3), None)?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
##INFO=<ID=DUMMY_R,Number=R,Type=Integer,Description=\"dummy R\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC\t.\t.\tDUMMY_R=100,50;AN=6;AC=3;AF=0.5;F_MISSING=0.25\tGT\t0/1\t./.\t1/1\t0/0\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}

#[test]
fn no_alt_recomputes_an_only() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\t.\t.\t.\t.\tGT\t0/0\t0/0\t./.\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, None, None)?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t100\t.\tA\t.\t.\t.\tAN=4\tGT\t0/0\t0/0\t./.\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}

#[test]
fn multiple_sites_some_dropped() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf");
    let output = dir.path().join("out.vcf");

    let in_vcf = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC,T\t.\t.\t.\tGT\t0/1\t0/2\t1/1\t0/0\n\
chr1\t200\t.\tG\tA,T\t.\t.\t.\tGT\t0/0\t0/0\t0/0\t0/1\n\
chr1\t300\t.\tT\tC\t.\t.\t.\tGT\t1/1\t1/1\t1/1\t1/1\n";

    fs::write(&input, in_vcf)?;
    run_filter(&input, &output, Some(3), None)?;

    let out = read_text(&output);
    let expected = "##fileformat=VCFv4.2\n\
##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count\">\n\
##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele frequency\">\n\
##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles\">\n\
##INFO=<ID=F_MISSING,Number=1,Type=Float,Description=\"Fraction of missing genotypes\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\tS4\n\
chr1\t100\t.\tA\tC\t.\t.\tAN=6;AC=3;AF=0.5;F_MISSING=0.25\tGT\t0/1\t./.\t1/1\t0/0\n\
chr1\t300\t.\tT\tC\t.\t.\tAN=8;AC=8;AF=1;F_MISSING=0\tGT\t1/1\t1/1\t1/1\t1/1\n";

    assert_eq!(out, expected, "output VCF mismatch");
    Ok(())
}

fn varanrs_binary() -> String {
    std::env::var("CARGO_BIN_EXE_varanrs").unwrap_or_else(|_| "target/debug/varanrs".into())
}

#[test]
fn stdin_stdout_filter_identity() -> Result<()> {
    use std::io::Write as _;
    use std::process::{Command, Stdio};

    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1
chr1\t200\t.\tG\tC\t.\t.\t.\tGT\t0/0
";

    let mut child = Command::new(varanrs_binary())
        .arg("allelefilter")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    child.stdin.take().unwrap().write_all(vcf.as_bytes())?;
    let output = child.wait_with_output()?;
    assert!(output.status.success());
    let out = String::from_utf8_lossy(&output.stdout);
    assert!(out.contains("chr1\t100"));
    assert!(out.contains("chr1\t200"));
    Ok(())
}

#[test]
fn stdin_stdout_filter_min_ac() -> Result<()> {
    use std::io::Write as _;
    use std::process::{Command, Stdio};

    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1\t0/0
chr1\t200\t.\tG\tC\t.\t.\t.\tGT\t1/1\t0/1
";

    let mut child = Command::new(varanrs_binary())
        .args(["allelefilter", "--min-ac", "3"])
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    child.stdin.take().unwrap().write_all(vcf.as_bytes())?;
    let output = child.wait_with_output()?;
    assert!(output.status.success());
    let out = String::from_utf8_lossy(&output.stdout);
    assert!(!out.contains("chr1\t100"), "low-AC should be filtered");
    assert!(out.contains("chr1\t200"), "high-AC should be kept");
    Ok(())
}

#[test]
fn stdin_filter_from_bcf_pipe() -> Result<()> {
    use std::process::{Command, Stdio};

    if !Command::new("bcftools").arg("--version").output().map(|o| o.status.success()).unwrap_or(false) {
        eprintln!("SKIP: bcftools not available");
        return Ok(());
    }

    let dir = tempfile::tempdir()?;
    let vcf_path = dir.path().join("in.vcf");
    std::fs::write(&vcf_path, "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/1
")?;

    let bcf_path = dir.path().join("in.bcf");
    assert!(Command::new("bcftools")
        .args(["convert", "-O", "b", "-o"])
        .arg(&bcf_path).arg(&vcf_path)
        .status()?.success());

    let producer = Command::new("bcftools")
        .args(["view", "-Ou"])
        .arg(&bcf_path)
        .stdout(Stdio::piped())
        .spawn()?;

    let consumer = Command::new(varanrs_binary())
        .arg("allelefilter")
        .stdin(producer.stdout.unwrap())
        .stdout(Stdio::piped())
        .spawn()?;

    let output = consumer.wait_with_output()?;
    assert!(output.status.success());
    let out = String::from_utf8_lossy(&output.stdout);
    assert!(out.contains("chr1\t100"), "record should survive filter on BCF pipe");
    Ok(())
}
