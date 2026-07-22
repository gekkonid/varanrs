//! End-to-end test for the `upper` subcommand.

use std::fs::File;
use std::io::Write as _;
use std::path::Path;

use anyhow::Result;
use noodles_util::variant;
use noodles_vcf as vcf;

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

#[test]
fn upper_command_uppercases_alleles_round_trip() -> Result<()> {
    let dir = tempfile::tempdir()?;
    let input = dir.path().join("in.vcf.gz");
    let output = dir.path().join("out.vcf.gz");

    let text = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=10000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\ta\tt\t.\t.\t.\tGT\t0/1
chr1\t200\t.\tc\tG,a\t.\t.\t.\tGT\t1/2
chr1\t300\t.\tn\tA\t.\t.\t.\tGT\t0/1
chr1\t400\t.\tACgt\tA\t.\t.\t.\tGT\t0/1
chr1\t500\t.\tA\tT\t.\t.\t.\tGT\t0/0
";
    write_indexed_vcf(&input, text);

    let args = varanrs::commands::upper::UpperArgs {
        indexed: varanrs::args::IndexedInput {
            input: Some(input.clone()),
            threads: Some(2),
            contig: vec![],
            fai: None,
        },
        output: Some(output.display().to_string()),
        window_size: Some(1000),
    };
    varanrs::commands::upper::run(args)?;

    let mut reader = variant::io::reader::Builder::default().build_from_path(&output)?;
    let header = reader.read_header()?;
    let mut count = 0;
    for r in reader.records(&header) {
        let r = r?;
        let buf = vcf::variant::RecordBuf::try_from_variant_record(&header, r.as_ref())?;
        let reference = buf.reference_bases();
        assert_eq!(reference, &reference.to_ascii_uppercase(), "REF not upper");
        for alt in buf.alternate_bases().as_ref() {
            assert_eq!(alt, &alt.to_ascii_uppercase(), "ALT not upper");
        }
        count += 1;
    }
    assert_eq!(count, 5);
    Ok(())
}

fn varanrs_binary() -> String {
    std::env::var("CARGO_BIN_EXE_varanrs").unwrap_or_else(|_| "target/debug/varanrs".into())
}

#[test]
fn stdin_stdout_upper_identity() -> Result<()> {
    use std::io::Write as _;
    use std::process::{Command, Stdio};

    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\ta\tt\t.\t.\t.\tGT\t0/1
chr1\t200\t.\tc\tG,a\t.\t.\t.\tGT\t1/2
";

    let mut child = Command::new(varanrs_binary())
        .arg("upper")
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .spawn()?;

    child.stdin.take().unwrap().write_all(vcf.as_bytes())?;
    let output = child.wait_with_output()?;
    assert!(output.status.success());
    let out = String::from_utf8_lossy(&output.stdout);
    assert!(out.contains("chr1\t100\t.\tA\tT"), "REF/ALT should be uppercased");
    assert!(out.contains("chr1\t200\t.\tC\tG,A"), "multi-allelic should be uppercased");
    Ok(())
}

#[test]
fn stdin_upper_from_bcf_pipe() -> Result<()> {
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
chr1\t100\t.\ta\tt\t.\t.\t.\tGT\t0/1
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
        .arg("upper")
        .stdin(producer.stdout.unwrap())
        .stdout(Stdio::piped())
        .spawn()?;

    let output = consumer.wait_with_output()?;
    assert!(output.status.success());
    let out = String::from_utf8_lossy(&output.stdout);
    assert!(out.contains("A\tT"), "uncompressed BCF pipe should uppercase alleles");
    Ok(())
}
