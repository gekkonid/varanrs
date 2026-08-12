use std::io::Cursor;

use indexmap::IndexMap;
use noodles_vcf as vcf;
use vcf::variant::RecordBuf;

use varanrs::snpsketch::SketchAccumulator;

fn parse_vcf_text(vcf_text: &str) -> (vcf::Header, Vec<RecordBuf>) {
    let mut reader = vcf::io::reader::Builder::default()
        .build_from_reader(Cursor::new(vcf_text.as_bytes().to_vec()))
        .unwrap();
    let header = reader.read_header().unwrap();
    let mut records = Vec::new();
    for result in reader.record_bufs(&header) {
        records.push(result.unwrap());
    }
    (header, records)
}

fn make_acc_from_vcf(vcf_text: &str) -> SketchAccumulator {
    let (header, records) = parse_vcf_text(vcf_text);
    let sample_ids: Vec<String> = header.sample_names().iter().cloned().collect();
    let contig_rank: IndexMap<String, usize> = header
        .contigs()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect();
    let mut acc = SketchAccumulator::new(sample_ids, contig_rank);
    for rec in &records {
        acc.process_record(rec);
    }
    acc
}

#[test]
fn smoke_empty_sample_set_errors_gracefully() {
    let acc = SketchAccumulator::new(vec![], IndexMap::new());
    assert_eq!(acc.n_samples(), 0);
    assert_eq!(acc.n_sites, 0);

    let mut csv = Cursor::new(Vec::new());
    acc.write_pairs_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();
    assert!(output.contains("sample_i,sample_j"));
    assert!(!output.contains("\n1"));
}

#[test]
fn single_sample_no_pairs() {
    let acc = SketchAccumulator::new(vec!["S1".to_string()], IndexMap::new());
    let mut csv = Cursor::new(Vec::new());
    acc.write_pairs_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();
    assert!(output.contains("sample_i"));
    assert!(!output.contains("\n1,"));
}

#[test]
fn genotype_classification_diploid() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	AA	HET	BB	MISS
chr1	100	.	A	G	.	.	.	GT	0/0	0/1	1/1	./.
";

    let (header, records) = parse_vcf_text(vcf);
    let ids: Vec<String> = header.sample_names().iter().cloned().collect();
    let contig_rank: IndexMap<String, usize> = header
        .contigs()
        .iter()
        .enumerate()
        .map(|(i, (name, _))| (name.to_string(), i))
        .collect();
    let mut acc = SketchAccumulator::new(ids, contig_rank);
    acc.process_record(&records[0]);

    assert_eq!(acc.n_sites, 1);

    assert!(acc.sample_called(0, 0));   // AA
    assert!(!acc.sample_any_alt(0, 0));
    assert!(!acc.sample_hom_alt(0, 0));

    assert!(acc.sample_called(1, 0));   // HET
    assert!(acc.sample_any_alt(1, 0));
    assert!(!acc.sample_hom_alt(1, 0));

    assert!(acc.sample_called(2, 0));   // BB
    assert!(acc.sample_any_alt(2, 0));
    assert!(acc.sample_hom_alt(2, 0));

    assert!(!acc.sample_called(3, 0));  // MISS
}

#[test]
fn genotype_classification_no_gt_in_format() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1
chr1	100	.	A	G	.	.	.	DP	15
";
    let acc = make_acc_from_vcf(vcf);
    assert_eq!(acc.n_sites, 1);
    assert!(!acc.sample_called(0, 0));
}

#[test]
fn pairwise_distance_all_same() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	0/0
";
    let acc = make_acc_from_vcf(vcf);
    let (n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 1);
    assert_eq!(n_diff, 0);
    assert_eq!(distance, 0.0);
}

#[test]
fn pairwise_distance_aa_vs_bb() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	1/1
";
    let acc = make_acc_from_vcf(vcf);
    let (n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 1);
    assert_eq!(n_diff, 2);
    assert_eq!(distance, 1.0);
}

#[test]
fn pairwise_distance_het_vs_aa() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	0/1
";
    let acc = make_acc_from_vcf(vcf);
    let (n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 1);
    assert_eq!(n_diff, 1);
    assert!((distance - 0.5).abs() < 0.001);
}

#[test]
fn pairwise_distance_het_vs_bb() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/1	1/1
";
    let acc = make_acc_from_vcf(vcf);
    let (n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 1);
    assert_eq!(n_diff, 1);
    assert!((distance - 0.5).abs() < 0.001);
}

#[test]
fn pairwise_distance_one_missing() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/1	./.
";
    let acc = make_acc_from_vcf(vcf);
    let (n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 0);
    assert_eq!(n_diff, 0);
    assert!(distance.is_nan());
}

#[test]
fn pairwise_distance_both_missing() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	./.	./.
";
    let acc = make_acc_from_vcf(vcf);
    let (_n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 0);
    assert!(distance.is_nan());
}

#[test]
fn pairwise_distance_multiple_sites() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	1/1
chr1	200	.	C	T	.	.	.	GT	0/1	0/1
chr1	300	.	G	A	.	.	.	GT	0/0	0/1
";
    let acc = make_acc_from_vcf(vcf);

    let (n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 3);
    assert_eq!(n_diff, 3); // 2 (aa-bb) + 0 (same het) + 1 (aa-het)
    let expected = 3.0 / 6.0; // n_diff / (2*n_common)
    assert!((distance - expected).abs() < 0.001);
}

#[test]
fn missingness_all_called() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	0/1
chr1	200	.	C	T	.	.	.	GT	1/1	0/0
";
    let acc = make_acc_from_vcf(vcf);
    assert_eq!(acc.sample_call_count(0), 2);
    assert_eq!(acc.sample_call_count(1), 2);
}

#[test]
fn missingness_partial() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	./.
chr1	200	.	C	T	.	.	.	GT	./.	0/1
";
    let acc = make_acc_from_vcf(vcf);
    assert_eq!(acc.sample_call_count(0), 1);
    assert_eq!(acc.sample_call_count(1), 1);
}

#[test]
fn sample_stats_csv_output() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	./.
chr1	200	.	C	T	.	.	.	GT	0/1	0/0
";
    let acc = make_acc_from_vcf(vcf);

    let mut csv = Cursor::new(Vec::new());
    acc.write_sample_stats_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();

    assert!(output.contains("sample_id,n_missing,n_total,miss_rate,n_het,het_rate,avg_dp"));
    assert!(output.contains("S1,0,2,0.000000,1,0.500000,"));
    assert!(output.contains("S2,1,2,0.500000,0,0.000000,"));
}

#[test]
fn sample_stats_csv_with_dp() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Depth\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT:DP	0/0:10	0/1:30
chr1	200	.	C	T	.	.	.	GT:DP	1/1:20	./.:.
";
    let acc = make_acc_from_vcf(vcf);

    let mut csv = Cursor::new(Vec::new());
    acc.write_sample_stats_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();

    assert!(output.contains("S1,0,2,0.000000,0,0.000000,15.0"));
    assert!(output.contains("S2,1,2,0.500000,1,1.000000,30.0"));
}

#[test]
fn pairs_csv_output_format() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	0/1
";
    let acc = make_acc_from_vcf(vcf);

    let mut csv = Cursor::new(Vec::new());
    acc.write_pairs_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();

    assert!(output.contains("sample_i,sample_j,sample_i_id,sample_j_id,n_diff,n_common,distance"));
    assert!(output.contains("0,1,S1,S2,"));
    assert!(output.contains(",1,1,"));
}

#[test]
fn genotypes_csv_single_site() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	1/1
";
    let acc = make_acc_from_vcf(vcf);

    let mut csv = Cursor::new(Vec::new());
    acc.write_genotypes_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();

    assert!(output.contains("site,S1,S2"));
    assert!(output.contains("0,0,2"));
}

#[test]
fn genotypes_csv_with_missing() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2	S3
chr1	100	.	A	G	.	.	.	GT	0/0	0/1	./.
";
    let acc = make_acc_from_vcf(vcf);

    let mut csv = Cursor::new(Vec::new());
    acc.write_genotypes_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();

    assert!(output.contains("site,S1,S2,S3"));
    assert!(output.contains("0,0,1,."));
}

#[test]
fn genotypes_csv_multiple_sites() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	0/0	0/1
chr1	200	.	C	T	.	.	.	GT	1/1	0/0
chr1	300	.	G	A	.	.	.	GT	./.	1/1
";
    let acc = make_acc_from_vcf(vcf);

    let mut csv = Cursor::new(Vec::new());
    acc.write_genotypes_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();

    assert!(output.contains("0,0,1"));
    assert!(output.contains("1,2,0"));
    assert!(output.contains("2,.,2"));
}

#[test]
fn computes_across_zero_sites() {
    let acc = SketchAccumulator::new(vec!["S1".to_string(), "S2".to_string()], IndexMap::new());
    let (_n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 0);
    assert!(distance.is_nan());
}

#[test]
fn distance_nan_when_no_common() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr1	100	.	A	G	.	.	.	GT	./.	0/1
chr1	200	.	C	T	.	.	.	GT	0/1	./.
";
    let acc = make_acc_from_vcf(vcf);
    let (_n_diff, n_common, distance) = acc.compute_pair_stats(0, 1);
    assert_eq!(n_common, 0);
    assert!(distance.is_nan());
}

#[test]
fn three_sample_gt_compare() {
    let vcf = "\
##fileformat=VCFv4.2
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	A	B	C
chr1	100	.	X	Y	.	.	.	GT	0/0	0/1	1/1
";
    let acc = make_acc_from_vcf(vcf);

    let (_diff, common, dist) = acc.compute_pair_stats(0, 1);
    assert_eq!(common, 1);
    assert!((dist - 0.5).abs() < 0.001); // aa vs het

    let (_diff, common, dist) = acc.compute_pair_stats(0, 2);
    assert_eq!(common, 1);
    assert!((dist - 1.0).abs() < 0.001); // aa vs bb

    let (_diff, common, dist) = acc.compute_pair_stats(1, 2);
    assert_eq!(common, 1);
    assert!((dist - 0.5).abs() < 0.001); // het vs bb
}

#[test]
fn genotypes_output_is_sorted_by_contig_then_pos() {
    let vcf = "\
##fileformat=VCFv4.2
##contig=<ID=chr1,length=5000>
##contig=<ID=chr2,length=5000>
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">
#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO	FORMAT	S1	S2
chr2	300	.	A	G	.	.	.	GT	1/1	0/0
chr1	200	.	C	T	.	.	.	GT	0/1	1/1
chr2	100	.	A	G	.	.	.	GT	0/0	0/1
chr1	100	.	A	G	.	.	.	GT	0/0	1/1
chr2	200	.	C	T	.	.	.	GT	1/1	0/1
chr1	300	.	T	C	.	.	.	GT	1/1	0/0
";
    let acc = make_acc_from_vcf(vcf);

    let mut csv = Cursor::new(Vec::new());
    acc.write_genotypes_csv(&mut csv).unwrap();
    let output = String::from_utf8(csv.into_inner()).unwrap();

    let data_lines: Vec<&str> = output
        .lines()
        .skip(1) // header
        .collect();

    // Extract the genotype columns (after "site_idx") from each row
    let rows: Vec<Vec<&str>> = data_lines
        .iter()
        .map(|line| {
            let cols: Vec<&str> = line.split(',').collect();
            cols[1..].to_vec() // skip site index
        })
        .collect();

    assert_eq!(rows.len(), 6);

    // Expected order: chr1:100, chr1:200, chr1:300, chr2:100, chr2:200, chr2:300
    // (chr1 comes before chr2 in header contig order)
    let expected = vec![
        vec!["0", "2"],   // chr1:100: 0/0, 1/1
        vec!["1", "2"],   // chr1:200: 0/1, 1/1
        vec!["2", "0"],   // chr1:300: 1/1, 0/0
        vec!["0", "1"],   // chr2:100: 0/0, 0/1
        vec!["2", "1"],   // chr2:200: 1/1, 0/1
        vec!["2", "0"],   // chr2:300: 1/1, 0/0
    ];

    assert_eq!(rows, expected);
}
