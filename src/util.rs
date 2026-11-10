//! Per-record utility callbacks suitable for `ParallelVariantWindowProcessor`.

use rust_htslib::bcf::Record;

/// Force REF and ALT alleles to ASCII-uppercase. Workaround for GLnexus, which
/// occasionally emits lowercase allele characters that downstream tools
/// reject. Always keeps the record.
pub fn uppercase_alleles(rec: &mut Record) -> bool {
    let alleles: Vec<Vec<u8>> = rec.alleles().iter().map(|a| a.to_ascii_uppercase()).collect();
    if alleles.iter().zip(rec.alleles()).any(|(u, a)| u.as_slice() != a) {
        let refs: Vec<&[u8]> = alleles.iter().map(Vec::as_slice).collect();
        rec.set_alleles(&refs).expect("setting alleles of an unpacked record");
    }
    true
}

pub fn human_bp(bp: u64) -> (f64, &'static str) {
    const KB: f64 = 1_000.0;
    const MB: f64 = 1_000_000.0;
    const GB: f64 = 1_000_000_000.0;
    if bp as f64 >= GB {
        (bp as f64 / GB, "Gbp")
    } else if bp as f64 >= MB {
        (bp as f64 / MB, "Mbp")
    } else if bp as f64 >= KB {
        (bp as f64 / KB, "kbp")
    } else {
        (bp as f64, "bp")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::hts::VcfReader;

    /// Parse one VCF record line (single sample-less site) via htslib.
    fn record(reference: &str, alts: &[&str]) -> (VcfReader, Record) {
        let dir = tempfile::tempdir().unwrap();
        let path = dir.path().join("t.vcf");
        let alt = if alts.is_empty() { ".".to_string() } else { alts.join(",") };
        std::fs::write(
            &path,
            format!("##fileformat=VCFv4.2\n##contig=<ID=chr1>\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\nchr1\t100\t.\t{reference}\t{alt}\t.\t.\t.\n"),
        )
        .unwrap();
        let mut r = VcfReader::open(path.to_str().unwrap()).unwrap();
        let mut rec = r.empty_record();
        assert!(r.read(&mut rec).unwrap());
        (r, rec)
    }

    fn alleles(rec: &Record) -> Vec<String> {
        rec.alleles().iter().map(|a| String::from_utf8_lossy(a).into_owned()).collect()
    }

    #[test]
    fn uppercases_lowercase_ref_and_alt() {
        let (_r, mut rec) = record("a", &["t"]);
        assert!(uppercase_alleles(&mut rec));
        assert_eq!(alleles(&rec), vec!["A", "T"]);
    }

    #[test]
    fn leaves_already_uppercase_unchanged() {
        let (_r, mut rec) = record("ACGT", &["G", "TT"]);
        uppercase_alleles(&mut rec);
        assert_eq!(alleles(&rec), vec!["ACGT", "G", "TT"]);
    }

    #[test]
    fn handles_multi_alt() {
        let (_r, mut rec) = record("a", &["c", "g", "t"]);
        uppercase_alleles(&mut rec);
        assert_eq!(alleles(&rec), vec!["A", "C", "G", "T"]);
    }

    #[test]
    fn uppercases_iupac_ambiguity_code() {
        let (_r, mut rec) = record("n", &["a"]);
        uppercase_alleles(&mut rec);
        assert_eq!(alleles(&rec), vec!["N", "A"]);
    }
}
