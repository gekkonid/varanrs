//! Per-record utility callbacks suitable for `ParallelVariantWindowProcessor`.

use noodles_vcf::variant::RecordBuf;

/// Force REF and ALT alleles to ASCII-uppercase. Workaround for GLnexus, which
/// occasionally emits lowercase allele characters that downstream tools reject.
pub fn uppercase_alleles(mut buf: RecordBuf) -> Option<RecordBuf> {
    buf.reference_bases_mut().make_ascii_uppercase();
    for alt in buf.alternate_bases_mut().as_mut().iter_mut() {
        alt.make_ascii_uppercase();
    }
    Some(buf)
}

#[cfg(test)]
mod tests {
    use super::*;
    use noodles_core::Position;

    fn record(reference: &str, alts: &[&str]) -> RecordBuf {
        let mut buf = RecordBuf::default();
        *buf.reference_sequence_name_mut() = String::from("chr1");
        *buf.variant_start_mut() = Some(Position::try_from(100).unwrap());
        *buf.reference_bases_mut() = reference.to_string();
        *buf.alternate_bases_mut().as_mut() = alts.iter().map(|a| a.to_string()).collect();
        buf
    }

    #[test]
    fn uppercases_lowercase_ref_and_alt() {
        let buf = uppercase_alleles(record("a", &["t"])).unwrap();
        assert_eq!(buf.reference_bases(), "A");
        assert_eq!(buf.alternate_bases().as_ref(), &["T".to_string()]);
    }

    #[test]
    fn leaves_already_uppercase_unchanged() {
        let buf = uppercase_alleles(record("ACGT", &["G", "TT"])).unwrap();
        assert_eq!(buf.reference_bases(), "ACGT");
        assert_eq!(
            buf.alternate_bases().as_ref(),
            &["G".to_string(), "TT".to_string()]
        );
    }

    #[test]
    fn handles_multi_alt() {
        let buf = uppercase_alleles(record("a", &["c", "g", "t"])).unwrap();
        assert_eq!(buf.reference_bases(), "A");
        assert_eq!(
            buf.alternate_bases().as_ref(),
            &["C".to_string(), "G".to_string(), "T".to_string()]
        );
    }

    #[test]
    fn uppercases_iupac_ambiguity_code() {
        let buf = uppercase_alleles(record("n", &["a"])).unwrap();
        assert_eq!(buf.reference_bases(), "N");
        assert_eq!(buf.alternate_bases().as_ref(), &["A".to_string()]);
    }
}
