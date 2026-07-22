//! Per-record utility callbacks suitable for `ParallelVariantWindowProcessor`.

use noodles_vcf as vcf;
use vcf::variant::RecordBuf;

/// Force REF and ALT alleles to ASCII-uppercase. Workaround for GLnexus, which
/// occasionally emits lowercase allele characters that downstream tools reject.
pub fn uppercase_alleles(mut buf: RecordBuf) -> Option<RecordBuf> {
    buf.reference_bases_mut().make_ascii_uppercase();
    for alt in buf.alternate_bases_mut().as_mut().iter_mut() {
        alt.make_ascii_uppercase();
    }
    Some(buf)
}

/// Work around a GLnexus header bug: it declares `RNC` as `Type=Character,
/// Number=2` but actually encodes it as a single string like `"II"`.  Patching
/// the header to `Type=String, Number=1` lets noodles parse the records.
pub fn normalize_header_for_noodles(header: &mut vcf::Header) {
    use vcf::header::record::value::{Map, map::Format};
    use vcf::header::record::value::map::format::{Number, Type};

    if header.formats().contains_key("RNC") {
        header.formats_mut().insert(
            String::from("RNC"),
            Map::<Format>::new(Number::Count(1), Type::String, "Reason for No Call"),
        );
    }
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

    #[test]
    fn normalize_header_fixes_rnc_type() {
        let mut header = vcf::Header::default();
        use vcf::header::record::value::{Map, map::Format};
        use vcf::header::record::value::map::format::{Number, Type};
        header.formats_mut().insert(
            String::from("RNC"),
            Map::<Format>::new(Number::Count(2), Type::Character, "Reason for No Call"),
        );

        normalize_header_for_noodles(&mut header);

        let rnc = header.formats().get("RNC").unwrap();
        assert_eq!(rnc.number(), Number::Count(1));
        assert_eq!(rnc.ty(), Type::String);
    }
}
