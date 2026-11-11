//! Thin wrappers over htslib (via rust-htslib's bindings) for the pieces of
//! VCF/BCF I/O varanrs needs, plus the vendored bcftools expression engine.
//!
//! Records are rust-htslib [`Record`]s, so parsing and serialisation are
//! exactly htslib's (`vcf_parse` / `vcf_format`), i.e. byte-identical to
//! `bcftools view`.

use std::ffi::{CStr, CString};
use std::os::raw::{c_char, c_int, c_void};
use std::path::Path;
use std::ptr;
use std::sync::Arc;

use anyhow::{Context, Result, anyhow, bail};
use rust_htslib::bcf::Record;
pub use rust_htslib::bcf::header::HeaderView;
use rust_htslib::htslib as hts;

// ---------------------------------------------------------------------------
// Missing-value conventions (htslib vcf.h macros)
// ---------------------------------------------------------------------------

pub const INT32_MISSING: i32 = i32::MIN;
pub const INT32_VECTOR_END: i32 = i32::MIN + 1;
const FLOAT_MISSING_BITS: u32 = 0x7F80_0001;
const FLOAT_VECTOR_END_BITS: u32 = 0x7F80_0002;

pub fn float_missing() -> f32 {
    f32::from_bits(FLOAT_MISSING_BITS)
}
pub fn float_vector_end() -> f32 {
    f32::from_bits(FLOAT_VECTOR_END_BITS)
}
pub fn is_float_missing(f: f32) -> bool {
    f.to_bits() == FLOAT_MISSING_BITS
}
pub fn is_float_vector_end(f: f32) -> bool {
    f.to_bits() == FLOAT_VECTOR_END_BITS
}

/// `bcf_gt_is_missing`
pub fn gt_is_missing(v: i32) -> bool {
    (v >> 1) == 0
}
/// `bcf_gt_allele`
pub fn gt_allele(v: i32) -> i32 {
    (v >> 1) - 1
}
/// `bcf_gt_missing` (unphased missing allele)
pub const GT_MISSING: i32 = 0;

// ---------------------------------------------------------------------------
// kstring
// ---------------------------------------------------------------------------

/// Owned htslib `kstring_t`.
pub struct KString(hts::kstring_t);

impl KString {
    pub fn new() -> Self {
        Self(hts::kstring_t { l: 0, m: 0, s: ptr::null_mut() })
    }
    pub fn clear(&mut self) {
        self.0.l = 0;
    }
    pub fn as_bytes(&self) -> &[u8] {
        if self.0.s.is_null() {
            &[]
        } else {
            unsafe { std::slice::from_raw_parts(self.0.s as *const u8, self.0.l) }
        }
    }
    pub fn as_mut_ptr(&mut self) -> *mut hts::kstring_t {
        &mut self.0
    }
}

impl Default for KString {
    fn default() -> Self {
        Self::new()
    }
}

impl Drop for KString {
    fn drop(&mut self) {
        if !self.0.s.is_null() {
            unsafe { libc_free(self.0.s as *mut c_void) };
        }
    }
}

unsafe extern "C" {
    #[link_name = "free"]
    fn libc_free(p: *mut c_void);
}

// SAFETY: a KString is a plain heap buffer owned by one value.
unsafe impl Send for KString {}

// ---------------------------------------------------------------------------
// Headers
// ---------------------------------------------------------------------------

fn hdr_ptr(h: &HeaderView) -> *mut hts::bcf_hdr_t {
    unsafe { h.as_ptr() }
}

/// Duplicate `base` and append the given `##` header lines (lines whose ID
/// already exists are left as they are, as with `bcf_hdr_append`).
pub fn header_with_lines(base: &HeaderView, lines: &[String]) -> Result<Arc<HeaderView>> {
    unsafe {
        let h = hts::bcf_hdr_dup(hdr_ptr(base));
        if h.is_null() {
            bail!("failed to duplicate VCF header");
        }
        for line in lines {
            let c = CString::new(line.as_str())?;
            if hts::bcf_hdr_append(h, c.as_ptr()) < 0 {
                hts::bcf_hdr_destroy(h);
                bail!("failed to append header line: {line}");
            }
        }
        if hts::bcf_hdr_sync(h) < 0 {
            hts::bcf_hdr_destroy(h);
            bail!("failed to sync VCF header");
        }
        Ok(Arc::new(HeaderView::from_ptr(h)))
    }
}

/// The VCF text of a header (as `bcf_hdr_write` would emit, without any
/// program/version lines).
pub fn format_header(h: &HeaderView) -> Result<Vec<u8>> {
    let mut ks = KString::new();
    if unsafe { hts::bcf_hdr_format(hdr_ptr(h), 0, ks.as_mut_ptr()) } < 0 {
        bail!("failed to format VCF header");
    }
    Ok(ks.as_bytes().to_vec())
}

/// Whether the header defines an INFO (`hl = BCF_HL_INFO`), FORMAT or FILTER
/// line with this ID.
pub fn header_has(h: &HeaderView, hl: u32, id: &str) -> bool {
    let Ok(c) = CString::new(id) else { return false };
    unsafe {
        let i = hts::bcf_hdr_id2int(hdr_ptr(h), hts::BCF_DT_ID as c_int, c.as_ptr());
        if i < 0 {
            return false;
        }
        let idinfo = (*hdr_ptr(h)).id[hts::BCF_DT_ID as usize].add(i as usize);
        let hrec = (*(*idinfo).val).hrec[hl as usize];
        !hrec.is_null()
    }
}

pub fn n_samples(h: &HeaderView) -> usize {
    h.sample_count() as usize
}

/// Contig names in header order, with declared lengths (`None` if absent).
pub fn header_contigs(h: &HeaderView) -> Vec<(String, Option<u64>)> {
    let mut out = Vec::new();
    for rec in h.header_records() {
        if let rust_htslib::bcf::header::HeaderRecord::Contig { values, .. } = rec
            && let Some(id) = values.get("ID")
        {
            let len = values.get("length").and_then(|l| l.parse::<u64>().ok());
            out.push((id.clone(), len));
        }
    }
    out
}

// ---------------------------------------------------------------------------
// Records
// ---------------------------------------------------------------------------

pub fn rec_ptr(r: &mut Record) -> *mut hts::bcf1_t {
    r.inner_mut() as *mut hts::bcf1_t
}

/// Append the VCF text line of `rec` (including the trailing newline).
pub fn format_record(rec: &Record, ks: &mut KString, out: &mut Vec<u8>) -> Result<()> {
    ks.clear();
    let r = unsafe { hts::vcf_format(hdr_ptr(rec.header()), rec.inner() as *const _, ks.as_mut_ptr()) };
    if r < 0 {
        bail!("failed to format VCF record");
    }
    out.extend_from_slice(ks.as_bytes());
    Ok(())
}

/// Contig name of a record.
pub fn rec_contig(rec: &Record) -> String {
    rec.rid()
        .and_then(|rid| rec.header().rid2name(rid).ok())
        .map(|n| String::from_utf8_lossy(n).into_owned())
        .unwrap_or_default()
}

/// `chrom:pos` for messages.
pub fn rec_locus(rec: &Record) -> String {
    format!("{}:{}", rec_contig(rec), rec.pos() + 1)
}

// ---------------------------------------------------------------------------
// Typed INFO/FORMAT access (bcf_get_*_values / bcf_update_*)
// ---------------------------------------------------------------------------

fn rec_raw(rec: &Record) -> *mut hts::bcf1_t {
    rec.inner() as *const hts::bcf1_t as *mut hts::bcf1_t
}

/// `bcf_get_info_values` / `bcf_get_format_values` into a Vec. Returns the
/// raw htslib return code on failure (-1: no header definition, -2: type
/// clash, -3: tag absent from the record).
unsafe fn get_values<T: Copy>(rec: &Record, tag: &str, ty: u32, info: bool) -> std::result::Result<Vec<T>, i32> {
    let c = CString::new(tag).map_err(|_| -1)?;
    let mut dst: *mut c_void = ptr::null_mut();
    let mut ndst: c_int = 0;
    let n = unsafe {
        if info {
            hts::bcf_get_info_values(hdr_ptr(rec.header()), rec_raw(rec), c.as_ptr(), &mut dst, &mut ndst, ty as c_int)
        } else {
            hts::bcf_get_format_values(hdr_ptr(rec.header()), rec_raw(rec), c.as_ptr(), &mut dst, &mut ndst, ty as c_int)
        }
    };
    let out = if n < 0 {
        Err(n)
    } else if ty == hts::BCF_HT_FLAG {
        Ok(Vec::new())
    } else {
        Ok(unsafe { std::slice::from_raw_parts(dst as *const T, n as usize) }.to_vec())
    };
    if !dst.is_null() {
        unsafe { libc_free(dst) };
    }
    out
}

/// INFO integer values (htslib missing/vector_end sentinels preserved), or
/// the htslib error code (-3 when the tag is absent).
pub fn info_i32(rec: &Record, tag: &str) -> std::result::Result<Vec<i32>, i32> {
    unsafe { get_values::<i32>(rec, tag, hts::BCF_HT_INT, true) }
}
pub fn info_f32(rec: &Record, tag: &str) -> std::result::Result<Vec<f32>, i32> {
    unsafe { get_values::<f32>(rec, tag, hts::BCF_HT_REAL, true) }
}
/// FORMAT integer values for all samples (n_samples × width, padded with
/// vector_end), as `bcf_get_format_int32`.
pub fn format_i32(rec: &Record, tag: &str) -> std::result::Result<Vec<i32>, i32> {
    unsafe { get_values::<i32>(rec, tag, hts::BCF_HT_INT, false) }
}

unsafe fn update(rec: &mut Record, tag: &str, data: *const c_void, n: usize, ty: u32, info: bool) -> Result<()> {
    let c = CString::new(tag)?;
    let h = hdr_ptr(rec.header());
    let r = unsafe {
        if info {
            hts::bcf_update_info(h, rec_ptr(rec), c.as_ptr(), data, n as c_int, ty as c_int)
        } else {
            hts::bcf_update_format(h, rec_ptr(rec), c.as_ptr(), data, n as c_int, ty as c_int)
        }
    };
    if r != 0 {
        bail!("failed to update {} {tag} at {}", if info { "INFO" } else { "FORMAT" }, rec_locus(rec));
    }
    Ok(())
}

/// `bcf_update_info_int32`; an empty slice removes the tag.
pub fn set_info_i32(rec: &mut Record, tag: &str, v: &[i32]) -> Result<()> {
    unsafe { update(rec, tag, v.as_ptr() as *const c_void, v.len(), hts::BCF_HT_INT, true) }
}
/// `bcf_update_info_float`; an empty slice removes the tag.
pub fn set_info_f32(rec: &mut Record, tag: &str, v: &[f32]) -> Result<()> {
    unsafe { update(rec, tag, v.as_ptr() as *const c_void, v.len(), hts::BCF_HT_REAL, true) }
}
/// `bcf_update_format_float`.
pub fn set_format_f32(rec: &mut Record, tag: &str, v: &[f32]) -> Result<()> {
    unsafe { update(rec, tag, v.as_ptr() as *const c_void, v.len(), hts::BCF_HT_REAL, false) }
}
/// `bcf_update_format_int32` (also used for GT, encoded as htslib does).
pub fn set_format_i32(rec: &mut Record, tag: &str, v: &[i32]) -> Result<()> {
    unsafe { update(rec, tag, v.as_ptr() as *const c_void, v.len(), hts::BCF_HT_INT, false) }
}

/// INFO definitions in the header: (ID, Number, Type) as written.
pub fn header_infos(h: &HeaderView) -> Vec<(String, String, String)> {
    let mut out = Vec::new();
    for rec in h.header_records() {
        if let rust_htslib::bcf::header::HeaderRecord::Info { values, .. } = rec {
            let get = |k: &str| values.get(k).cloned().unwrap_or_default();
            out.push((get("ID"), get("Number"), get("Type")));
        }
    }
    out
}

/// htslib's FILTER id for a name (`PASS` is always 0), if defined.
pub fn filter_id(h: &HeaderView, name: &str) -> Option<i32> {
    let c = CString::new(name).ok()?;
    let id = unsafe { hts::bcf_hdr_id2int(hdr_ptr(h), hts::BCF_DT_ID as c_int, c.as_ptr()) };
    (id >= 0 && header_has(h, hts::BCF_HL_FLT, name)).then_some(id)
}

/// `bcf_add_filter`
pub fn add_filter(rec: &mut Record, id: i32) {
    let h = hdr_ptr(rec.header());
    unsafe { hts::bcf_add_filter(h, rec_ptr(rec), id) };
}
/// `bcf_update_filter` with a single id.
pub fn set_filter(rec: &mut Record, id: i32) {
    let h = hdr_ptr(rec.header());
    let mut id = id;
    unsafe { hts::bcf_update_filter(h, rec_ptr(rec), &mut id, 1) };
}
/// Number of FILTER values on the record (0 = ".").
pub fn n_filters(rec: &mut Record) -> usize {
    unsafe {
        hts::bcf_unpack(rec_ptr(rec), hts::BCF_UN_FLT as c_int);
        (*rec_ptr(rec)).d.n_flt as usize
    }
}

// ---------------------------------------------------------------------------
// Readers
// ---------------------------------------------------------------------------

/// Subset samples as `bcftools view -S FILE` does: keep the samples listed
/// in `file` (one per line), in that order, or with a leading `^` drop them
/// (keeping header order). INFO/AC and INFO/AN are recomputed from the kept
/// genotypes, as `bcftools view` does by default.
#[derive(Clone, Debug)]
pub struct SampleSubset {
    pub file: String,
}

/// Per-reader sample subsetting state.
struct Subset {
    /// The full input header, used for parsing. The reader's public header is
    /// the subset header; every dictionary ID is shared between the two
    /// (htslib's IDX keys), and lines added to `full` are mirrored into it.
    full: Arc<HeaderView>,
    /// Number of `full` header records already mirrored.
    synced: c_int,
    /// For each output sample, its index in the input (`bcf_subset`).
    imap: Vec<c_int>,
    ac: Vec<c_int>,
}

/// A VCF/BCF reader (any compression), optionally with its index.
pub struct VcfReader {
    fp: *mut hts::htsFile,
    header: Arc<HeaderView>,
    subset: Option<Subset>,
    is_bcf: bool,
    tbx: *mut hts::tbx_t,
    idx: *mut hts::hts_idx_t,
    line: KString,
}

// SAFETY: the reader owns its htslib handles and is used by one thread at a
// time.
unsafe impl Send for VcfReader {}

impl VcfReader {
    /// Open a file, or stdin for "-".
    pub fn open(path: &str) -> Result<Self> {
        let c = CString::new(path)?;
        let mode = CString::new("r")?;
        unsafe {
            let fp = hts::hts_open(c.as_ptr(), mode.as_ptr());
            if fp.is_null() {
                bail!("failed to open {path}");
            }
            let fmt = hts::hts_get_format(fp);
            let exact = (*fmt).format;
            if exact != hts::htsExactFormat_vcf && exact != hts::htsExactFormat_bcf {
                hts::hts_close(fp);
                bail!("{path}: not a VCF or BCF file");
            }
            let h = hts::bcf_hdr_read(fp);
            if h.is_null() {
                hts::hts_close(fp);
                bail!("{path}: failed to read the VCF header");
            }
            Ok(Self {
                fp,
                header: Arc::new(HeaderView::from_ptr(h)),
                subset: None,
                is_bcf: exact == hts::htsExactFormat_bcf,
                tbx: ptr::null_mut(),
                idx: ptr::null_mut(),
                line: KString::new(),
            })
        }
    }

    /// Open a bgzipped VCF (tabix index) or BCF (CSI index).
    pub fn open_indexed(path: &Path) -> Result<Self> {
        let p = path.to_str().ok_or_else(|| anyhow!("non-UTF-8 path"))?;
        let mut r = Self::open(p)?;
        let c = CString::new(p)?;
        unsafe {
            if r.is_bcf {
                r.idx = hts::hts_idx_load(c.as_ptr(), hts::HTS_FMT_CSI as c_int);
                if r.idx.is_null() {
                    bail!("{p}: failed to load the CSI index (run `bcftools index`)");
                }
            } else {
                r.tbx = hts::tbx_index_load(c.as_ptr());
                if r.tbx.is_null() {
                    bail!("{p}: failed to load the index (run `tabix -p vcf` or `bcftools index -t`)");
                }
                r.idx = (*r.tbx).idx;
            }
        }
        Ok(r)
    }

    pub fn header(&self) -> &Arc<HeaderView> {
        &self.header
    }

    /// Append `##` header lines to this reader's header in place (lines whose
    /// ID already exists are left as they are). Must be called before any
    /// records are read. Records keep using this single header, so tags that
    /// htslib adds while parsing (undeclared in the file) and the appended
    /// ones share one consistent dictionary.
    pub fn append_header_lines(&mut self, lines: &[String]) -> Result<()> {
        if lines.is_empty() {
            return Ok(());
        }
        let h = self.parse_hdr();
        unsafe {
            for line in lines {
                let c = CString::new(line.as_str())?;
                if hts::bcf_hdr_append(h, c.as_ptr()) < 0 {
                    bail!("failed to append header line: {line}");
                }
            }
            if hts::bcf_hdr_sync(h) < 0 {
                bail!("failed to sync VCF header");
            }
        }
        self.sync_subset_header()
    }

    /// The header records are parsed with (the full input header).
    fn parse_hdr(&self) -> *mut hts::bcf_hdr_t {
        hdr_ptr(self.subset.as_ref().map_or(&self.header, |s| &s.full))
    }

    /// Mirror header records added to the full header (appended lines, or
    /// tags htslib declared while parsing) into the subset header.
    fn sync_subset_header(&mut self) -> Result<()> {
        let Some(sub) = self.subset.as_mut() else { return Ok(()) };
        let full = hdr_ptr(&sub.full);
        unsafe {
            let n = (*full).nhrec;
            if n == sub.synced {
                return Ok(());
            }
            let h = hdr_ptr(&self.header);
            for i in sub.synced..n {
                let hrec = hts::bcf_hrec_dup(*(*full).hrec.add(i as usize));
                if hrec.is_null() || hts::bcf_hdr_add_hrec(h, hrec) < 0 {
                    bail!("failed to update the subset VCF header");
                }
            }
            if hts::bcf_hdr_sync(h) < 0 {
                bail!("failed to sync VCF header");
            }
            sub.synced = n;
        }
        Ok(())
    }

    /// Restrict output to a subset of samples (see [`SampleSubset`]). Must be
    /// called before any header lines are appended or records are read.
    /// Returns warnings for the caller to report.
    pub fn subset_samples(&mut self, spec: &SampleSubset) -> Result<Vec<String>> {
        if self.subset.is_some() {
            bail!("samples already subset");
        }
        let (exclude, file) = match spec.file.strip_prefix('^') {
            Some(f) => (true, f),
            None => (false, spec.file.as_str()),
        };
        let listed = read_list_file(file)?;
        let in_header: Vec<String> = self
            .header
            .samples()
            .iter()
            .map(|s| String::from_utf8_lossy(s).into_owned())
            .collect();
        let known: std::collections::HashSet<&str> = in_header.iter().map(String::as_str).collect();
        for s in &listed {
            if !known.contains(s.as_str()) {
                bail!(
                    "{} called for sample that does not exist in header: \"{s}\"",
                    if exclude { "exclude" } else { "subset" }
                );
            }
        }
        let samples: Vec<String> = if exclude {
            let drop: std::collections::HashSet<&str> = listed.iter().map(String::as_str).collect();
            in_header.into_iter().filter(|s| !drop.contains(s.as_str())).collect()
        } else {
            listed
        };
        let mut warnings = Vec::new();
        if samples.is_empty() {
            warnings.push("subsetting has removed all samples".to_string());
        }

        // As bcftools view: declare AC/AN (no-op when present), then subset.
        let full = Arc::clone(&self.header);
        let fp = hdr_ptr(&full);
        let c_samples: Vec<CString> = samples.iter().map(|s| CString::new(s.as_str())).collect::<Result<_, _>>()?;
        let mut ptrs: Vec<*mut c_char> = c_samples.iter().map(|c| c.as_ptr() as *mut c_char).collect();
        let mut imap: Vec<c_int> = vec![0; samples.len()];
        unsafe {
            for line in [
                "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes\">",
                "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">",
            ] {
                let c = CString::new(line)?;
                if hts::bcf_hdr_append(fp, c.as_ptr()) < 0 {
                    bail!("failed to append header line: {line}");
                }
            }
            if hts::bcf_hdr_sync(fp) < 0 {
                bail!("failed to sync VCF header");
            }
            let sub = hts::bcf_hdr_subset(fp, samples.len() as c_int, ptrs.as_mut_ptr(), imap.as_mut_ptr());
            if sub.is_null() {
                bail!("error occurred while subsetting samples (duplicate sample names in {file}?)");
            }
            if samples.is_empty() {
                // Sites-only output, as bcftools view.
                hts::bcf_hdr_remove(sub, hts::BCF_HL_FMT as c_int, ptr::null());
                if hts::bcf_hdr_sync(sub) < 0 {
                    hts::bcf_hdr_destroy(sub);
                    bail!("failed to sync VCF header");
                }
            }
            self.header = Arc::new(HeaderView::from_ptr(sub));
            self.subset = Some(Subset { synced: (*fp).nhrec, full, imap, ac: Vec::new() });
        }
        Ok(warnings)
    }

    /// Post-parse step: subset samples (updating AC/AN) and unpack.
    fn finish_record(&mut self, rec: &mut Record) -> Result<()> {
        self.sync_subset_header()?;
        let r = rec_ptr(rec);
        if let Some(sub) = self.subset.as_mut() {
            let full = hdr_ptr(&sub.full);
            let h = hdr_ptr(&self.header);
            unsafe {
                let n_allele = (*r).n_allele() as usize;
                sub.ac.clear();
                sub.ac.resize(n_allele.max(1), 0);
                let n = sub.imap.len();
                // AC/AN from INFO if present, otherwise from all genotypes.
                hts::bcf_calc_ac(full, r, sub.ac.as_mut_ptr(), (hts::BCF_UN_INFO | hts::BCF_UN_FMT) as c_int);
                let mut update = true;
                if n == 0 {
                    hts::bcf_subset(full, r, 0, ptr::null_mut());
                } else {
                    hts::bcf_subset(full, r, n as c_int, sub.imap.as_mut_ptr());
                    let gt = CString::new("GT").unwrap();
                    if hts::bcf_get_fmt(h, r, gt.as_ptr()).is_null() {
                        update = false;
                    } else {
                        hts::bcf_calc_ac(h, r, sub.ac.as_mut_ptr(), hts::BCF_UN_FMT as c_int);
                    }
                }
                if update {
                    let ac: Vec<i32> = sub.ac[..n_allele].to_vec();
                    let an: i32 = ac.iter().sum();
                    set_info_i32(rec, "AC", if ac.len() > 1 { &ac[1..] } else { &[] })?;
                    set_info_i32(rec, "AN", &[an])?;
                }
            }
        }
        unsafe { hts::bcf_unpack(rec_ptr(rec), hts::BCF_UN_ALL as c_int) };
        Ok(())
    }

    pub fn empty_record(&self) -> Record {
        Record::new(Arc::clone(&self.header))
    }

    /// Read the next record sequentially. `Ok(false)` at end of input.
    pub fn read(&mut self, rec: &mut Record) -> Result<bool> {
        let r = unsafe { hts::bcf_read(self.fp, self.parse_hdr(), rec_ptr(rec)) };
        match r {
            0 => {
                self.finish_record(rec)?;
                Ok(true)
            }
            -1 => Ok(false),
            _ => bail!("failed to read/parse a VCF record (error {r})"),
        }
    }

    /// Index id of a contig, if the index knows it.
    pub fn tid(&self, name: &str) -> Option<i32> {
        let c = CString::new(name).ok()?;
        let id = unsafe {
            if !self.tbx.is_null() {
                hts::tbx_name2id(self.tbx, c.as_ptr())
            } else {
                hts::bcf_hdr_id2int(hdr_ptr(&self.header), hts::BCF_DT_CTG as c_int, c.as_ptr())
            }
        };
        (id >= 0).then_some(id)
    }

    /// Contig names known to a tabix index, in index (file) order.
    pub fn tabix_seqnames(&self) -> Option<Vec<String>> {
        if self.tbx.is_null() {
            return None;
        }
        let mut n: c_int = 0;
        unsafe {
            let names = hts::tbx_seqnames(self.tbx, &mut n);
            if names.is_null() {
                return Some(Vec::new());
            }
            let v = (0..n as usize)
                .map(|i| CStr::from_ptr(*names.add(i)).to_string_lossy().into_owned())
                .collect();
            libc_free(names as *mut c_void);
            Some(v)
        }
    }

    /// Index iterator over records overlapping the 0-based half-open
    /// interval `[beg, end)` of contig `tid`.
    pub fn query(&self, tid: i32, beg: i64, end: i64) -> Option<Itr> {
        if self.idx.is_null() {
            return None;
        }
        let readrec: hts::hts_readrec_func = if self.is_bcf { Some(hts::bcf_readrec) } else { Some(hts::tbx_readrec) };
        let it = unsafe { hts::hts_itr_query(self.idx, tid, beg, end, readrec) };
        (!it.is_null()).then_some(Itr(it))
    }

    /// Next record from an iterator. `Ok(false)` when exhausted.
    pub fn itr_next(&mut self, itr: &mut Itr, rec: &mut Record) -> Result<bool> {
        unsafe {
            let bgzf = hts::hts_get_bgzfp(self.fp);
            let r = if self.is_bcf {
                hts::hts_itr_next(bgzf, itr.0, rec_ptr(rec) as *mut c_void, ptr::null_mut())
            } else {
                let r = hts::hts_itr_next(bgzf, itr.0, self.line.as_mut_ptr() as *mut c_void, self.tbx as *mut c_void);
                if r >= 0 && hts::vcf_parse(self.line.as_mut_ptr(), self.parse_hdr(), rec_ptr(rec)) < 0 {
                    bail!("failed to parse VCF record: {}", String::from_utf8_lossy(self.line.as_bytes()));
                }
                r
            };
            if r >= 0 {
                self.finish_record(rec)?;
                Ok(true)
            } else if r == -1 {
                Ok(false)
            } else {
                bail!("failed to read an indexed VCF/BCF record (error {r})")
            }
        }
    }
}

impl Drop for VcfReader {
    fn drop(&mut self) {
        unsafe {
            if !self.tbx.is_null() {
                hts::tbx_destroy(self.tbx);
            } else if !self.idx.is_null() {
                hts::hts_idx_destroy(self.idx);
            }
            hts::hts_close(self.fp);
        }
    }
}

/// Non-empty lines of a list file (as htslib's `hts_readlist`, trailing
/// whitespace and CR removed).
fn read_list_file(path: &str) -> Result<Vec<String>> {
    let text = std::fs::read_to_string(path).with_context(|| format!("could not read the sample list: {path}"))?;
    Ok(text
        .lines()
        .map(|l| l.trim_end_matches(['\r', '\n']))
        .filter(|l| !l.is_empty())
        .map(str::to_string)
        .collect())
}

/// Read every record of a file ("-" for stdin) into memory. Intended for
/// tests and small inputs.
pub fn read_all(path: &str) -> Result<(Arc<HeaderView>, Vec<Record>)> {
    let mut r = VcfReader::open(path)?;
    let mut out = Vec::new();
    loop {
        let mut rec = r.empty_record();
        if !r.read(&mut rec)? {
            break;
        }
        out.push(rec);
    }
    Ok((Arc::clone(r.header()), out))
}

/// An htslib index iterator.
pub struct Itr(*mut hts::hts_itr_t);

impl Itr {
    /// True when the index proves no record overlaps the queried interval.
    pub fn is_empty(&self) -> bool {
        unsafe { (*self.0).n_off == 0 }
    }
    /// Smallest file offset of the chunks to read (file order of the data).
    pub fn min_offset(&self) -> Option<u64> {
        unsafe {
            let n = (*self.0).n_off;
            if n <= 0 || (*self.0).off.is_null() {
                return None;
            }
            (0..n as usize).map(|i| (*(*self.0).off.add(i)).u).min()
        }
    }
}

impl Drop for Itr {
    fn drop(&mut self) {
        unsafe { hts::hts_itr_destroy(self.0) };
    }
}

// ---------------------------------------------------------------------------
// bcftools expression engine (csrc/bcftools/filter.c)
// ---------------------------------------------------------------------------

#[repr(C)]
struct FilterT {
    _private: [u8; 0],
}

unsafe extern "C" {
    fn filter_parse(hdr: *mut hts::bcf_hdr_t, str: *const c_char) -> *mut FilterT;
    fn filter_status(filter: *mut FilterT) -> c_int;
    fn filter_test(filter: *mut FilterT, rec: *mut hts::bcf1_t, samples: *mut *const u8) -> c_int;
    fn filter_destroy(filter: *mut FilterT);
    fn filter_get_doubles(filter: *mut FilterT, nval: *mut c_int, nval1: *mut c_int) -> *const f64;
}

/// A compiled bcftools filtering expression, bound to one header. Not shared
/// between threads: create one per worker with [`Filter::new`].
pub struct Filter {
    ptr: *mut FilterT,
    // keep the header alive for as long as the filter references it
    _header: Arc<HeaderView>,
}

// SAFETY: a Filter is used by one thread at a time; filter.c keeps all
// evaluation state inside the filter_t.
unsafe impl Send for Filter {}

impl Filter {
    /// Compile `expr` against `header`. Invalid expressions are reported by
    /// bcftools' own error handling (message on stderr, exit status 1).
    pub fn new(header: &Arc<HeaderView>, expr: &str) -> Result<Self> {
        let c = CString::new(expr).context("expression contains a NUL byte")?;
        let ptr = unsafe { filter_parse(hdr_ptr(header), c.as_ptr()) };
        if ptr.is_null() {
            bail!("could not parse the expression: {expr}");
        }
        let status = unsafe { filter_status(ptr) };
        if status != 0 {
            unsafe { filter_destroy(ptr) };
            bail!("could not parse the expression (unknown tag or syntax error): {expr}");
        }
        Ok(Self { ptr, _header: Arc::clone(header) })
    }

    /// True when the expression is true for the site (bcftools `filter_test`).
    pub fn test(&mut self, rec: &mut Record) -> bool {
        unsafe { filter_test(self.ptr, rec_ptr(rec), ptr::null_mut()) != 0 }
    }

    /// Evaluate per sample: the site result plus, when the expression uses
    /// FORMAT fields, whether it is true for each sample (bcftools
    /// `filter_test`'s `samples` output). `None` for site-level expressions.
    pub fn test_samples(&mut self, rec: &mut Record) -> (bool, Option<Vec<bool>>) {
        let n = n_samples(rec.header());
        let mut smpl: *const u8 = ptr::null();
        let site = unsafe { filter_test(self.ptr, rec_ptr(rec), &mut smpl) != 0 };
        let per = (!smpl.is_null()).then(|| unsafe { std::slice::from_raw_parts(smpl, n) }.iter().map(|&b| b != 0).collect());
        (site, per)
    }

    /// Evaluate and return the numeric result values (used for F_MISSING,
    /// exactly as `bcftools +fill-tags` computes it).
    pub fn values(&mut self, rec: &mut Record) -> Vec<f64> {
        unsafe {
            filter_test(self.ptr, rec_ptr(rec), ptr::null_mut());
            let (mut n, mut n1) = (0, 0);
            let p = filter_get_doubles(self.ptr, &mut n, &mut n1);
            if p.is_null() || n <= 0 {
                Vec::new()
            } else {
                std::slice::from_raw_parts(p, n as usize).to_vec()
            }
        }
    }
}

impl Drop for Filter {
    fn drop(&mut self) {
        unsafe { filter_destroy(self.ptr) };
    }
}

// ---------------------------------------------------------------------------
// BGZF compression and on-the-fly CSI indexing
// ---------------------------------------------------------------------------

/// Uncompressed bytes per BGZF block, as htslib (`BGZF_BLOCK_SIZE`).
pub const BGZF_BLOCK_SIZE: usize = 0xff00;
const BGZF_MAX_BLOCK_SIZE: usize = 0x10000;
/// The standard empty BGZF block marking end of file.
pub const BGZF_EOF: [u8; 28] = [
    0x1f, 0x8b, 0x08, 0x04, 0, 0, 0, 0, 0, 0xff, 0x06, 0, 0x42, 0x43, 0x02, 0, 0x1b, 0, 0x03, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

/// Compress `text` into BGZF blocks of [`BGZF_BLOCK_SIZE`] uncompressed
/// bytes (htslib's `bgzf_compress`, default level), appending to `out`.
/// Returns each block's start offset relative to the start of `out`'s new
/// data.
pub fn bgzf_compress_blocks(text: &[u8], out: &mut Vec<u8>) -> Result<Vec<u64>> {
    let base = out.len();
    let mut starts = Vec::with_capacity(text.len().div_ceil(BGZF_BLOCK_SIZE));
    for chunk in text.chunks(BGZF_BLOCK_SIZE) {
        let at = out.len();
        starts.push((at - base) as u64);
        out.resize(at + BGZF_MAX_BLOCK_SIZE, 0);
        let mut dlen = BGZF_MAX_BLOCK_SIZE;
        let r = unsafe {
            hts::bgzf_compress(
                out.as_mut_ptr().add(at) as *mut c_void,
                &mut dlen,
                chunk.as_ptr() as *const c_void,
                chunk.len(),
                -1,
            )
        };
        if r != 0 {
            bail!("BGZF compression failed");
        }
        out.truncate(at + dlen);
    }
    Ok(starts)
}

/// Virtual offset, relative to the start of a run of compressed blocks
/// (`block_starts`, `total` compressed bytes, holding `len` uncompressed
/// bytes), of uncompressed position `u` taken as the *end* of the preceding
/// data — as `bgzf_tell` reports it after reading up to `u`: a position at
/// the end of a block's data is the start of the next block (possibly
/// `total`, i.e. just past this run).
pub fn end_voffset(block_starts: &[u64], total: u64, len: usize, u: usize) -> u64 {
    if u == len {
        return total << 16;
    }
    let b = u / BGZF_BLOCK_SIZE;
    (block_starts[b] << 16) | (u % BGZF_BLOCK_SIZE) as u64
}

const CSI_MIN_SHIFT: c_int = 14;
const DEFAULT_N_LVLS: c_int = (31 - CSI_MIN_SHIFT + 2) / 3; // TBX_MAX_SHIFT = 31

/// Builds a CSI index (min_shift 14, tabix VCF metadata) while a bgzipped
/// VCF is written, producing the same index as `bcftools index` /
/// `tabix --csi` would from the finished file.
pub struct CsiBuilder {
    /// Created at the first record, as tbx_index() does (an index of a file
    /// without records keeps the default number of levels).
    idx: *mut hts::hts_idx_t,
    offset0: u64,
    n_lvls: c_int,
    last_contig: Option<Arc<str>>,
    n_contigs: c_int,
}

// SAFETY: owned and used by one thread at a time.
unsafe impl Send for CsiBuilder {}

impl CsiBuilder {
    /// `header_text`: the VCF header as written; `offset0`: the virtual
    /// offset just after it.
    pub fn new(header_text: &[u8], offset0: u64) -> Result<Self> {
        // As tbx_index(): size the levels from the largest ##contig length
        // (default 100 Gbp).
        let mut max_len: i64 = 0;
        for line in header_text.split(|&b| b == b'\n') {
            if !line.starts_with(b"##contig") {
                continue;
            }
            let Some(p) = line[8..].windows(6).position(|w| w == b"length") else { continue };
            let rest = &line[8 + p + 6..];
            let rest = &rest[rest.iter().take_while(|&&c| c == b' ' || c == b'=').count()..];
            let neg = rest.first() == Some(&b'-');
            let digits = &rest[usize::from(neg || rest.first() == Some(&b'+'))..];
            let v = digits
                .iter()
                .take_while(|c| c.is_ascii_digit())
                .fold(0i64, |a, &c| a.saturating_mul(10).saturating_add((c - b'0') as i64));
            max_len = max_len.max(if neg { -v } else { v });
        }
        if max_len == 0 {
            max_len = 100 * 1024 * 1024 * 1024;
        }
        let mut n_lvls: c_int = DEFAULT_N_LVLS;
        let mut s: i64 = 1 << (CSI_MIN_SHIFT + n_lvls * 3);
        let max_len = max_len + 256;
        while max_len > s {
            n_lvls += 1;
            s <<= 3;
        }
        Ok(Self { idx: ptr::null_mut(), offset0, n_lvls, last_contig: None, n_contigs: 0 })
    }

    fn init(&mut self, n_lvls: c_int) -> Result<()> {
        unsafe {
            let idx = hts::hts_idx_init(0, hts::HTS_FMT_CSI as c_int, self.offset0, CSI_MIN_SHIFT, n_lvls);
            if idx.is_null() {
                bail!("failed to initialise the CSI index");
            }
            // Tabix VCF configuration: preset, name/beg/end columns, '#', skip,
            // names length (updated by hts_idx_tbi_name).
            let mut conf = [0u8; 28];
            for (i, v) in [hts::TBX_VCF as u32, 1, 2, 0, b'#' as u32, 0, 0].iter().enumerate() {
                conf[i * 4..i * 4 + 4].copy_from_slice(&v.to_le_bytes());
            }
            if hts::hts_idx_set_meta(idx, conf.len() as u32, conf.as_mut_ptr(), 1) < 0 {
                hts::hts_idx_destroy(idx);
                bail!("failed to set CSI index metadata");
            }
            self.idx = idx;
        }
        Ok(())
    }

    /// Add one record (0-based `beg`, exclusive `end`) whose data ends at
    /// virtual offset `end_voff`. Records must come in file order.
    pub fn push(&mut self, contig: &Arc<str>, beg: i64, end: i64, end_voff: u64) -> Result<()> {
        if self.idx.is_null() {
            self.init(self.n_lvls)?;
        }
        unsafe {
            if self.last_contig.as_ref() != Some(contig) {
                self.n_contigs += 1;
                self.last_contig = Some(Arc::clone(contig));
            }
            let name = CString::new(contig.as_bytes())?;
            let tid = hts::hts_idx_tbi_name(self.idx, self.n_contigs, name.as_ptr());
            if tid < 0 || hts::hts_idx_push(self.idx, tid, beg, end, end_voff, 1) < 0 {
                bail!("failed to index record at {contig}:{} (is the input sorted?)", beg + 1);
            }
        }
        Ok(())
    }

    /// Finish with the virtual offset of the end of the data and write
    /// `<vcf_path>.csi`.
    pub fn save(mut self, final_voff: u64, vcf_path: &Path) -> Result<()> {
        if self.idx.is_null() {
            self.init(DEFAULT_N_LVLS)?;
        }
        let p = CString::new(vcf_path.to_str().ok_or_else(|| anyhow!("non-UTF-8 path"))?)?;
        unsafe {
            if hts::hts_idx_finish(self.idx, final_voff) != 0 {
                bail!("failed to finish the CSI index");
            }
            if hts::hts_idx_save_as(self.idx, p.as_ptr(), ptr::null(), hts::HTS_FMT_CSI as c_int) != 0 {
                bail!("failed to write {}.csi", vcf_path.display());
            }
        }
        Ok(())
    }
}

impl Drop for CsiBuilder {
    fn drop(&mut self) {
        if !self.idx.is_null() {
            unsafe { hts::hts_idx_destroy(self.idx) };
        }
    }
}
