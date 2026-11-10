//! Compile the vendored bcftools expression engine (csrc/bcftools/filter.c)
//! against the htslib built by hts-sys.

fn main() {
    let hts_include = std::env::var("DEP_HTS_INCLUDE")
        .expect("hts-sys must export its include directory (DEP_HTS_INCLUDE)");
    cc::Build::new()
        .file("csrc/bcftools/filter.c")
        .file("csrc/shim.c")
        .include("csrc")
        .include("csrc/bcftools")
        .include(&hts_include)
        .flag_if_supported("-Wno-unused-parameter")
        .flag_if_supported("-Wno-sign-compare")
        .flag_if_supported("-Wno-unused-but-set-variable")
        .warnings(false)
        .compile("bcftools_filter");
    println!("cargo:rerun-if-changed=csrc");
}
