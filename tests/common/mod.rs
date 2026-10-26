//! Shared helpers for integration tests that shell out to external
//! bioinformatics tools (`bcftools`, `tabix`, `bgzip`).
//!
//! Tests that depend on these tools should start with:
//!
//! ```rust,ignore
//! if !common::check_tools(&["bcftools", "bgzip"]) {
//!     return; // or `return Ok(());`
//! }
//! ```
//!
//! `check_tools` prints `SKIP: missing required tools: ...` and returns
//! `false` when a tool is absent, unless `ERROR_ON_MISSING_BCFTOOLS=1` is
//! set in the environment, in which case it panics so the test fails with a
//! clear error (used in CI where the tools are expected to be installed).

use std::process::Command;

/// Check whether `tool` is available on PATH and runs successfully.
pub fn tool_available(tool: &str) -> bool {
    Command::new(tool)
        .arg("--version")
        .output()
        .map(|o| o.status.success())
        .unwrap_or(false)
}

/// Check that all `tools` are available on PATH.
///
/// Returns `true` if every tool is available.  If any tool is missing and
/// `ERROR_ON_MISSING_BCFTOOLS=1` is set, panics so the test fails loudly;
/// otherwise prints a skip message and returns `false`.
pub fn check_tools(tools: &[&str]) -> bool {
    let missing: Vec<&str> = tools
        .iter()
        .copied()
        .filter(|t| !tool_available(t))
        .collect();

    if missing.is_empty() {
        return true;
    }

    if std::env::var("ERROR_ON_MISSING_BCFTOOLS").as_deref() == Ok("1") {
        panic!(
            "required tools missing: {} (ERROR_ON_MISSING_BCFTOOLS=1)",
            missing.join(", ")
        );
    }

    eprintln!("SKIP: missing required tools: {}", missing.join(", "));
    false
}
