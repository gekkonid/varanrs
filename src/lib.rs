//! varanrs -- VCF/BCF maniupulation in Rust
//!
//! Exposes the parallel variant processor (`processor::ParallelVariantWindowProcessor`)

pub mod args;
pub mod commands;
pub mod filltags;
pub mod filter;
pub mod hts;
pub mod processor;
pub mod snpsketch;
pub mod util;
