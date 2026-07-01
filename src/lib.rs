//! Allele frequency and count tables from a VCF.
//!
//! Reimplements vcftools 0.1.17 `--freq`, `--freq2`, `--counts`, `--counts2`.
//!
//! * `--freq` / `--counts` — columns labelled `ALLELE:VALUE`
//! * `--freq2` / `--counts2` — columns are bare values, no allele prefix
//!
//! N_CHR counts called (non-missing) allele copies; missing GT alleles (`.`)
//! are excluded. Multiallelic sites emit all alleles: REF then ALT[0], ALT[1],
//! in order.

pub mod table;
pub mod vcf;

use std::path::Path;

use rsomics_common::Result;
pub use table::{Mode, SiteRow};

/// Compute the frequency/count table for `path`, streaming each record.
///
/// The returned `Vec<SiteRow>` is in VCF record order. Each row holds the raw
/// allele counts; `to_tsv` renders according to `mode`.
pub fn freq_table(path: Option<&Path>, mode: Mode) -> Result<Vec<SiteRow>> {
    vcf::read_sites(path, mode)
}
