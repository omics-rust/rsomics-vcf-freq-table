//! Streaming VCF parser for per-site allele counts.
//!
//! Reads plain VCF or gzipped VCF from a file path or stdin. For each data
//! record: locates the GT field in FORMAT, tallies allele indices across all
//! sample columns, and returns a `SiteRow` with raw counts.

use std::fs::File;
use std::io::{self, BufRead, BufReader, Read, stdin};
use std::path::Path;

use flate2::read::MultiGzDecoder;
use rsomics_common::{Result, RsomicsError};

use crate::table::{Mode, SiteRow};

fn open_reader(path: Option<&Path>) -> Result<Box<dyn BufRead>> {
    match path {
        None => Ok(Box::new(BufReader::new(stdin()))),
        Some(p) if p == Path::new("-") => Ok(Box::new(BufReader::new(stdin()))),
        Some(p) => {
            let file = File::open(p).map_err(|e| {
                RsomicsError::Io(io::Error::new(
                    e.kind(),
                    format!("cannot open {}: {e}", p.display()),
                ))
            })?;
            let is_gz = p.extension().is_some_and(|e| e.eq_ignore_ascii_case("gz"));
            let reader: Box<dyn Read> = if is_gz {
                Box::new(MultiGzDecoder::new(file))
            } else {
                Box::new(file)
            };
            Ok(Box::new(BufReader::new(reader)))
        }
    }
}

/// GT column index within the FORMAT field; returns `None` if GT absent.
fn gt_index(format: &str) -> Option<usize> {
    format.split(':').position(|f| f == "GT")
}

/// Tally all called allele indices in one sample's GT string. Missing alleles
/// (`.`) and phasing separators are handled; multi-digit allele indices are
/// supported.
fn tally_gt(gt: &str, counts: &mut [u32]) {
    // GT is `a[/|]b[/|]...` where each a/b is `.` or a decimal allele index.
    for part in gt.split(['/', '|']) {
        if part == "." {
            continue;
        }
        if let Ok(idx) = part.parse::<usize>() {
            if let Some(c) = counts.get_mut(idx) {
                *c += 1;
            }
            // allele index beyond the declared ALT list → skip (malformed VCF)
        }
    }
}

/// Parse one VCF data line into a `SiteRow`.
///
/// Returns `None` for lines that cannot be parsed (should not happen in valid
/// VCF, but guarding avoids crashing on INFO-only lines).
fn parse_line(line: &str) -> Option<SiteRow> {
    let mut cols = line.splitn(10, '\t');
    let chrom = cols.next()?;
    let pos_str = cols.next()?;
    let _id = cols.next()?;
    let ref_allele = cols.next()?;
    let alt_field = cols.next()?;
    let _qual = cols.next()?;
    let _filter = cols.next()?;
    let _info = cols.next()?;
    let format = cols.next()?;
    let samples_rest = cols.next()?; // everything after FORMAT

    let pos: u64 = pos_str.parse().ok()?;
    let gt_idx = gt_index(format)?;

    // Build allele list: REF first, then each comma-separated ALT.
    let mut alleles: Vec<String> = Vec::new();
    alleles.push(ref_allele.to_string());
    if alt_field != "." {
        for alt in alt_field.split(',') {
            alleles.push(alt.to_string());
        }
    }

    let mut counts = vec![0u32; alleles.len()];

    // Iterate sample columns (tab-separated within the trailing rest).
    let mut sample_str = samples_rest;
    loop {
        let (sample, rest) = match sample_str.split_once('\t') {
            Some((s, r)) => (s, r),
            None => (sample_str, ""),
        };
        if !sample.is_empty() {
            let gt = sample.split(':').nth(gt_idx).unwrap_or(".");
            tally_gt(gt, &mut counts);
        }
        if rest.is_empty() {
            break;
        }
        sample_str = rest;
    }

    let n_chr: u32 = counts.iter().sum();

    Some(SiteRow {
        chrom: chrom.to_string(),
        pos,
        alleles,
        counts,
        n_chr,
    })
}

/// Stream `path` (or stdin when `path` is `None`) and return one `SiteRow`
/// per VCF data record.
pub fn read_sites(path: Option<&Path>, _mode: Mode) -> Result<Vec<SiteRow>> {
    let reader = open_reader(path)?;
    let mut rows = Vec::new();

    for line in reader.lines() {
        let line = line
            .map_err(|e| RsomicsError::Io(io::Error::new(e.kind(), format!("reading VCF: {e}"))))?;
        if line.starts_with('#') {
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        if let Some(row) = parse_line(&line) {
            rows.push(row);
        }
    }

    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(line: &str) -> SiteRow {
        parse_line(line).expect("parse failed")
    }

    #[test]
    fn biallelic_diploid() {
        // chr1  100  .  A  T  50  PASS  .  GT  0/0  0/1  1/1
        let row = parse("chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1");
        assert_eq!(row.chrom, "chr1");
        assert_eq!(row.pos, 100);
        assert_eq!(row.alleles, ["A", "T"]);
        assert_eq!(row.counts, [3, 3]);
        assert_eq!(row.n_chr, 6);
    }

    #[test]
    fn multiallelic() {
        // chr1  300  .  T  A,G  40  PASS  .  GT  0/1  0/2  1/2
        let row = parse("chr1\t300\t.\tT\tA,G\t40\tPASS\t.\tGT\t0/1\t0/2\t1/2");
        assert_eq!(row.alleles, ["T", "A", "G"]);
        assert_eq!(row.counts, [2, 2, 2]);
        assert_eq!(row.n_chr, 6);
    }

    #[test]
    fn gt_index_in_format() {
        // FORMAT = DP:GT:GQ — GT is at index 1
        let row = parse("chr1\t100\t.\tA\tT\t.\t.\t.\tDP:GT:GQ\t10:0/1:30\t5:1/1:20");
        assert_eq!(row.counts, [1, 3]);
    }

    #[test]
    fn phased_gt() {
        let row = parse("chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0|1\t1|1");
        assert_eq!(row.counts, [1, 3]);
    }

    #[test]
    fn missing_allele_skipped() {
        // ./. contributes 0; 0/. contributes A:1; 0/1 contributes A:1 T:1
        let row = parse("chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t./.\t0/.\t0/1");
        assert_eq!(row.counts, [2, 1]);
        assert_eq!(row.n_chr, 3);
    }
}
