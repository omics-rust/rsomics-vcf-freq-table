//! Streaming VCF parser for per-site allele counts.
//!
//! Reads plain VCF or gzipped VCF from a file path or stdin. For each data
//! record: locates the GT field in FORMAT, tallies allele indices across all
//! sample columns, and returns a `SiteRow` with raw counts. A record whose
//! FORMAT lacks GT still yields a row — every genotype is treated as missing,
//! so N_CHR is 0, exactly as vcftools does.

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

/// Tally called alleles in one sample's GT string, returning the number of
/// called (non-missing) chromosomes — the site's N_CHR contribution.
///
/// A token is missing only when it is exactly `.`; anything else is a called
/// chromosome and counts toward N_CHR, matching vcftools. The allele bucket is
/// incremented from the token's leading decimal digits (`1X` reads as 1), and
/// only when that index falls inside the declared allele list — an out-of-range
/// or non-numeric token still counts as called but lands in no bucket.
fn tally_gt(gt: &str, counts: &mut [u32]) -> u32 {
    let mut called = 0;
    for part in gt.split(['/', '|']) {
        if part == "." || part.is_empty() {
            continue;
        }
        called += 1;
        let end = part
            .find(|c: char| !c.is_ascii_digit())
            .unwrap_or(part.len());
        if end > 0 {
            if let Ok(idx) = part[..end].parse::<usize>() {
                if let Some(c) = counts.get_mut(idx) {
                    *c += 1;
                }
            }
        }
    }
    called
}

/// Parse one VCF data line into a `SiteRow`.
///
/// `Ok(None)` marks a line with too few columns to be a data record. `Err`
/// signals a genotype vcftools refuses: a ploidy above 2 aborts the whole run,
/// matching `vcftools --freq`'s "Polyploidy found" bail.
fn parse_line(line: &str) -> Result<Option<SiteRow>> {
    let mut cols = line.splitn(10, '\t');
    let mut next = || cols.next();
    let (Some(chrom), Some(pos_str), Some(_id), Some(ref_allele), Some(alt_field)) =
        (next(), next(), next(), next(), next())
    else {
        return Ok(None);
    };
    let (Some(_qual), Some(_filter), Some(_info), Some(format), Some(samples_rest)) =
        (next(), next(), next(), next(), next())
    else {
        return Ok(None);
    };

    let Ok(pos) = pos_str.parse::<u64>() else {
        return Ok(None);
    };
    let gt_idx = gt_index(format);

    // Build allele list: REF first, then each comma-separated ALT. vcftools
    // upper-cases alleles before printing, so soft-masked bases collapse.
    let mut alleles: Vec<String> = Vec::new();
    alleles.push(ref_allele.to_ascii_uppercase());
    if alt_field != "." {
        for alt in alt_field.split(',') {
            alleles.push(alt.to_ascii_uppercase());
        }
    }

    let mut counts = vec![0u32; alleles.len()];
    let mut n_chr = 0u32;

    // With no GT field every genotype is missing: the row still prints, N_CHR 0.
    if let Some(gt_idx) = gt_idx {
        let mut sample_str = samples_rest;
        loop {
            let (sample, rest) = match sample_str.split_once('\t') {
                Some((s, r)) => (s, r),
                None => (sample_str, ""),
            };
            if !sample.is_empty() {
                let gt = sample.split(':').nth(gt_idx).unwrap_or(".");
                if gt.split(['/', '|']).count() > 2 {
                    return Err(RsomicsError::InvalidInput(format!(
                        "Polyploidy found, and not supported by vcftools: {chrom}:{pos}"
                    )));
                }
                n_chr += tally_gt(gt, &mut counts);
            }
            if rest.is_empty() {
                break;
            }
            sample_str = rest;
        }
    }

    Ok(Some(SiteRow {
        chrom: chrom.to_string(),
        pos,
        alleles,
        counts,
        n_chr,
    }))
}

/// Stream `path` (or stdin when `path` is `None`) and return one `SiteRow`
/// per VCF data record.
pub fn read_sites(path: Option<&Path>, mode: Mode) -> Result<Vec<SiteRow>> {
    let reader = open_reader(path)?;
    read_sites_from(reader, mode)
}

/// Parse VCF data records from an arbitrary reader into `SiteRow`s.
///
/// The `#CHROM` header line fixes the sample count: columns beyond the eight
/// mandatory fields plus FORMAT are individuals. A file with zero individuals
/// (a sites-only VCF, whether or not it carries a FORMAT column) cannot yield
/// frequency statistics, so vcftools bails with exit 1 and a fixed message; we
/// mirror that bail rather than emit an empty table.
pub fn read_sites_from(reader: impl BufRead, _mode: Mode) -> Result<Vec<SiteRow>> {
    let mut rows = Vec::new();

    for line in reader.lines() {
        let line = line
            .map_err(|e| RsomicsError::Io(io::Error::new(e.kind(), format!("reading VCF: {e}"))))?;
        if let Some(header) = line.strip_prefix("#CHROM") {
            let n_individuals = header.matches('\t').count().saturating_sub(8);
            if n_individuals == 0 {
                return Err(RsomicsError::InvalidInput(
                    "Require Genotypes in VCF file in order to output Frequency Statistics."
                        .to_string(),
                ));
            }
            continue;
        }
        if line.starts_with('#') {
            continue;
        }
        if line.trim().is_empty() {
            continue;
        }
        if let Some(row) = parse_line(&line)? {
            rows.push(row);
        }
    }

    Ok(rows)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn parse(line: &str) -> SiteRow {
        parse_line(line).expect("parse errored").expect("no row")
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

    #[test]
    fn gt_absent_from_format() {
        let row = parse("chr1\t100\t.\tA\tT\t.\t.\t.\tDP\t10");
        assert_eq!(row.counts, [0, 0]);
        assert_eq!(row.n_chr, 0);
    }

    #[test]
    fn out_of_range_allele_counts_toward_n_chr_only() {
        let row = parse("chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/2\t0/1");
        assert_eq!(row.counts, [2, 1]);
        assert_eq!(row.n_chr, 4);
    }

    #[test]
    fn non_numeric_token_is_called_but_unbucketed() {
        let row = parse("chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/X");
        assert_eq!(row.counts, [1, 0]);
        assert_eq!(row.n_chr, 2);
    }

    #[test]
    fn leading_digit_prefix_parses() {
        let row = parse("chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t1X/0");
        assert_eq!(row.counts, [1, 1]);
        assert_eq!(row.n_chr, 2);
    }

    #[test]
    fn alleles_upper_cased() {
        let row = parse("chr1\t100\t.\tacGt\taCg,n\t.\t.\t.\tGT\t0/1");
        assert_eq!(row.alleles, ["ACGT", "ACG", "N"]);
    }

    #[test]
    fn polyploid_aborts() {
        let err = parse_line("chr1\t200\t.\tA\tT\t.\t.\t.\tGT\t0/0/1").unwrap_err();
        assert!(err.to_string().contains("Polyploidy found"));
        assert!(err.to_string().contains("chr1:200"));
    }
}
