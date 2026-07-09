//! Per-site data model and TSV rendering for the four vcftools freq/count modes.

use rsomics_common::fmt::format_g6;
use serde::Serialize;

/// Which of the four vcftools output modes to produce.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "snake_case")]
pub enum Mode {
    /// `--freq`: `ALLELE:FREQ` columns, 6-sig-fig float.
    Freq,
    /// `--freq2`: bare FREQ columns, no allele labels.
    Freq2,
    /// `--counts`: `ALLELE:COUNT` columns, integer.
    Counts,
    /// `--counts2`: bare COUNT columns, no allele labels.
    Counts2,
}

impl Mode {
    pub fn is_freq(self) -> bool {
        matches!(self, Self::Freq | Self::Freq2)
    }

    pub fn labeled(self) -> bool {
        matches!(self, Self::Freq | Self::Counts)
    }

    /// Output file extension vcftools uses (informational only; we write stdout).
    pub fn extension(self) -> &'static str {
        match self {
            Self::Freq | Self::Freq2 => "frq",
            Self::Counts | Self::Counts2 => "frq.count",
        }
    }

    pub fn header(self) -> &'static str {
        match self {
            Self::Freq => "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}",
            Self::Freq2 => "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}",
            Self::Counts => "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}",
            Self::Counts2 => "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}",
        }
    }
}

impl std::str::FromStr for Mode {
    type Err = String;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        match s {
            "freq" => Ok(Self::Freq),
            "freq2" => Ok(Self::Freq2),
            "counts" => Ok(Self::Counts),
            "counts2" => Ok(Self::Counts2),
            other => Err(format!(
                "unknown mode {other:?}; expected freq, freq2, counts, or counts2"
            )),
        }
    }
}

/// Allele counts for one VCF record.
#[derive(Debug, Clone, Serialize)]
pub struct SiteRow {
    pub chrom: String,
    pub pos: u64,
    /// REF followed by each ALT allele, in VCF column order.
    pub alleles: Vec<String>,
    /// Observed copy count per allele (parallel to `alleles`).
    pub counts: Vec<u32>,
    /// Total non-missing allele copies (= sum of `counts`).
    pub n_chr: u32,
}

impl SiteRow {
    /// Format one value column for `mode`.
    fn col(&self, idx: usize, mode: Mode) -> String {
        let count = self.counts[idx];
        if mode.labeled() {
            let allele = &self.alleles[idx];
            if mode.is_freq() {
                let freq = count as f64 / self.n_chr as f64;
                format!("{}:{}", allele, format_g6(freq))
            } else {
                format!("{}:{}", allele, count)
            }
        } else if mode.is_freq() {
            let freq = count as f64 / self.n_chr as f64;
            format_g6(freq)
        } else {
            count.to_string()
        }
    }

    /// Render as one TSV line (no trailing newline).
    pub fn to_tsv_line(&self, mode: Mode) -> String {
        let n = self.alleles.len();
        let mut parts: Vec<String> = Vec::with_capacity(4 + n);
        parts.push(self.chrom.clone());
        parts.push(self.pos.to_string());
        parts.push(n.to_string());
        parts.push(self.n_chr.to_string());
        for i in 0..n {
            parts.push(self.col(i, mode));
        }
        parts.join("\t")
    }
}

/// Render the full table (header + one line per site) for `mode`.
pub fn to_tsv(rows: &[SiteRow], mode: Mode) -> String {
    let mut out = String::with_capacity(rows.len() * 64);
    out.push_str(mode.header());
    out.push('\n');
    for row in rows {
        out.push_str(&row.to_tsv_line(mode));
        out.push('\n');
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn format_g6_cases() {
        assert_eq!(format_g6(0.5), "0.5");
        assert_eq!(format_g6(1.0 / 3.0), "0.333333");
        assert_eq!(format_g6(2.0 / 3.0), "0.666667");
        assert_eq!(format_g6(5.0 / 6.0), "0.833333");
        assert_eq!(format_g6(1.0 / 6.0), "0.166667");
        assert_eq!(format_g6(0.0), "0");
        assert_eq!(format_g6(1.0), "1");
    }

    #[test]
    fn format_g6_six_sig_figs_below_tenth() {
        // %g keeps six significant figures, not six decimal places: 1/60.
        assert_eq!(format_g6(1.0 / 60.0), "0.0166667");
        assert_eq!(format_g6(59.0 / 60.0), "0.983333");
        assert_eq!(format_g6(11999.0 / 12000.0), "0.999917");
    }

    #[test]
    fn format_g6_scientific_below_1e_minus_4() {
        // Exponent < -4 switches to %e style with a two-digit signed exponent.
        assert_eq!(format_g6(1.0 / 12000.0), "8.33333e-05");
        assert_eq!(format_g6(5e-5), "5e-05");
        // Exactly 1e-4 stays fixed (exponent -4 is not < -4).
        assert_eq!(format_g6(1e-4), "0.0001");
    }

    #[test]
    fn format_g6_nan_is_lowercase() {
        assert_eq!(format_g6(f64::NAN), "nan");
        // The all-missing site computes 0 / N_CHR with N_CHR == 0.
        let n_chr = 0u32;
        assert_eq!(format_g6(0u32 as f64 / n_chr as f64), "nan");
    }

    #[test]
    fn mode_from_str() {
        assert_eq!("freq".parse::<Mode>().unwrap(), Mode::Freq);
        assert_eq!("freq2".parse::<Mode>().unwrap(), Mode::Freq2);
        assert_eq!("counts".parse::<Mode>().unwrap(), Mode::Counts);
        assert_eq!("counts2".parse::<Mode>().unwrap(), Mode::Counts2);
        assert!("bad".parse::<Mode>().is_err());
    }

    #[test]
    fn site_row_tsv_freq() {
        let row = SiteRow {
            chrom: "chr1".to_string(),
            pos: 100,
            alleles: vec!["A".to_string(), "T".to_string()],
            counts: vec![3, 3],
            n_chr: 6,
        };
        assert_eq!(row.to_tsv_line(Mode::Freq), "chr1\t100\t2\t6\tA:0.5\tT:0.5");
    }

    #[test]
    fn site_row_tsv_counts() {
        let row = SiteRow {
            chrom: "chr1".to_string(),
            pos: 100,
            alleles: vec!["A".to_string(), "T".to_string()],
            counts: vec![3, 3],
            n_chr: 6,
        };
        assert_eq!(row.to_tsv_line(Mode::Counts), "chr1\t100\t2\t6\tA:3\tT:3");
    }
}
