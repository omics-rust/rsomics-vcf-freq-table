//! Value-exact compatibility against vcftools 0.1.17 --freq/--freq2/--counts/--counts2.
//!
//! The reference VCF is embedded inline; expected outputs were derived from
//! black-box observation of vcftools 0.1.17. No vcftools binary is required
//! at test time.

use std::io::Write;
use std::path::PathBuf;
use std::process::Command;

fn bin() -> PathBuf {
    PathBuf::from(env!("CARGO_BIN_EXE_rsomics-vcf-freq-table"))
}

/// Write `content` to a temp file under the system temp dir (TMPDIR env var
/// or OS default); the path is returned. The file is not explicitly deleted —
/// the OS cleans it up on process exit.
fn write_temp_vcf(content: &str) -> PathBuf {
    let mut f = tempfile::NamedTempFile::new().expect("tempfile");
    f.write_all(content.as_bytes()).expect("write");
    f.keep().expect("keep").1
}

/// Black-box VCF: 3 samples, biallelic + multiallelic sites across 2 chroms.
const TEST_VCF: &str = "\
##fileformat=VCFv4.1\n\
##FILTER=<ID=PASS,Description=\"All filters passed\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3\n\
chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1\n\
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/1\n\
chr1\t300\t.\tT\tA,G\t40\tPASS\t.\tGT\t0/1\t0/2\t1/2\n\
chr2\t100\t.\tC\tT\t55\tPASS\t.\tGT\t0/0\t0/0\t0/1\n";

fn run(vcf: &std::path::Path, mode: &str) -> String {
    let out = Command::new(bin())
        .arg("--mode")
        .arg(mode)
        .arg(vcf)
        .output()
        .expect("run binary");
    assert!(
        out.status.success(),
        "binary exited non-zero for mode={mode}: {}",
        String::from_utf8_lossy(&out.stderr)
    );
    String::from_utf8(out.stdout).expect("utf8")
}

#[test]
fn freq_exact() {
    let vcf = write_temp_vcf(TEST_VCF);
    let got = run(&vcf, "freq");
    let want = "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n\
chr1\t100\t2\t6\tA:0.5\tT:0.5\n\
chr1\t200\t2\t6\tG:0.333333\tC:0.666667\n\
chr1\t300\t3\t6\tT:0.333333\tA:0.333333\tG:0.333333\n\
chr2\t100\t2\t6\tC:0.833333\tT:0.166667\n";
    assert_eq!(got, want, "freq output differs");
}

#[test]
fn freq2_exact() {
    let vcf = write_temp_vcf(TEST_VCF);
    let got = run(&vcf, "freq2");
    let want = "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\n\
chr1\t100\t2\t6\t0.5\t0.5\n\
chr1\t200\t2\t6\t0.333333\t0.666667\n\
chr1\t300\t3\t6\t0.333333\t0.333333\t0.333333\n\
chr2\t100\t2\t6\t0.833333\t0.166667\n";
    assert_eq!(got, want, "freq2 output differs");
}

#[test]
fn counts_exact() {
    let vcf = write_temp_vcf(TEST_VCF);
    let got = run(&vcf, "counts");
    let want = "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n\
chr1\t100\t2\t6\tA:3\tT:3\n\
chr1\t200\t2\t6\tG:2\tC:4\n\
chr1\t300\t3\t6\tT:2\tA:2\tG:2\n\
chr2\t100\t2\t6\tC:5\tT:1\n";
    assert_eq!(got, want, "counts output differs");
}

#[test]
fn counts2_exact() {
    let vcf = write_temp_vcf(TEST_VCF);
    let got = run(&vcf, "counts2");
    let want = "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\n\
chr1\t100\t2\t6\t3\t3\n\
chr1\t200\t2\t6\t2\t4\n\
chr1\t300\t3\t6\t2\t2\t2\n\
chr2\t100\t2\t6\t5\t1\n";
    assert_eq!(got, want, "counts2 output differs");
}

#[test]
fn missing_gt_excluded() {
    // A site with ./. and 0/. — only called alleles counted.
    const MISSING_VCF: &str = "\
##fileformat=VCFv4.1\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3\n\
chr1\t500\t.\tA\tT\t.\t.\t.\tGT\t./.\t0/.\t0/1\n";
    let vcf = write_temp_vcf(MISSING_VCF);
    let got = run(&vcf, "counts");
    // S1 contributes 0, S2 contributes A:1, S3 contributes A:1 T:1 → n_chr=3
    let want = "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n\
chr1\t500\t2\t3\tA:2\tT:1\n";
    assert_eq!(got, want, "missing-GT site differs");
}
