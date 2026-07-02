//! Value-exact compatibility against vcftools 0.1.17
//! --freq/--freq2/--counts/--counts2.
//!
//! Every expected output below was captured once from vcftools 0.1.17 on
//! macOS/arm64 and pasted as a constant; no vcftools (or any subprocess) runs
//! at test time. Inputs are parsed straight through the library reader, so the
//! tests exercise the real parse + tally + render path without touching the
//! filesystem.

use std::io::Cursor;

use rsomics_vcf_freq_table::table::to_tsv;
use rsomics_vcf_freq_table::vcf::read_sites_from;
use rsomics_vcf_freq_table::{Mode, SiteRow};

fn render(vcf: &str, mode: Mode) -> String {
    let rows: Vec<SiteRow> = read_sites_from(Cursor::new(vcf.as_bytes()), mode).expect("parse VCF");
    to_tsv(&rows, mode)
}

fn try_render(vcf: &str, mode: Mode) -> Result<String, String> {
    read_sites_from(Cursor::new(vcf.as_bytes()), mode)
        .map(|rows| to_tsv(&rows, mode))
        .map_err(|e| e.to_string())
}

fn check_all(name: &str, vcf: &str, freq: &str, freq2: &str, counts: &str, counts2: &str) {
    assert_eq!(render(vcf, Mode::Freq), freq, "{name}: freq");
    assert_eq!(render(vcf, Mode::Freq2), freq2, "{name}: freq2");
    assert_eq!(render(vcf, Mode::Counts), counts, "{name}: counts");
    assert_eq!(render(vcf, Mode::Counts2), counts2, "{name}: counts2");
}

/// Single-site VCF with one REF/ALT and an explicit per-sample GT list.
fn one_site_vcf(ref_a: &str, alt: &str, gts: &[&str]) -> String {
    let mut s =
        String::from("##fileformat=VCFv4.1\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT");
    for i in 0..gts.len() {
        s.push_str(&format!("\tS{i}"));
    }
    s.push_str("\nchr1\t100\t.\t");
    s.push_str(ref_a);
    s.push('\t');
    s.push_str(alt);
    s.push_str("\t.\t.\t.\tGT");
    for g in gts {
        s.push('\t');
        s.push_str(g);
    }
    s.push('\n');
    s
}

const BASIC_VCF: &str = "\
##fileformat=VCFv4.1
##FILTER=<ID=PASS,Description=\"All filters passed\">
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSample1\tSample2\tSample3
chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT\t0/0\t0/1\t1/1
chr1\t200\t.\tG\tC\t60\tPASS\t.\tGT\t0/1\t1/1\t0/1
chr1\t300\t.\tT\tA,G\t40\tPASS\t.\tGT\t0/1\t0/2\t1/2
chr2\t100\t.\tC\tT\t55\tPASS\t.\tGT\t0/0\t0/0\t0/1
";

#[test]
fn basic_biallelic_and_multiallelic() {
    check_all(
        "basic",
        BASIC_VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n\
chr1\t100\t2\t6\tA:0.5\tT:0.5\n\
chr1\t200\t2\t6\tG:0.333333\tC:0.666667\n\
chr1\t300\t3\t6\tT:0.333333\tA:0.333333\tG:0.333333\n\
chr2\t100\t2\t6\tC:0.833333\tT:0.166667\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\n\
chr1\t100\t2\t6\t0.5\t0.5\n\
chr1\t200\t2\t6\t0.333333\t0.666667\n\
chr1\t300\t3\t6\t0.333333\t0.333333\t0.333333\n\
chr2\t100\t2\t6\t0.833333\t0.166667\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n\
chr1\t100\t2\t6\tA:3\tT:3\n\
chr1\t200\t2\t6\tG:2\tC:4\n\
chr1\t300\t3\t6\tT:2\tA:2\tG:2\n\
chr2\t100\t2\t6\tC:5\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\n\
chr1\t100\t2\t6\t3\t3\n\
chr1\t200\t2\t6\t2\t4\n\
chr1\t300\t3\t6\t2\t2\t2\n\
chr2\t100\t2\t6\t5\t1\n",
    );
}

/// Minor-allele frequency 1/60 — a value below 0.1 where naive `{:.6}` fixed
/// decimals diverge from `%g` six-significant-figures.
#[test]
fn small_frequency_six_sig_figs() {
    let mut gts = vec!["0/0"; 29];
    gts.push("0/1");
    let vcf = one_site_vcf("A", "T", &gts);
    check_all(
        "small1",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t60\tA:0.983333\tT:0.0166667\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t60\t0.983333\t0.0166667\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t60\tA:59\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t60\t59\t1\n",
    );
}

/// Minor-allele frequency 1/12000 — `%g` switches to scientific notation
/// (`8.33333e-05`) once the exponent drops below -4.
#[test]
fn tiny_frequency_scientific_notation() {
    let mut gts = vec!["0/0"; 5999];
    gts.push("0/1");
    let vcf = one_site_vcf("A", "T", &gts);
    check_all(
        "small3",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t12000\tA:0.999917\tT:8.33333e-05\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t12000\t0.999917\t8.33333e-05\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t12000\tA:11999\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t12000\t11999\t1\n",
    );
}

/// All genotypes missing → N_CHR=0 and 0/0 frequencies print lowercase `nan`.
#[test]
fn all_missing_site_nan() {
    let vcf = one_site_vcf("A", "T", &["./.", "./.", "./."]);
    check_all(
        "allmiss",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t0\tA:nan\tT:nan\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t0\tnan\tnan\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t0\tA:0\tT:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t0\t0\t0\n",
    );
}

#[test]
fn single_sample() {
    let vcf = one_site_vcf("A", "T", &["0/1"]);
    check_all(
        "single",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t2\tA:0.5\tT:0.5\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t2\t0.5\t0.5\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t2\tA:1\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t2\t1\t1\n",
    );
}

/// Monomorphic site: ALT frequency is exactly 0, REF exactly 1.
#[test]
fn monomorphic() {
    let vcf = one_site_vcf("A", "T", &["0/0", "0/0", "0/0"]);
    check_all(
        "mono",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t6\tA:1\tT:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t6\t1\t0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t6\tA:6\tT:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t6\t6\t0\n",
    );
}

/// Half-calls (`0/.`) contribute only their called allele; `./.` contributes
/// nothing.
#[test]
fn half_calls_and_missing() {
    let vcf = one_site_vcf("A", "T", &["0/.", "0/1", "./."]);
    check_all(
        "half",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t3\tA:0.666667\tT:0.333333\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t3\t0.666667\t0.333333\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t3\tA:2\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t3\t2\t1\n",
    );
}

/// Phased separators (`|`) tally identically to unphased.
#[test]
fn phased_genotypes() {
    let vcf = one_site_vcf("A", "T", &["0|1", "1|1"]);
    check_all(
        "phased",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t4\tA:0.25\tT:0.75\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t4\t0.25\t0.75\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t4\tA:1\tT:3\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t4\t1\t3\n",
    );
}

/// A standalone triallelic site: REF then ALT[0], ALT[1] in column order.
#[test]
fn multiallelic() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t300\t.\tT\tA,G\t.\t.\t.\tGT\t0/1\t0/2\t1/2
";
    check_all(
        "multi",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t300\t3\t6\tT:0.333333\tA:0.333333\tG:0.333333\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t300\t3\t6\t0.333333\t0.333333\t0.333333\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t300\t3\t6\tT:2\tA:2\tG:2\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t300\t3\t6\t2\t2\t2\n",
    );
}

/// FORMAT with GT not first, plus a `./.` and a half-call.
#[test]
fn missing_mix() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2\tS3
chr1\t500\t.\tA\tT\t.\t.\t.\tGT\t./.\t0/.\t0/1
";
    check_all(
        "missmix",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t500\t2\t3\tA:0.666667\tT:0.333333\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t500\t2\t3\t0.666667\t0.333333\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t500\t2\t3\tA:2\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t500\t2\t3\t2\t1\n",
    );
}

/// FORMAT lacks a GT subfield: vcftools treats every genotype as missing and
/// still emits the row with N_CHR=0 (freqs `nan`, counts 0). We must not drop
/// the site.
#[test]
fn gt_absent_from_format() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tT\t.\t.\t.\tDP\t10
";
    check_all(
        "gtabsent",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t0\tA:nan\tT:nan\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t0\tnan\tnan\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t0\tA:0\tT:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t0\t0\t0\n",
    );
}

/// GT absent on a multiallelic site: all three alleles print with N_CHR=0.
#[test]
fn gt_absent_multiallelic() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\t.\tA\tG,C\t.\t.\t.\tDP\t5\t9
";
    check_all(
        "gtabsentmulti",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t3\t0\tA:nan\tG:nan\tC:nan\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t3\t0\tnan\tnan\tnan\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t3\t0\tA:0\tG:0\tC:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t3\t0\t0\t0\t0\n",
    );
}

/// Out-of-range allele index (`0/2` with a single ALT): the `2` copy is counted
/// in N_CHR (=4) but lands in no bucket, so freqs use 4 as denominator.
#[test]
fn out_of_range_allele() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\tS2
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/2\t0/1
";
    check_all(
        "oor",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t4\tA:0.5\tT:0.25\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t4\t0.5\t0.25\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t4\tA:2\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t4\t2\t1\n",
    );
}

/// Non-numeric allele token (`X`): counted as a called chromosome (N_CHR=2) but
/// bucketed nowhere, while the leading digit of `0` still lands in REF.
#[test]
fn non_numeric_allele_token() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tA\tT\t.\t.\t.\tGT\t0/X
";
    check_all(
        "nonnum",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t2\tA:0.5\tT:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t2\t0.5\t0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t2\tA:1\tT:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t2\t1\t0\n",
    );
}

/// Soft-masked (lowercase) alleles are upper-cased in the output, exactly as
/// vcftools does; `acGt`→`ACGT`, `n`→`N`.
#[test]
fn lowercase_alleles_upper_cased() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
chr1\t100\t.\tacGt\taCg,n\t.\t.\t.\tGT\t0/1
";
    check_all(
        "lowercase",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t3\t2\tACGT:0.5\tACG:0.5\tN:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t3\t2\t0.5\t0.5\t0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t3\t2\tACGT:1\tACG:1\tN:0\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t3\t2\t1\t1\t0\n",
    );
}

/// A whole-sample `.` counts as missing; the other diploid sample supplies the
/// two called chromosomes.
#[test]
fn bare_dot_sample() {
    let vcf = one_site_vcf("A", "T", &[".", "0/1"]);
    check_all(
        "baredot",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t2\tA:0.5\tT:0.5\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t2\t0.5\t0.5\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t2\tA:1\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t2\t1\t1\n",
    );
}

/// Haploid calls (`0`, `1`) each contribute a single chromosome.
#[test]
fn haploid_calls() {
    let vcf = one_site_vcf("A", "T", &["0", "1"]);
    check_all(
        "haploid",
        &vcf,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\nchr1\t100\t2\t2\tA:0.5\tT:0.5\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\nchr1\t100\t2\t2\t0.5\t0.5\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\nchr1\t100\t2\t2\tA:1\tT:1\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\nchr1\t100\t2\t2\t1\t1\n",
    );
}

/// Polyploid genotypes (ploidy > 2) are unsupported by vcftools, which aborts
/// with exit code 1 and a "Polyploidy found" message at the offending site. We
/// match by failing loud rather than emitting a table vcftools refuses to.
#[test]
fn polyploid_is_rejected() {
    let vcf = one_site_vcf("A", "T", &["0/0/1"]);
    let err = try_render(&vcf, Mode::Freq).expect_err("polyploid should abort");
    assert!(err.contains("Polyploidy found"), "got: {err}");
    assert!(err.contains("chr1:100"), "got: {err}");
}

/// Sites-only VCF (eight mandatory columns, no FORMAT/sample columns): vcftools
/// counts zero individuals and aborts with exit 1 before writing any output,
/// under all four modes. We match by failing loud with the same core message.
#[test]
fn sites_only_requires_genotypes() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO
chr1\t100\t.\tA\tT\t50\tPASS\t.
chr1\t200\t.\tG\tC\t60\tPASS\t.
";
    for mode in [Mode::Freq, Mode::Freq2, Mode::Counts, Mode::Counts2] {
        let err = try_render(VCF, mode).expect_err("sites-only should abort");
        assert!(
            err.contains("Require Genotypes in VCF file in order to output Frequency Statistics."),
            "mode {mode:?} got: {err}"
        );
    }
}

/// A FORMAT column with zero sample columns is likewise zero individuals: same
/// abort as the eight-column sites-only case.
#[test]
fn format_without_samples_requires_genotypes() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT
chr1\t100\t.\tA\tT\t50\tPASS\t.\tGT
";
    for mode in [Mode::Freq, Mode::Freq2, Mode::Counts, Mode::Counts2] {
        let err = try_render(VCF, mode).expect_err("zero-individual should abort");
        assert!(
            err.contains("Require Genotypes in VCF file in order to output Frequency Statistics."),
            "mode {mode:?} got: {err}"
        );
    }
}

/// Header only, zero data records: just the column header line.
#[test]
fn empty_no_records() {
    const VCF: &str = "\
##fileformat=VCFv4.1
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1
";
    check_all(
        "empty",
        VCF,
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:FREQ}\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{FREQ}\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{ALLELE:COUNT}\n",
        "CHROM\tPOS\tN_ALLELES\tN_CHR\t{COUNT}\n",
    );
}
