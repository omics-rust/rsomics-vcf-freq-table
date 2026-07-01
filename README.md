# rsomics-vcf-freq-table

Allele frequency and count tables from a VCF — a fast Rust reimplementation of
`vcftools 0.1.17` `--freq`, `--freq2`, `--counts`, and `--counts2`.

## Usage

```
rsomics-vcf-freq-table [OPTIONS] [VCF]

Arguments:
  [VCF]  Input VCF or VCF.gz; omit or pass `-` to read stdin

Options:
  --mode <MODE>    freq | freq2 | counts | counts2  [default: freq]
  -t, --threads    Worker threads
  -q, --quiet      Suppress progress
      --json       JSON envelope output
```

## Output formats

### `--mode freq` (default)
```
CHROM	POS	N_ALLELES	N_CHR	{ALLELE:FREQ}
chr1	100	2	6	A:0.5	T:0.5
chr1	300	3	6	T:0.333333	A:0.333333	G:0.333333
```

### `--mode freq2`
Same columns but without allele labels on the value fields; header says `{FREQ}`.

### `--mode counts`
```
CHROM	POS	N_ALLELES	N_CHR	{ALLELE:COUNT}
chr1	100	2	6	A:3	T:3
```

### `--mode counts2`
Same columns without allele labels; header says `{COUNT}`.

## Algorithm

- **N_CHR**: total non-missing allele copies across all samples (diploid → up to 2 per sample; missing `.` alleles excluded).
- **Allele order**: REF first, then each ALT in VCF column order.
- **Frequency format**: 6 significant figures, trailing zeros stripped — matching C `printf("%g", x)`.

## Origin

This crate is an independent Rust reimplementation of `vcftools 0.1.17`
`--freq`/`--freq2`/`--counts`/`--counts2` based on:

- Black-box behavioral testing against vcftools 0.1.17
- The VCFv4 format specification

No source code from the LGPL vcftools was read during implementation.

License: MIT OR Apache-2.0.  
Upstream credit: [vcftools](https://github.com/vcftools/vcftools) (LGPL).
