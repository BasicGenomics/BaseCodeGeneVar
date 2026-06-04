# BaseCode GeneVar

[![Version](https://img.shields.io/badge/version-0.1.1-blue.svg)](VERSION)

Inspect read support for known variant sites in RNA BaseCode data, per sample.

BaseCode GeneVar is a small, focused tool for BaseCode runs: you give it
a BAM and a list of known variants, and it tells you, for every sample
in the BAM, exactly how many reads support a variant at each position.

It is not a variant caller; it is a variant inspector
for sites you already care about (a returned ClinVar variant, a recurring
hotspot, a cohort marker, a positive control, you name it).

GeneVar is conversion-aware: it knows BaseCode chemistry deaminates G→A
on the + strand (C→T on the opposite strand) and by default flags
variants that collide with that conversion, so a base modification can't
silently masquerade as a true variant.

Samples are demultiplexed by the `SM:Z:` read-tag, so a single multiplexed
BAM produces one row per (variant, sample).

## Supported variant types

| Type  | Description                                                |
| ----- | ---------------------------------------------------------- |
| `snv` | Single-nucleotide variant.                                 |
| `ins` | Insertion, VCF-style anchored (REF = anchor base, ALT = anchor + inserted bases). |
| `del` | Deletion, VCF-style anchored (REF = anchor + deleted bases, ALT = anchor base). |

## Install

BaseCode GeneVar is a single script. Install its dependencies and run it directly:

```bash
pip install -r requirements.txt
python genevar.py --help
```

## Usage

```bash
python genevar.py \
    -i sample.bam \
    -v sample_variants.csv \
    -o calls.csv
```

The BAM must be coordinate-sorted and indexed.

### Which BAM to use

From a BaseCode run, two BAMs are typically usable:

- `*.stitched.molecules.sorted.bam` — **preferred.** One alignment per
  reconstructed molecule, so `depth` and `vaf` in the output are
  molecule-level (effectively PCR-deduplicated).
- `*.reads.aligned.trimmed.genetagged.sorted.markdup.reconstructed.sorted.bam`
  — also works. Counts are per read, so depth is higher but PCR families
  contribute multiple times to the totals.

### Options

| Flag           | Default | Description                                  |
| -------------- | ------- | -------------------------------------------- |
| `-i/--bam`     | —       | Input BAM.                         |
| `-v/--variants`| —       | Variants TSV/CSV. |
| `-o/--out`     | —       | Output TSV/CSV. |
| `--min-baseq`  | `0`     | Minimum base quality for SNV pileup.         |
| `--min-mapq`   | `0`     | Minimum mapping quality.                     |
| `--no-check-conversion` | (on by default) | Skip the BaseCode conversion-ambiguity check. |

### Conversion-aware mode (on by default)

The flag is driven by the **`strand`** column in the variants TSV/CSV (the
gene's strand on the reference, `+` or `−`):

- `strand = +` → flag SNVs whose REF/ALT is `G→A`.
- `strand = −` → flag SNVs whose REF/ALT is `C→T`.
- `strand` blank / column missing → flag both directions (conservative
  fallback, preserves the pre-strand behavior).

For each flagged variant GeneVar prints a stderr warning, still processes
it (the counts are still useful), and sets the `conversion_ambiguous`
column to `True` in the output so downstream filtering is easy.

Insertions and deletions don't add or remove a base from a single
deamination event, so they are flagged only in the unusual case where
the VCF anchor base on the REF side differs from the ALT side by one of
the ambiguous swaps.

## Variant TSV/CSV format

Tab- or comma-separated, 1-based VCF-style coordinates.

```
name                chrom  pos        ref   alt    type  strand
BRAF_V600E          7      140753336  A     T      snv   -
BRCA1_c.68_69dupAG  17     43044295   AG    AGAG   ins   -
CFTR_F508del        7      117559590  ATCT  A      del   +
```

The `strand` column is the gene's strand on the reference (`+` or `−`).
It is optional — variants without a strand still run, but the conversion
check falls back to flagging both `G→A` and `C→T` for them.

For **insertions**, `pos` is the anchor base (the base immediately
before the inserted sequence), `ref` is that anchor base, and `alt` is the
anchor followed by the inserted bases.

For **deletions**, `pos` is again the anchor base (immediately before the
deleted bases), `ref` is the anchor followed by the deleted bases, and
`alt` is the anchor base alone.

Ready-to-use examples are provided in both formats:
[`sample_variants.tsv`](sample_variants.tsv) and
[`sample_variants.csv`](sample_variants.csv).

## Output

One row per (variant, sample). Columns:

- `variant, chrom, pos, ref, alt, type` — copied from the input TSV/CSV.
- `sample` — value of the `SM` tag on the contributing reads.
- `depth` — total reads counted at the site (deduplicated by read name).
- `ref_count, alt_count` — reads supporting REF / ALT.
- `vaf` — `alt_count / depth`.
- For SNVs: `del_count`, `other_count`.
- For insertions: `ins_any_count`, `ins_other_count`, `top_inserted_seqs`
  (the three most common inserted sequences observed, for diagnostics).
- For deletions: `del_any_count`, `del_other_count`, `top_deleted_lengths`
  (the three most common deletion lengths observed).

  ---
© 2026 Basic Genomics AB · All rights reserved
