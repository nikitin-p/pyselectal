# pyselectal
Python script for filtering BAM alignments by 5′-end soft-clipping length and sequence, supporting single-end and paired-end reads with exact length and range-based selection. Designed to be easily integrated into NGS processing pipelines.

## Contents

- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Modes](#modes)
- [Options](#options)
- [Examples](#examples)
- [Paired-end behaviour](#paired-end-behaviour)
- [Notes](#notes)
- [Test data](#test-data)

## Concept and motivation

This tool is conceptually inspired by the read-level filtering strategy introduced in [Oguchi *et al.*, Science (2024)](https://www.science.org/doi/10.1126/science.add8394), where transcription start sites (TSSs) were inferred from precise 5′-end positions of 5′ single-cell RNA-seq reads. Specifically, Oguchi and colleagues distinguished transcription initiation from other events by the presence of the characteristic 5′ soft-clipped cap-dependent unencoded G base added by the reverse transcriptase during template switching.

Building on this approach, our tool enables general alignment filtering based on 5′-end soft-clipping patterns, mapped 5′-ends, and optional sequence constraints. While method-agnostic, it is particularly useful for CAGE, nanoCAGE, CAGEscan and other 5′-end-focused transcriptomics experiments, including bulk and single-cell protocols, where precise control over the 5′ read structure is critical for downstream analyses.

## Requirements

### Python

- **Python ≥ v3.6**
  - Required due to the use of f-strings.
  - **Python ≥ v3.8** is recommended for the best compatibility with modern `pysam` [link](https://pysam.readthedocs.io/en/latest/release.html).

### Python libraries

- **pysam ≥ 0.15.0**

`pysam 0.15.0` is the earliest version that supports the `threads=` argument in
`pysam.AlignmentFile`, which is used for parallel BGZF compression/decompression.
```bash
pip install pysam
```

## Installation

Clone the repository and make the script executable:

```bash
git clone https://github.com/nikitin-p/pyselectal.git
cd pyselectal
chmod +x pyselectal.py
```

You can then run it directly:

`python pyselectal.py --help`

## Usage
`python pyselectal.py [options] in.bam out.bam`
You can also use the script through the pipe:
```bash
script_1.sh \
  | python pyselectal.py [options] - - \
  | script_2.sh > output.file
```

### Modes of operation

The behaviour is determined by the relationship between `-n` and `-m`.

#### Mode 1 — Exact 5′ soft-clip (n = m > 0)

Require an exact n-bp 5′ soft-clip. Optional `--prefix` of length n must match the soft-clipped sequence.

Forward reads: soft-clipped sequence == prefix.
Reverse reads: soft-clipped sequence == revcomp(prefix).

#### Mode 2 — Mapped 5′-end (n = m = 0)

Require the 5′-end to be mapped (5′-end CIGAR operation = M). There are two ways to use this mode:

A. Prefix-based matching (`--prefix`)
Forward: read sequence starts with the full prefix
Reverse: read sequence ends with revcomp(prefix)

`-k` interaction:
If k <= len(prefix) (including k = 0), `-k` is ignored and selection is only by full prefix (you may emit a warning).
If k > len(prefix), require both:
full prefix match (as above), and at least k aligned MATCH bases at the 5′ end (CIGAR 5′ M length ≥ k).

B. Length-only matching (`-k` without `--prefix`)
Require at least k aligned MATCH bases at the 5′-end (CIGAR 5′ M length ≥ k). No sequence content is checked.
If k = 0, this reduces to requiring only a mapped 5′-end (5′ CIGAR op = M).

#### Mode 3 — Soft-clip range (n < m)

Require a 5′ soft-clip of length x, where n ≤ x ≤ m. 
Optional `--prefix` must be a single base (A/C/G/T/N).
Forward: all soft-clipped bases equal that base.
Reverse: all soft-clipped bases equal its complement.

## Options

`-n, --min-softclip`	Minimum 5′ soft-clip length (required)
`-m, --max-softclip`	Maximum 5′ soft-clip length (required)
`-p, --prefix`	Prefix / base (meaning depends on mode)
`-k`	Mapped mode only: minimum 5′ MATCH length or prefix match length
`-s, --sort`	Internally name-sort BAM before processing
`-t, --threads`	Number of BGZF threads (default: 1)
`--paired`	Enable paired-end mode
`-h, --help`	Show full manual

## Examples
**1. Single-end:** exact 3-bp soft-clip with prefix ATG.
```bash
python pyselectal.py \
    -n 3 -m 3 \
    -p ATG \
    in.bam out.bam
```

**2. Single-end:** mapped 5′-end with prefix ATG.
```bash
python pyselectal.py \
    -n 0 -m 0 \
    -p ATG \
    in.bam out.bam
```

**3. Single-end:** mapped 5′-end with at least 10 aligned bases (no sequence check).
```bash
python pyselectal.py \
    -n 0 -m 0 \
    -k 10 \
    in.bam out.bam
```

**4. Single-end:** 2–5 bp soft-clip, G homopolymer.
```bash
python pyselectal.py \
    -n 2 -m 5 \
    -p G \
    in.bam out.bam
```

**5. Paired-end:** exact 3-bp soft-clip on R1, emit matching R2 mates.
```bash
python pyselectal.py \
    -n 3 -m 3 \
    -p ATG \
    --paired \
    --sort \
    in.bam out.bam
```

**6. Single-end:** Prefix given, k ignored because k <= len(prefix).
```bash
python pyselectal.py \
    -n 0 -m 0 \
    -p ATGC \
    -k 3 \
    in.bam out.bam
```

**7. Single-end:** Prefix given, and additionally require at least 20 mapped bases at 5′-end.
```bash
python pyselectal.py \
    -n 0 -m 0 \
    -p ATGC \
    -k 20 \
    in.bam out.bam
```

## Paired-end behaviour

When --paired is enabled, Alignments are grouped by query_name. Only read1 (R1) is evaluated using the selected mode. 
For each passing R1:
- R1 is written and
- any R2 whose (reference_id, reference_start) matches R1’s recorded mate
coordinates is also written.

Input must be name-sorted, or you must use --sort.

## Notes
1. Unmapped reads are always rejected
2. 3′-soft-clipping is ignored (only the 5′-end matters)
3. Internal sorting uses samtools sort -n via pysam.sort

## Test data

The repository includes small, synthetic BAM files under `testdata/` that are
designed to **test all major modes of `pyselectal`**.  
These files are intended for **functional testing, debugging, and examples**,
not for benchmarking or performance evaluation.

### `test_softclip_se.bam`

Single-end test BAM containing reads with diverse 5′-end configurations:

- Exact 5′ soft-clips of varying lengths (`1S`, `2S`, `3S`, `4S`)
- Reads with **mapped 5′-ends** and **3′ soft-clips**
- Forward and reverse strand alignments
- Homopolymer soft-clips (`G` / `C`) suitable for range mode testing
- Multi-mapping reads (`NH:i:2`, `HI:i:*`) to verify that filtering is purely
  CIGAR- and sequence-based

### `test_softclip_pe.bam`

Paired-end test BAM designed specifically for paired-end mode (--paired):
- Read pairs are grouped by query_name
- Only read1 (R1) carries the relevant 5′ soft-clip or mapped pattern
- Corresponding read2 (R2) alignments are fully mapped

Multiple scenarios include:
- exact 5′ soft-clips on R1
- range soft-clips on R1
- mapped 5′ ends on R1
- cases where R1 has multiple alignments
- cases with multiple R2 candidates

Test BAMs are deliberately small and manually inspectable with:
```bash
samtools view testdata/test_softclip_se.bam
samtools view testdata/test_softclip_pe.bam
```