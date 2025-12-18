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
- [License](#license)

## Concept and motivation

This script is conceptually inspired by the read-level filtering strategy introduced in [Oguchi *et al.*, *Science* (2024)](https://www.science.org/doi/10.1126/science.add8394), where enhancer and transcription start site (TSS) activity is inferred from precise 5′-end signatures of capped RNAs in 5′ single-cell RNA-seq data.

In that work, bona fide transcription initiation events are distinguished from background by leveraging characteristic 5′ soft-clipping and unencoded nucleotides generated during template switching at capped RNA ends. This principle is central to accurately identifying transcription start sites and bidirectionally transcribed enhancers.

Building on this idea, the present script provides a general, alignment-level filtering tool for selecting reads based on 5′-end soft-clipping patterns, mapped 5′ ends, and optional sequence constraints. While method-agnostic, it is particularly useful for CAGE, nanoCAGE, CAGEscan and other 5′-end–focused transcriptomics experiments, including bulk and single-cell protocols, where precise control over 5′ read structure is critical for downstream analyses.

## Requirements

### Python

- **Python ≥ v3.6**
  - Required due to use of f-strings.
  - Python ≥ v3.8 recommended for best compatibility with modern `pysam`.

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
git clone https://github.com/nikitin-p/softclip5.git
cd softclip5
chmod +x softclip5.py
```

You can then run it directly:

`./softclip5.py --help`

## Usage
`softclip5.py [options] in.bam out.bam`

### Modes of operation

The behaviour is determined by the relationship between -n and -m.

#### Mode 1 — Exact 5′ soft-clip (n = m > 0)

Require an exact n-bp 5′ soft-clip

Optional --prefix of length n must match the soft-clipped sequence

Forward reads:

soft-clipped sequence == prefix


Reverse reads:

soft-clipped sequence == revcomp(prefix)

#### Mode 2 — Mapped 5′ end (n = m = 0)

Require the 5′ end to be mapped (CIGAR op = M)

Two independent alternatives:

A. Prefix-based matching (--prefix)

Forward: read sequence starts with prefix

Reverse: read sequence ends with revcomp(prefix)

-k optionally restricts matching to the first k bases of the prefix
(-k 0 = full prefix length)

B. Length-only matching (-k without --prefix)

Require at least k aligned MATCH bases at the 5′ end (CIGAR only)

No sequence content is checked

#### Mode 3 — Soft-clip range (n < m)

Require a 5′ soft-clip of length x, where n ≤ x ≤ m

Optional --prefix must be a single base (A/C/G/T/N)

Forward: all soft-clipped bases equal that base

Reverse: all soft-clipped bases equal its complement

## Options

-n, --min-softclip	Minimum 5′ soft-clip length (required)
-m, --max-softclip	Maximum 5′ soft-clip length (required)
-p, --prefix	Prefix / base (meaning depends on mode)
-k	Mapped mode only: minimum 5′ MATCH length or prefix match length
-s, --sort	Internally name-sort BAM before processing
-t, --threads	Number of BGZF threads (default: 1)
--paired	Enable paired-end mode
-h, --help	Show full manual

## Examples
1. Single-end: exact 3-bp soft-clip with prefix ATG
softclip5.py \
  -n 3 -m 3 \
  -p ATG \
  in.bam out.bam

2. Single-end: mapped 5′ end with prefix ATG
softclip5.py \
  -n 0 -m 0 \
  -p ATG \
  in.bam out.bam

3. Single-end: mapped 5′ end with at least 10 aligned bases (no sequence check)
softclip5.py \
  -n 0 -m 0 \
  -k 10 \
  in.bam out.bam

4. Single-end: 2–5 bp soft-clip, G homopolymer
softclip5.py \
  -n 2 -m 5 \
  -p G \
  in.bam out.bam

5. Paired-end: exact 3-bp soft-clip on R1, emit matching R2 mates
softclip5.py \
  -n 3 -m 3 \
  -p ATG \
  --paired \
  --sort \
  in.bam out.bam
