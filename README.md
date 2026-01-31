# pyselectal
Pyselectal (Python selection of alignments) is a Python script for filtering BAM alignments by 5′-end soft-clipping length and sequence, supporting single-end and paired-end reads with exact length and range-based selection. Designed to be easily integrated into NGS processing pipelines.

## Contents

- [Concept and motivation](#concept-and-motivation)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Examples](#examples)
- [Test data](#test-data)

## Concept and motivation

This tool is conceptually inspired by the alignment-level filtering strategy introduced in [Oguchi *et al.*, Science (2024)](https://www.science.org/doi/10.1126/science.add8394), where transcription start sites (TSSs) were inferred from precise 5′-end positions of 5′ single-cell RNA-seq reads. Specifically, Oguchi and colleagues distinguished transcription initiation from other events by the presence of the characteristic 5′ soft-clipped cap-dependent unencoded G base added by the reverse transcriptase during template switching.

Building on this approach, our tool enables general alignment filtering based on 5′-end soft-clipping patterns, mapped 5′-ends, and optional sequence constraints. While method-agnostic, it is particularly useful for CAGE, nanoCAGE, CAGEscan and other 5′-end-focused transcriptomics experiments, including bulk and single-cell protocols, where precise control over the 5′ alignment structure is critical for downstream analyses.

## Requirements

### Python

- **Python ≥ v3.6**
  - Required due to the use of f-strings.
  - **Python ≥ v3.8** is recommended for the best compatibility with modern `pysam` ([link](https://pysam.readthedocs.io/en/latest/release.html)).

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

`pyselectal.py --help`

## Usage

This script is designed to work with BAM files generated using local alignment modes, specifically [STAR](https://github.com/alexdobin/STAR) with the</br>
`--alignEndsType Local` option, or [HISAT2](https://github.com/DaehwanKimLab/hisat2) / [Bowtie2](https://github.com/BenLangmead/bowtie2) with the`--local` option. SAM files are not processed. Only the 5′ end of each alignment is evaluated; any 3′ soft-clipping is ignored.

```bash
pyselectal.py [options] in.bam out.bam
```  
You can also use the script through the pipe:  
```bash
tool_1 | \
  pyselectal.py [options] - - | \
  tool_2 > output.file
```

## Options

`-n, --min-softclip`	Select alignments whose 5′ end has at least `n` soft-clipped bases (required).<br> 
`-m, --max-softclip`	Select alignments whose 5′ end has at most `m` soft-clipped bases (required).<br> 
`-x, --prefix`	Restrict selection to alignments whose 5′ end matches either a specific prefix sequence or a single-base homopolymer, with the exact interpretation depending on the selected mode (see details below).<br> 
`-k, --match`	In mapped-end mode (`-n 0 -m 0`) only, require a minimum number of 5′ matched bases in the CIGAR string or require a minimum prefix length, depending on whether `--prefix` is provided.<br> 
`-s, --sort`	Internally name-sort the input BAM file before processing via pysam.sort.<br> 
`-t, --threads`	Use the specified number of BGZF compression/decompression threads (default: 1).<br> 
`-p, --paired`	Indicate that the input alignments are paired-end reads; by default, the program assumes single-end read alignments.<br> 
`-h, --help`	Display the full manual and exit.<br> 
`-v, --version` Print the program version and exit.<br> 

### Modes of operation

The behaviour is determined by the relationship between `-n` and `-m`. In all modes, unmapped reads are ignored and are never selected.

#### Mode 1 – Exact 5′ soft-clip (n = m > 0)

Select alignments that have an exact n-bp soft-clip at the 5′ end.
If `--prefix` is provided, select only alignments such that the soft-clipped sequence matches the specified prefix of length n.

Plus strand alignments: the 5′ soft-clipped sequence equals prefix.
Minus strand alignments: the 5′ soft-clipped sequence equals reverse complement of prefix.

#### Mode 2 — Mapped 5′-end (n = m = 0)

Select alignments that have a mapped 5′ end (i.e. the 5′ CIGAR operation is M).
This mode supports two selection strategies.

**A.** Length-only matching (`-k` without `--prefix`)
Select alignments that have at least k aligned MATCH bases at the 5′ end (5′ CIGAR M length ≥ k), without checking sequence content.

If k = 0, this reduces to selecting alignments that have a mapped 5′ end (5′ CIGAR operation = M).

**B.** Prefix-based matching (`--prefix`)
Select alignments whose mapped 5′ end matches the full prefix sequence.

Plus strand alignments: the read sequence starts with prefix.
Minus strand alignments: the read sequence ends with reverse complement of prefix.

Interaction with `-k`:

If k ≤ len(prefix) (including k = 0), selection is based only on the full prefix match and `-k` is ignored. A warning message will be shown.

If k > len(prefix), select alignments that satisfy both:
a full prefix match (as above), and
at least k aligned MATCH bases at the 5′ end (5′ CIGAR M length ≥ k).

#### Mode 3 — Soft-clip range (n < m)

Select alignments that have a 5′ soft-clip of length x, where n ≤ x ≤ m.
If `--prefix` is provided, it must be a single base (A, C, G, T, or N).

Plus strand alignments: all 5′ soft-clipped bases equal the specified base.
Minus strand alignments: all 5′ soft-clipped bases equal the complement of that base.

### Paired-end behaviour

When `--paired` is enabled, input alignments are treated as paired-end reads and are grouped by query_name.
Selection is applied only to forward reads (read 1, R1) according to the chosen mode.

For each R1 that passes the selection criteria:

the R1 alignment is written to the output, and
any corresponding reverse read (R2) whose (reference_id, reference_start) matches the mate coordinates recorded in R1 is also written.

Input alignments must be name-sorted; otherwise, internal name sorting must be enabled using --sort.

## Examples
**1. Single-end:** exact 3-bp soft-clip with prefix ATG.
```bash
pyselectal.py \
    -n 3 -m 3 \
    -x ATG \
    in.bam out.bam
```

**2. Single-end:** mapped 5′-end with prefix ATG.
```bash
pyselectal.py \
    -n 0 -m 0 \
    -x ATG \
    in.bam out.bam
```

**3. Single-end:** mapped 5′-end with at least 10 aligned bases (no sequence check).
```bash
pyselectal.py \
    -n 0 -m 0 \
    -k 10 \
    in.bam out.bam
```

**4. Single-end:** 2–5 bp soft-clip, G homopolymer.
```bash
pyselectal.py \
    -n 2 -m 5 \
    -x G \
    in.bam out.bam
```

**5. Paired-end:** exact 3-bp soft-clip on R1, emit matching R2 mates.
```bash
pyselectal.py \
    -n 3 -m 3 \
    -x ATG \
    --paired \
    --sort \
    in.bam out.bam
```

**6. Single-end:** Prefix given, k ignored because k <= len(prefix).
```bash
pyselectal.py \
    -n 0 -m 0 \
    -x ATGC \
    -k 3 \
    in.bam out.bam
```

**7. Single-end:** Prefix given, and additionally require at least 20 mapped bases at 5′-end.
```bash
pyselectal.py \
    -n 0 -m 0 \
    -x ATGC \
    -k 20 \
    in.bam out.bam
```

## Test data

The repository includes small, synthetic BAM files under `testdata/` that are
designed to **test all major modes of `pyselectal`**.  
These files are intended for **functional testing, debugging, and examples**,
not for benchmarking or performance evaluation.

### `test_softclip_se.bam`

Single-end test BAM containing alignments with diverse 5′-end configurations:

- Exact 5′ soft-clips of varying lengths (`1S`, `2S`, `3S`, `4S`)
- Alignments with **mapped 5′-ends** and **3′ soft-clips**
- Plus and minus strand alignments
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