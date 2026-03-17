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

**Important:** At least one of `-n` or `-m` must be specified.<br> 

`-n, --min-softclip`	Select alignments whose 5′ end has at least `n` soft-clipped bases (optional). If missing, n=0 by default.<br> 
`-m, --max-softclip`	Select alignments whose 5′ end has at most `m` soft-clipped bases (optional). If missing, there is no upper bound (i.e., any number of soft-clipped bases above the minimum is allowed).<br> 
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

When `--paired` is enabled, input alignments are treated as paired-end reads and are grouped by `query_name`.
Selection is applied only to forward reads (read 1, R1) according to the chosen mode.

For each R1 that passes the selection criteria:

the R1 alignment is written to the output, and
any corresponding reverse read (R2) whose (`reference_id`, `reference_start`) matches the mate coordinates recorded in R1 is also written.

Input alignments must be name-sorted; otherwise, internal name sorting must be enabled using `--sort`.

## Examples
**1. Single-end CAGE ([Murata *et al.*, 2014](https://link.springer.com/protocol/10.1007/978-1-4939-0805-9_7)).** Select alignments that have an exact 1-bp soft-clip at the 5′-end with prefix G.
```bash
pyselectal.py \
    -n 1 -m 1 \
    -x G \
    in.bam out.bam
```

**2. CAGEscan ([Bertin *et al.*, 2011](https://onlinelibrary.wiley.com/doi/abs/10.1002/9783527644582.ch3)).** Select paired-end alignments such that R1 has an exact 3-bp soft-clip at the 5′-end with prefix GGG, and include the corresponding R2 mates.
```bash
pyselectal.py \
    -n 3 -m 3 \
    -x GGG \
    --paired \
    --sort \
    in.bam out.bam
```

**3. Paired-end CAGEscan.** Select paired-end alignments such that R1 has a 5′ soft-clip of length 1–3 bp of a G homopolymer, and include the corresponding R2 mates.
```bash
pyselectal.py \
    -n 1 -m 3 \
    -x G \
    --paired \
    --sort \
    in.bam out.bam
```


**4. Single-end CAGE.** Select alignments that have an exact 1-bp soft-clip at the 5′-end with a prefix A that may occasionally correspond to the canonical m<sup>7</sup>G cap ([Ohtake *et al.*, 2004](https://academic.oup.com/dnaresearch/article/11/4/305/336335)).
```bash
pyselectal.py \
    -n 1 -m 1 \
    -x A \
    in.bam out.bam
```

**5. Single-end.** Select alignments that have a mapped 5′-end with prefix GG.
```bash
pyselectal.py \
    -n 0 -m 0 \
    -x GG \
    in.bam out.bam
```

**6. Single-end.** Select alignments that have at least 10 aligned bases at the 5′-end, without applying any sequence constraint.
```bash
pyselectal.py \
    -n 0 -m 0 \
    -k 10 \
    in.bam out.bam
```

**7. Single-end.** Select alignments that have a 5′ soft-clip of length 2–5 bp composed of a G homopolymer.
```bash
pyselectal.py \
    -n 2 -m 5 \
    -x G \
    in.bam out.bam
```

**8. Single-end.** Select alignments that have the specified prefix at the 5′-end and at least 20 aligned bases at the 5′-end.
```bash
pyselectal.py \
    -n 0 -m 0 \
    -x ATGC \
    -k 20 \
    in.bam out.bam
```

**9. Single-end.** Select alignments that have the specified prefix at the 5′ end; the -k parameter is ignored because k is shorter than the prefix length.
```bash
pyselectal.py \
    -n 0 -m 0 \
    -x ATGC \
    -k 3 \
    in.bam out.bam
```

**10. Paired-end.** Select paired-end alignments such that R1 has an exact 1-bp soft-clip at the 5′ end with prefix G, and include the corresponding R2 mates.
```bash
pyselectal.py \
    -n 1 -m 1 \
    -x G \
    --paired \
    --sort \
    in.bam out.bam
```

## Test data

The repository includes small, synthetic BAM files under `testdata/` that are
designed for **functional testing of `pyselectal`**.  
These files are intended for testing, debugging, and illustrating tool behaviour, not for benchmarking or performance evaluation.

### `test_softclip_se.bam`

Single-end test BAM containing a curated set of alignments with diverse 5′-end configurations:

- Alignments with exact 5′ soft-clips of varying lengths (`1S`, `2S`, `3S`, `4S`)
- Alignments with **mapped 5′-ends** and **soft-clipping at the 3′ end**
- Alignments on both plus and minus strands
- Homopolymer soft-clips (`G` / `C`), suitable for testing range-based selection
- Multiple alignments per query, to ensure selection depends only on CIGAR structure and sequence content, and not on mapping multiplicity

### `test_softclip_pe.bam`

Paired-end test BAM designed specifically for paired-end mode (`--paired`):
- Alignments are grouped by query name, representing paired-end fragments
- For each fragment, selection criteria are evaluated exclusively on read1 (R1)
- Corresponding read2 (R2) alignments are fully mapped and included only if their R1 counterpart satisfies the selection criteria

The file includes explicit paired-end scenarios such as:
- Fragments where R1 has an exact 5′ soft-clip
- Fragments where R1 has a 5′ soft-clip within a specified length range
- Fragments where R1 has a mapped 5′ end matching a given prefix
- Fragments where R1 has multiple alternative alignments
- Fragments where a selected R1 alignment is associated with more than one valid R2 alignment, testing correct mate inclusion

Test BAMs are deliberately small and can be manually inspected with:
```bash
samtools view testdata/test_softclip_se.bam
samtools view testdata/test_softclip_pe.bam
```

We also provide test data to help users verify the installation and functionality of the pipeline. The repository includes three small FASTQ files that can be used to test the mapping step and the overall workflow of the script, including the filtering of the resulting BAM file. The subset data are derived from the human LUHMES cell line at 1, 3, and 6 days of neuronal differentiation ([Yoshihara *et al.*, 2025](https://link.springer.com/article/10.1038/s44319-025-00372-1)).

Example usage:
```bash
STAR \
  --runThreadN 20 \
  --genomeDir star_index \
  --readFilesIn testdata/luhmesDay1Rep1_NETCAGE_subsample_L001_R1_001.fastq.gz \
  --outSAMtype BAM Unsorted \
  --readFilesCommand gunzip -c \
  --alignEndsType Local \
  --outSAMunmapped Within | \
    pyselectal.py [options] - \
    luhmesDay1Rep1_NETCAGE_subsample_filtered_L001_R1_001.fastq.gz
```