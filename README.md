# pyselectal

Pyselectal (Python selection of alignments) is a Python script for filtering alignments in the BAM format by the length and sequence of soft-clipped or mapped 5′-end of reads. It supports single-end and paired-end reads and is designed to be easily integrated into NGS pipelines.

## Contents

- [Concept and motivation](#concept-and-motivation)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Input](#input)
- [Modes](#modes)
- [Optional arguments](#optional-arguments)
- [Examples](#examples)
- [Test data](#test-data)

## Concept and motivation

Pyselectal is conceptually inspired by the alignment filtering strategy described in [Oguchi *et al.*, 2024](https://www.science.org/doi/10.1126/science.add8394), where transcription start sites (TSSs) were inferred from the precise 5′-end positions of 5′ single-cell RNA-seq reads. Specifically, Oguchi and colleagues distinguished transcription initiation from other events by the presence of the characteristic 5′ soft-clipped cap-dependent unencoded base `G` added by the reverse transcriptase during template switching.

Building on this approach, our tool enables general alignment filtering based on 5′-end soft-clipping patterns, mapped 5′ ends and optional sequence constraints. While sequencing method-agnostic, it is particularly useful for CAGE, nanoCAGE, CAGEscan and other 5′-end-focused RNA sequencing experiments, including bulk and single-cell protocols, where precise control over the structure of the 5′-end alignment is critical for downstream analyses.

## Requirements

### Python

- **Python ≥ v3.8** is recommended ([for the best compatibility with modern `pysam`](https://pysam.readthedocs.io/en/latest/release.html)).
  - **Python ≥ v3.6** is required due to the use of [f-strings](https://peps.python.org/pep-0498/) in the `pyselectal` implementation.

### Python libraries

- **pysam ≥ 0.15.0**

`pysam 0.15.0` is the earliest version that supports the `threads` argument in `pysam.AlignmentFile` which is used for parallel BGZF compression/decompression. You can install `pysam` using `pip`:

```bash
pip install pysam
```

## Installation

Clone the repository and make the `pyselectal.py` script executable:

```bash
git clone https://github.com/nikitin-p/pyselectal.git
cd pyselectal
chmod +x pyselectal.py
```

Then, you can run it directly:

`pyselectal.py --help`

## Usage

`pyselectal` takes one or more alignment files (in the BAM, SAM or [CRAM](https://samtools.github.io/hts-specs/CRAMv3.pdf) format) containing local alignments (i.e., with possible soft-clipping). They can be obtained by running, for example, [STAR](https://github.com/alexdobin/STAR) with `--alignEndsType Local` or [HISAT2](https://github.com/DaehwanKimLab/hisat2) / [Bowtie2](https://github.com/BenLangmead/bowtie2) with `--local`. The 3′-end mapping pattern is ignored. Unmapped reads are never selected.

```bash
pyselectal.py -i FILE[,FILE,...] <mode> [optional arguments]
```

Exactly one mode must be specified: `--select`, `--count`, or `--all`.

For paired-end input, reads must be name-sorted; use `--name` to sort internally if needed.

You can also use `pyselectal` with a pipe:

```bash
tool_1 | \
  pyselectal.py -i - <mode> [optional arguments] | \
  tool_2 > output.file
```

The dash (`-`) as the argument to `-i` designates reading input alignments from `stdin` (the BAM, SAM or CRAM format is auto-detected by magic bytes; CRAM requires `-r`, see below). When a single input and a single `--select` spec are given without `-o` (see below), output goes to `stdout`.

## Input

`-i, --input FILE[,FILE,...]` — Input alignment file(s). The format is auto-detected from the extension (`.bam`, `.sam` or `.cram`). Use `-` for `stdin`. Required.

## Modes

### Definitions

Modes are mutually exclusive, and exactly one is required:

`-s, --select SPEC[,SPEC,...]` — Select alignments whose 5′ end matches one or more specs (see the [Spec grammar](#spec-grammar) below). With a single input file and a single spec, output goes to `stdout`; otherwise, output files are named `{stem}_{spec}.bam`. Use `--merge` to combine multiple specs into one output file (`{stem}_merged.bam`).

`-c, --count` — Scan all alignments and print a histogram of all types of 5′ ends, present in input, in a TSV format with columns `type` and `count`. Rows are sorted by descending count. Categories strictly below `--collapse-threshold` are summed up into the `other` category. The output histogram goes to `stdout` in the case of a single input file, or into `{stem}_5prime_counts.tsv` per input file, for multiple inputs.

`-a, --all` — Write each alignment to a respective 5'-end type-specific file (`{stem}_{type}.bam`). With `-o DIR`, files are placed inside the directory `DIR`. Unmapped reads are silently dropped. Use `--collapse-threshold` to route rare types into a single `{stem}_other` file, instead of individual per-type files.

### Spec grammar

A spec describes a 5′-end type as `[n[.m]]<S|M>[seq]` (case-insensitive):

| Part | Meaning |
| --- | --- |
| `S` | Soft-clipped 5′ end |
| `M` | Mapped 5′ end |
| `n` (before S) | Exactly n soft-clipped bases |
| `n` (before M) | At least n matched bases |
| `n.m` | Between n and m bases (inclusive) |
| `n.` | At least n bases |
| `.m` | At most m bases |
| `seq` | Required sequence at the 5′ end of the read |

Reverse-strand reads are handled automatically (sequence is reverse-complemented before matching).

Examples:

| Spec | Meaning |
| --- | --- |
| `1Sg` | Exactly 1 soft-clipped G |
| `2.Sgg` | 2 or more soft-clipped Gs |
| `1.3S` | 1–3 soft-clipped bases (any sequence) |
| `S` | Any soft-clipped 5′ end |
| `Mg` | Any mapped 5′ end starting with G |
| `2.3MA` | 2–3 matched bases starting with A |
| `M` | Any mapped 5′ end |

## Optional arguments

`-o, --output PATH` — Output file path (for `--select`) or directory (for `--all`). Overrides automatic naming.

`-m, --merge` — With `--select` and multiple specs, write all matches to one output file instead of one file per spec.

`-n, --name` — Internally name-sort the input before processing (required for `--paired` if input is not already name-sorted).

`-t, --threads N` — BGZF compression/decompression threads (default: 1).

`-p, --paired` — Paired-end mode: selection is applied to R1; matching R2 mates are included automatically.

`-S, --sam` / `-B, --bam` / `-C, --cram` — Force output format. Default: matches input format. `-C` requires `-r`.

`-r, --reference FASTA` — Reference FASTA for CRAM input or output.

`--mapped-prefix N` — Number of 5′ matched bases to show in `--count` output (default: 5; 0 = length only).

`--collapse-threshold PCT` — Collapse 5′-end type categories strictly below PCT% of total into an `other` row (`--count`) or `{stem}_other` file (`--all`) (default: 1; 0 = off).

`-h, --help` — Display a full manual and exit.

`-v, --version` — Print the program version and exit.

## Examples

**1. Single-end CAGE ([Murata *et al.*, 2014](https://link.springer.com/protocol/10.1007/978-1-4939-0805-9_7)).** Select alignments with exactly 1 soft-clipped G at the 5′ end.

```bash
pyselectal.py -i in.bam --select 1Sg -o out.bam
```

**2. CAGEscan ([Bertin *et al.*, 2011](https://onlinelibrary.wiley.com/doi/abs/10.1002/9783527644582.ch3)).** Select paired-end alignments where R1 has exactly 3 soft-clipped Gs at the 5′ end; include matching R2 mates.

```bash
pyselectal.py -i in.bam --select 3Sggg --paired --name -o out.bam
```

**3. Paired-end CAGEscan.** Select paired-end alignments where R1 has 1–3 soft-clipped Gs at the 5′ end.

```bash
pyselectal.py -i in.bam --select 1.3Sg --paired --name -o out.bam
```

**4. Single-end CAGE.** Select alignments with exactly 1 soft-clipped A at the 5′ end (may correspond to the canonical m7G cap ([Ohtake *et al.*, 2004](https://academic.oup.com/dnaresearch/article/11/4/305/336335))).

```bash
pyselectal.py -i in.bam --select 1Sa -o out.bam
```

**5. Single-end.** Select alignments with a mapped 5′ end starting with GG.

```bash
pyselectal.py -i in.bam --select Mgg -o out.bam
```

**6. Single-end.** Select alignments with at least 10 matched bases at the 5′ end, without sequence constraint.

```bash
pyselectal.py -i in.bam --select 10.M -o out.bam
```

**7. Single-end.** Select alignments with 2–5 soft-clipped Gs at the 5′ end.

```bash
pyselectal.py -i in.bam --select 2.5Sg -o out.bam
```

**8. Merge multiple specs.** Select alignments matching either `1Sg` or `2Sgg` and write all to one file.

```bash
pyselectal.py -i in.bam --select 1Sg,2Sgg --merge -o out.bam
```

**9. Count 5′-end types.** Output a TSV histogram of all 5′-end types.

```bash
pyselectal.py -i in.bam --count -o counts.tsv
```

**10. Split by type.** Write each alignment to a separate file by 5′-end type, placed in `out_dir/`.

```bash
pyselectal.py -i in.bam --all -o out_dir/
```

**11. Split by type with rare-type collapsing.** As above, but types accounting for less than 5% of reads are written to `out_dir/in_other.bam` instead of individual files.

```bash
pyselectal.py -i in.bam --all --collapse-threshold 5 -o out_dir/
```

**12. CRAM from stdin.** Pipe a CRAM stream into `pyselectal` with a reference.

```bash
cat in.cram | pyselectal.py -i - --select 1Sg -r ref.fa -o out.bam
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

You can also verify the full mapping-to-filtering workflow using your own data. For example, piping STAR output directly into `pyselectal`:

```bash
STAR \
  --runThreadN 20 \
  --genomeDir star_index \
  --readFilesIn reads_R1.fastq.gz \
  --outSAMtype BAM Unsorted \
  --readFilesCommand gunzip -c \
  --alignEndsType Local \
  --outSAMunmapped Within | \
    pyselectal.py -i - --select 1Sg -o filtered.bam
```
