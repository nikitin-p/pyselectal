# pyselectal

Pyselectal (Python selection of alignments) is a Python script for filtering alignments in the BAM format by the length and sequence of soft-clipped or mapped 5′-end of alignments. It supports single-end and paired-end reads and is designed to be easily integrated into NGS pipelines.

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

`pyselectal` takes as input one or more alignment files (in the BAM, SAM or [CRAM](https://samtools.github.io/hts-specs/CRAMv3.pdf) format) containing local alignments (i.e., with possible soft-clipping). They can be obtained by running, for example, [STAR](https://github.com/alexdobin/STAR) with `--alignEndsType Local` or [HISAT2](https://github.com/DaehwanKimLab/hisat2) / [Bowtie2](https://github.com/BenLangmead/bowtie2) with `--local`.

```bash
pyselectal.py -i FILE[,FILE,...] <mode> [optional arguments]
```

Input is indicated by the mandatory `-i` option and consists of one alignment file or several files separated by commas. `<mode>` is also mandatory and specifies the mode of action (`--select`, `--count`, or `--all`; see below). Exactly one mode must be set. Additionally, you can provide optional arguments, including the name of the output file or directory. See the description of optional arguments below.

You can also use `pyselectal` with a pipe:

```bash
tool_1 | \
  pyselectal.py -i - <mode> [optional arguments] | \
  tool_2 > output.file
```

The dash (`-`) as the argument to `-i` designates reading input alignments from `stdin` (the BAM, SAM or CRAM format is auto-detected by magic bytes; CRAM requires `-r`, see below). When a single input and a single `--select` spec are given without `-o` (see below), output goes to `stdout`.

### Important notes

1. The 3′-end mapping pattern is ignored.
2. Unmapped alignments are never selected.
3. For paired-end input, reads must be name-sorted; use `--name` to sort internally if needed.
4. Multi-mapped reads produce multiple alignments, each processed independently. Consider filtering for primary alignments beforehand (e.g., `samtools view -F 2304` to exclude secondary and supplementary) unless studying multi-mappers.

## Input

`-i, --input FILE[,FILE,...]` — Comma-separated alignment file(s). The file format is auto-detected from the extension (`.bam`, `.sam` or `.cram`). Use `-` for `stdin`. At least one input file is required.

## Modes

### Definitions

Modes are mutually exclusive, and setting exactly one is required:

`-s, --select SPEC[,SPEC,...]` — Select alignments whose 5′ end matches one or more specs (see the [Spec grammar](#spec-grammar) below). With a single input file and a single spec, output goes to `stdout`; otherwise, output files are named `{stem}_{spec}.bam`. Use `--merge` to combine multiple specs into one output file (`{stem}_merged.bam`); see [Optional arguments](#optional-arguments) below.

`-c, --count` — Scan all alignments and print a histogram of all types of 5′ ends, present in input, in a TSV format with columns `type` and `count`. Rows (5' end types) are sorted by decreasing count. Types strictly below `--collapse-threshold` percent are summed up into an `other` type. The output histogram goes to `stdout`, in the case of a single input file, or into `{stem}_5prime_counts.tsv` per input file, for multiple inputs.

`-a, --all` — Write each alignment to a respective 5'-end type-specific file (`{stem}_{type}.bam`). With `-o DIR`, files are placed inside the directory `DIR`. Multiple input files are supported; each input file generates its own set of per-type output files (e.g., `foo_1sg.bam`, `bar_1sg.bam`). Unmapped alignments are silently dropped. Use `--collapse-threshold` to route rare types into a single `{stem}_other` file, instead of individual per-type files.

### Spec grammar

A spec describes a 5′ end type as `[n[..m]]<S|M>[regex]` (case-insensitive):

| Part | Meaning |
| --- | --- |
| `S` | Soft-clipped 5′ end |
| `M` | Mapped 5′ end |
| `n` (before S) | Exactly n soft-clipped bases |
| `n` (before M) | At least n matched bases |
| `n..m` | Between n and m bases (inclusive) |
| `n..` | At least n bases |
| `..m` | At most m bases |
| `regex` | Python regex pattern matched via `re.fullmatch()` |

Reverse-strand alignments are handled automatically (sequence is reverse-complemented before matching).

Examples:

| Spec | Meaning |
| --- | --- |
| `1Sg` | Exactly 1 soft-clipped G |
| `3Sggg` | Exactly 3 soft-clipped GGG |
| `3..Sg+` | 3+ soft-clipped G-homopolymer (via `g+` regex) |
| `2..Sgg.*` | 2+ soft-clipped bases starting with GG (via `gg.*` regex) |
| `1..3S` | 1–3 soft-clipped bases (any sequence) |
| `S` | Any soft-clipped 5′ end |
| `3Mggg` | Exactly 3 matched GGG |
| `3..Mg+` | 3+ matched G-homopolymer |
| `2..3Ma.*` | 2–3 matched bases starting with A (via `a.*` regex) |
| `2..M` | 2+ matched bases |
| `Sg.c` | Soft-clipped gac, ggc, gcc, gtc |
| `M` | Any mapped 5′ end |

## Optional arguments

`-o, --output PATH` — Output file path (for `--select`) or directory (for `--all`). Overrides automatic naming.

`-m, --merge` — With `--select` and multiple specs, write all matches to one output file instead of one file per spec.

`-n, --name` — Internally name-sort the input before processing (required for `--paired` if input is not already name-sorted).

`-t, --threads N` — BGZF compression/decompression threads (default: 1).

`-p, --paired` — Paired-end mode: selection is applied to R1; R2 mates corresponding to selected R1 mates are included automatically if mapped.

`-S, --sam` / `-B, --bam` / `-C, --cram` — Force output format. Default: matches input format. `-C` requires `-r`.

`-r, --reference FASTA` — Reference FASTA for CRAM input or output.

`--mapped-prefix N` — Number of 5′ matched bases to show in `--count` output (default: 5; 0 = length only).

`--collapse-threshold PCT` — Collapse 5′ end types strictly below PCT% of the total number of alignments into an `other` type (`--count`) or `{stem}_other` file (`--all`) (default: 1; 0 = off).

`-h, --help` — Display a full manual and exit.

`-v, --version` — Print the program version and exit.

## Examples

1. Select paired-end alignments where R1 has exactly 1 soft-clipped `G` at the 5′ end; include alignments of the R2 mates corresponding to selected R1s:

```bash
pyselectal.py -i in.bam --select 1Sg --paired --name -o out.bam
```

2. Select single-end alignments with 1–3 soft-clipped bases starting with `G` at the 5′ end:

```bash
pyselectal.py -i in.bam --select 1..3Sg.* -o out.bam
```

3. Select single-end alignments with exactly 1 soft-clipped `A`, `C`, or `T`, or exactly 2 soft-clipped `G`s, or with a mapped 5' end, at the 5′ end and put the selected alignments in the respective output BAM files (one per 5' end type):

```bash
pyselectal.py -i in.bam --select 1Sa,1Sc,1St,2Sgg,M
```

4. Select single-end alignments with either one or two `G` bases soft-clipped at the 5' end and write all selected alignments into one output file:

```bash
pyselectal.py -i in.bam --select 1Sg,2Sgg --merge -o out.bam
```

5. Select single-end alignments with a mapped 5′ end starting with `GG` (single input + single spec outputs to stdout):

```bash
pyselectal.py -i in.bam --select Mgg > out.bam
```

6. Select single-end alignments with at least 10 matched bases at the 5′ end, without a sequence constraint (stdout piped to another tool):

```bash
pyselectal.py -i in.bam --select 10..M | samtools view -c
```

7. Count all 5′ end types present in the input file and generate the corresponding TSV histogram. By default, types accounting for less than 1% of alignments are collapsed into a single "other" row:

```bash
pyselectal.py -i in.bam --count -o counts.tsv
```

8. Write each alignment to a separate file by its 5′ end type and place the output files into `out_dir/`. By default, types accounting for less than 1% of alignments are routed to `in_other.bam`:

```bash
pyselectal.py -i in.bam --all -o out_dir/
```

9. As above, but 5' end types accounting for less than 5% of alignments are written to `out_dir/in_other.bam`, instead of individual files.

```bash
pyselectal.py -i in.bam --all --collapse-threshold 5 -o out_dir/
```

10. Split multiple input files by 5' end type into a shared output directory; each input file generates its own set of type-specific files (e.g., `sample1_1sg.bam`, `sample2_1sg.bam`):

```bash
pyselectal.py -i sample1.bam,sample2.bam --all -o out_dir/
```

11. Pipe a CRAM stream into `pyselectal` (convert BAM to CRAM, then filter):

```bash
samtools view -C -T ref.fa in.bam | pyselectal.py -i - --select 1Sg -r ref.fa -o out.bam
```

## Test data

The repository includes small, synthetic BAM files under `testdata/` that are
designed for **functional testing of `pyselectal`**.  
These files are intended for testing, debugging, and illustrating tool behaviour, not for benchmarking or performance evaluation.

### `test_softclip_se.bam`

Single-end test BAM containing a curated set of alignments with diverse 5′-end configurations:

- Alignments with **5′ soft-clips** of varying lengths (`1S`, `2S`, `3S`, `4S`).
- Alignments with **mapped 5′-ends** and **soft-clipping at the 3′ end**.
- Alignments on both plus and minus strands.
- Homopolymer soft-clips (`G` or `C`), suitable for testing range-based selection.
- Multiple alignments per query, to ensure selection depends only on CIGAR structure and sequence content, and not on mapping multiplicity.

### `test_softclip_pe.bam`

Paired-end test BAM containing alignments grouped by query name and matching various conditions:

- R1 has a 5′ soft-clip.
- R1 has a mapped 5′ end.
- R1 has multiple alternative alignments.
- R1 alignment is associated with more than one R2 alignment.

Test BAM files are small and can be inspected manually with:

```bash
samtools view testdata/test_softclip_se.bam
samtools view testdata/test_softclip_pe.bam
```
