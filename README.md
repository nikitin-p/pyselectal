# pyselectal

Pyselectal (Python selection of alignments) is a Python script for filtering alignments in BAM, SAM and CRAM formats by the length and sequence of soft-clipped or mapped 5′-end of alignments. It supports single-end and paired-end reads and is designed to be easily integrated into NGS pipelines. You can then **select** reads matching a 5'-end pattern, **count** the distribution of 5'-end types, or **split** an alignment file into one file per type.

## Contents

- [Concept and motivation](#concept-and-motivation)
- [Requirements](#requirements)
- [Installation](#installation)
- [Usage](#usage)
- [Options](#options)
- [Examples](#examples)
- [Testing](#testing)
- [Resource usage](#resource-usage)
- [Test data](#test-data)
- [Citation](#citation)

## Concept and motivation

Pyselectal is conceptually inspired by the alignment filtering strategy described in [Oguchi *et al.*, 2024](https://www.science.org/doi/10.1126/science.add8394), where transcription start sites (TSSs) were inferred from the precise 5′-end positions of 5′ single-cell RNA-seq reads. Specifically, Oguchi and colleagues distinguished transcription initiation from other events by the presence of the characteristic 5′ soft-clipped cap-dependent unencoded base `G` added by the reverse transcriptase during template switching.

Building on this approach, our tool enables general alignment filtering based on 5′-end soft-clipping patterns, mapped 5′ ends and optional sequence constraints. While sequencing method-agnostic, it is particularly useful for CAGE ([Shiraki *et al.*, 2003](https://pubmed.ncbi.nlm.nih.gov/14663149/); [Murata *et al.*, 2014](https://pubmed.ncbi.nlm.nih.gov/24927836/)), nanoCAGE ([Salimullah *et al.*, 2011](https://pmc.ncbi.nlm.nih.gov/articles/PMC4181851/)), CAGEscan ([Bertin *et al.*, 2017](https://pubmed.ncbi.nlm.nih.gov/28972578/)) and other 5′-end-focused RNA sequencing experiments, including bulk and single-cell protocols, where precise control over the structure of the 5′-end alignment is critical for downstream analyses.

## Requirements

### Python

- **Python ≥ v3.8** is recommended ([for the best compatibility with modern *pysam*](https://pysam.readthedocs.io/en/latest/release.html)).
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

You can also read from `stdin`:

```bash
tool_1 | pyselectal.py -i - <mode> [optional arguments]
```

The dash (`-`) as the argument to `-i` designates reading input alignments from `stdin` (the BAM, SAM or CRAM format is auto-detected by magic bytes ([SAM/BAM](https://samtools.github.io/hts-specs/SAMv1.pdf), [CRAM](https://samtools.github.io/hts-specs/CRAMv3.pdf)); CRAM requires `-r`, see below). Output always goes to named files on disk.

```
-i, --input FILE[,FILE,...]   input file(s)
-s, --select SPEC[,SPEC,...]  select by 5′-end spec
-c, --count                   count 5′-end types
-a, --all                     split by 5′-end type
-o, --output PATH             output file or directory
-m, --merge                   merge specs into one output
-u, --unmatched               write unmatched alignments
-f, --unmatched-file FILE     custom unmatched filename
-x, --exclude                 invert selection
-p, --paired                  paired-end mode
-n, --no-name-sort            skip name sorting
-t, --threads N               BGZF threads
-S, --sam                     force SAM output
-B, --bam                     force BAM output
-C, --cram                    force CRAM output
-r, --reference FASTA         reference for CRAM
--mapped-prefix N             matched bases in --count
--collapse-threshold PCT      collapse rare types
-h, --help                    show help
-v, --version                 show version
```

### Important notes

1. The 3′-end mapping pattern is ignored.
2. Unmapped alignments are never selected.
3. For paired-end input, reads must be name-sorted; pyselectal name-sorts internally by default (use `--no-name-sort` to skip if input is already sorted).
4. Multi-mapped reads produce multiple alignments, each processed independently. Consider filtering for primary alignments beforehand (e.g., `samtools view -F 2304` to exclude secondary and supplementary) unless studying multi-mappers.

## Options

### Definitions

A **5′ type** is a string that describes the exact 5′-end structure of an alignment, such as `Sg` (1 soft-clipped G), `Sgg` (2 soft-clipped Gs), or `Mgaggg` (5 matched bases GAGGG). A **5′ spec** is a pattern used to select alignments by their 5′ type (see the [Spec grammar](#spec-grammar) below).

### Input

`-i, --input FILE[,FILE,...]` — Comma-separated alignment file(s). The file format is auto-detected from the extension (`.bam`, `.sam` or `.cram`). Use `-` for `stdin`. At least one input file is required.

### Mode (exactly one required)

`-s, --select SPEC[,SPEC,...]` — Select alignments whose 5′ end matches one or more specs (see the [Spec grammar](#spec-grammar) below). Output files are named `{stem}_{spec}{ext}`. An alignment matching multiple specs is written to each corresponding output file (e.g., a `3Sg` alignment matches both `2..5Sg` and `2..6Sg` and appears in both output files). Use `--merge` to combine multiple specs into one output file per input file (`{stem}_merged{ext}`); each alignment is written only once, even if it matches multiple specs. Use `--unmatched` to additionally write alignments that do not match any spec to `{stem}_unmatched{ext}`. Use `--exclude` to invert the selection: output only alignments that do not match any spec (the complement of `--select`). `--exclude` and `--unmatched` are mutually exclusive; `--exclude` and `--merge` are mutually exclusive.

`-c, --count` — Scan all alignments and write a histogram of all types of 5′ ends, present in input, in a TSV format with columns `type` and `count` to `{stem}_5prime_counts.tsv`. Rows (5′ end types) are sorted by decreasing count. Types strictly below `--collapse-threshold` percent are summed up into `other_soft_clipped` and `other_mapped` rows.

`-a, --all` — Write each alignment to a respective 5'-end type-specific file (`{stem}_{type}.bam`). With `-o DIR`, files are placed inside the directory `DIR`. Multiple input files are supported; each input file generates its own set of per-type output files (e.g., `foo_1sg.bam`, `bar_1sg.bam`). Unmapped alignments are silently dropped. Use `--collapse-threshold` to route rare types into `{stem}_other_soft_clipped.bam` (soft-clipped) and `{stem}_other_mapped.bam` (mapped) instead of individual per-type files.

### Output

By default, the output format matches the input format; with multiple input files in different formats, each output inherits its corresponding input's format. Use `-S` (SAM), `-B` (BAM), or `-C` (CRAM) to force a single format for all outputs; `-C` requires `-r/--reference`. The output file extension is adjusted accordingly (e.g., forcing SAM output produces `{stem}_{spec}.sam` instead of `.bam`).

In the output naming patterns below, `{stem}` refers to the input filename without its extension (e.g., `sample.bam` → `sample`), and `{ext}` refers to the output file extension (`.bam`, `.sam`, or `.cram`).

Output file naming depends on the mode:

- `--select`: `{stem}_{spec}{ext}` per spec;
- `--select --unmatched`: additionally writes `{stem}_unmatched{ext}`;
- `--select --merge`: `{stem}_merged{ext}` per input file;
- `--select --merge --unmatched`: additionally writes `{stem}_unmatched{ext}` per input file;
- `--select --exclude`: `{stem}_excluded{ext}` per input file;
- `--count`: `{stem}_5prime_counts.tsv` per input file;
- `--all`: `{stem}_{type}{ext}` per 5′ end type, per input file.

`-o, --output PATH` — Output file path (for `--select` and `--count`) or directory (for `--all`). Overrides automatic naming.

`-S, --sam` / `-B, --bam` / `-C, --cram` — Force output format. Default: matches input format. `-C` requires `-r`.

`-r, --reference FASTA` — Reference FASTA for CRAM input or output.

### Selection modifiers

`-m, --merge` — With `--select` and multiple specs, write all matches to one output file instead of one file per spec.

`-u, --unmatched` — Write alignments that do not match any spec to a separate auto-named file (`{stem}_unmatched{ext}`). Only valid with `--select`. Not affected by `-o`. Incompatible with `-x/--exclude`.

`-f, --unmatched-file FILE` — Write non-matching alignments to FILE instead of auto-naming. Implies `-u/--unmatched`.

`-x, --exclude` — Invert the `--select` logic: output every alignment that matches *none* of the given specs (analogous to `grep -v`). Output is written to `{stem}_excluded{ext}` per input file. Only valid with `--select`. Incompatible with `-m/--merge` and `-u/--unmatched`.

### Processing

`-n, --no-name-sort` — Skip internal name sorting (use if input is already name-sorted). By default, pyselectal name-sorts the input before processing.

`-t, --threads N` — BGZF compression/decompression threads (default: 1).

`-p, --paired` — Paired-end mode: selection is applied to R1; R2 mates corresponding to selected R1 mates are included automatically if mapped.

`--mapped-prefix N` — Number of 5′ matched bases to show in `--count` output (default: 5; 0 = length only).

`--collapse-threshold PCT` — Collapse 5′ end types strictly below PCT% of the total number of alignments into `other_soft_clipped` and `other_mapped` rows (`--count`) or into `{stem}_other_soft_clipped.bam` / `{stem}_other_mapped.bam` (`--all`) (default: 1; 0 = off).

### Informational

`-h, --help` — Display a full manual and exit.

`-v, --version` — Print the program version and exit.

### Spec grammar

A spec describes a 5′ end type as `[n|n..|..m|n..m]<S|M>[pattern]` (case-insensitive; `<>` = required):

| Part | Meaning |
| --- | --- |
| `S` | Soft-clipped 5′ end |
| `M` | Mapped 5′ end |
| `n` (before S or M) | Exactly n bases |
| `n..m` | Between n and m bases (inclusive) |
| `n..` | At least n bases |
| `..m` | At most m bases |
| `pattern` | Sequence pattern (literal or [Python regex](https://docs.python.org/3/library/re.html#regular-expression-syntax), matched via `re.fullmatch()`) |

The pattern can be a literal sequence (`g`, `ttc`) or use regex syntax (`g+`, `g.c`, `[acgt]+`). The entire extracted 5′ sequence must match the pattern (fullmatch semantics). Reverse-strand alignments are handled automatically (the sequence is reverse-complemented before matching).

If the pattern is a literal (no regex metacharacters), the length is inferred from the pattern itself — `Sg` equals `1Sg`, `Sttc` equals `3Sttc`. When the pattern contains regex metacharacters (`+`, `*`, `.`, `[]`, etc.), an explicit length quantifier controls the match range.

Examples:

| Spec | Meaning |
| --- | --- |
| `S` | Any soft-clipped 5′ end |
| `M` | Any mapped 5′ end |
| `2S` | Any 2 soft-clipped bases |
| `3M` | Any 3 matched bases |
| `2..5S` | 2–5 soft-clipped bases |
| `3..M` | 3+ matched bases |
| `..4S` | Up to 4 soft-clipped bases (inclusive) |
| `Sg` | Exactly 1 soft-clipped G (equals `1Sg`) |
| `Ma` | Exactly 1 matched A (equals `1Ma`) |
| `Sttc` | Exactly 3 soft-clipped TTC (equals `3Sttc`) |
| `Maat` | Exactly 3 matched AAT (equals `3Maat`) |
| `2Ma` | Empty — pattern `a` is 1 base but 2 required |
| `2..3Sa` | Empty — pattern `a` is 1 base but 2–3 required |
| `2..3Saa` | Equals `Saa` — quantifier redundant with fixed-length literal |
| `Sg+` | Soft-clipped G-homopolymer of any length |
| `..5Sg+` | Soft-clipped G-homopolymer of length 1–5 |
| `Mac+t` | Matched 5′ end matching `ac+t`: act, acct, accct, … |
| `..4Mac+t` | Only act and acct (length ≤ 4) |
| `Sg.c` | Soft-clipped gac, ggc, gcc, gtc (`.` = any char) |
| `4Sg.c` | Empty — regex `g.c` matches length 3, but 4 required |
| `Sg[ga]c` | Soft-clipped ggc or gac (`[ga]` = character class) |
| `3..4Sca+g` | Soft-clipped cag or caag (length 3–4 from `ca+g`) |
| `4Sca+g` | Only caag (exactly 4 bases matching `ca+g`) |

## Examples

1. **SE CAGE:** Select single-end alignments with exactly 1 soft-clipped `G` at the 5′ end (cap-dependent unencoded G):

```bash
pyselectal.py -i in.bam --select Sg -o out.bam
```

2. **PE CAGEscan:** Select paired-end alignments where R1 has exactly 3 soft-clipped `G`s at the 5′ end; include R2 mates:

```bash
pyselectal.py -i in.bam --select Sggg --paired -o out.bam
```

3. Select single-end alignments with exactly 1 soft-clipped `A`, `C`, or `T`, or exactly 2 soft-clipped `G`s, or with a mapped 5′ end, and put the selected alignments in the respective output BAM files (one per 5′ end type):

```bash
pyselectal.py -i in.bam --select Sa,Sc,St,Sgg,M
```

4. Select single-end alignments with either one or two `G` bases soft-clipped at the 5' end and write all selected alignments into one output file:

```bash
pyselectal.py -i in.bam --select Sg,Sgg --merge -o out.bam
```

5. Select single-end alignments with a mapped 5′ end starting with `GG`:

```bash
pyselectal.py -i in.bam --select Mgg
# output: in_mgg.bam
```

6. Select single-end alignments with at least 10 matched bases at the 5′ end, without a sequence constraint, and count them:

```bash
pyselectal.py -i in.bam --select 10..M -o out.bam && samtools view -c out.bam
```

7. Select alignments with 1 soft-clipped G and save unmatched to a separate file:

```bash
pyselectal.py -i in.bam -s Sg -u
# output: in_sg.bam (matched), in_unmatched.bam (unmatched)
```

8. Same as above, but specify a custom filename for unmatched alignments:

```bash
pyselectal.py -i in.bam -s Sg -f rejected.bam
# output: in_sg.bam (matched), rejected.bam (unmatched)
```

9. Exclude single-end alignments with exactly 1 soft-clipped `G` at the 5′ end:

```bash
pyselectal.py -i in.bam -s Sg -x
# output: in_excluded.bam
```

10. Exclude all soft-clipped reads across multiple specs (the complement of all G-cap variants) and write the result to a named file:

```bash
pyselectal.py -i in.bam -s Sg,Sgg,Sggg -x -o non_gcap.bam
```

11. Count all 5′ end types present in the input file and generate the corresponding TSV histogram. By default, types accounting for less than 1% of alignments are collapsed into a single "other" row:

```bash
pyselectal.py -i in.bam --count -o counts.tsv
```

12. Write each alignment to a separate file by its 5′ end type and place the output files into `out_dir/`. By default, types accounting for less than 1% of alignments are routed to `in_other_soft_clipped.bam` or `in_other_mapped.bam`:

```bash
pyselectal.py -i in.bam --all -o out_dir/
```

13. As above, but 5' end types accounting for less than 5% of alignments are written to `out_dir/in_other_soft_clipped.bam` or `out_dir/in_other_mapped.bam`, instead of individual files.

```bash
pyselectal.py -i in.bam --all --collapse-threshold 5 -o out_dir/
```

14. Split multiple input files by 5' end type into a shared output directory; each input file generates its own set of type-specific files (e.g., `sample1_1sg.bam`, `sample2_1sg.bam`):

```bash
pyselectal.py -i sample1.bam,sample2.bam --all -o out_dir/
```

15. Pipe a CRAM stream into `pyselectal` (convert BAM to CRAM, then filter):

```bash
samtools view -C -T ref.fa in.bam | pyselectal.py -i - --select Sg -r ref.fa -o out.bam
```

## Testing

The project uses [pytest](https://docs.pytest.org/) for automated testing. To run the full test suite:

```bash
python -m pytest test_pyselectal.py -v
```

All tests are contained in `test_pyselectal.py` and cover spec parsing, alignment matching, SE/PE processing, and end-to-end mode execution.

## Resource usage

Benchmarks on three BAM files with 5 million alignments each (VM with 20 CPUs, 40 GB RAM):

| Mode | Time (avg) | Memory |
| --- | --- | --- |
| `--select` | ~37 s | ~950 MB |
| `--count` | ~43 s | ~950 MB |
| `--all` | ~68 s | ~950 MB |

Notes:

- **SELECT** is the fastest mode, processing ~135,000 alignments per second;
- **COUNT** is slightly slower due to per-alignment type classification;
- **ALL** takes roughly twice as long because it performs two passes over the input (first to count types for `--collapse-threshold`, then to route alignments);
- Memory usage is stable regardless of mode and input size, since pyselectal processes alignments in a streaming fashion (PE mode buffers one read group at a time).

## Test data

The repository includes small, synthetic BAM files under `testdata/` that are
designed for **functional testing of `pyselectal`**.  
These files are intended for testing, debugging, and illustrating tool behaviour, not for benchmarking or performance evaluation.

### `test_softclip_se.bam`

Single-end test BAM containing a curated set of alignments with diverse 5′-end configurations:

- Alignments with **5′ soft-clips** of varying lengths (`1S`, `2S`, `3S`, `4S`);
- Alignments with **mapped 5′-ends** and **soft-clipping at the 3′ end**;
- Alignments on both plus and minus strands;
- Homopolymer soft-clips (`G` or `C`), suitable for testing range-based selection;
- Multiple alignments per query, to ensure selection depends only on CIGAR structure and sequence content, and not on mapping multiplicity.

### `test_softclip_pe.bam`

Paired-end test BAM containing alignments grouped by query name and matching various conditions:

- R1 has a 5′ soft-clip;
- R1 has a mapped 5′ end;
- R1 has multiple alternative alignments;
- R1 alignment is associated with more than one R2 alignment.

Test BAM files are small and can be inspected manually with:

```bash
samtools view testdata/test_softclip_se.bam
samtools view testdata/test_softclip_pe.bam
```

## Citation

If you use pyselectal in your research, please cite:

> Nikitin P, Sidorov S. pyselectal: Python selection of alignments by 5′-end type. 2026. https://github.com/nikitin-p/pyselectal

BibTeX:

```bibtex
@software{pyselectal,
  author       = {Nikitin, Pavel and Sidorov, Sviatoslav},
  title        = {pyselectal: Python selection of alignments by 5′-end type},
  year         = {2026},
  url          = {https://github.com/nikitin-p/pyselectal},
  note         = {Version 1.0}
}
```
