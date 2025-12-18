# pyselectal
Python script for filtering BAM alignments by 5′-end soft-clipping length and sequence, supporting single-end and paired-end reads with exact length and range-based selection. Designed to be easily integrated into NGS processing pipelines.

## Requirements

### Python

- **Python ≥ v3.6**
  - Required due to use of f-strings.
  - Python ≥ v3.8 recommended for best compatibility with modern `pysam`.

### Python libraries

- **pysam ≥ 0.15.0**

`pysam 0.15.0` is the earliest version that supports the `threads=` argument in
`pysam.AlignmentFile`, which is used for parallel BGZF compression/decompression.

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

