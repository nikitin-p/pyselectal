#!/usr/bin/env python3

VERSION = "pyselectal v2.0"

MANUAL = """\
pyselectal — filter and analyze alignments by 5'-end type

Complete manual: https://github.com/nikitin-p/pyselectal

Modes (mutually exclusive, exactly one required):
  -s, --select SPEC[,SPEC,...]   Select alignments by 5' end type
  -c, --count                    Output TSV histogram of 5' end types
  -a, --all                      Split alignments into per-type BAM files

Options:
  -i, --input FILE[,FILE,...]    Input file(s), comma-separated; format auto-detected from extension (required)
  -o, --output PATH              Output file path / directory (overrides automatic naming)
  -m, --merge                    Merge multiple --select specs into one output
  -n, --name                     Name-sort input internally
  -t, --threads N                BGZF threads (default: 1)
  -p, --paired                   Paired-end mode
  -S, --sam                      Force SAM output
  -B, --bam                      Force BAM output
  -C, --cram                     Force CRAM output (requires -r)
  -r, --reference FASTA          Reference for CRAM input/output
  -h, --help                     Show this manual and exit
  -v, --version                  Print version and exit

5' end type notation (for --select):
  [n[.m]]<S|M>[seq]    case-insensitive
  Examples: 1Sg, 2.4S, .5M, 3Saac, M, St

"""

import argparse
import os
import re
import sys
import tempfile
import pysam

SOFT = 4   # 'S' soft-clip in pysam CIGAR codes
MATCH = 0  # 'M' (alignment match) in pysam CIGAR codes

_WARNED = set()

class HelpfulArgumentParser(argparse.ArgumentParser):
    def error(self, message):
        sys.stderr.write(f"Error: {message}\nUse -h for help.\n")
        raise SystemExit(2)

def die(message):
    sys.stderr.write(f"Error: {message}\nUse -h for help.\n")
    raise SystemExit(2)

def warn_once(key: str, message: str):
    if key in _WARNED:
        return
    _WARNED.add(key)
    sys.stderr.write(f"Warning: {message}\n")

def _detect_format(path):
    """Detect alignment format from file extension. Returns 'bam', 'sam', or 'cram'."""
    ext = os.path.splitext(path)[1].lower()
    mapping = {'.bam': 'bam', '.sam': 'sam', '.cram': 'cram'}
    fmt = mapping.get(ext)
    if fmt is None:
        die(f"Cannot determine format for '{path}': unknown extension '{ext}'. "
            "Rename with .bam/.sam/.cram or use -B/-S/-C.")
    return fmt


def _pysam_read_mode(fmt):
    return {'bam': 'rb', 'sam': 'r', 'cram': 'rc'}[fmt]


def _pysam_write_mode(fmt):
    return {'bam': 'wb', 'sam': 'w', 'cram': 'wc'}[fmt]


def _fmt_ext(fmt):
    return {'bam': '.bam', 'sam': '.sam', 'cram': '.cram'}[fmt]


def _resolve_out_fmt(in_fmt, args):
    """Determine output format: explicit flag overrides; default matches input."""
    if args.sam:  return 'sam'
    if args.bam:  return 'bam'
    if args.cram: return 'cram'
    return in_fmt


def _spool_and_detect_stdin():
    """Spool stdin to a temp file; detect BAM/CRAM/SAM by magic bytes.

    Returns (tmp_path, fmt) where fmt is 'bam', 'cram', or 'sam'.
    Caller is responsible for removing tmp_path.
    """
    fd, tmp_path = tempfile.mkstemp(prefix="pyselectal.stdin.", suffix=".tmp")
    os.close(fd)
    try:
        with open(tmp_path, "wb") as f:
            while True:
                chunk = sys.stdin.buffer.read(1024 * 1024)
                if not chunk:
                    break
                f.write(chunk)
    except Exception as e:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        die(f"Failed to spool stdin: {e}")
    # CRAM starts with b'CRAM'; BAM (BGZF) starts with gzip magic \x1f\x8b; SAM is plain text
    fmt = 'sam'
    try:
        with open(tmp_path, "rb") as f:
            magic = f.read(4)
        if magic[:4] == b'CRAM':
            fmt = 'cram'
        elif magic[:2] == b'\x1f\x8b':
            fmt = 'bam'
    except Exception:
        pass
    return tmp_path, fmt


def _cram_to_temp_bam(in_path, threads, reference):
    """Stream a CRAM file to a temp BAM (needed before pysam.sort which cannot pass a CRAM reference).

    Returns tmp_path; caller is responsible for removal.
    """
    fd, tmp_path = tempfile.mkstemp(prefix="pyselectal.cram2bam.", suffix=".bam")
    os.close(fd)
    kw_in  = {'reference_filename': reference} if reference else {}
    kw_out = {}
    if threads > 1:
        kw_in['threads']  = threads
        kw_out['threads'] = threads
    try:
        with pysam.AlignmentFile(in_path, 'rc', **kw_in) as src:
            with pysam.AlignmentFile(tmp_path, 'wb', template=src, **kw_out) as dst:
                for aln in src.fetch(until_eof=True):
                    dst.write(aln)
    except Exception as e:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        die(f"Failed to convert CRAM to BAM: {e}")
    return tmp_path


def revcomp(seq: str) -> str:
    """
    Reverse-complement a DNA sequence (ACGTN, case-insensitive).
    """
    table = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(table)[::-1]


def comp_base(base: str) -> str:
    """
    Complement of a single DNA base (A/C/G/T/N), case-insensitive.
    """
    base = base.upper()
    table = {"A": "T", "C": "G", "G": "C", "T": "A", "N": "N"}
    return table.get(base, "N")


def get_5prime_cigar(aln):
    """
    Get the 5'-end CIGAR operation and length, orientation-aware.

    Forward reads: first CIGAR element.
    Reverse reads: last  CIGAR element.
    """
    cig = aln.cigartuples
    if not cig:
        return None, None
    return (cig[-1] if aln.is_reverse else cig[0])


def classify_5prime_type(aln, mapped_prefix=5):
    """
    Classify an alignment's 5' end type as a string.

    Returns e.g. '1Sg', '2Sgg', '3Saac', '75Mgaggg', or None for unmapped.
    For soft-clipped reads: {length}S{sequence} (forward-strand sequence,
    reverse reads report the reverse-complement of the 3'-end soft-clip).
    For mapped 5' ends: {match_length}M{first_N_bases} where N=mapped_prefix.
    If mapped_prefix=0: returns {match_length}M with no sequence.
    """
    if aln.is_unmapped:
        return None

    op, length = get_5prime_cigar(aln)
    if op == SOFT:
        seq = aln.query_sequence
        if not seq:
            return f"{length}S"
        if aln.is_reverse:
            sc_seq = revcomp(seq[-length:])
        else:
            sc_seq = seq[:length]
        return f"{length}S{sc_seq.lower()}"
    elif op == MATCH:
        if mapped_prefix == 0:
            return f"{length}M"
        seq = aln.query_sequence
        if not seq:
            return f"{length}M"
        n = min(mapped_prefix, length)
        if aln.is_reverse:
            pfx = revcomp(seq[-n:])
        else:
            pfx = seq[:n]
        return f"{length}M{pfx.lower()}"
    return None


def has_5prime_mapped_exact(aln, prefix, rc_prefix, k, warn_k_ignored=False):
    """
    Exact match mode (n = m = 0):
        - require that the 5'-end is mapped;
        - 3'-end soft-clips are allowed;
    If --prefix is NOT given:
        - keep reads whose 5'-end CIGAR operation is MATCH (M) and M-length >= k
           - if k==0: keep reads with 5'-end operation == M (i.e., no 5' soft-clip)

    If --prefix IS given:
        1) if k <= len(prefix) (including k==0):
            - ignore k; select ONLY by full prefix match (orientation-aware)
            - emit a warning to stderr once
        2) if k > len(prefix):
            - require BOTH:
                - 5'-end operation == M with M-length >= k
                - full prefix match (orientation-aware)
    """
    if aln.is_unmapped:
        return False

    op_5p, len_5p = get_5prime_cigar(aln)

    # Need a MATCH at 5'-end to enforce a prefix
    if op_5p != MATCH:
        return False

    if len_5p is None:
        return False
    

    # If no prefix, only enforce minimum 5' MATCH length (>=k); if k==0, any 5' M passes.
    if not prefix:
        if k <= 0:
            return True
        return int(len_5p) >= k
    
    # If prefix provided:
    # Case 1: k <= len(prefix) => ignore k, select only by FULL prefix match
    if k <= len(prefix):
        if warn_k_ignored:
            warn_once("k_ignored", f"-k={k} is <= len(--prefix)={len(prefix)} in mapped mode; ignoring -k and selecting only by full prefix.")
        seq = aln.query_sequence
        if not seq:
            return False
        if aln.is_reverse:
            return seq[-len(prefix):].upper() == rc_prefix.upper()
        return seq[:len(prefix)].upper() == prefix.upper()

    # Case 2: k > len(prefix) => require BOTH: 5' MATCH length >= k AND full prefix match
    if len_5p < k:
        return False

    seq = aln.query_sequence
    if not seq:
        return False

    if aln.is_reverse:
        return seq[-len(prefix):].upper() == rc_prefix.upper()
    return seq[:len(prefix)].upper() == prefix.upper()

def has_5prime_softclip_exact(aln, n, prefix, rc_prefix):
    """
    Exact mode: n = m > 0 (exact n-bp 5' soft-clip):
      - require 5'-end is (S, n)
      - if prefix provided (length n):
          forward: first n bases == prefix
          reverse: last  n bases == revcomp(prefix)
    """
    if aln.is_unmapped:
        return False

    op_5p, len_5p = get_5prime_cigar(aln)
    if op_5p != SOFT or len_5p != n:
        return False

    if not prefix:
        return True

    seq = aln.query_sequence
    if not seq:
        return False

    if aln.is_reverse:
        return seq[-n:].upper() == rc_prefix.upper()
    return seq[:n].upper() == prefix.upper()

def has_5prime_softclip_range(aln, n, m, seq=None, rc_seq=None):
    """
    Range mode: n <= length <= m (m=None means unbounded).

      - require a 5' soft-clip of length x where:
            x >= n
            and (if m is not None) x <= m
      - if seq is None: only length is checked;
      - if seq is provided: soft-clip must start with seq (prefix match)
          forward: starts with seq
          reverse: starts with rc_seq (revcomp of seq)
    """
    if aln.is_unmapped:
        return False

    op_5p, len_5p = get_5prime_cigar(aln)

    if op_5p != SOFT:
        return False

    x = int(len_5p)
    if x < n:
        return False
    if m is not None and x > m:
        return False

    if not seq:
        return True

    query = aln.query_sequence
    if not query:
        return False

    if aln.is_reverse:
        sc = query[-x:]
        target = rc_seq.upper()
    else:
        sc = query[:x]
        target = seq.upper()

    return sc.upper().startswith(target)


def alignment_matches(aln, n, m, prefix, k, warn_k_ignored=False):
    """
    Wrapper that dispatches to the correct mode based on n_min, n_max.

      if n == m > 0 and m is not None,   then: exact soft-clip regime with optional sequence.
      if n == m == 0,                    then: mapped 5'-end regime with optional prefix.
      if n < m (including m is None),    then: range soft-clip regime with optional prefix.
    """
    # mapped 5'-end regime
    if n == 0 and m == 0:
        rc_prefix = revcomp(prefix) if prefix else None
        return has_5prime_mapped_exact(aln, prefix, rc_prefix, k, warn_k_ignored=warn_k_ignored)

    # exact soft-clip regime only if m is provided and equals n (>0)
    if m is not None and n == m and n > 0:
        rc_prefix = revcomp(prefix) if prefix else None
        return has_5prime_softclip_exact(aln, n, prefix, rc_prefix)

    # range regime (includes unbounded top if m is None)
    rc_prefix = revcomp(prefix) if prefix else None
    return has_5prime_softclip_range(aln, n, m, prefix, rc_prefix)


def parse_spec(spec_str):
    """
    Parse a 5' end type spec into a structured dict.

    Grammar: [n[.m]]<S|M>[seq]   (case-insensitive)

    Returns dict with keys: type, n, m, seq, raw.
    """
    raw = spec_str.strip()
    if not raw:
        die("Empty spec.")
    s = raw.upper()

    pattern = r'^(\d+)?(\.(\d*))?([SM])([ACGT]*)$'
    match = re.match(pattern, s)
    if not match:
        die(f"Invalid spec: '{raw}'")

    n_str, dot_group, m_str, mode, seq = match.groups()
    has_dot = dot_group is not None
    seq = seq if seq else None

    if has_dot:
        n_val = int(n_str) if n_str else 0
        m_val = int(m_str) if m_str else None
    elif n_str is not None:
        n_val = int(n_str)
        m_val = n_val  # nS or nM without dot = exact
    else:
        # Bare S or M with no numbers
        n_val = None
        m_val = None

    if mode == "S":
        if n_val is None:
            # Bare "S[seq]": any softclip
            n_val = 1
            m_val = None
        if n_val == m_val and n_val > 0:
            # Exact softclip — validate / expand seq
            if seq is not None:
                if len(seq) == 1:
                    seq = seq * n_val  # expand homopolymer
                elif len(seq) != n_val:
                    die(f"Spec '{raw}': sequence length {len(seq)} != n={n_val}.")
        else:
            # Range softclip — seq is a prefix to match
            pass
    else:
        # M mode
        if n_val is None:
            n_val = 0
            m_val = None

    # Normalized lowercase raw for file naming
    norm = raw.lower()

    return {
        "type": mode,
        "n": n_val,
        "m": m_val,
        "seq": seq,
        "raw": norm,
    }


def spec_matches_aln(aln, spec):
    """Check if an alignment matches a parsed spec."""
    if aln.is_unmapped:
        return False

    if spec["type"] == "S":
        n, m, seq = spec["n"], spec["m"], spec["seq"]
        if m is not None and n == m and n > 0:
            # Exact softclip
            prefix = seq
            rc_prefix = revcomp(prefix) if prefix else None
            return has_5prime_softclip_exact(aln, n, prefix, rc_prefix)
        else:
            # Range softclip — seq is prefix to match
            rc_seq = revcomp(seq) if seq else None
            return has_5prime_softclip_range(aln, n, m, seq, rc_seq)
    else:
        # M mode
        prefix = spec["seq"]
        rc_prefix = revcomp(prefix) if prefix else None
        k = spec["n"]
        if not has_5prime_mapped_exact(aln, prefix, rc_prefix, k):
            return False
        # Enforce max match length if specified
        if spec["m"] is not None:
            op, length = get_5prime_cigar(aln)
            if op == MATCH and length > spec["m"]:
                return False
        return True


def process_group_pe(records, n, m, prefix, k, out_bam, warn_k_ignored=False):
    """
    Paired-end mode: records = list of AlignedSegment with the same query_name.

    - Select R1 alignments passing alignment_matches.
    - For each selected R1:
        - collect its (next_reference_id, next_reference_start).
    - Write:
        - all passing R1s;
        - all R2s with (reference_id, reference_start) collected above.
    """
    if not records:
        return

    r1_selected = []
    for r in records:
        if r.is_read1 and alignment_matches(r, n, m, prefix, k, warn_k_ignored=warn_k_ignored):
            r1_selected.append(r)

    if not r1_selected:
        return

    mate_coords = []
    for r in r1_selected:
        if r.next_reference_id is not None and r.next_reference_start is not None:
            if r.next_reference_id >= 0 and r.next_reference_start >= 0:
                mate_coords.append((r.next_reference_id, r.next_reference_start))

    r2_selected = []
    for r in records:
        if r.is_read2 and (r.reference_id, r.reference_start) in mate_coords:
            r2_selected.append(r)

    for r in r1_selected + r2_selected:
        out_bam.write(r)


def process_stream_se(in_bam, out_bam, n, m, prefix, k, warn_k_ignored=False):
    """
    Single-end mode: filter each alignment independently.
    """
    for aln in in_bam.fetch(until_eof=True):
        if alignment_matches(aln, n, m, prefix, k, warn_k_ignored=warn_k_ignored):
            out_bam.write(aln)


def process_stream_pe(in_bam, out_bam, n, m, prefix, k, warn_k_ignored=False):
    """
    Paired-end mode: assumes name-sorted BAM (all alignments with the same
    query_name are contiguous).
    """
    current_qname = None
    group = []

    for aln in in_bam.fetch(until_eof=True):
        qn = aln.query_name
        if current_qname is None:
            current_qname = qn
            group = [aln]
        elif qn == current_qname:
            group.append(aln)
        else:
            process_group_pe(group, n, m, prefix, k, out_bam, warn_k_ignored=warn_k_ignored)
            current_qname = qn
            group = [aln]

    if group:
        process_group_pe(group, n, m, prefix, k, out_bam, warn_k_ignored=warn_k_ignored)


def parse_args(argv):

    pre = argparse.ArgumentParser(add_help=False)
    pre.add_argument("-h", "--help", action="store_true")
    pre.add_argument("-v", "--version", action="store_true")
    pre_args, _ = pre.parse_known_args(argv)

    if pre_args.version:
        sys.stdout.write(f"{VERSION}\n")
        raise SystemExit(0)

    parser = HelpfulArgumentParser(
        add_help=False,
        description="Filter and analyze BAM alignments by 5'-end type."
    )

    # Help/version
    parser.add_argument("-h", "--help", action="store_true",
                        help="Show the manual and exit.")
    parser.add_argument("-v", "--version", action="store_true",
                        help="Print version information and exit.")

    # Input
    parser.add_argument("-i", "--input", type=str, default=None,
                        help="Input BAM file(s), comma-separated.")

    # Modes (mutually exclusive)
    parser.add_argument("-s", "--select", type=str, default=None,
                        help="Select alignments by 5' end type spec(s), comma-separated.")
    parser.add_argument("-c", "--count", action="store_true",
                        help="Output TSV histogram of 5' end types.")
    parser.add_argument("-a", "--all", action="store_true",
                        help="Split alignments into per-type BAM files.")

    # Options
    parser.add_argument("-m", "--merge", action="store_true",
                        help="Merge multiple --select specs into one output file.")
    parser.add_argument("-n", "--name", action="store_true",
                        help="Name-sort input BAM internally (pysam.sort -n).")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Number of BGZF threads for pysam I/O.")
    parser.add_argument("-p", "--paired", action="store_true",
                        help="Paired-end mode: select read1 and emit matching read2 mates.")
    parser.add_argument("-o", "--output", type=str, default=None,
                        help="Output file path (overrides automatic naming).")

    # Output format overrides (mutually exclusive)
    fmt_group = parser.add_mutually_exclusive_group()
    fmt_group.add_argument("-S", "--sam", action="store_true",
                           help="Force SAM output.")
    fmt_group.add_argument("-B", "--bam", action="store_true",
                           help="Force BAM output.")
    fmt_group.add_argument("-C", "--cram", action="store_true",
                           help="Force CRAM output (requires -r/--reference).")

    parser.add_argument("-r", "--reference", type=str, default=None,
                        help="Reference FASTA for CRAM input/output.")

    parser.add_argument("--collapse-threshold", type=float, default=1.0,
                        metavar="PCT",
                        help="Collapse categories below PCT%% into 'other' row in --count output (default: 1; 0 = off).")
    parser.add_argument("--mapped-prefix", type=int, default=5,
                        metavar="N",
                        help="Bases of 5' matched sequence to show in --count output (default: 5; 0 = length only).")

    if pre_args.help or len(argv) == 0:
        sys.stdout.write(MANUAL)
        parser.print_help(sys.stdout)
        raise SystemExit(0)

    args = parser.parse_args(argv)

    if args.help:
        sys.stdout.write(MANUAL)
        parser.print_help(sys.stdout)
        raise SystemExit(0)

    # Require -i/--input
    if not args.input:
        parser.error("-i/--input is required.")

    # Parse comma-separated input files
    args.input_files = [f.strip() for f in args.input.split(",") if f.strip()]
    if not args.input_files:
        parser.error("-i/--input must specify at least one file.")

    # Exactly one mode required
    modes = sum([args.select is not None, args.count, getattr(args, 'all')])
    if modes == 0:
        parser.error("Exactly one mode is required: -s/--select, -c/--count, or -a/--all.")
    if modes > 1:
        parser.error("Modes -s/--select, -c/--count, and -a/--all are mutually exclusive.")

    # --merge only with --select
    if args.merge and args.select is None:
        parser.error("-m/--merge can only be used with -s/--select.")

    if args.threads < 1:
        parser.error("--threads must be >= 1.")

    if args.collapse_threshold < 0 or args.collapse_threshold >= 100:
        parser.error("--collapse-threshold must be in [0, 100).")

    if args.mapped_prefix < 0:
        parser.error("--mapped-prefix must be >= 0.")

    if args.cram and not args.reference:
        parser.error("-C/--cram requires -r/--reference.")

    return args

def name_sort_bam(in_bam_path, threads):
    """
    Name-sort BAM using pysam.sort (-n). Returns path to sorted BAM.
    Caller is responsible for cleanup of the returned temp file.
    """
    # Put temp next to output for large files. Using system temp by default.
    fd, tmp_path = tempfile.mkstemp(prefix="softclip5.namesort.", suffix=".bam")
    os.close(fd)

    try:
        # pysam exposes samtools sort as pysam.sort()
        if threads > 1:
            # samtools sort thread option is -@; passthrough in pysam.sort
            pysam.sort("-n", "-@", str(threads), "-o", tmp_path, in_bam_path, catch_stdout=False)
        else:
            pysam.sort("-n", "-o", tmp_path, in_bam_path, catch_stdout=False)
    except Exception as e:
        try:
            os.remove(tmp_path)
        except OSError:
            pass
        die(f"Failed to name-sort BAM: {e}")

    return tmp_path

def open_alignment_files(in_path, out_path, threads,
                         in_fmt='bam', out_fmt='bam', reference=None):
    """Open input and output alignment files with format and CRAM-reference support.

    in_path/out_path may be '-' for stdin/stdout.
    """
    read_mode  = _pysam_read_mode(in_fmt)
    write_mode = _pysam_write_mode(out_fmt)

    kw_in = {}
    if threads > 1:
        kw_in['threads'] = threads
    if in_fmt == 'cram' and reference:
        kw_in['reference_filename'] = reference

    try:
        if in_path == "-":
            in_bam = pysam.AlignmentFile(sys.stdin.buffer, read_mode, **kw_in)
        else:
            in_bam = pysam.AlignmentFile(in_path, read_mode, **kw_in)

        kw_out = {'template': in_bam}
        if threads > 1:
            kw_out['threads'] = threads
        if out_fmt == 'cram' and reference:
            kw_out['reference_filename'] = reference

        if out_path == "-":
            out_bam = pysam.AlignmentFile(sys.stdout.buffer, write_mode, **kw_out)
        else:
            out_bam = pysam.AlignmentFile(out_path, write_mode, **kw_out)
    except Exception as e:
        die(f"Failed to open alignment files: {e}")
    return in_bam, out_bam
                                                        
def process_select_se(in_bam, out_bam, spec):
    """Single-end --select: write alignments matching spec."""
    for aln in in_bam.fetch(until_eof=True):
        if spec_matches_aln(aln, spec):
            out_bam.write(aln)


def process_select_pe(in_bam, out_bam, spec):
    """Paired-end --select: group by query_name, select R1, emit R1+R2."""
    current_qname = None
    group = []

    for aln in in_bam.fetch(until_eof=True):
        qn = aln.query_name
        if current_qname is None:
            current_qname = qn
            group = [aln]
        elif qn == current_qname:
            group.append(aln)
        else:
            _flush_pe_group(group, spec, out_bam)
            current_qname = qn
            group = [aln]

    if group:
        _flush_pe_group(group, spec, out_bam)


def _flush_pe_group(records, spec, out_bam):
    """Select R1s matching spec, then emit their R2 mates."""
    if not records:
        return

    r1_selected = [r for r in records
                   if r.is_read1 and spec_matches_aln(r, spec)]
    if not r1_selected:
        return

    mate_coords = set()
    for r in r1_selected:
        if (r.next_reference_id is not None and r.next_reference_start is not None
                and r.next_reference_id >= 0 and r.next_reference_start >= 0):
            mate_coords.add((r.next_reference_id, r.next_reference_start))

    for r in r1_selected:
        out_bam.write(r)
    for r in records:
        if r.is_read2 and (r.reference_id, r.reference_start) in mate_coords:
            out_bam.write(r)


def _select_output_path(in_path, spec, args, out_fmt):
    """Determine output path for a --select run."""
    if args.output:
        return args.output
    specs = [s.strip() for s in args.select.split(",") if s.strip()]
    if len(args.input_files) == 1 and len(specs) == 1:
        return "-"  # stdout
    stem = os.path.splitext(os.path.basename(in_path))[0]
    return f"{stem}_{spec['raw']}{_fmt_ext(out_fmt)}"


def _merge_output_path(in_path, args, out_fmt):
    """Determine output path for a --select --merge run."""
    if args.output:
        return args.output
    stem = os.path.splitext(os.path.basename(in_path))[0]
    return f"{stem}_merged{_fmt_ext(out_fmt)}"


def process_select_merge_se(in_bam, out_bam, specs):
    """Single-end --select --merge: write alignments matching any spec, once each."""
    for aln in in_bam.fetch(until_eof=True):
        for spec in specs:
            if spec_matches_aln(aln, spec):
                out_bam.write(aln)
                break


def process_select_merge_pe(in_bam, out_bam, specs):
    """Paired-end --select --merge: group by query_name, select R1 matching any spec."""
    current_qname = None
    group = []

    for aln in in_bam.fetch(until_eof=True):
        qn = aln.query_name
        if current_qname is None:
            current_qname = qn
            group = [aln]
        elif qn == current_qname:
            group.append(aln)
        else:
            _flush_pe_group_merge(group, specs, out_bam)
            current_qname = qn
            group = [aln]

    if group:
        _flush_pe_group_merge(group, specs, out_bam)


def _flush_pe_group_merge(records, specs, out_bam):
    """Select R1s matching any spec (deduplicated), then emit their R2 mates."""
    if not records:
        return

    r1_selected = []
    for r in records:
        if not r.is_read1:
            continue
        for spec in specs:
            if spec_matches_aln(r, spec):
                r1_selected.append(r)
                break

    if not r1_selected:
        return

    mate_coords = set()
    for r in r1_selected:
        if (r.next_reference_id is not None and r.next_reference_start is not None
                and r.next_reference_id >= 0 and r.next_reference_start >= 0):
            mate_coords.add((r.next_reference_id, r.next_reference_start))

    for r in r1_selected:
        out_bam.write(r)
    for r in records:
        if r.is_read2 and (r.reference_id, r.reference_start) in mate_coords:
            out_bam.write(r)


def run_select(args):
    """Execute --select mode."""
    specs_raw = [s.strip() for s in args.select.split(",") if s.strip()]
    specs = [parse_spec(s) for s in specs_raw]

    for in_path in args.input_files:
        actual_in = in_path
        tmp_stdin   = None
        tmp_cram2bam = None
        tmp_sorted  = None

        try:
            # Detect input format
            if in_path == "-":
                tmp_stdin, in_fmt = _spool_and_detect_stdin()
                actual_in = tmp_stdin
            else:
                in_fmt = _detect_format(in_path)

            if in_fmt == 'cram' and not args.reference:
                die(f"CRAM input '{in_path}' requires -r/--reference.")

            out_fmt = _resolve_out_fmt(in_fmt, args)
            if out_fmt == 'cram' and not args.reference:
                die("CRAM output requires -r/--reference.")

            # Name sort: CRAM must first be converted to BAM
            if args.name:
                if in_fmt == 'cram':
                    tmp_cram2bam = _cram_to_temp_bam(actual_in, args.threads, args.reference)
                    actual_in = tmp_cram2bam
                tmp_sorted = name_sort_bam(actual_in, args.threads)
                actual_in = tmp_sorted
                in_fmt = 'bam'

            if args.merge:
                out_path = _merge_output_path(in_path, args, out_fmt)
                in_bam, out_bam = open_alignment_files(
                    actual_in, out_path, args.threads, in_fmt, out_fmt, args.reference)
                try:
                    if args.paired:
                        process_select_merge_pe(in_bam, out_bam, specs)
                    else:
                        process_select_merge_se(in_bam, out_bam, specs)
                finally:
                    out_bam.close()
                    in_bam.close()
            else:
                for spec in specs:
                    out_path = _select_output_path(in_path, spec, args, out_fmt)
                    in_bam, out_bam = open_alignment_files(
                        actual_in, out_path, args.threads, in_fmt, out_fmt, args.reference)
                    try:
                        if args.paired:
                            process_select_pe(in_bam, out_bam, spec)
                        else:
                            process_select_se(in_bam, out_bam, spec)
                    finally:
                        out_bam.close()
                        in_bam.close()
        finally:
            for tmp in [tmp_stdin, tmp_cram2bam, tmp_sorted]:
                if tmp:
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass


def _count_output_path(in_path, args):
    """Determine output path for a --count run."""
    if args.output:
        return args.output
    if len(args.input_files) == 1:
        return "-"  # stdout
    stem = os.path.splitext(os.path.basename(in_path))[0]
    return f"{stem}_5prime_counts.tsv"


def process_count_se(in_bam, mapped_prefix=5):
    """Count 5' end types for all alignments. Returns dict {type: count}."""
    counts = {}
    for aln in in_bam.fetch(until_eof=True):
        t = classify_5prime_type(aln, mapped_prefix=mapped_prefix)
        if t is not None:
            counts[t] = counts.get(t, 0) + 1
    return counts


def process_count_pe(in_bam, mapped_prefix=5):
    """Count 5' end types for R1 alignments only. Returns dict {type: count}."""
    counts = {}
    for aln in in_bam.fetch(until_eof=True):
        if not aln.is_read1:
            continue
        t = classify_5prime_type(aln, mapped_prefix=mapped_prefix)
        if t is not None:
            counts[t] = counts.get(t, 0) + 1
    return counts


def write_count_tsv(counts, out_path, collapse_threshold=1.0):
    """Write type\\tcount TSV sorted by descending count.

    Categories whose percentage of total is strictly below collapse_threshold
    are collapsed into a single 'other (<N%)' row appended at the end.
    collapse_threshold=0 disables collapsing.
    """
    total = sum(counts.values())
    if collapse_threshold == 0 or total == 0:
        main_rows = sorted(counts.items(), key=lambda x: (-x[1], x[0]))
        other_count = 0
    else:
        main_items = []
        other_count = 0
        for t, c in counts.items():
            if c / total * 100 < collapse_threshold:
                other_count += c
            else:
                main_items.append((t, c))
        main_rows = sorted(main_items, key=lambda x: (-x[1], x[0]))

    if out_path == "-":
        f = sys.stdout
    else:
        f = open(out_path, "w")
    try:
        f.write("type\tcount\n")
        for t, c in main_rows:
            f.write(f"{t}\t{c}\n")
        if other_count > 0:
            label = f"other (<{collapse_threshold:g}%)"
            f.write(f"{label}\t{other_count}\n")
    finally:
        if out_path != "-":
            f.close()


def run_count(args):
    """Execute --count mode."""
    for in_path in args.input_files:
        actual_in = in_path
        tmp_stdin   = None
        tmp_cram2bam = None
        tmp_sorted  = None

        try:
            # Detect input format
            if in_path == "-":
                tmp_stdin, in_fmt = _spool_and_detect_stdin()
                actual_in = tmp_stdin
            else:
                in_fmt = _detect_format(in_path)

            if in_fmt == 'cram' and not args.reference:
                die(f"CRAM input '{in_path}' requires -r/--reference.")

            # Name sort
            if args.name:
                if in_fmt == 'cram':
                    tmp_cram2bam = _cram_to_temp_bam(actual_in, args.threads, args.reference)
                    actual_in = tmp_cram2bam
                tmp_sorted = name_sort_bam(actual_in, args.threads)
                actual_in = tmp_sorted
                in_fmt = 'bam'

            kw = {}
            if args.threads > 1:
                kw['threads'] = args.threads
            if in_fmt == 'cram' and args.reference:
                kw['reference_filename'] = args.reference
            in_bam = pysam.AlignmentFile(actual_in, _pysam_read_mode(in_fmt), **kw)

            try:
                if args.paired:
                    counts = process_count_pe(in_bam, mapped_prefix=args.mapped_prefix)
                else:
                    counts = process_count_se(in_bam, mapped_prefix=args.mapped_prefix)
            finally:
                in_bam.close()

            out_path = _count_output_path(in_path, args)
            write_count_tsv(counts, out_path, collapse_threshold=args.collapse_threshold)
        finally:
            for tmp in [tmp_stdin, tmp_cram2bam, tmp_sorted]:
                if tmp:
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass


def _all_output_path(in_path, type_str, args, out_fmt):
    """Determine output path for a --all type file."""
    stem = os.path.splitext(os.path.basename(in_path))[0]
    fname = f"{stem}_{type_str.lower()}{_fmt_ext(out_fmt)}"
    if args.output:
        return os.path.join(args.output, fname)
    return fname


def _get_or_open_all_file(type_str, in_bam, in_path, args, open_files):
    """Return open output file for type_str, opening lazily if needed.

    Expects args.out_fmt and args.reference to be set by run_all.
    """
    if type_str not in open_files:
        out_path = _all_output_path(in_path, type_str, args, args.out_fmt)
        kw = {'template': in_bam}
        if args.threads > 1:
            kw['threads'] = args.threads
        if args.out_fmt == 'cram' and args.reference:
            kw['reference_filename'] = args.reference
        open_files[type_str] = pysam.AlignmentFile(
            out_path, _pysam_write_mode(args.out_fmt), **kw)
    return open_files[type_str]


def process_all_se(in_bam, in_path, args, collapse_types=None):
    """Single-end --all: write each alignment to its type-specific BAM."""
    open_files = {}
    try:
        for aln in in_bam.fetch(until_eof=True):
            t = classify_5prime_type(aln)
            if t is None:
                continue
            effective_t = "other" if (collapse_types and t in collapse_types) else t
            out = _get_or_open_all_file(effective_t, in_bam, in_path, args, open_files)
            out.write(aln)
    finally:
        for f in open_files.values():
            f.close()


def _flush_pe_group_all(records, in_bam, in_path, args, open_files, collapse_types=None):
    """Route each R1 (by type) and its R2 mates into the matching output BAM."""
    if not records:
        return

    mate_type = {}
    for r in records:
        if not r.is_read1:
            continue
        t = classify_5prime_type(r)
        if t is None:
            continue
        effective_t = "other" if (collapse_types and t in collapse_types) else t
        out = _get_or_open_all_file(effective_t, in_bam, in_path, args, open_files)
        out.write(r)
        if (r.next_reference_id is not None and r.next_reference_start is not None
                and r.next_reference_id >= 0 and r.next_reference_start >= 0):
            mate_type[(r.next_reference_id, r.next_reference_start)] = effective_t

    for r in records:
        if not r.is_read2:
            continue
        t = mate_type.get((r.reference_id, r.reference_start))
        if t is not None:
            open_files[t].write(r)


def process_all_pe(in_bam, in_path, args, collapse_types=None):
    """Paired-end --all: group by query_name, route R1+R2 pairs by R1 type."""
    open_files = {}
    current_qname = None
    group = []

    try:
        for aln in in_bam.fetch(until_eof=True):
            qn = aln.query_name
            if current_qname is None:
                current_qname = qn
                group = [aln]
            elif qn == current_qname:
                group.append(aln)
            else:
                _flush_pe_group_all(group, in_bam, in_path, args, open_files, collapse_types)
                current_qname = qn
                group = [aln]

        if group:
            _flush_pe_group_all(group, in_bam, in_path, args, open_files, collapse_types)
    finally:
        for f in open_files.values():
            f.close()


def run_all(args):
    """Execute --all mode."""
    for in_path in args.input_files:
        actual_in = in_path
        tmp_stdin   = None
        tmp_cram2bam = None
        tmp_sorted  = None

        try:
            # Detect input format
            if in_path == "-":
                tmp_stdin, in_fmt = _spool_and_detect_stdin()
                actual_in = tmp_stdin
            else:
                in_fmt = _detect_format(in_path)

            if in_fmt == 'cram' and not args.reference:
                die(f"CRAM input '{in_path}' requires -r/--reference.")

            out_fmt = _resolve_out_fmt(in_fmt, args)
            if out_fmt == 'cram' and not args.reference:
                die("CRAM output requires -r/--reference.")

            # Attach resolved output format to args so helpers can access it
            args.out_fmt = out_fmt

            # Name sort
            if args.name:
                if in_fmt == 'cram':
                    tmp_cram2bam = _cram_to_temp_bam(actual_in, args.threads, args.reference)
                    actual_in = tmp_cram2bam
                tmp_sorted = name_sort_bam(actual_in, args.threads)
                actual_in = tmp_sorted
                in_fmt = 'bam'

            kw = {}
            if args.threads > 1:
                kw['threads'] = args.threads
            if in_fmt == 'cram' and args.reference:
                kw['reference_filename'] = args.reference

            # Pass 1: count types to determine which to collapse
            collapse_types = set()
            if args.collapse_threshold > 0:
                in_bam = pysam.AlignmentFile(actual_in, _pysam_read_mode(in_fmt), **kw)
                try:
                    if args.paired:
                        counts = process_count_pe(in_bam, mapped_prefix=args.mapped_prefix)
                    else:
                        counts = process_count_se(in_bam, mapped_prefix=args.mapped_prefix)
                finally:
                    in_bam.close()
                total = sum(counts.values())
                if total > 0:
                    collapse_types = {t for t, c in counts.items()
                                      if c / total * 100 < args.collapse_threshold}

            if args.output:
                os.makedirs(args.output, exist_ok=True)

            # Pass 2: route alignments, collapsing rare types into "other"
            in_bam = pysam.AlignmentFile(actual_in, _pysam_read_mode(in_fmt), **kw)
            try:
                if args.paired:
                    process_all_pe(in_bam, in_path, args, collapse_types=collapse_types)
                else:
                    process_all_se(in_bam, in_path, args, collapse_types=collapse_types)
            finally:
                in_bam.close()
        finally:
            for tmp in [tmp_stdin, tmp_cram2bam, tmp_sorted]:
                if tmp:
                    try:
                        os.remove(tmp)
                    except OSError:
                        pass


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    args = parse_args(argv)

    if args.select is not None:
        run_select(args)
    elif args.count:
        run_count(args)
    elif getattr(args, 'all'):
        run_all(args)

if __name__ == "__main__":
    main()
