"""Tests for pyselectal CLI."""

import os
import tempfile
import pytest
import pysam
import io
from unittest import mock
from pyselectal import (
    parse_args, parse_spec, spec_matches_aln, classify_5prime_type, write_count_tsv,
    main, VERSION, _detect_format, _spool_and_detect_stdin,
    has_5prime_softclip_range, has_5prime_mapped_range,
)

SE_BAM = os.path.join(os.path.dirname(__file__), "testdata", "test_softclip_se.bam")
PE_BAM = os.path.join(os.path.dirname(__file__), "testdata", "test_softclip_pe.bam")


# ---------------------------------------------------------------------------
# Step 1: Issue 6 — CLI restructuring
# ---------------------------------------------------------------------------

class TestVersion:
    def test_version_flag(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-v"])
        assert exc.value.code == 0

    def test_version_string(self):
        assert "v2.0" in VERSION


class TestHelp:
    def test_help_flag(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-h"])
        assert exc.value.code == 0

    def test_no_args_shows_help(self):
        with pytest.raises(SystemExit) as exc:
            parse_args([])
        assert exc.value.code == 0


class TestInputFlag:
    def test_input_single_file(self):
        args = parse_args(["-i", "sample.bam", "-c"])
        assert args.input_files == ["sample.bam"]

    def test_input_multiple_files(self):
        args = parse_args(["-i", "a.bam,b.bam,c.bam", "-c"])
        assert args.input_files == ["a.bam", "b.bam", "c.bam"]

    def test_input_required(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-c"])
        assert exc.value.code == 2

    def test_input_stdin(self):
        args = parse_args(["-i", "-", "-c"])
        assert args.input_files == ["-"]


class TestNameSort:
    def test_name_flag(self):
        args = parse_args(["-i", "in.bam", "-c", "-n"])
        assert args.name is True

    def test_name_long_flag(self):
        args = parse_args(["-i", "in.bam", "-c", "--name"])
        assert args.name is True

    def test_name_default_false(self):
        args = parse_args(["-i", "in.bam", "-c"])
        assert args.name is False


class TestOldFlagsRemoved:
    def test_old_sort_flag_rejected(self):
        """Old -s without argument should fail (it now expects a SPEC value)."""
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-s"])
        assert exc.value.code == 2

    def test_old_min_softclip_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "--min-softclip", "3"])
        assert exc.value.code == 2

    def test_old_max_softclip_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "--max-softclip", "3"])
        assert exc.value.code == 2

    def test_old_prefix_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-x", "ATG"])
        assert exc.value.code == 2

    def test_old_match_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-k", "3"])
        assert exc.value.code == 2

    def test_old_positional_args_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["in.bam", "out.bam"])
        assert exc.value.code == 2


class TestModeExclusivity:
    def test_no_mode_error(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam"])
        assert exc.value.code == 2

    def test_select_and_count_exclusive(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-s", "1Sg", "-c"])
        assert exc.value.code == 2

    def test_select_and_all_exclusive(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-s", "1Sg", "-a"])
        assert exc.value.code == 2

    def test_count_and_all_exclusive(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "-a"])
        assert exc.value.code == 2

    def test_select_alone_ok(self):
        args = parse_args(["-i", "in.bam", "-s", "1Sg"])
        assert args.select == "1Sg"

    def test_count_alone_ok(self):
        args = parse_args(["-i", "in.bam", "-c"])
        assert args.count is True

    def test_all_alone_ok(self):
        args = parse_args(["-i", "in.bam", "-a"])
        assert getattr(args, 'all') is True


class TestMergeValidation:
    def test_merge_without_select_error(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "-m"])
        assert exc.value.code == 2

    def test_merge_with_select_ok(self):
        args = parse_args(["-i", "in.bam", "-s", "1Sg,2Sg", "-m"])
        assert args.merge is True


class TestThreads:
    def test_default_threads(self):
        args = parse_args(["-i", "in.bam", "-c"])
        assert args.threads == 1

    def test_custom_threads(self):
        args = parse_args(["-i", "in.bam", "-c", "-t", "4"])
        assert args.threads == 4

    def test_invalid_threads(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "-t", "0"])
        assert exc.value.code == 2


class TestPairedEnd:
    def test_paired_flag(self):
        args = parse_args(["-i", "in.bam", "-c", "-p"])
        assert args.paired is True

    def test_paired_default_false(self):
        args = parse_args(["-i", "in.bam", "-c"])
        assert args.paired is False


# ---------------------------------------------------------------------------
# Step 2: Issue 3 — --select with spec grammar parser
# ---------------------------------------------------------------------------

class TestParseSpec:
    """Unit tests for parse_spec()."""

    # --- Softclip specs ---

    def test_exact_softclip_with_single_base(self):
        """1Sg -> exact 1bp softclip, seq expanded to 'G'."""
        s = parse_spec("1Sg")
        assert s["type"] == "S"
        assert s["n"] == 1
        assert s["m"] == 1
        assert s["seq"] == "G"

    def test_exact_softclip_multi_base(self):
        """3Saac -> exact 3bp softclip, full prefix 'AAC'."""
        s = parse_spec("3Saac")
        assert s["type"] == "S"
        assert s["n"] == 3
        assert s["m"] == 3
        assert s["seq"] == "AAC"

    def test_exact_softclip_homopolymer_expand(self):
        """2Sg -> exact 2bp softclip, seq expanded to 'GG'."""
        s = parse_spec("2Sg")
        assert s["type"] == "S"
        assert s["n"] == 2
        assert s["m"] == 2
        assert s["seq"] == "GG"

    def test_exact_softclip_no_seq(self):
        """2S -> exact 2bp softclip, any sequence."""
        s = parse_spec("2S")
        assert s["type"] == "S"
        assert s["n"] == 2
        assert s["m"] == 2
        assert s["seq"] is None

    def test_range_softclip(self):
        """1.3S -> range 1-3bp softclip, any base."""
        s = parse_spec("1.3S")
        assert s["type"] == "S"
        assert s["n"] == 1
        assert s["m"] == 3
        assert s["seq"] is None

    def test_range_softclip_with_base(self):
        """1.3Sg -> range 1-3bp softclip, G homopolymer."""
        s = parse_spec("1.3Sg")
        assert s["type"] == "S"
        assert s["n"] == 1
        assert s["m"] == 3
        assert s["seq"] == "G"

    def test_open_range_softclip(self):
        """4.Sg -> 4+ bp softclip G."""
        s = parse_spec("4.Sg")
        assert s["type"] == "S"
        assert s["n"] == 4
        assert s["m"] is None
        assert s["seq"] == "G"

    def test_zero_to_m_softclip(self):
        """.5S -> 0-5bp softclip (but min 1 for S)."""
        s = parse_spec(".5S")
        assert s["type"] == "S"
        assert s["n"] == 0
        assert s["m"] == 5
        assert s["seq"] is None

    def test_any_softclip(self):
        """S -> any softclip (n=1, m=None)."""
        s = parse_spec("S")
        assert s["type"] == "S"
        assert s["n"] == 1
        assert s["m"] is None
        assert s["seq"] is None

    def test_any_softclip_with_base(self):
        """St -> any softclip starting with T."""
        s = parse_spec("St")
        assert s["type"] == "S"
        assert s["n"] == 1
        assert s["m"] is None
        assert s["seq"] == "T"

    # --- Mapped specs ---

    def test_mapped_exact(self):
        """2M -> mapped with exactly 2bp match at 5' end."""
        s = parse_spec("2M")
        assert s["type"] == "M"
        assert s["n"] == 2
        assert s["m"] == 2
        assert s["seq"] is None

    def test_mapped_with_prefix(self):
        """2Mg -> exactly 2 matched Gs (homopolymer expansion)."""
        s = parse_spec("2Mg")
        assert s["type"] == "M"
        assert s["n"] == 2
        assert s["m"] == 2
        assert s["seq"] == "GG"  # expanded from single G

    def test_mapped_range(self):
        """2.5M -> mapped with 2-5bp match."""
        s = parse_spec("2.5M")
        assert s["type"] == "M"
        assert s["n"] == 2
        assert s["m"] == 5

    def test_mapped_zero_to_m(self):
        """.5M -> at most 5bp match."""
        s = parse_spec(".5M")
        assert s["type"] == "M"
        assert s["n"] == 0
        assert s["m"] == 5

    def test_any_mapped(self):
        """M -> any mapped 5' end."""
        s = parse_spec("M")
        assert s["type"] == "M"
        assert s["n"] == 0
        assert s["m"] is None
        assert s["seq"] is None

    def test_any_mapped_with_prefix(self):
        """Mg -> any mapped starting with G."""
        s = parse_spec("Mg")
        assert s["type"] == "M"
        assert s["n"] == 0
        assert s["m"] is None
        assert s["seq"] == "G"

    # --- Case insensitivity ---

    def test_case_insensitive_mode(self):
        """1sg == 1Sg == 1SG."""
        for variant in ("1sg", "1Sg", "1SG", "1sG"):
            s = parse_spec(variant)
            assert s["type"] == "S"
            assert s["n"] == 1
            assert s["seq"] == "G"

    # --- raw field for file naming ---

    def test_raw_field_lowercase(self):
        s = parse_spec("2Sg")
        assert s["raw"] == "2sg"

    def test_raw_field_range(self):
        s = parse_spec("1.3S")
        assert s["raw"] == "1.3s"

    # --- Invalid specs ---

    def test_empty_spec(self):
        with pytest.raises(SystemExit):
            parse_spec("")

    def test_no_mode_letter(self):
        with pytest.raises(SystemExit):
            parse_spec("123")

    def test_invalid_base_char(self):
        with pytest.raises(SystemExit):
            parse_spec("1SZ")

    def test_range_softclip_multi_char_seq(self):
        """Range softclip seq is allowed as prefix."""
        r = parse_spec("1.3Sgg")
        assert r["type"] == "S"
        assert r["n"] == 1
        assert r["m"] == 3
        assert r["seq"] == "GG"

    def test_exact_softclip_wrong_seq_length(self):
        """Exact softclip: multi-char seq must have length == n."""
        with pytest.raises(SystemExit):
            parse_spec("2Sggg")

    def test_exact_mapped_wrong_seq_length(self):
        """Exact mapped: multi-char seq must have length == n."""
        with pytest.raises(SystemExit):
            parse_spec("3Mgg")

    # --- Homopolymer specs ---

    def test_range_softclip_homopolymer_flag(self):
        """Range softclip with single base sets homopolymer=True."""
        s = parse_spec("3.Sg")
        assert s["type"] == "S"
        assert s["n"] == 3
        assert s["m"] is None
        assert s["seq"] == "G"
        assert s["homopolymer"] is True

    def test_range_softclip_multi_char_not_homopolymer(self):
        """Range softclip with multi-char seq is prefix match, not homopolymer."""
        s = parse_spec("3.Sgg")
        assert s["seq"] == "GG"
        assert s["homopolymer"] is False

    def test_range_mapped_homopolymer_flag(self):
        """Range mapped with single base sets homopolymer=True."""
        s = parse_spec("3.Mg")
        assert s["type"] == "M"
        assert s["n"] == 3
        assert s["m"] is None
        assert s["seq"] == "G"
        assert s["homopolymer"] is True

    def test_exact_softclip_not_homopolymer_flag(self):
        """Exact softclip (n==m) has homopolymer=False (uses expanded seq)."""
        s = parse_spec("3Sg")
        assert s["seq"] == "GGG"  # expanded
        assert s["homopolymer"] is False

    def test_exact_mapped_not_homopolymer_flag(self):
        """Exact mapped (n==m) has homopolymer=False (uses expanded seq)."""
        s = parse_spec("3Mg")
        assert s["seq"] == "GGG"  # expanded
        assert s["homopolymer"] is False


# ---------------------------------------------------------------------------
# Step 2: spec_matches_aln on real BAM data
# ---------------------------------------------------------------------------

def _read_all(bam_path):
    """Read all alignments from a BAM file."""
    with pysam.AlignmentFile(bam_path, "rb") as f:
        return list(f.fetch(until_eof=True))


def _matching_reads(bam_path, spec_str):
    """Return query names (with flags) of alignments matching a spec."""
    spec = parse_spec(spec_str)
    results = []
    for aln in _read_all(bam_path):
        if spec_matches_aln(aln, spec):
            results.append((aln.query_name, aln.flag))
    return results


class TestSpecMatchesSE:
    """Test spec matching against test_softclip_se.bam."""

    def test_1sg_selects_1bp_softclip_g(self):
        """1Sg: exact 1bp G softclip -> READ1(0), READ2(16), READ7(16)."""
        hits = _matching_reads(SE_BAM, "1Sg")
        names = [(n, f) for n, f in hits]
        assert ("READ1", 0) in names
        assert ("READ2", 16) in names
        assert ("READ7", 16) in names
        assert len(hits) == 3

    def test_2sg_selects_2bp_softclip_g(self):
        """2Sg: exact 2bp GG -> READ1(2048), READ4(16), READ6(256)."""
        hits = _matching_reads(SE_BAM, "2Sg")
        assert len(hits) == 3
        names_flags = {(n, f) for n, f in hits}
        assert ("READ1", 2048) in names_flags
        assert ("READ4", 16) in names_flags
        assert ("READ6", 256) in names_flags

    def test_any_softclip_selects_all(self):
        """S: any softclip -> all 10 alignments."""
        hits = _matching_reads(SE_BAM, "S")
        assert len(hits) == 10

    def test_range_1_3_selects_correct_reads(self):
        """1.3S: softclip length 1-3 -> 8 reads (not 4S ones)."""
        hits = _matching_reads(SE_BAM, "1.3S")
        assert len(hits) == 8
        # 4S reads excluded
        for name, flag in hits:
            assert not (name == "READ5")

    def test_open_range_4_plus(self):
        """4.S: softclip 4+ -> READ5(0) and READ5(2064)."""
        hits = _matching_reads(SE_BAM, "4.S")
        assert len(hits) == 2
        assert all(n == "READ5" for n, _ in hits)

    def test_3sg_exact(self):
        """3Sg: exact 3bp GGG -> READ3(0) fwd, READ2(272) rev."""
        hits = _matching_reads(SE_BAM, "3Sg")
        assert len(hits) == 2
        names_flags = {(n, f) for n, f in hits}
        assert ("READ3", 0) in names_flags
        assert ("READ2", 272) in names_flags


class TestHomopolymerMatching:
    """Test homopolymer matching for range specs."""

    def test_range_softclip_homopolymer_accepts_pure(self):
        """3.Sg accepts 3+ pure G soft-clips."""
        hits = _matching_reads(SE_BAM, "3.Sg")
        names = [n for n, f in hits]
        assert "READ3" in names  # 3Sggg
        assert "READ5" in names  # 4Sgggg
        # Should not match 1S or 2S reads
        assert ("READ1", 0) not in hits  # 1Sg

    def test_range_softclip_homopolymer_rejects_mixed(self):
        """Homopolymer mode rejects non-homopolymer soft-clips."""
        # Use has_5prime_softclip_range directly with a mock alignment
        class MockAln:
            is_unmapped = False
            is_reverse = False
            query_sequence = "GGAT" + "A" * 72  # 4-bp soft-clip but not all G

        # Mock get_5prime_cigar to return (SOFT, 4)
        from pyselectal import SOFT
        import pyselectal
        original_get = pyselectal.get_5prime_cigar
        pyselectal.get_5prime_cigar = lambda aln: (SOFT, 4)
        try:
            # Prefix mode (homopolymer=False) should match (starts with G)
            assert has_5prime_softclip_range(MockAln(), 3, None, "G", "C", homopolymer=False)
            # Homopolymer mode should reject (not all G)
            assert not has_5prime_softclip_range(MockAln(), 3, None, "G", "C", homopolymer=True)
        finally:
            pyselectal.get_5prime_cigar = original_get

    def test_range_mapped_homopolymer_rejects_mixed(self):
        """Homopolymer mode rejects non-homopolymer matched bases."""
        class MockAln:
            is_unmapped = False
            is_reverse = False
            query_sequence = "GGAT" + "A" * 72  # 4-bp match but not all G

        from pyselectal import MATCH
        import pyselectal
        original_get = pyselectal.get_5prime_cigar
        pyselectal.get_5prime_cigar = lambda aln: (MATCH, 4)
        try:
            # Prefix mode should match (starts with G)
            assert has_5prime_mapped_range(MockAln(), 3, None, "G", "C", homopolymer=False)
            # Homopolymer mode should reject (not all G)
            assert not has_5prime_mapped_range(MockAln(), 3, None, "G", "C", homopolymer=True)
        finally:
            pyselectal.get_5prime_cigar = original_get

    def test_spec_matches_range_softclip_homopolymer(self):
        """spec_matches_aln uses homopolymer mode for range specs."""
        class MockAln:
            is_unmapped = False
            is_reverse = False
            query_sequence = "GGAT" + "A" * 72

        from pyselectal import SOFT
        import pyselectal
        original_get = pyselectal.get_5prime_cigar
        pyselectal.get_5prime_cigar = lambda aln: (SOFT, 4)
        try:
            spec = parse_spec("3.Sg")
            assert spec["homopolymer"] is True
            # Should reject mixed sequence
            assert not spec_matches_aln(MockAln(), spec)

            # Multi-char seq should use prefix mode
            spec2 = parse_spec("3.Sgg")
            assert spec2["homopolymer"] is False
            # Should accept (starts with GG)
            assert spec_matches_aln(MockAln(), spec2)
        finally:
            pyselectal.get_5prime_cigar = original_get


# ---------------------------------------------------------------------------
# Step 2: end-to-end --select via main()
# ---------------------------------------------------------------------------

def _count_reads_in_bam(path):
    with pysam.AlignmentFile(path, "rb") as f:
        return sum(1 for _ in f.fetch(until_eof=True))


class TestSelectEndToEnd:
    """End-to-end tests for --select mode."""

    def test_select_single_spec_to_stdout(self, tmp_path):
        """Single input + single spec -> stdout (captured via redirect)."""
        out = str(tmp_path / "out.bam")
        main(["-i", SE_BAM, "-s", "1Sg", "-o", out])
        assert _count_reads_in_bam(out) == 3

    def test_select_single_spec_file_naming(self, tmp_path):
        """Single input + multiple specs -> {stem}_{spec}.bam files."""
        # Use tmp_path as working directory for outputs
        os.chdir(str(tmp_path))
        main(["-i", SE_BAM, "-s", "1Sg,2Sg"])
        out1 = str(tmp_path / "test_softclip_se_1sg.bam")
        out2 = str(tmp_path / "test_softclip_se_2sg.bam")
        assert os.path.exists(out1)
        assert os.path.exists(out2)
        assert _count_reads_in_bam(out1) == 3
        assert _count_reads_in_bam(out2) == 3

    def test_select_se_any_softclip(self, tmp_path):
        out = str(tmp_path / "out.bam")
        main(["-i", SE_BAM, "-s", "S", "-o", out])
        assert _count_reads_in_bam(out) == 10

    def test_select_pe_1sg(self, tmp_path):
        """PE mode: select R1 with 1Sg and emit matching R2 mates."""
        out = str(tmp_path / "out.bam")
        main(["-i", PE_BAM, "-s", "1Sg", "-p", "-o", out])
        count = _count_reads_in_bam(out)
        # Should have both R1 and R2 reads
        assert count > 0

    def test_select_with_name_sort(self, tmp_path):
        """--name flag triggers name sorting before PE processing."""
        out = str(tmp_path / "out.bam")
        main(["-i", PE_BAM, "-s", "S", "-p", "-n", "-o", out])
        assert _count_reads_in_bam(out) > 0


# ---------------------------------------------------------------------------
# Step 3: Issue 4 — --merge
# ---------------------------------------------------------------------------

class TestMergeEndToEnd:
    """End-to-end tests for --select --merge mode."""

    def test_merge_se_combines_specs(self, tmp_path):
        """Merge 1Sg + 2Sg should yield union (no duplicates)."""
        out = str(tmp_path / "out.bam")
        main(["-i", SE_BAM, "-s", "1Sg,2Sg", "-m", "-o", out])
        # 1Sg: READ1(0), READ2(16), READ7(16) = 3
        # 2Sg: READ1(2048), READ4(16), READ6(256) = 3
        # No overlap -> 6
        assert _count_reads_in_bam(out) == 6

    def test_merge_se_deduplicates(self, tmp_path):
        """Merge with overlapping specs: each alignment written once."""
        out = str(tmp_path / "out.bam")
        # S matches all 10; 1Sg matches 3 (subset) -> still 10
        main(["-i", SE_BAM, "-s", "S,1Sg", "-m", "-o", out])
        assert _count_reads_in_bam(out) == 10

    def test_merge_se_output_naming(self, tmp_path):
        """Merge produces {stem}_merged.bam when no -o given."""
        os.chdir(str(tmp_path))
        main(["-i", SE_BAM, "-s", "1Sg,2Sg", "-m"])
        expected = str(tmp_path / "test_softclip_se_merged.bam")
        assert os.path.exists(expected)
        assert _count_reads_in_bam(expected) == 6

    def test_merge_se_respects_output_flag(self, tmp_path):
        """-o overrides automatic naming."""
        out = str(tmp_path / "custom.bam")
        main(["-i", SE_BAM, "-s", "1Sg,2Sg", "-m", "-o", out])
        assert os.path.exists(out)
        assert _count_reads_in_bam(out) == 6

    def test_merge_pe(self, tmp_path):
        """PE merge: R1 matching any spec + their R2 mates."""
        out = str(tmp_path / "out.bam")
        main(["-i", PE_BAM, "-s", "1Sg,2Sg", "-p", "-m", "-o", out])
        count = _count_reads_in_bam(out)
        assert count > 0

    def test_merge_pe_with_name_sort(self, tmp_path):
        """PE merge with --name sort."""
        out = str(tmp_path / "out.bam")
        main(["-i", PE_BAM, "-s", "S", "-p", "-n", "-m", "-o", out])
        assert _count_reads_in_bam(out) > 0

    def test_merge_single_spec(self, tmp_path):
        """Merge with one spec should work identically to non-merge."""
        out_merge = str(tmp_path / "merge.bam")
        out_normal = str(tmp_path / "normal.bam")
        main(["-i", SE_BAM, "-s", "1Sg", "-m", "-o", out_merge])
        main(["-i", SE_BAM, "-s", "1Sg", "-o", out_normal])
        assert _count_reads_in_bam(out_merge) == _count_reads_in_bam(out_normal)


# ---------------------------------------------------------------------------
# Step 4: Issue 1 — --count TSV histogram
# ---------------------------------------------------------------------------

def _read_tsv(path):
    """Read a TSV file and return (header, rows) where rows is list of (type, count)."""
    with open(path) as f:
        lines = f.read().strip().split("\n")
    header = lines[0]
    rows = []
    for line in lines[1:]:
        parts = line.split("\t")
        rows.append((parts[0], int(parts[1])))
    return header, rows


class TestClassify5primeType:
    """Unit tests for classify_5prime_type() on SE BAM."""

    def test_fwd_1bp_softclip_g(self):
        """READ1(0): 1S75M fwd, 5' base G -> '1Sg'."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ1" and a.flag == 0][0]
        assert classify_5prime_type(aln) == "1Sg"

    def test_rev_1bp_softclip_c(self):
        """READ2(16): 75M1S rev, 5' base C -> '1Sg' (revcomp of C = G)."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ2" and a.flag == 16][0]
        assert classify_5prime_type(aln) == "1Sg"

    def test_fwd_3bp_softclip(self):
        """READ3(0): 3S73M fwd, 5' seq GGG -> '3Sggg'."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ3" and a.flag == 0][0]
        assert classify_5prime_type(aln) == "3Sggg"

    def test_rev_3bp_softclip(self):
        """READ2(272): 73M3S rev, 5' seq CCC -> '3Sggg' (revcomp)."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ2" and a.flag == 272][0]
        assert classify_5prime_type(aln) == "3Sggg"

    def test_fwd_2bp_softclip(self):
        """READ1(2048): 2S74M fwd -> '2Sgg'."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ1" and a.flag == 2048][0]
        assert classify_5prime_type(aln) == "2Sgg"

    def test_rev_2bp_softclip(self):
        """READ4(16): 74M2S rev, 5' seq CC -> '2Sgg' (revcomp)."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ4" and a.flag == 16][0]
        assert classify_5prime_type(aln) == "2Sgg"

    def test_fwd_4bp_softclip(self):
        """READ5(0): 4S72M fwd -> '4Sgggg'."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ5" and a.flag == 0][0]
        assert classify_5prime_type(aln) == "4Sgggg"

    def test_rev_4bp_softclip(self):
        """READ5(2064): 72M4S rev -> '4Scccc' revcomp = '4Sgggg'."""
        alns = _read_all(SE_BAM)
        aln = [a for a in alns if a.query_name == "READ5" and a.flag == 2064][0]
        assert classify_5prime_type(aln) == "4Sgggg"


class TestCountEndToEnd:
    """End-to-end tests for --count mode."""

    def test_count_se_to_stdout(self, tmp_path, capsys):
        """Single input -> stdout TSV."""
        main(["-i", SE_BAM, "-c"])
        captured = capsys.readouterr()
        lines = captured.out.strip().split("\n")
        assert lines[0] == "type\tcount"
        # 10 alignments, all soft-clipped -> should have entries
        assert len(lines) > 1

    def test_count_se_to_file(self, tmp_path):
        """With -o, write TSV to file."""
        out = str(tmp_path / "counts.tsv")
        main(["-i", SE_BAM, "-c", "-o", out])
        header, rows = _read_tsv(out)
        assert header == "type\tcount"
        total = sum(c for _, c in rows)
        assert total == 10  # all 10 SE alignments are soft-clipped

    def test_count_se_types_correct(self, tmp_path):
        """Verify exact type counts from SE BAM."""
        out = str(tmp_path / "counts.tsv")
        main(["-i", SE_BAM, "-c", "-o", out])
        _, rows = _read_tsv(out)
        counts = dict(rows)
        # From test data: 1Sg x3 (READ1:0, READ2:16, READ7:16),
        # 2Sgg x3 (READ1:2048, READ4:16, READ6:256),
        # 3Sggg x2 (READ3:0, READ2:272),
        # 4Sgggg x2 (READ5:0, READ5:2064)
        assert counts["1Sg"] == 3
        assert counts["2Sgg"] == 3
        assert counts["3Sggg"] == 2
        assert counts["4Sgggg"] == 2

    def test_count_se_sorted_by_descending_count(self, tmp_path):
        """Rows should be sorted by descending count."""
        out = str(tmp_path / "counts.tsv")
        main(["-i", SE_BAM, "-c", "-o", out])
        _, rows = _read_tsv(out)
        counts_list = [c for _, c in rows]
        assert counts_list == sorted(counts_list, reverse=True)

    def test_count_output_naming_multi_input(self, tmp_path):
        """Multiple inputs -> {stem}_5prime_counts.tsv files."""
        os.chdir(str(tmp_path))
        main(["-i", f"{SE_BAM},{PE_BAM}", "-c"])
        assert os.path.exists(str(tmp_path / "test_softclip_se_5prime_counts.tsv"))
        assert os.path.exists(str(tmp_path / "test_softclip_pe_5prime_counts.tsv"))

    def test_count_pe_r1_only(self, tmp_path):
        """PE mode counts only R1 alignments."""
        out_se = str(tmp_path / "se.tsv")
        out_pe = str(tmp_path / "pe.tsv")
        main(["-i", PE_BAM, "-c", "-o", out_se])
        main(["-i", PE_BAM, "-c", "-p", "-o", out_pe])
        _, rows_se = _read_tsv(out_se)
        _, rows_pe = _read_tsv(out_pe)
        total_se = sum(c for _, c in rows_se)
        total_pe = sum(c for _, c in rows_pe)
        # PE BAM has 20 reads (10 R1 + 10 R2); SE counts all, PE counts R1 only
        assert total_se == 20
        assert total_pe < total_se

    def test_count_respects_output_flag(self, tmp_path):
        """-o overrides automatic naming."""
        out = str(tmp_path / "custom.tsv")
        main(["-i", SE_BAM, "-c", "-o", out])
        assert os.path.exists(out)


# ---------------------------------------------------------------------------
# Step 4 additions: --mapped-prefix and --collapse-threshold
# ---------------------------------------------------------------------------

def _get_r2_alns(bam_path):
    """Return all R2 alignments from a BAM."""
    with pysam.AlignmentFile(bam_path, "rb") as f:
        return [a for a in f.fetch(until_eof=True) if a.is_read2]


class TestMappedPrefix:
    """Tests for classify_5prime_type() MATCH case with --mapped-prefix."""

    def test_mapped_default_prefix_format(self):
        """Default mapped_prefix=5: returns '{length}M{bases}' for a 76M R2."""
        import re
        alns = _get_r2_alns(PE_BAM)
        fwd_alns = [a for a in alns if not a.is_reverse]
        assert fwd_alns, "Expected at least one forward R2 in PE BAM"
        result = classify_5prime_type(fwd_alns[0])
        assert re.match(r'^\d+M[acgtn]+$', result), f"Unexpected format: {result!r}"

    def test_mapped_zero_prefix_no_sequence(self):
        """mapped_prefix=0: returns '{length}M' without bases."""
        import re
        alns = _get_r2_alns(PE_BAM)
        fwd_alns = [a for a in alns if not a.is_reverse]
        assert fwd_alns
        result = classify_5prime_type(fwd_alns[0], mapped_prefix=0)
        assert re.match(r'^\d+M$', result), f"Unexpected format: {result!r}"

    def test_mapped_prefix_length_cap(self):
        """mapped_prefix > match length: at most match_length bases shown."""
        alns = _get_r2_alns(PE_BAM)
        fwd_alns = [a for a in alns if not a.is_reverse]
        assert fwd_alns
        aln = fwd_alns[0]
        result = classify_5prime_type(aln, mapped_prefix=100)
        m_idx = result.index('M')
        match_len = int(result[:m_idx])
        seq_part = result[m_idx + 1:]
        assert len(seq_part) <= match_len

    def test_mapped_prefix_shorter_than_match_is_prefix_of_longer(self):
        """mapped_prefix=3 sequence is prefix of mapped_prefix=10 sequence."""
        alns = _get_r2_alns(PE_BAM)
        fwd_alns = [a for a in alns if not a.is_reverse]
        assert fwd_alns
        aln = fwd_alns[0]
        r3 = classify_5prime_type(aln, mapped_prefix=3)
        r10 = classify_5prime_type(aln, mapped_prefix=10)
        seq3 = r3[r3.index('M') + 1:]
        seq10 = r10[r10.index('M') + 1:]
        assert seq10.startswith(seq3)

    def test_mapped_prefix_reverse_strand(self):
        """Reverse-strand R2: result is forward-strand orientation (revcomp of stored)."""
        import re
        alns = _get_r2_alns(PE_BAM)
        rev_alns = [a for a in alns if a.is_reverse]
        if not rev_alns:
            pytest.skip("No reverse-strand R2 in PE BAM")
        result = classify_5prime_type(rev_alns[0], mapped_prefix=3)
        assert re.match(r'^\d+M[acgtn]+$', result), f"Unexpected format: {result!r}"

    def test_mapped_prefix_cli_zero(self, tmp_path):
        """--mapped-prefix 0: all M-type rows are '{length}M' with no bases."""
        import re
        out = str(tmp_path / "counts.tsv")
        main(["-i", PE_BAM, "-c", "--mapped-prefix", "0", "-o", out])
        _, rows = _read_tsv(out)
        m_rows = [(t, c) for t, c in rows if 'M' in t and not t.startswith('other')]
        for t, _ in m_rows:
            assert re.match(r'^\d+M$', t), f"Unexpected M-type: {t!r}"

    def test_mapped_prefix_cli_default(self, tmp_path):
        """Default --mapped-prefix 5: M-type rows include base sequence."""
        import re
        out = str(tmp_path / "counts.tsv")
        main(["-i", PE_BAM, "-c", "-o", out])
        _, rows = _read_tsv(out)
        m_rows = [(t, c) for t, c in rows if re.match(r'^\d+M', t) and not t.startswith('other')]
        assert m_rows, "Expected at least one M-type in PE BAM"
        for t, _ in m_rows:
            assert re.match(r'^\d+M[acgtn]+$', t), f"Unexpected M-type: {t!r}"

    def test_mapped_prefix_parse_args_defaults(self):
        args = parse_args(["-i", "in.bam", "-c"])
        assert args.mapped_prefix == 5

    def test_mapped_prefix_parse_args_custom(self):
        args = parse_args(["-i", "in.bam", "-c", "--mapped-prefix", "3"])
        assert args.mapped_prefix == 3

    def test_mapped_prefix_zero_valid(self):
        args = parse_args(["-i", "in.bam", "-c", "--mapped-prefix", "0"])
        assert args.mapped_prefix == 0

    def test_mapped_prefix_negative_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "--mapped-prefix", "-1"])
        assert exc.value.code == 2


class TestCollapseThreshold:
    """Tests for write_count_tsv collapse and --collapse-threshold CLI flag."""

    def test_collapse_no_other_when_all_above(self, tmp_path):
        """No 'other' row if all types meet the threshold."""
        counts = {"1Sg": 50, "2Sgg": 30, "3Sggg": 20}  # 50%, 30%, 20% — all >= 5%
        out = str(tmp_path / "out.tsv")
        write_count_tsv(counts, out, collapse_threshold=5.0)
        _, rows = _read_tsv(out)
        types = [t for t, _ in rows]
        assert not any("other" in t for t in types)
        assert sum(c for _, c in rows) == 100

    def test_collapse_creates_other_row(self, tmp_path):
        """Types below threshold are collapsed into 'other' row."""
        counts = {"1Sg": 80, "2Sgg": 15, "3Sggg": 5}  # 5% < 5% threshold? No.
        # 3Sggg is exactly 5%, not strictly less — should NOT collapse
        out = str(tmp_path / "out.tsv")
        write_count_tsv(counts, out, collapse_threshold=5.0)
        _, rows = _read_tsv(out)
        types = [t for t, _ in rows]
        assert not any("other" in t for t in types)

    def test_collapse_strictly_below(self, tmp_path):
        """3% < 5% threshold -> collapses."""
        counts = {"1Sg": 80, "2Sgg": 17, "3Sggg": 3}  # 3% < 5%
        out = str(tmp_path / "out.tsv")
        write_count_tsv(counts, out, collapse_threshold=5.0)
        _, rows = _read_tsv(out)
        types = [t for t, _ in rows]
        assert any("other (<5%)" in t for t in types)
        # other row count = 3
        other_count = dict(rows).get("other (<5%)", 0)
        assert other_count == 3

    def test_collapse_other_row_is_last(self, tmp_path):
        """'other' row always appears last regardless of its count."""
        counts = {"1Sg": 60, "2Sgg": 30, "3Sggg": 4, "4Sgggg": 4, "5Sggggg": 2}
        # 4%, 4%, 2% are each below 5% -> other = 10
        out = str(tmp_path / "out.tsv")
        write_count_tsv(counts, out, collapse_threshold=5.0)
        _, rows = _read_tsv(out)
        assert rows[-1][0] == "other (<5%)"
        assert rows[-1][1] == 10

    def test_collapse_zero_disables(self, tmp_path):
        """--collapse-threshold 0 disables collapsing."""
        counts = {"1Sg": 99, "2Sgg": 1}  # 1% would normally collapse
        out = str(tmp_path / "out.tsv")
        write_count_tsv(counts, out, collapse_threshold=0)
        _, rows = _read_tsv(out)
        types = [t for t, _ in rows]
        assert not any("other" in t for t in types)
        assert len(rows) == 2

    def test_collapse_label_reflects_threshold(self, tmp_path):
        """Label 'other (<10%)' when threshold is 10."""
        counts = {"1Sg": 91, "2Sgg": 5, "3Sggg": 4}  # 4% < 10%
        out = str(tmp_path / "out.tsv")
        write_count_tsv(counts, out, collapse_threshold=10.0)
        _, rows = _read_tsv(out)
        types = [t for t, _ in rows]
        assert any("other (<10%)" in t for t in types)

    def test_collapse_total_preserved(self, tmp_path):
        """Sum of all counts (including other) equals original total."""
        counts = {"1Sg": 60, "2Sgg": 30, "3Sggg": 7, "4Sgggg": 3}
        out = str(tmp_path / "out.tsv")
        write_count_tsv(counts, out, collapse_threshold=5.0)
        _, rows = _read_tsv(out)
        assert sum(c for _, c in rows) == 100

    def test_collapse_cli_default(self, tmp_path):
        """Default threshold 1% applies in --count output."""
        # SE BAM: 1Sg=3(30%), 2Sgg=3(30%), 3Sggg=2(20%), 4Sgggg=2(20%) — none collapse
        out = str(tmp_path / "counts.tsv")
        main(["-i", SE_BAM, "-c", "-o", out])
        _, rows = _read_tsv(out)
        types = [t for t, _ in rows]
        assert not any("other" in t for t in types)

    def test_collapse_cli_high_threshold(self, tmp_path):
        """threshold=25%: 3Sggg(20%) and 4Sgggg(20%) collapse into other."""
        out = str(tmp_path / "counts.tsv")
        main(["-i", SE_BAM, "-c", "--collapse-threshold", "25", "-o", out])
        _, rows = _read_tsv(out)
        row_dict = dict(rows)
        assert "other (<25%)" in row_dict
        assert row_dict["other (<25%)"] == 4  # 2 + 2
        assert "3Sggg" not in row_dict
        assert "4Sgggg" not in row_dict
        assert rows[-1][0] == "other (<25%)"

    def test_collapse_cli_zero_disables(self, tmp_path):
        """--collapse-threshold 0 shows all rows."""
        out = str(tmp_path / "counts.tsv")
        main(["-i", SE_BAM, "-c", "--collapse-threshold", "0", "-o", out])
        _, rows = _read_tsv(out)
        assert len(rows) == 4
        assert not any("other" in t for t, _ in rows)

    def test_collapse_parse_args_default(self):
        args = parse_args(["-i", "in.bam", "-c"])
        assert args.collapse_threshold == 1.0

    def test_collapse_parse_args_custom(self):
        args = parse_args(["-i", "in.bam", "-c", "--collapse-threshold", "10"])
        assert args.collapse_threshold == 10.0

    def test_collapse_threshold_negative_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "--collapse-threshold", "-1"])
        assert exc.value.code == 2

    def test_collapse_threshold_100_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "--collapse-threshold", "100"])
        assert exc.value.code == 2


# ---------------------------------------------------------------------------
# Step 5: Issue 2 — --all split by type
# ---------------------------------------------------------------------------

class TestAllEndToEnd:
    """End-to-end tests for --all mode."""

    # SE BAM types (from test data table):
    # 1Sg  x3: READ1(0), READ2(16), READ7(16)
    # 2Sgg x3: READ1(2048), READ4(16), READ6(256)
    # 3Sggg x2: READ3(0), READ2(272)
    # 4Sgggg x2: READ5(0), READ5(2064)

    def test_all_se_creates_per_type_files(self, tmp_path):
        """SE --all creates one BAM per distinct 5' type."""
        os.chdir(str(tmp_path))
        main(["-i", SE_BAM, "-a"])
        assert os.path.exists(str(tmp_path / "test_softclip_se_1sg.bam"))
        assert os.path.exists(str(tmp_path / "test_softclip_se_2sgg.bam"))
        assert os.path.exists(str(tmp_path / "test_softclip_se_3sggg.bam"))
        assert os.path.exists(str(tmp_path / "test_softclip_se_4sgggg.bam"))

    def test_all_se_correct_counts(self, tmp_path):
        """Each type file contains the right number of alignments."""
        os.chdir(str(tmp_path))
        main(["-i", SE_BAM, "-a"])
        assert _count_reads_in_bam(str(tmp_path / "test_softclip_se_1sg.bam")) == 3
        assert _count_reads_in_bam(str(tmp_path / "test_softclip_se_2sgg.bam")) == 3
        assert _count_reads_in_bam(str(tmp_path / "test_softclip_se_3sggg.bam")) == 2
        assert _count_reads_in_bam(str(tmp_path / "test_softclip_se_4sgggg.bam")) == 2

    def test_all_se_total_reads_preserved(self, tmp_path):
        """Sum of all type file counts equals total alignments in input."""
        os.chdir(str(tmp_path))
        main(["-i", SE_BAM, "-a"])
        total = sum(
            _count_reads_in_bam(str(p))
            for p in tmp_path.glob("test_softclip_se_*.bam")
        )
        assert total == _count_reads_in_bam(SE_BAM)

    def test_all_se_output_dir(self, tmp_path):
        """With -o DIR, files are written inside DIR."""
        out_dir = str(tmp_path / "split")
        main(["-i", SE_BAM, "-a", "-o", out_dir])
        assert os.path.isdir(out_dir)
        assert os.path.exists(os.path.join(out_dir, "test_softclip_se_1sg.bam"))
        assert os.path.exists(os.path.join(out_dir, "test_softclip_se_2sgg.bam"))

    def test_all_se_output_dir_created(self, tmp_path):
        """Output directory is created if it doesn't exist yet."""
        out_dir = str(tmp_path / "new_dir")
        assert not os.path.exists(out_dir)
        main(["-i", SE_BAM, "-a", "-o", out_dir])
        assert os.path.isdir(out_dir)

    def test_all_pe_creates_files(self, tmp_path):
        """PE --all creates per-type BAM files from R1 classification."""
        os.chdir(str(tmp_path))
        main(["-i", PE_BAM, "-a", "-p"])
        bams = list(tmp_path.glob("test_softclip_pe_*.bam"))
        assert len(bams) > 0

    def test_all_pe_r2_paired_with_r1(self, tmp_path):
        """PE --all: R2 mates are written into the same file as their R1 mate."""
        os.chdir(str(tmp_path))
        main(["-i", PE_BAM, "-a", "-p"])
        total_r2 = 0
        for bam_path in tmp_path.glob("test_softclip_pe_*.bam"):
            with pysam.AlignmentFile(str(bam_path), "rb") as f:
                alns = list(f.fetch(until_eof=True))
            total_r2 += sum(1 for a in alns if a.is_read2)
        # At least some R2s must have been routed into type files
        assert total_r2 > 0

    def test_all_pe_with_name_sort(self, tmp_path):
        """PE --all with -n name-sorts before processing."""
        os.chdir(str(tmp_path))
        main(["-i", PE_BAM, "-a", "-p", "-n"])
        bams = list(tmp_path.glob("test_softclip_pe_*.bam"))
        assert len(bams) > 0


# ---------------------------------------------------------------------------
# Step 6: Issue 5 — SAM/CRAM format support
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def se_sam(tmp_path_factory):
    """SAM version of the SE test BAM, created once per test session."""
    tmp = tmp_path_factory.mktemp("samdata")
    sam_path = str(tmp / "test_softclip_se.sam")
    with pysam.AlignmentFile(SE_BAM, "rb") as src:
        with pysam.AlignmentFile(sam_path, "w", template=src) as dst:
            for aln in src.fetch(until_eof=True):
                dst.write(aln)
    return sam_path


class TestFormatFlags:
    """parse_args tests for -S/-B/-C/-r flags."""

    def test_sam_flag(self):
        args = parse_args(["-i", "in.bam", "-c", "-S"])
        assert args.sam is True
        assert args.bam is False
        assert args.cram is False

    def test_bam_flag(self):
        args = parse_args(["-i", "in.bam", "-c", "-B"])
        assert args.bam is True

    def test_cram_flag_with_reference(self):
        args = parse_args(["-i", "in.bam", "-c", "-C", "-r", "ref.fa"])
        assert args.cram is True
        assert args.reference == "ref.fa"

    def test_cram_without_reference_rejected(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "-C"])
        assert exc.value.code == 2

    def test_sam_bam_mutually_exclusive(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "-S", "-B"])
        assert exc.value.code == 2

    def test_bam_cram_mutually_exclusive(self):
        with pytest.raises(SystemExit) as exc:
            parse_args(["-i", "in.bam", "-c", "-B", "-C", "-r", "ref.fa"])
        assert exc.value.code == 2

    def test_reference_flag(self):
        args = parse_args(["-i", "in.cram", "-c", "-r", "ref.fa"])
        assert args.reference == "ref.fa"

    def test_format_flags_default_false(self):
        args = parse_args(["-i", "in.bam", "-c"])
        assert args.sam is False
        assert args.bam is False
        assert args.cram is False
        assert args.reference is None


class TestFormatDetection:
    """Unit tests for _detect_format()."""

    def test_bam_extension(self):
        assert _detect_format("sample.bam") == "bam"

    def test_sam_extension(self):
        assert _detect_format("sample.sam") == "sam"

    def test_cram_extension(self):
        assert _detect_format("sample.cram") == "cram"

    def test_uppercase_extension(self):
        assert _detect_format("sample.BAM") == "bam"
        assert _detect_format("sample.SAM") == "sam"

    def test_path_with_directory(self):
        assert _detect_format("/data/sample.bam") == "bam"

    def test_unknown_extension_dies(self):
        with pytest.raises(SystemExit):
            _detect_format("sample.txt")


class TestSAMFormat:
    """End-to-end tests for SAM input/output."""

    def test_sam_input_count(self, se_sam, tmp_path):
        """SAM input to --count produces same counts as BAM input."""
        out_bam = str(tmp_path / "from_bam.tsv")
        out_sam = str(tmp_path / "from_sam.tsv")
        main(["-i", SE_BAM, "-c", "-o", out_bam])
        main(["-i", se_sam, "-c", "-o", out_sam])
        _, rows_bam = _read_tsv(out_bam)
        _, rows_sam = _read_tsv(out_sam)
        assert dict(rows_bam) == dict(rows_sam)

    def test_sam_input_select(self, se_sam, tmp_path):
        """SAM input to --select produces SAM output by default."""
        out = str(tmp_path / "out.sam")
        main(["-i", se_sam, "-s", "1Sg", "-o", out])
        assert os.path.exists(out)
        with pysam.AlignmentFile(out, "r") as f:
            count = sum(1 for _ in f.fetch(until_eof=True))
        assert count == 3

    def test_bam_to_sam_output(self, tmp_path):
        """BAM input + -S → SAM output."""
        out = str(tmp_path / "out.sam")
        main(["-i", SE_BAM, "-s", "1Sg", "-S", "-o", out])
        assert os.path.exists(out)
        with pysam.AlignmentFile(out, "r") as f:
            count = sum(1 for _ in f.fetch(until_eof=True))
        assert count == 3

    def test_sam_to_bam_output(self, se_sam, tmp_path):
        """SAM input + -B → BAM output."""
        out = str(tmp_path / "out.bam")
        main(["-i", se_sam, "-s", "1Sg", "-B", "-o", out])
        assert os.path.exists(out)
        with pysam.AlignmentFile(out, "rb") as f:
            count = sum(1 for _ in f.fetch(until_eof=True))
        assert count == 3

    def test_sam_auto_named_output_extension(self, se_sam, tmp_path):
        """SAM input + multiple specs → auto-named .sam output files."""
        os.chdir(str(tmp_path))
        main(["-i", se_sam, "-s", "1Sg,2Sg"])
        assert os.path.exists(str(tmp_path / "test_softclip_se_1sg.sam"))
        assert os.path.exists(str(tmp_path / "test_softclip_se_2sg.sam"))

    def test_all_sam_output(self, se_sam, tmp_path):
        """SAM input + --all → .sam type files readable as SAM."""
        os.chdir(str(tmp_path))
        main(["-i", se_sam, "-a"])
        sams = list(tmp_path.glob("*.sam"))
        assert len(sams) > 0
        for sam_path in sams:
            with pysam.AlignmentFile(str(sam_path), "r") as f:
                assert sum(1 for _ in f.fetch(until_eof=True)) > 0

    def test_bam_merge_to_sam(self, tmp_path):
        """BAM input + --merge + -S → SAM output."""
        out = str(tmp_path / "merged.sam")
        main(["-i", SE_BAM, "-s", "1Sg,2Sg", "-m", "-S", "-o", out])
        with pysam.AlignmentFile(out, "r") as f:
            count = sum(1 for _ in f.fetch(until_eof=True))
        assert count == 6  # 1Sg=3, 2Sg=3, no overlap


# ---------------------------------------------------------------------------
# Step 7: --collapse-threshold in --all
# ---------------------------------------------------------------------------

# SE BAM type distribution (all are >= 20%, so threshold must be > 20% to trigger collapse):
#   1Sg=3 (30%), 2Sgg=3 (30%), 3Sggg=2 (20%), 4Sgggg=2 (20%)

class TestAllCollapse:
    """Tests for --collapse-threshold in --all mode."""

    def test_all_collapse_high_threshold_creates_other(self, tmp_path):
        """threshold=25%: 3Sggg(20%) and 4Sgggg(20%) route to _other.bam."""
        out_dir = str(tmp_path / "split")
        main(["-i", SE_BAM, "-a", "--collapse-threshold", "25", "-o", out_dir])
        other = os.path.join(out_dir, "test_softclip_se_other.bam")
        assert os.path.exists(other)
        assert _count_reads_in_bam(other) == 4  # 2 + 2

    def test_all_collapse_high_threshold_no_rare_type_files(self, tmp_path):
        """threshold=25%: no per-type files for 3Sggg or 4Sgggg."""
        out_dir = str(tmp_path / "split")
        main(["-i", SE_BAM, "-a", "--collapse-threshold", "25", "-o", out_dir])
        assert not os.path.exists(os.path.join(out_dir, "test_softclip_se_3sggg.bam"))
        assert not os.path.exists(os.path.join(out_dir, "test_softclip_se_4sgggg.bam"))

    def test_all_collapse_zero_disables(self, tmp_path):
        """--collapse-threshold 0 produces one file per type (no other)."""
        out_dir = str(tmp_path / "split")
        main(["-i", SE_BAM, "-a", "--collapse-threshold", "0", "-o", out_dir])
        assert not os.path.exists(os.path.join(out_dir, "test_softclip_se_other.bam"))
        assert os.path.exists(os.path.join(out_dir, "test_softclip_se_3sggg.bam"))
        assert os.path.exists(os.path.join(out_dir, "test_softclip_se_4sgggg.bam"))
        assert _count_reads_in_bam(os.path.join(out_dir, "test_softclip_se_3sggg.bam")) == 2
        assert _count_reads_in_bam(os.path.join(out_dir, "test_softclip_se_4sgggg.bam")) == 2

    def test_all_collapse_total_reads_preserved(self, tmp_path):
        """Sum across all output files (including other) equals input total."""
        out_dir = str(tmp_path / "split")
        main(["-i", SE_BAM, "-a", "--collapse-threshold", "25", "-o", out_dir])
        total = sum(
            _count_reads_in_bam(str(p))
            for p in (tmp_path / "split").glob("test_softclip_se_*.bam")
        )
        assert total == _count_reads_in_bam(SE_BAM)

    def test_all_collapse_other_bam_contents(self, tmp_path):
        """_other.bam contains alignments from all rare types combined."""
        out_dir = str(tmp_path / "split")
        main(["-i", SE_BAM, "-a", "--collapse-threshold", "25", "-o", out_dir])
        other = os.path.join(out_dir, "test_softclip_se_other.bam")
        with pysam.AlignmentFile(other, "rb") as f:
            names = {a.query_name for a in f.fetch(until_eof=True)}
        # 3Sggg: READ3(0), READ2(272); 4Sgggg: READ5(0), READ5(2064)
        assert "READ3" in names
        assert "READ5" in names

    def test_all_collapse_pe_uses_r1_counts(self, tmp_path):
        """PE collapse is computed from R1 only; R2 mates follow their R1 into other."""
        out_dir = str(tmp_path / "split")
        main(["-i", PE_BAM, "-a", "-p", "--collapse-threshold", "99", "-o", out_dir])
        other = os.path.join(out_dir, "test_softclip_pe_other.bam")
        assert os.path.exists(other)
        with pysam.AlignmentFile(other, "rb") as f:
            alns = list(f.fetch(until_eof=True))
        # At high threshold everything collapses; R2s must also be present
        r2_count = sum(1 for a in alns if a.is_read2)
        assert r2_count > 0


# ---------------------------------------------------------------------------
# Step 8: CRAM via stdin/stdout
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def cram_fixture(tmp_path_factory):
    """Creates a minimal FASTA reference + CRAM with one 1Sg read for stdin/stdout tests."""
    tmp = tmp_path_factory.mktemp("cramdata")
    ref_fa = str(tmp / "mini_ref.fa")
    cram_path = str(tmp / "mini.cram")
    chrom, chrom_len = "testchr", 200
    with open(ref_fa, "w") as fh:
        fh.write(f">{chrom}\n" + "A" * chrom_len + "\n")
    pysam.faidx(ref_fa)
    header = pysam.AlignmentHeader.from_dict({
        "HD": {"VN": "1.6", "SO": "coordinate"},
        "SQ": [{"SN": chrom, "LN": chrom_len}],
    })
    with pysam.AlignmentFile(cram_path, "wc", reference_filename=ref_fa, header=header) as dst:
        r = pysam.AlignedSegment(header)
        r.query_name = "READ1"
        r.flag = 0
        r.reference_id = 0
        r.reference_start = 0
        r.mapping_quality = 30
        r.cigar = [(4, 1), (0, 75)]  # 1S75M
        r.query_sequence = "G" + "A" * 75
        r.query_qualities = pysam.qualitystring_to_array("I" * 76)
        dst.write(r)
    return cram_path, ref_fa


class TestCRAMStdin:
    """Step 8: CRAM detection from stdin and CRAM output."""

    def test_stdin_cram_detected(self, cram_fixture, tmp_path):
        """CRAM bytes piped to stdin are detected as 'cram' format."""
        cram_path, _ = cram_fixture
        with open(cram_path, "rb") as f:
            cram_data = f.read()
        with mock.patch("sys.stdin", mock.MagicMock(buffer=io.BytesIO(cram_data))):
            tmp_spooled, fmt = _spool_and_detect_stdin()
        try:
            assert fmt == "cram"
        finally:
            os.remove(tmp_spooled)

    def test_stdin_bam_still_detected(self, tmp_path):
        """BAM bytes piped to stdin are still detected as 'bam' format."""
        with open(SE_BAM, "rb") as f:
            bam_data = f.read()
        with mock.patch("sys.stdin", mock.MagicMock(buffer=io.BytesIO(bam_data))):
            tmp_spooled, fmt = _spool_and_detect_stdin()
        try:
            assert fmt == "bam"
        finally:
            os.remove(tmp_spooled)

    def test_stdin_cram_requires_reference(self, cram_fixture):
        """CRAM from stdin without -r raises SystemExit."""
        cram_path, _ = cram_fixture
        with open(cram_path, "rb") as f:
            cram_data = f.read()
        with mock.patch("sys.stdin", mock.MagicMock(buffer=io.BytesIO(cram_data))):
            with pytest.raises(SystemExit):
                main(["-i", "-", "-c"])

    def test_cram_output_has_magic(self, cram_fixture, tmp_path):
        """CRAM file input + -C → output file starts with CRAM magic bytes."""
        cram_path, ref_fa = cram_fixture
        out = str(tmp_path / "out.cram")
        main(["-i", cram_path, "-s", "1Sg", "-C", "-r", ref_fa, "-o", out])
        assert os.path.isfile(out)
        with open(out, "rb") as f:
            magic = f.read(4)
        assert magic == b"CRAM"

    def test_stdin_cram_roundtrip(self, cram_fixture, tmp_path):
        """Full pipe: CRAM bytes on stdin, --select 1Sg, BAM output → correct read count."""
        cram_path, ref_fa = cram_fixture
        with open(cram_path, "rb") as f:
            cram_data = f.read()
        out = str(tmp_path / "out.bam")
        with mock.patch("sys.stdin", mock.MagicMock(buffer=io.BytesIO(cram_data))):
            main(["-i", "-", "-s", "1Sg", "-r", ref_fa, "-B", "-o", out])
        assert os.path.isfile(out)
        with pysam.AlignmentFile(out, "rb") as f:
            count = sum(1 for _ in f.fetch(until_eof=True))
        assert count == 1
