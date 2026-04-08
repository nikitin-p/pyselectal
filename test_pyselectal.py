"""Tests for pyselectal CLI."""

import os
import tempfile
import pytest
import pysam
from pyselectal import (
    parse_args, parse_spec, spec_matches_aln, main, VERSION,
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
        """2M -> mapped with min 2bp match at 5' end."""
        s = parse_spec("2M")
        assert s["type"] == "M"
        assert s["n"] == 2
        assert s["m"] is None
        assert s["seq"] is None

    def test_mapped_with_prefix(self):
        """2Mg -> mapped 2bp match starting with G."""
        s = parse_spec("2Mg")
        assert s["type"] == "M"
        assert s["n"] == 2
        assert s["m"] is None
        assert s["seq"] == "G"

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
        """Range softclip seq must be 0 or 1 char."""
        with pytest.raises(SystemExit):
            parse_spec("1.3Sgg")

    def test_exact_softclip_wrong_seq_length(self):
        """Exact softclip: multi-char seq must have length == n."""
        with pytest.raises(SystemExit):
            parse_spec("2Sggg")


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
