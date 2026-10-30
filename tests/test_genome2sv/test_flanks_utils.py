"""genome2sv 侧翼序列纯函数测试|genome2sv flank helper tests"""

import pytest

from biopytools.genome2sv.utils import (
    FLANK_FORMAT_FIELDS,
    FaidxReader,
    flank_interval,
    format_flank_tsv,
    parse_sample_fields,
)


class TestFaidxReaderLengths:
    """FaidxReader 暴露染色体长度|FaidxReader exposes contig lengths"""

    def test_lengths(self, tmp_path):
        fa = tmp_path / "ref.fa"
        fa.write_text(">c1\n" + "ACGT" * 100 + "\n>c2\nAAAAAAAAAA\n")
        fa.with_suffix(".fa.fai").write_text(
            "c1\t400\t5\t400\t401\nc2\t10\t410\t10\t11\n")
        reader = FaidxReader(str(fa))
        try:
            assert reader.lengths == {"c1": 400, "c2": 10}
        finally:
            reader.close()


class TestFlankInterval:
    """flank_interval 坐标计算|flank interval coordinates"""

    def test_ins_centered(self):
        """INS(end=pos)插入位点在序列中部|INS centered (end=pos)"""
        assert flank_interval(5796, 5796, 50_000_000, 300) == (5496, 6096, False)

    def test_start_clamped(self):
        """上游越界截到 1 并标记|start clamped to 1, trimmed flagged"""
        assert flank_interval(100, 209, 5000, 300) == (1, 509, True)

    def test_end_clamped(self):
        """下游越界截到染色体长并标记|end clamped to contig length, flagged"""
        assert flank_interval(4800, 4900, 5000, 300) == (4500, 5000, True)

    def test_end_lt_pos_uses_min_max(self):
        """部分 INV/DUP 记录 END<POS,按 min/max 求覆盖区间|END<POS handled via min/max"""
        assert flank_interval(600, 100, 5000, 50) == (50, 650, False)

    def test_no_flank(self):
        """flank=0 时区间即 SV 本体|flank=0 keeps the SV interval only"""
        assert flank_interval(100, 209, 5000, 0) == (100, 209, False)


class TestParseSampleFields:
    """FORMAT 键驱动的样本字段解析|FORMAT-key-driven sample field parsing"""

    def test_simple(self):
        assert parse_sample_fields("1/1:57:60", ["GT", "LN", "QV"]) == (
            "1/1", "57", "60")

    def test_format_driven_not_positional(self):
        """SURVIVOR FORMAT 键序 GT:PSV:LN:DR:ST:QV...,按键取值|lookup by FORMAT keys"""
        fmt = "GT:PSV:LN:DR:ST:QV:TY:ID:RAL:AAL:CO".split(":")
        field = "0|1:0:52:.:+:61:DEL:1:A:T:c"
        assert parse_sample_fields(field, fmt) == ("0|1", "52", "61")

    def test_missing_trailing(self):
        """值少于 FORMAT 键,缺失补 .|fewer values than keys; pad with ."""
        assert parse_sample_fields("1/1:57", ["GT", "LN", "QV"]) == ("1/1", "57", ".")

    def test_single_token(self):
        assert parse_sample_fields("./.", ["GT", "LN", "QV"]) == ("./.", ".", ".")

    def test_field_absent_from_format(self):
        """FORMAT 无 QV 键时补 .|wanted key missing from FORMAT pads ."""
        assert parse_sample_fields("1/1:57", ["GT", "LN"]) == ("1/1", "57", ".")

    def test_default_want_gt_ln_qv(self):
        assert FLANK_FORMAT_FIELDS == ("GT", "LN", "QV")


class TestFormatFlankTsv:
    """flank TSV 格式化|flank TSV formatting"""

    @staticmethod
    def make_rows():
        """两条记录:INS 双样本 + 单样本 DEL|two records for assertions"""
        ins = ("pan_sv.INS.00001", "c1", 200, 200, "INS", "10", 150, 250,
               [("1/1", "10", "60"), ("0/0", ".", ".")], "ACGTACGTAC")
        dele = ("pan_sv.DEL.00001", "c1", 300, 309, "DEL", "-10", 250, 359,
                [("./.", "10", ".")], "TTTTAAAA")
        return [ins, dele]

    def test_header(self):
        tsv = format_flank_tsv(self.make_rows(), ["S1", "S2"])
        header = tsv.splitlines()[0]
        assert header.split("\t") == [
            "sv_id", "chrom", "pos", "end", "svtype", "svlen",
            "flank_lo", "flank_hi",
            "S1_GT", "S1_LN", "S1_QV", "S2_GT", "S2_LN", "S2_QV",
            "sequence"]

    def test_row_cells(self):
        lines = format_flank_tsv(self.make_rows(), ["S1", "S2"]).splitlines()
        assert lines[1].split("\t") == [
            "pan_sv.INS.00001", "c1", "200", "200", "INS", "10", "150", "250",
            "1/1", "10", "60", "0/0", ".", ".", "ACGTACGTAC"]

    def test_row_sample_count_flexible(self):
        """单样本记录列数与样本数一致|column count follows sample count"""
        lines = format_flank_tsv(self.make_rows(), ["S1", "S2"]).splitlines()
        assert lines[2].split("\t") == [
            "pan_sv.DEL.00001", "c1", "300", "309", "DEL", "-10", "250", "359",
            "./.", "10", ".", "TTTTAAAA"]

    def test_trailing_newline(self):
        assert format_flank_tsv(self.make_rows(), ["S1"]).endswith("\n")
