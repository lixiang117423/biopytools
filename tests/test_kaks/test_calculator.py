"""计算引擎输入准备测试|Calculator input preparation tests

回归|Regression: 原实现对不等长序列静默 [:min_len] 截断,indel配对全部错位。
新行为: align模式(默认)走翻译→比对→回译;严格模式长度不等直接报错。
"""

import pytest

from biopytools.kaks.calculator import KaKsCalculator
from biopytools.kaks.utils import KaKsLogger


@pytest.fixture
def calc(tmp_path):
    """不触发安装检查的计算器|Calculator without installation check"""
    logger = KaKsLogger(tmp_path)
    return KaKsCalculator(logger, kaks_path="KaKs_Calculator", check_installation=False)


def _read_axt(path):
    """读取AXT记录列表|Read AXT records as list of (name, seq1, seq2)"""
    records = []
    with open(path) as f:
        lines = [line.rstrip("\n") for line in f]
    i = 0
    while i < len(lines):
        if lines[i].strip():
            records.append((lines[i], lines[i + 1], lines[i + 2]))
            i += 4
        else:
            i += 1
    return records


class TestPrepareAlignMode:
    """align模式(默认)|Align mode (default)"""

    def test_unequal_pair_written_with_gaps_not_truncated(self, calc, tmp_path):
        """不等长配对产出含---的等长比对,而非截断|Unequal pair yields gapped equal-length alignment, not truncation"""
        seqs = {"s1": "ATGGCTAAATTGTAA", "s2": "ATGGCTTTGTAA"}
        pairs = [("s1", "s2", "p1")]
        axt_path, records = calc.prepare_input_file(seqs, seqs, pairs, str(tmp_path))
        recs = _read_axt(axt_path)
        assert len(recs) == 1
        _, aln1, aln2 = recs[0]
        assert len(aln1) == len(aln2)
        assert "---" in aln1 or "---" in aln2
        # 截断会把s1削成9bp;比对保留全部残基|Truncation would cut s1 to 9bp; alignment keeps all residues
        assert aln1.replace("-", "").startswith("ATGGCTAAATTG")

    def test_returns_alignment_records_with_stats(self, calc, tmp_path):
        """返回逐对比对记录含统计|Returns per-pair alignment records with stats"""
        seqs = {"s1": "ATGGCTAAATTGTAA", "s2": "ATGGCTTTGTAA"}
        pairs = [("s1", "s2", "p1")]
        _, records = calc.prepare_input_file(seqs, seqs, pairs, str(tmp_path))
        assert len(records) == 1
        record = records[0]
        assert record.pair_name == "p1"
        assert record.skip_reason == ""
        # 纯indel差异:非gap重叠区identity=1.0|Indel-only difference: identity 1.0 over non-gap overlap
        assert record.protein_identity == 1.0
        assert record.gap_codons >= 1

    def test_equal_length_pair_also_aligned(self, calc, tmp_path):
        """等长配对同样过比对(防补偿性indel错位)|Equal-length pairs also aligned (guards compensating indels)"""
        seqs = {
            "s1": "ATGGCTAAATTGTAA",
            "s2": "ATGGCTAAATTGTAA",
        }
        pairs = [("s1", "s2", "p1")]
        axt_path, records = calc.prepare_input_file(seqs, seqs, pairs, str(tmp_path))
        assert records[0].skip_reason == ""
        assert records[0].protein_identity == 1.0

    def test_skipped_pair_excluded_from_axt(self, calc, tmp_path):
        """无法翻译的配对不入AXT但保留记录|Untranslatable pair excluded from AXT but kept in records"""
        seqs = {
            "good1": "ATGGCTAAATTGTAA",
            "good2": "ATGGCTTTGTAA",
            "bad": "ATGTAAAAGGCTTAA",  # 内部终止子|internal stop
        }
        pairs = [("good1", "good2", "p_ok"), ("good1", "bad", "p_bad")]
        axt_path, records = calc.prepare_input_file(seqs, seqs, pairs, str(tmp_path))
        recs = _read_axt(axt_path)
        assert len(recs) == 1
        assert recs[0][0] == "p_ok"
        assert len(records) == 2
        bad = [r for r in records if r.pair_name == "p_bad"][0]
        assert bad.skip_reason == "internal_stop"


class TestPrepareStrictMode:
    """严格模式(align=False)|Strict mode (align=False)"""

    def test_unequal_lengths_raise_with_count(self, calc, tmp_path):
        """长度不等直接报错并给出总数|Unequal lengths raise with total count"""
        seqs = {"s1": "ATGGCTAAATTGTAA", "s2": "ATGGCTTTGTAA"}
        pairs = [("s1", "s2", "p1")]
        with pytest.raises(ValueError, match="长度不等|unequal"):
            calc.prepare_input_file(seqs, seqs, pairs, str(tmp_path), align=False)

    def test_multiple_unequal_collected_once(self, calc, tmp_path):
        """多个不等长配对一次性收集报错|Multiple unequal pairs collected into one error"""
        seqs = {
            "a1": "ATGGCTAAATTGTAA", "a2": "ATGGCTTTGTAA",
            "b1": "ATGGCTAAATTGCGTTAA", "b2": "ATGGCTCGTTAA",
        }
        pairs = [("a1", "a2", "p1"), ("b1", "b2", "p2")]
        with pytest.raises(ValueError) as exc_info:
            calc.prepare_input_file(seqs, seqs, pairs, str(tmp_path), align=False)
        assert "2" in str(exc_info.value)

    def test_equal_lengths_pass_through_without_gaps(self, calc, tmp_path):
        """等长配对原样写入,无gap|Equal-length pairs pass through untouched, no gaps"""
        seqs = {"s1": "ATGGCTAAATTG", "s2": "ATGGCGAAATTG"}
        pairs = [("s1", "s2", "p1")]
        axt_path, records = calc.prepare_input_file(seqs, seqs, pairs, str(tmp_path), align=False)
        recs = _read_axt(axt_path)
        assert recs[0][1] == "ATGGCTAAATTG"
        assert recs[0][2] == "ATGGCGAAATTG"
        assert "-" not in recs[0][1] + recs[0][2]
