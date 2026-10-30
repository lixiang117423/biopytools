"""detect_suffixes 旧命名识别测试|Legacy naming detection tests

同根因|Same root cause: 命名迁移(95d775d)后 SUFFIX_CANDIDATES 漏掉旧规范
_1.clean.fq.gz/_2.clean.fq.gz,旧命名数据目录将检测失败返回 None。
"""

from biopytools.insert2locus.config import detect_suffixes


def test_legacy_dot_suffix_detected(tmp_path):
    """旧命名 _1.clean.fq.gz 应被检测|Legacy dot naming must be detected"""
    (tmp_path / "sampleL_1.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "sampleL_2.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    assert detect_suffixes(tmp_path) == ("_1.clean.fq.gz", "_2.clean.fq.gz")


def test_new_underscore_suffix_still_first(tmp_path):
    """新旧并存时新命名优先|New naming wins when both present"""
    for suf in ("_1_clean.fq.gz", "_2_clean.fq.gz", "_1.clean.fq.gz", "_2.clean.fq.gz"):
        (tmp_path / f"sampleM{suf}").write_text("@r\nACGT\n+\nII\n")

    assert detect_suffixes(tmp_path) == ("_1_clean.fq.gz", "_2_clean.fq.gz")
