"""discover_samples 旧命名后缀测试|Legacy naming suffix tests

同根因|Same root cause: 命名迁移(95d775d)后 DEFAULT_R1/R2_SUFFIXES 漏掉旧规范
_1.clean.fq.gz/_2.clean.fq.gz,历史数据(如 2nd/ 软链接目录)将匹配 0 样本。
"""

from biopytools.rnaseq2vcf.utils import discover_samples


def test_legacy_dot_suffix_discovered(tmp_path):
    """旧命名 _1.clean.fq.gz 应被发现|Legacy dot naming must be discovered"""
    (tmp_path / "sampleX_1.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "sampleX_2.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    samples = discover_samples(str(tmp_path))

    assert len(samples) == 1
    assert samples[0][0] == "sampleX"


def test_new_underscore_suffix_still_works(tmp_path):
    """新命名 _1_clean.fq.gz 默认直配|New underscore naming matched by default"""
    (tmp_path / "sampleY_1_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "sampleY_2_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    samples = discover_samples(str(tmp_path))

    assert len(samples) == 1
    assert samples[0][0] == "sampleY"


def test_symlinked_pair_discovered(tmp_path):
    """软链接配对必须被发现(回归守卫)|Symlinked pair must be found (regression guard)"""
    import os
    target = tmp_path / "target.fq.gz"
    target.write_text("@r\nACGT\n+\nII\n")
    os.symlink(target, tmp_path / "sampleZ_1.clean.fq.gz")
    os.symlink(target, tmp_path / "sampleZ_2.clean.fq.gz")

    samples = discover_samples(str(tmp_path))

    assert len(samples) == 1
    assert samples[0][0] == "sampleZ"
