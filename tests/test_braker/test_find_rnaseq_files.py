"""find_rnaseq_files_in_directory 旧命名回退测试|Legacy naming fallback tests

复现事故|Reproduced incident: 2026-09-05 joblogs/2026-09-05-21-32-genome-anno-41.out
2nd/ 目录下全是旧命名软链接(*_1.clean.fq.gz),命名迁移(95d775d)后默认
pattern 为 _1_clean.fq.gz,alternative_patterns 漏掉旧名 _1.clean.fq.gz,
导致 "Found 0 pairs" 退出。软链接本身glob完全支持,与软链接无关。
"""

import os

from biopytools.braker.utils import find_rnaseq_files_in_directory


def _make_fastq_pair(directory, sample, suffix_r1, suffix_r2, symlink=False):
    """造一对FASTQ文件(可选软链接)|Create a FASTQ pair (optionally as symlinks)"""
    r1 = directory / f"{sample}{suffix_r1}"
    r2 = directory / f"{sample}{suffix_r2}"
    target = directory / "target.fq.gz"
    target.write_text("@read\nACGT\n+\nIIII\n")
    for f in (r1, r2):
        if symlink:
            os.symlink(target, f)
        else:
            f.write_text("@read\nACGT\n+\nIIII\n")
    return r1, r2


def test_legacy_dot_pattern_fallback(tmp_path):
    """旧命名 _1.clean.fq.gz 应被回退识别|Legacy dot naming must be recognized via fallback"""
    _make_fastq_pair(tmp_path, "sampleA", "_1.clean.fq.gz", "_2.clean.fq.gz")

    pairs = find_rnaseq_files_in_directory(str(tmp_path))

    assert len(pairs) == 1
    assert pairs[0][0].endswith("sampleA_1.clean.fq.gz")
    assert pairs[0][1].endswith("sampleA_2.clean.fq.gz")


def test_symlinked_fastq_found(tmp_path):
    """软链接FASTQ必须被识别(回归守卫)|Symlinked FASTQ must be found (regression guard)"""
    _make_fastq_pair(tmp_path, "sampleS", "_1.clean.fq.gz", "_2.clean.fq.gz", symlink=True)

    pairs = find_rnaseq_files_in_directory(str(tmp_path))

    assert len(pairs) == 1
    assert os.path.islink(pairs[0][0])


def test_new_underscore_pattern_still_works(tmp_path):
    """新命名 _1_clean.fq.gz 默认直配|New underscore naming matched by default pattern"""
    _make_fastq_pair(tmp_path, "sampleB", "_1_clean.fq.gz", "_2_clean.fq.gz")

    pairs = find_rnaseq_files_in_directory(str(tmp_path))

    assert len(pairs) == 1


def test_missing_r2_not_paired(tmp_path):
    """R2缺失时不入配对|R1 without R2 must not be paired"""
    _make_fastq_pair(tmp_path, "sampleC", "_1.clean.fq.gz", "_2.clean.fq.gz")
    os.remove(tmp_path / "sampleC_2.clean.fq.gz")

    pairs = find_rnaseq_files_in_directory(str(tmp_path))

    assert pairs == []
