"""vcf2deepbsa 测试夹具|vcf2deepbsa test fixtures"""

import gzip
import logging
from pathlib import Path

import pytest


def write_lines(path: Path, lines):
    """写文本文件(可选 gz 压缩)|Write text file (optional gz)"""
    if str(path).endswith(".gz"):
        with gzip.open(path, "wt") as fh:
            fh.write("\n".join(lines) + "\n")
    else:
        path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def logger():
    """静默 logger|Silent logger"""
    return logging.getLogger("test_vcf2deepbsa")


@pytest.fixture
def make_vcf(tmp_path):
    """构造带 AD 的测试 VCF|Build a test VCF with AD fields

    用法|Usage:
        make_vcf()                         # 标准 2 样本 VCF
        make_vcf(records=[...])            # 自定义数据行
        make_vcf(name="s.vcf.gz")          # gz 输入
    """

    def _make(records=None, samples=("pool_R", "pool_S"),
              fmt="GT:AD", name="x.vcf"):
        if records is None:
            records = [
                ("chr1", 100, "A", "T", ["0/0:10,5", "0/1:20,8"]),
                ("chr1", 200, "G", "C", ["1/1:0,30", "0/0:15,3"]),
                ("chr2", 300, "AT", "A", ["0/1:7,9", "1/1:2,22"]),
            ]
        lines = [
            "##fileformat=VCFv4.2",
            "##contig=<ID=chr1>",
            "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" +
            "\t".join(samples),
        ]
        for chrom, pos, ref, alt, gts in records:
            lines.append("\t".join(
                [chrom, str(pos), ".", ref, alt, "50", "PASS", ".", fmt] +
                list(gts)))
        return str(write_lines(tmp_path / name, lines))

    return _make
