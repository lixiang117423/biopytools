"""phy2fa utils/calculator 测试|phy2fa utils & calculator tests"""

import gzip
import os
from pathlib import Path

import pytest

from biopytools.phy2fa.calculator import Phy2FaCalculator
from biopytools.phy2fa.config import Phy2FaConfig
from biopytools.phy2fa.utils import (ModuleLogger, parse_header, parse_phy,
                                     write_fasta)


def make_calculator(tmp_path, phy_text=None, **kwargs):
    """构造转换器(带假输入)|Build a calculator with a fake input"""
    if phy_text is None:
        phy_text = make_phy_text()
    src = tmp_path / "in.phy"
    src.write_text(phy_text)
    kwargs.setdefault("output_dir", str(tmp_path / "out"))
    config = Phy2FaConfig(input=str(src), **kwargs)
    logger = ModuleLogger(log_level="ERROR").get_logger()
    return Phy2FaCalculator(config, logger), config


def make_phy_text(taxa=("A", "B", "C"), n_char=6):
    """sequential phylip|sequential layout"""
    rows = [f"{len(taxa)} {n_char}"]
    for t in taxa:
        rows.append(f"{t} {'ACGTAC'[:n_char]}")
    return "\n".join(rows) + "\n"


def make_interleaved(n_tax=3, blocks=2, width_per_row=3):
    """interleaved phylip|interleaved layout"""
    seqs = {}
    for i in range(n_tax):
        seqs[f"S{i+1}"] = "".join("ACGT"[(i + j) % 4] for j in range(blocks * width_per_row))
    lines = [f"{n_tax} {blocks * width_per_row}", ""]
    for b in range(blocks):
        for i in range(n_tax):
            s = seqs[f"S{i+1}"][b*width_per_row:(b+1)*width_per_row]
            # 标准 interleaved:每行都带名称|standard interleaved keeps names
            lines.append(f"S{i+1} {s}")
        lines.append("")
    return "\n".join(lines)


class TestHeader:
    """头行解析|header parsing"""

    def test_valid(self):
        assert parse_header("365 147433") == (365, 147433)

    def test_extra_fields(self):
        with pytest.raises(ValueError, match="头行格式错误"):
            parse_header("3 60 extra")

    def test_non_integer(self):
        with pytest.raises(ValueError, match="非整数"):
            parse_header("three sixty")

    def test_non_positive(self):
        with pytest.raises(ValueError, match="正"):
            parse_header("0 100")


class TestSequential:
    """sequential 形态|sequential layout"""

    def test_parse_and_lengths(self, tmp_path):
        calc, cfg = make_calculator(tmp_path)
        assert calc.run() == 0
        fasta = cfg.output_fasta.read_text()
        assert ">A\nACGTAC\n>B\nACGTAC\n>C\nACGTAC\n" in fasta

    def test_resume_skips_existing(self, tmp_path):
        calc, cfg = make_calculator(tmp_path)
        calc.run()
        size_before = cfg.output_fasta.stat().st_size
        # 改坏输入再跑——不应重写(断点续传)|stale input must not rewrite output
        (tmp_path / "in.phy").write_text("garbage")
        assert calc.run() == 0
        assert cfg.output_fasta.stat().st_size == size_before


class TestInterleaved:
    """interleaved 形态|interleaved layout"""

    def test_full_concat(self, tmp_path):
        text = make_interleaved(n_tax=3, blocks=2, width_per_row=3)
        calc, cfg = make_calculator(tmp_path, phy_text=text)
        assert calc.run() == 0
        content = cfg.output_fasta.read_text()
        # 每样本完整序列 = 两块拼接(每块各 3 字符)
        assert ">S1\n" in content

    def test_mismatched_length_reported(self, tmp_path):
        """声明 nchar 与实际不符 → 转换失败|declared length mismatch fails"""
        text = make_interleaved(n_tax=3, blocks=2, width_per_row=3)
        # 把头行 n_char 从 6 改成 99
        lines = text.splitlines()
        lines[0] = lines[0].replace("6", "99")
        d = tmp_path / "in_bad.phy"
        d.write_text("\n".join(lines) + "\n")
        config = Phy2FaConfig(input=str(d), output_dir=str(tmp_path / "out"))
        logger = ModuleLogger(log_level="ERROR").get_logger()
        calculator = Phy2FaCalculator(config, logger)
        assert calculator.run() == 1


class TestWriteFasta:
    """FASTA 写出|FASTA writing"""

    def test_line_width_zero_single_line(self, tmp_path):
        out = tmp_path / "x.fa"
        write_fasta({"S": "ACGTACGT"}, str(out), line_width=0)
        assert out.read_text() == ">S\nACGTACGT\n"

    def test_line_width_2(self, tmp_path):
        out = tmp_path / "x.fa"
        write_fasta({"S": "ACGTACGT"}, str(out), line_width=2)
        assert out.read_text() == ">S\nAC\nGT\nAC\nGT\n"


class TestGzip:
    """gz 输入|gzipped input"""

    def test_gz_input(self, tmp_path):
        import io
        calc, cfg = make_calculator(tmp_path)
        gz_src = tmp_path / "in2.phy.gz"
        with gzip.open(gz_src, "wt") as fh:
            fh.write(make_phy_text(taxa=("X", "Y"), n_char=4))
        config2 = Phy2FaConfig(input=str(gz_src),
                               output_dir=str(tmp_path / "out2"),
                               line_width=0)
        logger = ModuleLogger(log_level="ERROR").get_logger()
        calc2 = Phy2FaCalculator(config2, logger)
        assert calc2.run() == 0
        fa = config2.output_fasta.read_text()
        assert ">X\nACGT\n>Y\nACGT\n" in fa
