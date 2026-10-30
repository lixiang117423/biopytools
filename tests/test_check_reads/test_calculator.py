"""check_reads calculator 测试|check_reads calculator tests"""

import gzip
import io

from biopytools.check_reads.calculator import CheckReadsCalculator
from biopytools.check_reads.config import CheckReadsConfig
from biopytools.check_reads.utils import ModuleLogger


def make_calculator(tmp_path, files=None, **kwargs):
    """构造检查器(带假 fastq 目录)|Build a calculator with a fake fastq dir"""
    d = tmp_path / "reads"
    d.mkdir(exist_ok=True)
    for name, content in (files or {}).items():
        (d / name).write_bytes(content)
    kwargs.setdefault("output_dir", str(tmp_path / "out"))
    config = CheckReadsConfig(input_dir=str(d), **kwargs)
    config.validate()
    logger = ModuleLogger(log_file=str(config.log_file),
                          log_level="ERROR").get_logger()
    return CheckReadsCalculator(config, logger), d


def gz_bytes(text: str) -> bytes:
    """文本转 gz 字节|Compress text to gz bytes"""
    buf = io.BytesIO()
    with gzip.GzipFile(fileobj=buf, mode="wb") as fh:
        fh.write(text.encode())
    return buf.getvalue()


def make_pair(d, sample, ok=True):
    """造一对 R1/R2|Create an R1/R2 pair"""
    content = gz_bytes(f"@{sample}\nACGT\n+\n!!!!\n") if ok else b"PK\x03\x04truncated"
    (d / f"{sample}_1_clean.fq.gz").write_bytes(content)
    (d / f"{sample}_2_clean.fq.gz").write_bytes(content)


class TestCheckGzip:
    """gz 完整性|gzip integrity"""

    def test_valid_gz(self, tmp_path):
        calc, _ = make_calculator(tmp_path)
        p = tmp_path / "ok.fq.gz"
        p.write_bytes(gz_bytes("data"))
        assert calc.check_gzip(str(p)) is True

    def test_corrupt_gz(self, tmp_path):
        calc, _ = make_calculator(tmp_path)
        p = tmp_path / "bad.fq.gz"
        p.write_bytes(b"PK\x03\x04truncated")
        assert calc.check_gzip(str(p)) is False

    def test_empty_file(self, tmp_path):
        calc, _ = make_calculator(tmp_path)
        p = tmp_path / "empty.fq.gz"
        p.write_bytes(b"")
        assert calc.check_file(str(p)) == "EMPTY"

    def test_plain_file(self, tmp_path):
        calc, _ = make_calculator(tmp_path)
        p = tmp_path / "plain.fq"
        p.write_text("@x\nACGT\n+\n!!!!\n")
        assert calc.check_file(str(p)) == "PLAIN"


class TestPairing:
    """配对检查|pairing check"""

    def test_complete_pairs(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A")
        make_pair(d, "B")
        assert calc.check_pairing(
            [str(p) for p in d.iterdir()]) == []

    def test_missing_mate_detected(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A")
        (d / "B_1_clean.fq.gz").write_bytes(gz_bytes("x"))
        issues = calc.check_pairing([str(p) for p in d.iterdir()])
        assert ("B", "R2") in issues

    def test_single_end_ignored(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A")
        (d / "C.fq.gz").write_bytes(gz_bytes("x"))
        issues = calc.check_pairing([str(p) for p in d.iterdir()])
        assert issues == []


class TestRun:
    """主流程|main run"""

    def test_all_clean_returns_zero(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A")
        make_pair(d, "B")
        assert calc.run() == 0
        report = calc.config.report_file.read_text()
        assert "全部通过" in report
        assert "OK" in report

    def test_corrupt_returns_one(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A", ok=False)
        make_pair(d, "B")
        assert calc.run() == 1
        report = calc.config.report_file.read_text()
        assert "CORRUPT" in report and "发现问题" in report

    def test_empty_returns_one(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A")
        (d / "B_1_clean.fq.gz").write_bytes(b"")
        (d / "B_2_clean.fq.gz").write_bytes(gz_bytes("x"))
        assert calc.run() == 1

    def test_missing_mate_returns_one(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A")
        (d / "B_1_clean.fq.gz").write_bytes(gz_bytes("x"))
        assert calc.run() == 1
        report = calc.config.report_file.read_text()
        assert "缺R2" in report

    def test_no_fastq_returns_one(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        (d / "notes.txt").write_text("x")
        assert calc.run() == 1

    def test_report_counts(self, tmp_path):
        calc, d = make_calculator(tmp_path)
        make_pair(d, "A")
        (d / "C.fq.gz").write_bytes(gz_bytes("single"))
        assert calc.run() == 0
        report = calc.config.report_file.read_text()
        assert "单端|single-end: 1" in report

    def test_multi_dir_input(self, tmp_path):
        """逗号分隔多目录|comma-separated multiple dirs"""
        d1 = tmp_path / "d1"
        d2 = tmp_path / "d2"
        d1.mkdir()
        d2.mkdir()
        make_pair(d1, "A")
        make_pair(d2, "B")
        config = CheckReadsConfig(
            input_dir=f"{d1},{d2}", output_dir=str(tmp_path / "out"))
        config.validate()
        logger = ModuleLogger(log_file=str(config.log_file),
                              log_level="ERROR").get_logger()
        calc = CheckReadsCalculator(config, logger)
        assert calc.run() == 0
