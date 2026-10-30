"""reads2tree utils 测试|reads2tree utils tests"""

import gzip

import pytest

from biopytools.reads2tree.utils import (
    concat_reads, merge_paired_reads, scan_fastq_dir, split_fastq_filename,
)


class TestSplitFastq:
    """文件名解析|filename parsing"""

    def test_plain_fastq(self):
        assert split_fastq_filename("s1.fq") == ("s1", False)
        assert split_fastq_filename("s1.fastq") == ("s1", False)

    def test_gz_double_suffix(self):
        assert split_fastq_filename("s1_R1.fq.gz") == ("s1_R1", True)

    def test_case_insensitive(self):
        assert split_fastq_filename("S1.FQ") == ("S1", False)

    def test_non_fastq_rejected(self):
        assert split_fastq_filename("s1.fa") is None
        assert split_fastq_filename("s1.txt.gz") is None


class TestScanDir:
    """目录扫描与配对|dir scanning and pairing"""

    def test_paired_and_single(self, tmp_path):
        for name in ["A_R1.fq.gz", "A_R2.fq.gz", "B_1.fastq", "B_2.fastq",
                     "C.fq", "notes.txt"]:
            (tmp_path / name).write_text("x")
        fastq, ignored, paired, singles = scan_fastq_dir(str(tmp_path))
        assert len(fastq) == 5
        assert ignored == ["notes.txt"]
        assert [(s, r1, r2) for s, r1, r2 in paired] == [
            ("A", [str(tmp_path / "A_R1.fq.gz")], [str(tmp_path / "A_R2.fq.gz")]),
            ("B", [str(tmp_path / "B_1.fastq")], [str(tmp_path / "B_2.fastq")]),
        ]
        assert singles == [str(tmp_path / "C.fq")]

    def test_clean_suffix_paired(self, tmp_path):
        """_1_clean.fq.gz/_2_clean.fq.gz 形态|the _1_clean.fq.gz naming"""
        for name in ["SRR445529_1_clean.fq.gz", "SRR445529_2_clean.fq.gz"]:
            (tmp_path / name).write_text("x")
        _, _, paired, singles = scan_fastq_dir(str(tmp_path))
        assert len(paired) == 1
        sample, r1s, r2s = paired[0]
        assert sample == "SRR445529"
        assert len(r1s) == 1 and len(r2s) == 1

    def test_recursive_subdirs(self, tmp_path):
        raw = tmp_path / "raw"
        raw.mkdir()
        (raw / "A_R1.fq.gz").write_text("x")
        (raw / "A_R2.fq.gz").write_text("x")
        (tmp_path / "run.sh").write_text("echo hi")
        fastq, ignored, paired, singles = scan_fastq_dir(str(tmp_path))
        assert len(fastq) == 2
        assert ignored == ["run.sh"]
        assert len(paired) == 1

    def test_same_sample_across_subdirs_conflicts(self, tmp_path):
        raw = tmp_path / "raw"
        clean = tmp_path / "clean"
        raw.mkdir()
        clean.mkdir()
        (raw / "SRR1_1.fq.gz").write_text("x")
        (raw / "SRR1_2.fq.gz").write_text("x")
        (clean / "SRR1_1_clean.fq.gz").write_text("x")
        (clean / "SRR1_2_clean.fq.gz").write_text("x")
        with pytest.raises(ValueError, match="多个子目录"):
            scan_fastq_dir(str(tmp_path))

    def test_multi_lane_grouped(self, tmp_path):
        for name in ["A_R1_001.fq.gz", "A_R2_001.fq.gz",
                     "A_R1_002.fq.gz", "A_R2_002.fq.gz"]:
            (tmp_path / name).write_text("x")
        _, _, paired, singles = scan_fastq_dir(str(tmp_path))
        assert len(paired) == 1
        sample, r1s, r2s = paired[0]
        assert sample == "A"
        assert len(r1s) == 2 and len(r2s) == 2

    def test_one_sided_pair_falls_to_single(self, tmp_path):
        (tmp_path / "A_R1.fq").write_text("x")
        _, _, paired, singles = scan_fastq_dir(str(tmp_path))
        assert paired == []
        assert len(singles) == 1


class TestConcatReads:
    """reads 拼接|read concatenation"""

    def test_plain_and_gz_mixed(self, tmp_path):
        src_dir = tmp_path / "src"
        src_dir.mkdir()
        (src_dir / "r1.fq").write_text("@r1\nACGT\n+\n!!!!\n")
        with gzip.open(src_dir / "r2.fq.gz", "wt") as fh:
            fh.write("@r2\nTGCA\n+\n!!!!\n")
        dest = tmp_path / "out" / "A.fq"
        assert concat_reads([str(src_dir / "r1.fq"), str(src_dir / "r2.fq.gz")],
                            str(dest)) is True
        content = dest.read_text()
        assert "@r1" in content and "@r2" in content

    def test_checkpoint_skips_existing(self, tmp_path):
        src = tmp_path / "r.fq"
        src.write_text("data")
        dest = tmp_path / "A.fq"
        dest.write_text("old")
        assert concat_reads([str(src)], str(dest)) is True
        assert dest.read_text() == "old"  # 不覆盖|not overwritten

    def test_missing_source_fails(self, tmp_path):
        assert concat_reads([str(tmp_path / "nope.fq")],
                            str(tmp_path / "A.fq")) is False


class TestMergePaired:
    """BBMerge 合并|BBMerge merging"""

    def test_command_and_concat(self, tmp_path, monkeypatch):
        import subprocess

        r1 = tmp_path / "A_R1.fq"
        r2 = tmp_path / "A_R2.fq"
        r1.write_text("@r1\nACGT\n+\n!!!!\n")
        r2.write_text("@r2\nTGCA\n+\n!!!!\n")

        calls = {}

        def fake_run(cmd, **kwargs):
            calls["cmd"] = cmd
            # 模拟 BBMerge 产出|simulate BBMerge output
            part = cmd[cmd.index("-out") + 1]
            with open(part, "w") as fh:
                fh.write("@merged\nACGT\n+\n!!!!\n")
            return type("R", (), {"returncode": 0, "stderr": "", "stdout": ""})()

        # merge_paired_reads 内 `import subprocess` 拿同一全局模块|same module object
        monkeypatch.setattr(subprocess, "run", fake_run)
        dest = merge_paired_reads("A", [str(r1)], [str(r2)],
                                  str(tmp_path / "out"), "bbmerge.sh")
        assert dest is not None
        assert calls["cmd"][0] == "bbmerge.sh"
        assert calls["cmd"][1:3] == ["-in1", str(r1)]
        assert "@merged" in open(dest).read()
