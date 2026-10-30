"""IsoQuant 执行器测试|Tests for IsoQuant executor"""

import pytest

import biopytools.rna_iso.steps_isoquant as steps
from biopytools.rna_iso.config import RnaIsoConfig


class FakeLogger:
    def __init__(self):
        self.infos = []
        self.warnings = []

    def info(self, msg):
        self.infos.append(msg)

    def warning(self, msg):
        self.warnings.append(msg)

    def error(self, msg):
        pass


def _cfg(tmp_path, **kw):
    (tmp_path / "s.fq.gz").write_text(">x\nACGT\n")
    (tmp_path / "ref.fa").write_text(">x\nACGT\n")
    args = dict(
        reads=[str(tmp_path / "s.fq.gz")], data_type="ont",
        reference=str(tmp_path / "ref.fa"), output_dir=str(tmp_path / "out"),
        prefix="S1",
    )
    args.update(kw)
    return RnaIsoConfig(**args)


@pytest.fixture
def captured(monkeypatch):
    import os

    calls = []

    def fake_run(cmd, logger, work_dir=None, capture_output=False):
        calls.append(cmd)
        for arg in cmd[1:]:
            if arg.startswith("/") and "/out/" in arg and not arg.endswith((".gz", ".fa", ".gtf")):
                pass
        return True

    def fake_run_and_touch(cmd, logger, work_dir=None, capture_output=False):
        calls.append(cmd)
        # 模拟IsoQuant产出目录树|Simulate IsoQuant output tree
        if "--prefix" in cmd:
            i = cmd.index("--prefix")
            sample_dir = os.path.join(cmd[cmd.index("-o") + 1], cmd[i + 1])
            os.makedirs(sample_dir, exist_ok=True)
            with open(os.path.join(sample_dir, f"{cmd[i + 1]}.transcript_models.gtf"), "w") as f:
                f.write("gtf")
        return True

    monkeypatch.setattr(steps, "run_command", fake_run_and_touch)
    return calls


class TestRunIsoquant:
    def test_command_construction(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        ok = steps.run_isoquant(cfg, FakeLogger(),
                                reads_inputs=[str(tmp_path / "s.fq.gz")])
        assert ok is True
        cmd = captured[0]
        assert cmd[0] == cfg.isoquant_path
        assert "--reference" in cmd
        i = cmd.index("--reference")
        assert cmd[i + 1].endswith("ref.fa")
        assert "--data_type" in cmd
        i = cmd.index("--data_type")
        assert cmd[i + 1] == "nanopore"
        assert "-t" in cmd and "12" in cmd
        i = cmd.index("--prefix")
        assert cmd[i + 1] == "S1"

    def test_genedb_passed_when_given(self, tmp_path, captured):
        (tmp_path / "anno.gtf").write_text("x")
        cfg = _cfg(tmp_path, genedb=str(tmp_path / "anno.gtf"))
        steps.run_isoquant(cfg, FakeLogger(),
                           reads_inputs=[str(tmp_path / "s.fq.gz")])
        cmd = captured[0]
        assert "--genedb" in cmd

    def test_genedb_absent_when_none(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        steps.run_isoquant(cfg, FakeLogger(),
                           reads_inputs=[str(tmp_path / "s.fq.gz")])
        assert "--genedb" not in captured[0]

    def test_pacbio_data_type(self, tmp_path, captured):
        cfg = _cfg(tmp_path, data_type=None)
        cfg.isoquant_data_type = "pacbio_ccs"
        steps.run_isoquant(cfg, FakeLogger(),
                           reads_inputs=[str(tmp_path / "out" / "02_refine" / "m_flnc.fa")])
        cmd = captured[0]
        i = cmd.index("--data_type")
        assert cmd[i + 1] == "pacbio_ccs"

    def test_multi_reads_all_passed(self, tmp_path, captured):
        (tmp_path / "s2.fq.gz").write_text(">y\nACGT\n")
        cfg = _cfg(tmp_path, reads=[str(tmp_path / "s.fq.gz"), str(tmp_path / "s2.fq.gz")])
        steps.run_isoquant(cfg, FakeLogger(),
                           reads_inputs=cfg.reads)
        cmd = captured[0]
        i = cmd.index("--fastq")
        # --fastq 后跟两个文件直到下一个选项|two files after --fastq
        assert cmd[i + 1].endswith("s.fq.gz") and cmd[i + 2].endswith("s2.fq.gz")

    def test_resume_skips(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        sample_dir = tmp_path / "out" / "03_isoquant" / "S1"
        sample_dir.mkdir(parents=True)
        (sample_dir / "S1.transcript_models.gtf").write_text("gtf")
        steps.run_isoquant(cfg, FakeLogger(),
                           reads_inputs=[str(tmp_path / "s.fq.gz")])
        assert len(captured) == 0
