"""ccs+refine 前处理测试|Tests for PacBio preprocessing (ccs/refine/fasta)"""

import pytest

import biopytools.rna_iso.steps_ccs as steps
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
    (tmp_path / "m1.subreads.bam").write_text("x")
    args = dict(
        reads=[str(tmp_path / "m1.subreads.bam")], data_type=None,
        output_dir=str(tmp_path / "out"),
    )
    args.update(kw)
    return RnaIsoConfig(**args)


@pytest.fixture
def captured(monkeypatch):
    """捕获run_command调用并模拟产出输出文件|Capture run_command calls and fake outputs"""
    import os

    calls = []

    def fake_run(cmd, logger, work_dir=None, capture_output=False):
        calls.append(cmd)
        # 模拟工具产出输出文件(命令中以/out/开头目录下的bam/fa目标)|Fake output creation
        for arg in cmd[1:]:
            if "/out/" in arg and arg.endswith((".bam", ".fa")):
                os.makedirs(os.path.dirname(arg), exist_ok=True)
                with open(arg, "w") as f:
                    f.write("fake")
        return True

    monkeypatch.setattr(steps, "run_command", fake_run)
    return calls


class TestRunCcs:
    def test_command_construction(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        logger = FakeLogger()
        outs = steps.run_ccs(cfg, logger)
        assert len(captured) == 1
        cmd = captured[0]
        assert cmd[0] == cfg.ccs_path
        assert cmd[1].endswith("m1.subreads.bam")
        assert cmd[2].endswith("01_ccs/m1_ccs.bam")
        assert "--min-passes" in cmd and "1" in cmd
        assert "--num-threads" in cmd and "12" in cmd

    def test_pbi_missing_warning(self, tmp_path, captured):
        cfg = _cfg(tmp_path)  # 未创建 .pbi
        logger = FakeLogger()
        steps.run_ccs(cfg, logger)
        assert any("pbi" in w for w in logger.warnings)

    def test_resume_skips_completed(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        out_dir = tmp_path / "out" / "01_ccs"
        out_dir.mkdir(parents=True)
        (out_dir / "m1_ccs.bam").write_text("done")
        steps.run_ccs(cfg, FakeLogger())
        assert len(captured) == 0

    def test_multiple_movies(self, tmp_path, captured):
        (tmp_path / "m2.subreads.bam").write_text("x")
        cfg = _cfg(tmp_path, reads=[str(tmp_path / "m1.subreads.bam"),
                                    str(tmp_path / "m2.subreads.bam")])
        steps.run_ccs(cfg, FakeLogger())
        assert len(captured) == 2


class TestRunRefine:
    def test_command_construction(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        ccs_bams = [str(tmp_path / "out" / "01_ccs" / "m1_ccs.bam")]
        flncs = steps.run_refine(cfg, FakeLogger(), ccs_bams)
        cmd = captured[0]
        assert cmd[0] == cfg.isoseq3_path
        assert cmd[1] == "refine"
        assert cmd[2].endswith("m1_ccs.bam")
        assert cmd[4].endswith("02_refine/m1_flnc.bam")
        assert "-j" in cmd and "12" in cmd
        # 引物为tmp写出的fasta|primers fasta written to tmp
        assert cmd[3] == cfg.tmp_dir + "/primers.fasta"

    def test_resume_skips(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        out_dir = tmp_path / "out" / "02_refine"
        out_dir.mkdir(parents=True)
        (out_dir / "m1_flnc.bam").write_text("done")
        steps.run_refine(cfg, FakeLogger(),
                         [str(tmp_path / "out" / "01_ccs" / "m1_ccs.bam")])
        assert len(captured) == 0


class TestFlncToFasta:
    def test_samtools_fasta(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        flncs = [str(tmp_path / "out" / "02_refine" / "m1_flnc.bam")]
        steps.flnc_to_fasta(cfg, FakeLogger(), flncs)
        cmd = captured[0]
        assert cmd[0] == cfg.samtools_path
        assert cmd[1] == "fasta"
        assert cmd[2] == "-@"

    def test_resume_skips(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        out_dir = tmp_path / "out" / "02_refine"
        out_dir.mkdir(parents=True)
        (out_dir / "m1_flnc.fa").write_text("done")
        steps.flnc_to_fasta(cfg, FakeLogger(),
                            [str(out_dir / "m1_flnc.bam")])
        assert len(captured) == 0
