"""cluster2 转录本聚类测试|Tests for isoseq3 cluster2 engine"""

import pytest

import biopytools.rna_iso.steps_isoseq3 as steps
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
    (tmp_path / "m.subreads.bam").write_text("x")
    args = dict(
        reads=[str(tmp_path / "m.subreads.bam")], data_type=None,
        output_dir=str(tmp_path / "out"), prefix="S1",
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
            if "/out/" in arg and arg.endswith((".bam", ".fa")):
                os.makedirs(os.path.dirname(arg), exist_ok=True)
                with open(arg, "w") as f:
                    f.write("fake")
        return True

    monkeypatch.setattr(steps, "run_command", fake_run)
    return calls


class TestCluster2:
    def test_command_construction(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        flnc = str(tmp_path / "out" / "02_refine" / "m_flnc.bam")
        ok = steps.run_cluster2(cfg, FakeLogger(), [flnc])
        assert ok is True
        cmd = captured[0]
        assert cmd[0] == cfg.isoseq3_path
        assert cmd[1] == "cluster2"
        assert cmd[2] == flnc  # 单文件直接传|single file passed directly
        assert cmd[3].endswith("04_isoseq3/S1_transcripts.bam")
        assert "-j" in cmd and "12" in cmd

    def test_multi_files_use_fofn(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        flncs = [str(tmp_path / "out" / "02_refine" / "m1_flnc.bam"),
                 str(tmp_path / "out" / "02_refine" / "m2_flnc.bam")]
        steps.run_cluster2(cfg, FakeLogger(), flncs)
        cmd = captured[0]
        # 多文件走fofn(cluster2原生支持)|multiple files via fofn
        assert cmd[2].endswith("flnc.fofn")
        fofn_content = open(cmd[2]).read()
        assert "m1_flnc.bam" in fofn_content and "m2_flnc.bam" in fofn_content

    def test_resume_skips(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        out_dir = tmp_path / "out" / "04_isoseq3"
        out_dir.mkdir(parents=True)
        (out_dir / "S1_transcripts.bam").write_text("done")
        steps.run_cluster2(cfg, FakeLogger(),
                           [str(tmp_path / "out" / "02_refine" / "m_flnc.bam")])
        assert len(captured) == 0

    def test_summary_logged(self, tmp_path, captured):
        cfg = _cfg(tmp_path)
        flnc = str(tmp_path / "out" / "02_refine" / "m_flnc.bam")
        # cluster2 后 fasta 保底转换也会被捕获|fasta fallback also captured
        logger = FakeLogger()
        steps.run_cluster2(cfg, logger, [flnc])
        assert any("transcripts" in i or "转录本" in i for i in logger.infos)


class TestFastaStats:
    def test_counts_sequences(self, tmp_path):
        fa = tmp_path / "t.fa"
        fa.write_text(">a\nACGT\n>b\nACGTACGT\n")
        assert steps.count_fasta_records(str(fa)) == 2

    def test_empty(self, tmp_path):
        assert steps.count_fasta_records(str(tmp_path / "none.fa")) == 0
