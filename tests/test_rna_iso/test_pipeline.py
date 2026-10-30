"""流程编排测试|Tests for pipeline orchestration"""

import pytest

import biopytools.rna_iso.pipeline as pipeline_mod
from biopytools.rna_iso.config import RnaIsoConfig


class FakeLogger:
    def info(self, msg):
        pass

    def warning(self, msg):
        pass

    def error(self, msg):
        pass


def _ont_cfg(tmp_path, **kw):
    (tmp_path / "s.fq.gz").write_text(">x\nACGT\n")
    (tmp_path / "ref.fa").write_text(">x\nACGT\n")
    args = dict(
        reads=[str(tmp_path / "s.fq.gz")], data_type="ont",
        reference=str(tmp_path / "ref.fa"), output_dir=str(tmp_path / "out"),
    )
    args.update(kw)
    return RnaIsoConfig(**args)


def _subreads_cfg(tmp_path, **kw):
    (tmp_path / "m.subreads.bam").write_text("x")
    (tmp_path / "ref.fa").write_text(">x\nACGT\n")
    args = dict(
        reads=[str(tmp_path / "m.subreads.bam")], data_type=None,
        reference=str(tmp_path / "ref.fa"), output_dir=str(tmp_path / "out"),
    )
    args.update(kw)
    return RnaIsoConfig(**args)


@pytest.fixture
def steps_calls(monkeypatch):
    """mock三段执行器记录调用顺序|Mock step executors, record call order"""
    order = []

    monkeypatch.setattr(pipeline_mod, "run_ccs",
                        lambda cfg, logger: order.append("ccs") or ["/fake/m.ccs.bam"])
    monkeypatch.setattr(pipeline_mod, "run_refine",
                        lambda cfg, logger, ccs: order.append("refine") or ["/fake/m.flnc.bam"])
    monkeypatch.setattr(pipeline_mod, "flnc_to_fasta",
                        lambda cfg, logger, flnc: order.append("fasta") or ["/fake/m.flnc.fa"])
    monkeypatch.setattr(pipeline_mod, "run_isoquant",
                        lambda cfg, logger, reads_inputs: order.append("isoquant") or True)
    monkeypatch.setattr(pipeline_mod, "run_cluster2",
                        lambda cfg, logger, flnc: order.append("cluster2") or True)
    monkeypatch.setattr(pipeline_mod, "generate_software_versions_yml",
                        lambda cfg, out, logger: order.append("versions") or None)
    return order


class TestDispatch:
    def test_ont_only_isoquant(self, tmp_path, steps_calls):
        cfg = _ont_cfg(tmp_path)
        p = pipeline_mod.RnaIsoPipeline(cfg)
        with_fake_logger = FakeLogger()
        assert p.run(with_fake_logger) is True
        assert steps_calls == ["versions", "isoquant"]

    def test_subreads_both_full_chain(self, tmp_path, steps_calls):
        cfg = _subreads_cfg(tmp_path, engine="both")
        assert pipeline_mod.RnaIsoPipeline(cfg).run(FakeLogger()) is True
        assert steps_calls == ["versions", "ccs", "refine", "fasta",
                               "isoquant", "cluster2"]

    def test_subreads_isoseq3_only(self, tmp_path, steps_calls):
        cfg = _subreads_cfg(tmp_path, engine="isoseq3", reference=None)
        assert pipeline_mod.RnaIsoPipeline(cfg).run(FakeLogger()) is True
        assert steps_calls == ["versions", "ccs", "refine", "cluster2"]

    def test_ccs_bam_skips_ccs(self, tmp_path, steps_calls):
        (tmp_path / "m.ccs.bam").write_text("x")
        cfg = _ont_cfg(tmp_path, reads=[str(tmp_path / "m.ccs.bam")],
                       data_type=None, engine="both")
        assert pipeline_mod.RnaIsoPipeline(cfg).run(FakeLogger()) is True
        assert steps_calls == ["versions", "refine", "fasta", "isoquant", "cluster2"]

    def test_reads_fed_directly_to_isoquant(self, tmp_path, steps_calls):
        """ONT输入原样传给isoquant(不经refine)|ONT reads passed through"""
        seen = {}

        def fake_isoquant(cfg, logger, reads_inputs):
            seen["reads"] = list(reads_inputs)
            return True

        monkeypatch = pytest.MonkeyPatch()
        monkeypatch.setattr(pipeline_mod, "run_isoquant", fake_isoquant)
        monkeypatch.setattr(pipeline_mod, "generate_software_versions_yml",
                            lambda cfg, out, logger: None)
        cfg = _ont_cfg(tmp_path)
        pipeline_mod.RnaIsoPipeline(cfg).run(FakeLogger())
        monkeypatch.undo()
        assert seen["reads"] == cfg.reads

    def test_step_failure_aborts(self, tmp_path, monkeypatch):
        monkeypatch.setattr(pipeline_mod, "generate_software_versions_yml",
                            lambda cfg, out, logger: None)
        monkeypatch.setattr(pipeline_mod, "run_isoquant",
                            lambda cfg, logger, reads: False)
        cfg = _ont_cfg(tmp_path)
        assert pipeline_mod.RnaIsoPipeline(cfg).run(FakeLogger()) is False


class TestDirectories:
    def test_conditional_dirs_created(self, tmp_path, steps_calls):
        cfg = _ont_cfg(tmp_path)
        pipeline_mod.RnaIsoPipeline(cfg).run(FakeLogger())
        assert (tmp_path / "out" / "03_isoquant").is_dir()
        assert (tmp_path / "out" / "99_logs").is_dir()
        assert not (tmp_path / "out" / "01_ccs").exists()  # ONT不建ccs目录|no ccs dir for ONT

    def test_pacbio_dirs_created(self, tmp_path, steps_calls):
        cfg = _subreads_cfg(tmp_path, engine="isoseq3", reference=None)
        pipeline_mod.RnaIsoPipeline(cfg).run(FakeLogger())
        assert (tmp_path / "out" / "01_ccs").is_dir()
        assert (tmp_path / "out" / "02_refine").is_dir()
        assert (tmp_path / "out" / "04_isoseq3").is_dir()
