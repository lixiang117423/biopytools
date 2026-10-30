"""单倍型去冗余(purge haplotypes)测试|Haplotype purging tests

v1.1.0 新能力|hifi_hic 默认对 {prefix}_hap{i}.fa 也执行 purge_dups:
- n_hap>=2 且未 --no-purge-haplotypes 时,hap1..hapN 逐个去冗余
- 每个输入独立输出目录 04_purge_dups/hap{i}/(cutoffs/dups.bed 是固定名,分目录不互踩)
- hap 固定用全量HiFi reads(NGS筛选子集按primary高覆盖contig筛出,会带偏hap覆盖度阈值)
- 缺失hap降级跳过;purge失败仅WARNING不中断(primary失败仍exit 1)
- 断点续传:每输入独立判存在跳过;--no-resume强制重跑
"""

import os
import sys
from pathlib import Path
from unittest.mock import MagicMock

import pytest

from biopytools.hifi_hic.config import AssemblyConfig


def _write_fasta(path, name="c1", seq="ACGTACGTAC"):
    """写一个最小FASTA|Write a minimal FASTA"""
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    Path(path).write_text(f">{name}\n{seq}\n")


@pytest.fixture
def assembly_config(tmp_path):
    """最小AssemblyConfig(hap1/hap2齐备)|Minimal AssemblyConfig with hap1/hap2"""
    hifi = tmp_path / "hifi.fq"
    hifi.write_text(">r\nACGT\n")
    cfg = AssemblyConfig(
        hifi_data=str(hifi),
        prefix="sample1",
        base_dir=str(tmp_path / "out"),
        threads=4,
    )
    _write_fasta(os.path.join(cfg.fasta_dir, "sample1_primary.fa"))
    _write_fasta(os.path.join(cfg.fasta_dir, "sample1_hap1.fa"))
    _write_fasta(os.path.join(cfg.fasta_dir, "sample1_hap2.fa"))
    return cfg


# ---------------------------------------------------------------- config


class TestPurgeHaplotypesConfig:
    def test_default_enabled(self, assembly_config):
        """默认对单倍型去冗余|Haplotype purging enabled by default"""
        assert assembly_config.purge_haplotypes is True
        assert assembly_config.has_purge_haplotypes is True

    def test_disable_by_flag(self, tmp_path):
        """--no-purge-haplotypes 关闭|Disabled by flag"""
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        cfg = AssemblyConfig(
            hifi_data=str(hifi), prefix="s", base_dir=str(tmp_path / "o"), purge_haplotypes=False
        )
        assert cfg.has_purge_haplotypes is False

    def test_disabled_when_purge_dups_off(self, tmp_path):
        """purge_dups整体关闭时单倍型也关闭|Off when purge_dups disabled"""
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        cfg = AssemblyConfig(
            hifi_data=str(hifi), prefix="s", base_dir=str(tmp_path / "o"),
            enable_purge_dups=False,
        )
        assert cfg.has_purge_haplotypes is False


# ------------------------------------------------------------ 目标构建


class TestBuildPurgeTargets:
    def _build(self, cfg):
        from biopytools.hifi_hic.main import build_purge_targets
        return build_purge_targets(cfg)

    def test_primary_plus_haps(self, assembly_config):
        """默认含primary+hap1+hap2,hap带独立子目录|primary+hap1+hap2 with subdirs"""
        targets = self._build(assembly_config)
        labels = [t[0] for t in targets]
        assert labels == ["primary", "hap1", "hap2"]
        # primary写顶层(兼容历史输出),hap写hap{i}/子目录
        # |primary at top level (legacy-compatible), haps in hap{i}/ subdirs
        assert targets[0][2] is None
        assert targets[1][2] == "hap1"
        assert targets[2][2] == "hap2"
        assert targets[1][1].endswith("sample1_hap1.fa")

    def test_disabled_targets_primary_only(self, assembly_config):
        """关闭后仅primary|Primary only when disabled"""
        assembly_config.purge_haplotypes = False
        assembly_config.has_purge_haplotypes = False
        targets = self._build(assembly_config)
        assert [t[0] for t in targets] == ["primary"]

    def test_missing_hap_skipped(self, assembly_config):
        """缺失hap2仅告警降级|Missing hap2 degrades with warning"""
        os.remove(os.path.join(assembly_config.fasta_dir, "sample1_hap2.fa"))
        from biopytools.hifi_hic.main import build_purge_targets
        warnings = []
        logger = MagicMock()
        logger.warning.side_effect = lambda msg: warnings.append(msg)
        targets = build_purge_targets(assembly_config, logger=logger)
        assert [t[0] for t in targets] == ["primary", "hap1"]
        assert any("hap2" in w for w in warnings)

    def test_no_hap_targets_when_n_hap_1(self, tmp_path):
        """单倍性组装不对hap去冗余|No hap purge for haploid assembly"""
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        cfg = AssemblyConfig(
            hifi_data=str(hifi), prefix="s", base_dir=str(tmp_path / "o"), n_hap=1
        )
        _write_fasta(os.path.join(cfg.fasta_dir, "s_hap1.fa"))
        targets = self._build(cfg)
        assert [t[0] for t in targets] == ["primary"]

    def test_polished_preferred_when_ngs(self, tmp_path):
        """有NGS时primary目标用polished产物|Polished genome preferred with NGS"""
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        ngs_dir = tmp_path / "ngs"
        ngs_dir.mkdir()
        cfg = AssemblyConfig(
            hifi_data=str(hifi), prefix="s", base_dir=str(tmp_path / "o"),
            ngs_data=str(ngs_dir), purge_haplotypes=False,
        )
        _write_fasta(os.path.join(cfg.fasta_dir, "s_primary.fa"))
        _write_fasta(os.path.join(cfg.ngs_polish_dir, "s_polished.fa"))
        targets = self._build(cfg)
        assert targets[0][1].endswith("s_polished.fa")


# ------------------------------------------------------------ 断点检测


class TestGetCompletedStepsPurge:
    def test_requires_hap_outputs(self, assembly_config):
        """hap未purge完成不算完成|Not done until hap outputs exist"""
        seqs = Path(assembly_config.purge_dups_dir) / "seqs"
        _write_fasta(seqs / "sample1_primary_purged.purged.fa")
        steps = assembly_config.get_completed_steps()
        assert steps["purge_dups_seqs"] is False

    def test_all_present_done(self, assembly_config):
        """primary+hap1+hap2全齐才算完成|Done when primary+haps all present"""
        seqs = Path(assembly_config.purge_dups_dir) / "seqs"
        _write_fasta(seqs / "sample1_primary_purged.purged.fa")
        for i in (1, 2):
            hap_seqs = Path(assembly_config.purge_dups_dir) / f"hap{i}" / "seqs"
            _write_fasta(hap_seqs / f"sample1_hap{i}_purged.purged.fa")
        steps = assembly_config.get_completed_steps()
        assert steps["purge_dups_seqs"] is True

    def test_polished_output_recognized(self, tmp_path):
        """NGS时识别polished命名|Recognize polished output name"""
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        ngs_dir = tmp_path / "ngs"
        ngs_dir.mkdir()
        cfg = AssemblyConfig(
            hifi_data=str(hifi), prefix="s", base_dir=str(tmp_path / "o"),
            ngs_data=str(ngs_dir), purge_haplotypes=False,
        )
        seqs = Path(cfg.purge_dups_dir) / "seqs"
        _write_fasta(seqs / "s_polished_purged.purged.fa")
        steps = cfg.get_completed_steps()
        assert steps["purge_dups_seqs"] is True


# ------------------------------------------------------------ wrapper


class _FakeRunner:
    """假PurgeDupsRunner:记录kwargs并产出purged.fa|Fake runner: record kwargs, emit purged.fa"""

    calls = []

    def __init__(self, **kwargs):
        self.kwargs = kwargs
        _FakeRunner.calls.append(kwargs)

    def run_full_pipeline(self):
        stem = Path(self.kwargs["input"]).stem
        out_dir = Path(self.kwargs["output_dir"]) / "seqs"
        _write_fasta(out_dir / f"{stem}_purged.purged.fa")
        return True


class _FakeConfig:
    def __init__(self, **kwargs):
        self.kwargs = kwargs


def _make_wrapper(cfg):
    from biopytools.hifi_hic.purge_dups_wrapper import PurgeDupsWrapper
    wrapper = PurgeDupsWrapper(cfg, MagicMock())
    wrapper.PurgeDupsRunner = _FakeRunner
    wrapper.PurgeDupsConfig = _FakeConfig
    _FakeRunner.calls = []
    return wrapper


class TestPurgeDupsWrapperHap:
    def test_hap_output_subdir(self, assembly_config):
        """hap写04_purge_dups/hap{i}/|Haps write into hap{i}/ subdir"""
        wrapper = _make_wrapper(assembly_config)
        hap1 = os.path.join(assembly_config.fasta_dir, "sample1_hap1.fa")
        result = wrapper.run_purge_dups(hap1, output_subdir="hap1", force_hifi_reads=True)
        expected_dir = os.path.join(assembly_config.purge_dups_dir, "hap1")
        assert _FakeRunner.calls[0]["output_dir"] == expected_dir
        assert result == os.path.join(expected_dir, "seqs", "sample1_hap1_purged.purged.fa")

    def test_hap_uses_full_hifi_reads(self, assembly_config):
        """hap固定全量HiFi reads,不用NGS筛选子集|Haps force full HiFi reads"""
        filtered = os.path.join(assembly_config.ngs_polish_dir, "03_filtered_reads",
                                "sample1_high_quality_reads.fq.gz")
        _write_fasta(filtered)  # 内容无关,存在即可|existence is what matters
        assembly_config.has_ngs = True
        wrapper = _make_wrapper(assembly_config)
        hap1 = os.path.join(assembly_config.fasta_dir, "sample1_hap1.fa")
        wrapper.run_purge_dups(hap1, output_subdir="hap1", force_hifi_reads=True)
        assert _FakeRunner.calls[0]["reads"] == assembly_config.hifi_data

    def test_primary_still_prefers_ngs_reads(self, assembly_config):
        """primary保持原reads选择(NGS筛选优先)|Primary keeps NGS-filtered preference"""
        filtered = os.path.join(assembly_config.ngs_polish_dir, "03_filtered_reads",
                                "sample1_high_quality_reads.fq.gz")
        _write_fasta(filtered)
        assembly_config.has_ngs = True
        wrapper = _make_wrapper(assembly_config)
        primary = os.path.join(assembly_config.fasta_dir, "sample1_primary.fa")
        wrapper.run_purge_dups(primary)
        assert _FakeRunner.calls[0]["reads"] == filtered

    def test_resume_skips_existing(self, assembly_config):
        """断点续传:输出已存在不再跑|Resume: existing output skips run"""
        out = Path(assembly_config.purge_dups_dir) / "hap1" / "seqs" / "sample1_hap1_purged.purged.fa"
        _write_fasta(out)

        class _Boom:
            def __init__(self, **kw):
                raise AssertionError("不应实例化runner|runner must not be instantiated")

        wrapper = _make_wrapper(assembly_config)
        wrapper.PurgeDupsRunner = _Boom
        hap1 = os.path.join(assembly_config.fasta_dir, "sample1_hap1.fa")
        result = wrapper.run_purge_dups(hap1, output_subdir="hap1", force_hifi_reads=True)
        assert result == str(out)

    def test_no_resume_forces_rerun(self, tmp_path):
        """--no-resume强制重跑|Force rerun with --no-resume"""
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        cfg = AssemblyConfig(
            hifi_data=str(hifi), prefix="s", base_dir=str(tmp_path / "o"),
            threads=4, resume=False,
        )
        out = Path(cfg.purge_dups_dir) / "seqs" / "s_primary_purged.purged.fa"
        _write_fasta(out)
        wrapper = _make_wrapper(cfg)
        primary = os.path.join(cfg.fasta_dir, "s_primary.fa")
        _write_fasta(primary)
        result = wrapper.run_purge_dups(primary)
        assert len(_FakeRunner.calls) == 1
        assert result == str(out)

    def test_glob_fallback_filters_by_stem(self, assembly_config):
        """glob兜底按输入stem过滤,多产物不拿错|Glob fallback filters by input stem"""
        seqs = Path(assembly_config.purge_dups_dir) / "seqs"
        _write_fasta(seqs / "sample1_hap2_purged.purged.fa")  # 先创建,无stem过滤时会拿错
        _write_fasta(seqs / "sample1_hap1_purged.purged.fa")
        wrapper = _make_wrapper(assembly_config)
        hap1 = os.path.join(assembly_config.fasta_dir, "sample1_hap1.fa")
        result = wrapper._find_purged_output(hap1)
        assert result == str(seqs / "sample1_hap1_purged.purged.fa")


# ------------------------------------------------------------ 阶段执行


def _make_assembler(tmp_path, **kwargs):
    from biopytools.hifi_hic.main import GenomeAssembler
    hifi = tmp_path / "hifi.fq"
    hifi.write_text(">r\nACGT\n")
    assembler = GenomeAssembler(
        hifi_data=str(hifi), prefix="sample1", base_dir=str(tmp_path / "out"),
        threads=4, **kwargs,
    )
    _write_fasta(os.path.join(assembler.config.fasta_dir, "sample1_primary.fa"))
    _write_fasta(os.path.join(assembler.config.fasta_dir, "sample1_hap1.fa"))
    _write_fasta(os.path.join(assembler.config.fasta_dir, "sample1_hap2.fa"))
    return assembler


class TestRunPurgeStage:
    def test_hap_failure_degrades(self, tmp_path):
        """hap失败WARNING继续|Hap failure degrades with WARNING"""
        assembler = _make_assembler(tmp_path)
        purged_primary = tmp_path / "primary_purged.fa"
        _write_fasta(purged_primary)

        def fake_run(input_fa, output_subdir=None, force_hifi_reads=False):
            return str(purged_primary) if not force_hifi_reads else None

        assembler.purge_dups_wrapper = MagicMock()
        assembler.purge_dups_wrapper.run_purge_dups.side_effect = fake_run
        results = assembler._run_purge_stage()
        assert results["primary"] == str(purged_primary)
        assert results["hap1"] is None
        assert results["hap2"] is None

    def test_primary_failure_exits(self, tmp_path):
        """primary失败exit 1|Primary failure exits"""
        assembler = _make_assembler(tmp_path)
        assembler.purge_dups_wrapper = MagicMock()
        assembler.purge_dups_wrapper.run_purge_dups.return_value = None
        with pytest.raises(SystemExit):
            assembler._run_purge_stage()

    def test_all_success(self, tmp_path):
        """全部成功返回各自路径|All success returns per-target paths"""
        assembler = _make_assembler(tmp_path)
        outputs = {}
        for label in ("primary", "hap1", "hap2"):
            p = tmp_path / f"{label}_purged.fa"
            _write_fasta(p)
            outputs[label] = str(p)

        def fake_run(input_fa, output_subdir=None, force_hifi_reads=False):
            label = output_subdir or "primary"
            return outputs[label]

        assembler.purge_dups_wrapper = MagicMock()
        assembler.purge_dups_wrapper.run_purge_dups.side_effect = fake_run
        results = assembler._run_purge_stage()
        assert results == outputs
        # hap调用必须带force_hifi_reads与子目录|Hap calls carry force flag and subdir
        calls = assembler.purge_dups_wrapper.run_purge_dups.call_args_list
        assert calls[1].kwargs["output_subdir"] == "hap1"
        assert calls[1].kwargs["force_hifi_reads"] is True


# ------------------------------------------------------------ CLI接线


class TestCliAndArgparse:
    def test_click_option_exists(self):
        """Click包装暴露flag且默认False|Click exposes flag, default False"""
        import biopytools.cli.commands.hifi_hic as hifi_cli
        param = next(p for p in hifi_cli.hifi_hic.params if p.name == "no_purge_haplotypes")
        assert param.default is False

    def test_click_flag_forwarded(self, monkeypatch):
        """flag透传到模块argv|Flag forwarded to module argv"""
        import biopytools.cli.commands.hifi_hic as hifi_cli
        from click.testing import CliRunner
        captured = {}

        def fake_main():
            captured["argv"] = list(sys.argv)

        monkeypatch.setattr(hifi_cli, "_lazy_import_assembler_main", lambda: fake_main)
        runner = CliRunner()
        with runner.isolated_filesystem():
            with open("hifi.fq", "w") as f:
                f.write(">r\nACGT\n")
            runner.invoke(hifi_cli.hifi_hic, ["-i", "hifi.fq", "--no-purge-haplotypes"])
        assert "--no-purge-haplotypes" in captured["argv"]

    def test_argparse_wires_flag(self, tmp_path, monkeypatch):
        """argparse接purge_haplotypes=False|argparse wires purge_haplotypes=False"""
        import biopytools.hifi_hic.main as m
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        captured = {}

        class _FakeAssembler:
            def __init__(self, **kw):
                captured.update(kw)

            def run_assembly(self):
                return None

        monkeypatch.setattr(m, "GenomeAssembler", _FakeAssembler)
        monkeypatch.setattr(sys, "argv",
                            ["hifi_hic.py", "--hifi", str(hifi), "--no-purge-haplotypes"])
        m.main()
        assert captured["purge_haplotypes"] is False

    def test_argparse_default_true(self, tmp_path, monkeypatch):
        """不传flag默认True|Default True without flag"""
        import biopytools.hifi_hic.main as m
        hifi = tmp_path / "hifi.fq"
        hifi.write_text(">r\nACGT\n")
        captured = {}

        class _FakeAssembler:
            def __init__(self, **kw):
                captured.update(kw)

            def run_assembly(self):
                return None

        monkeypatch.setattr(m, "GenomeAssembler", _FakeAssembler)
        monkeypatch.setattr(sys, "argv", ["hifi_hic.py", "--hifi", str(hifi)])
        m.main()
        assert captured["purge_haplotypes"] is True
