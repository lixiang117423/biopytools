"""genome_size 默认值测试|Tests for genome_size default

需求|Requirement: 移除写死的 1.45g 默认——用户不传 --genome-size 时
不向 hifiasm 传 --hg-size(由 hifiasm 自动估计 read coverage,官方默认 auto);
显式传入时才透传。事故:某样品实际基因组远小于 1.45Gb,静默按 1.45g 跑偏了覆盖度推断。
"""

import logging

import pytest

import biopytools.hifi_hic.assembler as asm_mod
from biopytools.hifi_hic.assembler import HifiasmAssembler
from biopytools.hifi_hic.config import AssemblyConfig


def _make_config(tmp_path, **kwargs):
    """构造最小AssemblyConfig(不执行外部命令)|Build minimal AssemblyConfig"""
    hifi = tmp_path / "hifi.fq"
    hifi.write_text(">r\nACGT\n")
    cfg_kwargs = dict(hifi_data=str(hifi), base_dir=str(tmp_path / "out"))
    cfg_kwargs.update(kwargs)
    return AssemblyConfig(**cfg_kwargs)


def _capture_assembly_cmd(monkeypatch, config):
    """捕获hifiasm组装命令行|Capture hifiasm assembly command line"""
    captured = {}

    def fake_run_shell(shell_cmd, tool_paths, logger, work_dir=None):
        captured["cmd"] = shell_cmd
        return True

    monkeypatch.setattr(asm_mod, "run_shell_command", fake_run_shell)
    assembler = HifiasmAssembler(config, logging.getLogger("test_hifi_hic"))
    assert assembler.run_assembly()
    return captured["cmd"]


class TestGenomeSizeDefault:
    def test_default_is_none(self, tmp_path):
        """默认 genome_size=None(不预设 1.45g)|default is None, not 1.45g"""
        assert _make_config(tmp_path).genome_size is None

    def test_validate_ok_without_genome_size(self, tmp_path):
        """不传 genome_size 校验通过|validation passes when unset"""
        assert _make_config(tmp_path).validate()

    def test_validate_rejects_bad_format(self, tmp_path):
        """显式传入非法格式仍报错|bad format still rejected when set"""
        with pytest.raises(ValueError, match="基因组大小格式错误"):
            _make_config(tmp_path, genome_size="1.45").validate()

    def test_validate_accepts_good_format(self, tmp_path):
        assert _make_config(tmp_path, genome_size="800m").validate()


class TestHgSizePassing:
    def test_no_hg_size_when_unset(self, tmp_path, monkeypatch):
        """未设置时不传 --hg-size|hifiasm auto 估覆盖度|omit --hg-size when unset"""
        cmd = _capture_assembly_cmd(monkeypatch, _make_config(tmp_path))
        assert "--hg-size" not in cmd

    def test_hg_size_passed_when_set(self, tmp_path, monkeypatch):
        """设置时按值透传|forward value when set"""
        cmd = _capture_assembly_cmd(
            monkeypatch, _make_config(tmp_path, genome_size="1.2g")
        )
        assert "--hg-size" in cmd
        assert "1.2g" in cmd


def test_workflow_config_default_is_none():
    """hifi_hic_workflow 的 genome_size 默认同样为 None|workflow default is None too"""
    from biopytools.hifi_hic_workflow.config import HifiHicWorkflowConfig

    field = HifiHicWorkflowConfig.__dataclass_fields__["genome_size"]
    assert field.default is None
