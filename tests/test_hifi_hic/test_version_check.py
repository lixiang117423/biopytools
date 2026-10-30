"""hifiasm版本检查参数测试|Tests for hifiasm version check flag

事故|Incident: hifiasm 0.25.0-r726 不支持 -V(实测退出码139),
utils.check_dependencies 误报"缺少必需软件: hifiasm"致组装流程终止
(2026-09-03, joblogs/2026-09-03-10-41-hifi-hic.err)。正确参数是 --version。
"""

import subprocess
from types import SimpleNamespace

import pytest

import biopytools.hifi_hic.utils as hifi_utils
from biopytools.hifi_hic.config import AssemblyConfig


@pytest.fixture
def captured_run(monkeypatch):
    """捕获subprocess.run的命令并返回伪结果|Capture subprocess.run cmd, return fake result"""
    captured = {}

    def fake_run(cmd, **kwargs):
        captured.setdefault("cmds", []).append(cmd)
        return SimpleNamespace(returncode=0, stdout="0.25.0-r726", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    return captured


def test_check_dependencies_uses_long_flag(captured_run):
    """check_dependencies 用 --version 而非 -V|uses --version, not -V"""
    assert hifi_utils.check_dependencies(
        __import__("logging").getLogger("t"), hifiasm_bin="/nonexistent/bin/hifiasm"
    )
    assert len(captured_run["cmds"]) == 1
    cmd = captured_run["cmds"][0]
    assert "--version" in cmd
    assert "-V" not in cmd


def test_get_software_info_uses_long_flag(captured_run, tmp_path):
    """get_software_info 对 hifiasm 也用 --version|hifiasm also uses --version"""
    hifi = tmp_path / "hifi.fq"
    hifi.write_text(">r\nACGT\n")
    config = AssemblyConfig(hifi_data=str(hifi), base_dir=str(tmp_path / "out"))
    config.get_software_info()
    hifiasm_cmd = next(c for c in captured_run["cmds"] if c[0] == config.hifiasm_path)
    assert "--version" in hifiasm_cmd
    assert "-V" not in hifiasm_cmd
