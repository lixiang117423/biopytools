"""filter_snp_indel conda调用走公共层测试|Conda invocation delegation to common layer tests

根因(2026-09-13 gtx-71gene作业)|Root cause: 模块utils.py持有本地build_conda_command
副本, 裸调 'conda run -n align'——作业环境PATH上的系统anaconda解析不到align环境,
依赖检查误判"BCFtools不可用", step5在GTX联合检测成功后崩溃。
修复|Fix: 删除本地副本, 统一走 common/conda_runner(同源conda绝对路径 + run -p, §13)。
"""

import os
import sys
from pathlib import Path
from types import SimpleNamespace

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from biopytools.common import conda_runner
from biopytools.filter_snp_indel import utils as filter_utils


def _make_fake_install(tmp_path: Path):
    """伪造一套conda安装(bin/conda + envs/align/bin/bcftools)|Fake a conda installation"""
    base = tmp_path / "miniforge3"
    (base / "bin").mkdir(parents=True)
    conda_exe = base / "bin" / "conda"
    conda_exe.write_text("#!/bin/sh\nexit 0\n")
    conda_exe.chmod(0o755)
    env_bin = base / "envs" / "align" / "bin"
    env_bin.mkdir(parents=True)
    (env_bin / "bcftools").write_text("")
    return conda_exe, base / "envs" / "align"


def test_build_conda_command_uses_absolute_conda_and_prefix(tmp_path, monkeypatch):
    """模块的build_conda_command必须走公共层契约: 绝对conda + -p, 严禁裸调conda"""
    conda_exe, env_prefix = _make_fake_install(tmp_path)
    monkeypatch.setenv("CONDA_EXE", str(conda_exe))
    monkeypatch.setattr(conda_runner.shutil, "which", lambda *_: None)

    cmd = filter_utils.build_conda_command("bcftools", ["--version"])

    assert cmd[0] == str(conda_exe), "必须用与环境目录同源的绝对conda可执行文件"
    assert cmd[0] != "conda", "严禁裸调conda(作业PATH上可能是系统anaconda)"
    assert cmd[1:3] == ["run", "-p"]
    assert os.path.realpath(cmd[3]) == os.path.realpath(str(env_prefix))
    assert "--no-capture-output" in cmd
    assert cmd[5:] == ["bcftools", "--version"]


def test_check_dependencies_probes_through_wrapped_command(tmp_path, monkeypatch):
    """依赖检查的探测命令也必须是包装后的形态(裸调时作业环境PATH上没有bcftools)"""
    conda_exe, env_prefix = _make_fake_install(tmp_path)
    monkeypatch.setenv("CONDA_EXE", str(conda_exe))
    monkeypatch.setattr(conda_runner.shutil, "which", lambda *_: None)

    captured = {}

    def fake_run(cmd, **kwargs):
        captured["cmd"] = cmd
        return SimpleNamespace(returncode=0, stdout="bcftools 1.17", stderr="")

    monkeypatch.setattr(filter_utils.subprocess, "run", fake_run)
    cfg = SimpleNamespace(bcftools_path="bcftools")

    assert filter_utils.check_dependencies(cfg, _silent_logger()) is True
    cmd = captured["cmd"]
    assert os.path.basename(cmd[0]) == "conda" and cmd[0] != "conda"
    assert cmd[1:3] == ["run", "-p"]
    assert cmd[3].rstrip("/").endswith("/envs/align")


def _silent_logger():
    import logging
    return logging.getLogger("test_filter_snp_indel")
