"""conda解析修复守卫测试|Conda resolution fix guard tests

根因|Root cause: build_conda_command 检测到环境名后裸调 'conda run -n <env>',
作业环境 PATH 上的系统 anaconda(/share/apps/anaconda3, 只读)解析不到
miniforge3 的环境, 报 EnvironmentLocationNotFound。
修复|Fix: 用与环境目录同源的绝对 conda 可执行文件 + 'run -p <环境绝对前缀>'。
"""

import os
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from biopytools.common import conda_runner
from biopytools.common.conda_runner import (
    _extract_actual_command,
    build_conda_command,
)


def _make_fake_install(tmp_path: Path):
    """伪造一套conda安装(bin/conda + envs/<env>/bin/fastp)|Fake a conda installation"""
    base = tmp_path / "miniforge3"
    (base / "bin").mkdir(parents=True)
    conda_exe = base / "bin" / "conda"
    conda_exe.write_text("#!/bin/sh\nexit 0\n")
    conda_exe.chmod(0o755)  # _resolve_env_prefix 要求可执行|Must be executable
    env_bin = base / "envs" / "rna" / "bin"
    env_bin.mkdir(parents=True)
    (env_bin / "fastp").write_text("")
    return conda_exe, base / "envs" / "rna"


def test_build_uses_absolute_conda_and_env_prefix(tmp_path, monkeypatch):
    """环境解析到前缀后, 必须用绝对conda + -p, 不得裸调 'conda'"""
    conda_exe, env_prefix = _make_fake_install(tmp_path)
    monkeypatch.setenv("CONDA_EXE", str(conda_exe))
    # 隔离宿主机 PATH, 保证环境名来自伪安装的扫描|Isolate from host PATH
    monkeypatch.setattr(conda_runner.shutil, "which", lambda *_: None)

    cmd = build_conda_command("fastp", ["-i", "a.fq"])

    assert cmd[0] == str(conda_exe), "必须用与环境目录同源的绝对conda可执行文件"
    assert cmd[0] != "conda", "严禁裸调conda(作业PATH上可能是系统anaconda)"
    assert cmd[1] == "run"
    assert cmd[2] == "-p"
    assert os.path.realpath(cmd[3]) == os.path.realpath(str(env_prefix))
    assert "--no-capture-output" in cmd
    assert cmd[5:] == ["fastp", "-i", "a.fq"]


def test_build_falls_back_to_name_when_prefix_unresolvable(monkeypatch):
    """前缀解析不到时回退旧式 -n 形式(保持向后兼容)"""
    monkeypatch.setattr(conda_runner, "get_conda_env", lambda *a, **k: "some_env")
    monkeypatch.setattr(conda_runner, "_conda_installations", lambda: [])

    cmd = build_conda_command("fastp", ["-i", "a.fq"])

    assert cmd[:2] == ["conda", "run"]
    assert "some_env" in cmd
    idx = cmd.index("-n")
    assert cmd[idx + 1] == "some_env"
    assert "--no-capture-output" in cmd


def test_extract_actual_command_handles_prefix_form():
    """方案B提取: 新式 [conda_exe, run, -p, prefix, --no-capture-output, cmd...]"""
    wrapped = ["/x/miniforge3/bin/conda", "run", "-p", "/x/envs/rna",
               "--no-capture-output", "fastp", "-i", "a.fq"]
    assert _extract_actual_command(wrapped) == ["fastp", "-i", "a.fq"]


def test_extract_actual_command_handles_legacy_name_form():
    """方案B提取: 旧式 ['conda', 'run', '-n', env, '--no-capture-output', cmd...]"""
    wrapped = ["conda", "run", "-n", "rna", "--no-capture-output", "fastp", "-i", "a.fq"]
    assert _extract_actual_command(wrapped) == ["fastp", "-i", "a.fq"]


def test_extract_actual_command_handles_flag_order_variant():
    """方案B提取: --no-capture-output 在 -n 之前的变体"""
    wrapped = ["conda", "run", "--no-capture-output", "-n", "rna", "fastp", "-i", "a.fq"]
    assert _extract_actual_command(wrapped) == ["fastp", "-i", "a.fq"]


def test_extract_actual_command_passthrough_non_conda():
    """非conda包装的命令原样返回|Non-conda commands pass through unchanged"""
    assert _extract_actual_command(["samtools", "view", "-b"]) == ["samtools", "view", "-b"]
