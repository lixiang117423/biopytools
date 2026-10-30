"""conda_run_prefix守卫测试|conda_run_prefix guard tests

整条shell命令(含管道)包装场景需要字符串前缀; 必须用同源conda绝对路径 + -p,
严禁裸调 'conda'(作业PATH上可能是系统anaconda, §13)。
"""

import os
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from biopytools.common import conda_runner
from biopytools.common.conda_runner import conda_run_prefix


def _make_fake_install(tmp_path: Path):
    base = tmp_path / "miniforge3"
    (base / "bin").mkdir(parents=True)
    conda_exe = base / "bin" / "conda"
    conda_exe.write_text("#!/bin/sh\nexit 0\n")
    conda_exe.chmod(0o755)
    env_bin = base / "envs" / "align" / "bin"
    env_bin.mkdir(parents=True)
    (env_bin / "samtools").write_text("")
    return conda_exe, base / "envs" / "align"


def test_prefix_uses_absolute_conda_and_env_prefix(tmp_path, monkeypatch):
    conda_exe, env_prefix = _make_fake_install(tmp_path)
    monkeypatch.setenv("CONDA_EXE", str(conda_exe))
    monkeypatch.setattr(conda_runner.shutil, "which", lambda *_: None)

    prefix = conda_run_prefix("samtools")

    assert prefix is not None
    assert prefix.startswith(str(conda_exe)), "前缀必须以同源conda绝对路径开头"
    assert " run -p " in prefix
    assert os.path.realpath(prefix.split(" run -p ")[1].split(" ")[0]) == os.path.realpath(str(env_prefix))
    assert prefix.endswith("--no-capture-output")


def test_prefix_none_for_non_conda_tool(tmp_path, monkeypatch):
    _make_fake_install(tmp_path)
    monkeypatch.setattr(conda_runner.shutil, "which", lambda *_: None)
    assert conda_run_prefix("definitely_not_a_real_tool_xyz") is None
