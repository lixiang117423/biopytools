"""conda_env_run_prefix守卫测试|conda_env_run_prefix guard tests

显式环境名的整命令包装(rnaseq_val/vcf2gwas CLI等)需要不校验命令的前缀;
必须用同源conda绝对路径 + -p(§13)。
"""

import os
import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from biopytools.common import conda_runner
from biopytools.common.conda_runner import conda_env_run_prefix


def _make_fake_install(tmp_path: Path):
    base = tmp_path / "miniforge3"
    (base / "bin").mkdir(parents=True)
    conda_exe = base / "bin" / "conda"
    conda_exe.write_text("#!/bin/sh\nexit 0\n")
    conda_exe.chmod(0o755)
    (base / "envs" / "myenv" / "bin").mkdir(parents=True)
    return conda_exe, base / "envs" / "myenv"


def test_env_prefix_uses_absolute_conda(tmp_path, monkeypatch):
    conda_exe, env_dir = _make_fake_install(tmp_path)
    monkeypatch.setenv("CONDA_EXE", str(conda_exe))
    monkeypatch.setattr(conda_runner.shutil, "which", lambda *_: None)

    prefix = conda_env_run_prefix("myenv")

    assert prefix.startswith(str(conda_exe))
    assert " run -p " in prefix
    assert os.path.realpath(prefix.split(" run -p ")[1].split(" ")[0]) == os.path.realpath(str(env_dir))
    assert prefix.endswith("--no-capture-output")


def test_env_prefix_falls_back_when_env_missing(tmp_path, monkeypatch):
    conda_exe, _ = _make_fake_install(tmp_path)
    monkeypatch.setenv("CONDA_EXE", str(conda_exe))
    monkeypatch.setattr(conda_runner.shutil, "which", lambda *_: None)

    prefix = conda_env_run_prefix("no_such_env_xyz")

    assert prefix == "conda run -n no_such_env_xyz --no-capture-output"
