"""mixrace nohost 输入嗅探测试|nohost-input sniffing tests

根因(2026-09-14 两作业 GTX 0 R1)|Root cause: 用户 --clean-fastq-dir 指向已是
寄主剔除产物的目录时,mixrace 仍无条件重跑寄主剔除并把 GTX 输入钉死在
<output>/02_host_filter——用户指定的目录被旁路。修复: 嗅探输入目录是否已含
nohost 产物,是则跳过寄主剔除、GTX 直接用该目录。
"""

import sys
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

from biopytools.mixrace.samples import dir_has_nohost_input


def test_underscore_nohost_detected(tmp_path):
    """下划线风格 nohost(mixrace 自身写出的命名)须被识别"""
    (tmp_path / "S1_1_nohost.fq.gz").write_bytes(b"x")
    (tmp_path / "S1_2_nohost.fq.gz").write_bytes(b"x")
    assert dir_has_nohost_input(str(tmp_path)) is True


def test_dot_nohost_detected(tmp_path):
    """旧点号风格 nohost 也须被识别|legacy dot-style nohost also detected"""
    (tmp_path / "S1_1.nohost.fq.gz").write_bytes(b"x")
    (tmp_path / "S1_2.nohost.fq.gz").write_bytes(b"x")
    assert dir_has_nohost_input(str(tmp_path)) is True


def test_plain_clean_dir_not_detected(tmp_path):
    """普通 clean/raw 目录不得误判(否则会跳过必要的寄主剔除)"""
    (tmp_path / "S1_1_clean.fq.gz").write_bytes(b"x")
    (tmp_path / "S1_2_clean.fq.gz").write_bytes(b"x")
    (tmp_path / "S2_1.fq.gz").write_bytes(b"x")
    assert dir_has_nohost_input(str(tmp_path)) is False


def test_missing_or_empty_dir_not_detected(tmp_path):
    """目录不存在或为空 → False(不嗅探崩溃)"""
    assert dir_has_nohost_input(str(tmp_path / "no_such_dir")) is False
    assert dir_has_nohost_input(str(tmp_path)) is False
