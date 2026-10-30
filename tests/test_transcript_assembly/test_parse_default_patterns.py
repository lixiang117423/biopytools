"""SampleParser 默认模式旧命名测试|Legacy naming tests for default patterns

同根因|Same root cause: 命名迁移(95d775d)后 _parse_with_default_patterns
的模式列表补了新名 *_1_clean.fq.gz 却漏旧名 *_1.clean.fq.gz,
旧命名输入目录解析出 0 样本。
"""

import importlib.util
import logging
import sys
import types
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]


def _load_assembly():
    """伪造包上下文直载 assembly(绕过 __init__ 重导入)|Load module bypassing heavy __init__"""
    pkg_name = "biopytools.transcript_assembly"
    if pkg_name not in sys.modules:
        pkg = types.ModuleType(pkg_name)
        pkg.__path__ = [str(_REPO / "biopytools" / "transcript_assembly")]
        sys.modules[pkg_name] = pkg
    mod_name = f"{pkg_name}.assembly"
    if mod_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            mod_name, _REPO / "biopytools" / "transcript_assembly" / "assembly.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[mod_name] = mod
        spec.loader.exec_module(mod)
    return sys.modules[mod_name]


def _parse(tmp_path):
    mod = _load_assembly()
    parser = mod.SampleParser(logging.getLogger("test_ta"))
    return parser.parse_input_samples(str(tmp_path))


def test_legacy_dot_pair_parsed_by_default(tmp_path):
    """旧命名 _1.clean.fq.gz 不传pattern也应解析|Legacy dot naming parsed without explicit pattern"""
    (tmp_path / "smpA_1.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "smpA_2.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    samples = _parse(tmp_path)

    assert len(samples) == 1
    assert samples[0]["name"] == "smpA"
    assert samples[0]["fastq1"].endswith("smpA_1.clean.fq.gz")
    assert samples[0]["fastq2"].endswith("smpA_2.clean.fq.gz")


def test_new_underscore_pair_parsed_by_default(tmp_path):
    """新命名 _1_clean.fq.gz 默认解析|New underscore naming parsed by default"""
    (tmp_path / "smpB_1_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "smpB_2_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    samples = _parse(tmp_path)

    assert len(samples) == 1
    assert samples[0]["name"] == "smpB"
