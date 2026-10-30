"""FileProcessor R2推导测试|R2 derivation tests for FileProcessor

命名迁移(95d775d)配套遗漏|Migration (95d775d) follow-up gap:
pattern_suffix 推导R2用 replace('_1.', '_2.'),对下划线新命名
_1_clean.fq.gz 无效(不含'_1.'子串),R2推导不变→双端数据被静默降级单端。
"""

import importlib.util
import logging
import sys
import types
from pathlib import Path
from types import SimpleNamespace

_REPO = Path(__file__).resolve().parents[2]


def _load_file_processor():
    """伪造包上下文直载(绕过 __init__ 拉 pyfastx)|Load module bypassing heavy __init__"""
    for pkg, sub in (("biopytools.kmertools", "kmertools"),
                     ("biopytools.kmertools.kmer_count", "kmertools/kmer_count")):
        if pkg not in sys.modules:
            m = types.ModuleType(pkg)
            m.__path__ = [str(_REPO / "biopytools" / sub)]
            sys.modules[pkg] = m
    mod_name = "biopytools.kmertools.kmer_count.file_processor"
    if mod_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            mod_name, _REPO / "biopytools" / "kmertools" / "kmer_count" / "file_processor.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[mod_name] = mod
        spec.loader.exec_module(mod)
    return sys.modules[mod_name]


def _find(tmp_path, pattern):
    mod = _load_file_processor()
    cfg = SimpleNamespace(input_dir=str(tmp_path), pattern=pattern)
    processor = mod.FileProcessor(cfg, logging.getLogger("test_kmc"), None)
    return processor.find_fastq_files()


def test_new_underscore_pattern_pairs_r2(tmp_path):
    """新命名pattern必须推导出正确R2(防R1自配)|New-style pattern must derive real R2 (no self-pairing)"""
    (tmp_path / "smpN_1_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "smpN_2_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    samples = _find(tmp_path, "*_1_clean.fq.gz")

    assert len(samples) == 1
    assert samples[0][0] == "smpN"
    assert samples[0][2] is not None, "R2应配对成功,不应降级单端|R2 must pair, not degrade to single-end"
    assert samples[0][2].endswith("smpN_2_clean.fq.gz"), "R2不得配到R1自身|R2 must not resolve to the R1 file itself"


def test_legacy_dot_pattern_pairs_r2(tmp_path):
    """旧命名pattern推导R2(原有行为守卫)|Legacy dot pattern R2 derivation (existing behavior guard)"""
    (tmp_path / "smpO_1.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "smpO_2.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    samples = _find(tmp_path, "*_1.clean.fq.gz")

    assert len(samples) == 1
    assert samples[0][2] is not None
    assert samples[0][2].endswith("smpO_2.clean.fq.gz")
