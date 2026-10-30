"""generate_fof_file 旧命名识别测试|Legacy naming tests for generate_fof_file

同根因|Same root cause: 命名迁移(95d775d)后 pipeline_processor 的 R1 glob 链
漏掉旧规范 _1.clean.fq.gz,旧命名输入目录将报 "未找到FASTQ文件"。
"""

import importlib.util
import logging
import sys
import types
from pathlib import Path
from types import SimpleNamespace

_REPO = Path(__file__).resolve().parents[2]


def _load_pipeline_processor():
    """伪造包上下文直载(绕过 __init__ 拉 pyfastx)|Load module bypassing heavy __init__"""
    pkg_name = "biopytools.kmertools"
    if pkg_name not in sys.modules:
        pkg = types.ModuleType(pkg_name)
        pkg.__path__ = [str(_REPO / "biopytools" / "kmertools")]
        sys.modules[pkg_name] = pkg
    rdb_name = f"{pkg_name}.rocksdb_importer"
    if rdb_name not in sys.modules:
        sys.modules[rdb_name] = types.ModuleType(rdb_name)
    mod_name = f"{pkg_name}.pipeline_processor"
    if mod_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            mod_name, _REPO / "biopytools" / "kmertools" / "pipeline_processor.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[mod_name] = mod
        spec.loader.exec_module(mod)
    return sys.modules[mod_name]


def _make_processor(tmp_path):
    mod = _load_pipeline_processor()
    cfg = SimpleNamespace(input_dir=str(tmp_path), output_path=tmp_path / "out")
    cfg.output_path.mkdir()
    return mod.PipelineProcessor(cfg, logging.getLogger("test_kmertools"), None)


def test_legacy_dot_pair_generates_fof(tmp_path):
    """旧命名 _1.clean.fq.gz 应产出FOF|Legacy dot naming must yield FOF"""
    (tmp_path / "sampleK_1.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "sampleK_2.clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    processor = _make_processor(tmp_path)
    result = processor.generate_fof_file()

    assert result is not None
    fof = tmp_path / "out" / "input.fof"
    assert fof.exists()
    assert "sampleK_1.clean.fq.gz" in fof.read_text()


def test_new_underscore_pair_still_works(tmp_path):
    """新命名 _1_clean.fq.gz 默认直配|New underscore naming matched by default"""
    (tmp_path / "sampleJ_1_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (tmp_path / "sampleJ_2_clean.fq.gz").write_text("@r\nACGT\n+\nII\n")

    processor = _make_processor(tmp_path)
    result = processor.generate_fof_file()

    assert result is not None
    assert (tmp_path / "out" / "input.fof").exists()
