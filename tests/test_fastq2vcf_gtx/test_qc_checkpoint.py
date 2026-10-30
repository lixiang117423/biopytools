"""fastq2vcf_gtx 质控检查点加固测试|QC checkpoint hardening tests

根因|Root cause: biopytools fastp 全失败仍退出码 0 时, step1 误判成功并写下
quality_control.done, 重跑时质控被检查点直接跳过。加固两道防线:
1. 质控"成功"但清洁目录 0 文件 → 判失败, 不写检查点;
2. 检查点存在但清洁目录 0 文件(陈旧检查点) → 作废重跑, 不得直接跳过。
"""

import importlib.util
import logging
import sys
import types
from pathlib import Path
from types import SimpleNamespace

_REPO = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(_REPO))

_LOGGER = logging.getLogger("test_fastq2vcf_gtx")


def _load_data_processing():
    """伪造包上下文直载 data_processing(绕过 __init__ 重导入)|Load module bypassing heavy __init__"""
    pkg_name = "biopytools.fastq2vcf_gtx"
    if pkg_name not in sys.modules:
        pkg = types.ModuleType(pkg_name)
        pkg.__path__ = [str(_REPO / "biopytools" / "fastq2vcf_gtx")]
        sys.modules[pkg_name] = pkg
    mod_name = f"{pkg_name}.data_processing"
    if mod_name not in sys.modules:
        spec = importlib.util.spec_from_file_location(
            mod_name, _REPO / "biopytools" / "fastq2vcf_gtx" / "data_processing.py")
        mod = importlib.util.module_from_spec(spec)
        sys.modules[mod_name] = mod
        spec.loader.exec_module(mod)
    return sys.modules[mod_name]


class _FakeRunner:
    """记录命令并返回预设结果的执行器|Recording runner with preset result"""

    def __init__(self, success=True):
        self.calls = []
        self.success = success

    def run(self, command, description=""):
        self.calls.append(command)
        return self.success


def _make_controller(tmp_path: Path, runner: _FakeRunner):
    mod = _load_data_processing()
    raw = tmp_path / "raw"
    raw.mkdir()
    (raw / "a_1.fq.gz").write_bytes(b"x")
    (raw / "a_2.fq.gz").write_bytes(b"x")
    clean = tmp_path / "clean"
    clean.mkdir()
    cfg = SimpleNamespace(
        enable_checkpoint=True, skip_qc=False, raw_fastq_dir=str(raw),
        clean_fastq_dir=str(clean), threads=4,
        read1_pattern_fastp="_1.fq.gz", read2_pattern_fastp="_2.fq.gz",
        output_dir=str(tmp_path),
    )
    return mod.QualityController(cfg, _LOGGER, runner), cfg


def _checkpoint_file(cfg) -> Path:
    return Path(cfg.output_dir) / "00_pipeline_info" / "checkpoints" / "quality_control.done"


def test_success_with_clean_files_creates_checkpoint(tmp_path):
    runner = _FakeRunner(success=True)
    qc, cfg = _make_controller(tmp_path, runner)
    (Path(cfg.clean_fastq_dir) / "a_1_clean.fq.gz").write_bytes(b"x")

    assert qc.run_quality_control() is True
    assert _checkpoint_file(cfg).exists()


def test_success_with_zero_clean_files_fails_without_checkpoint(tmp_path):
    """质控命令退出码0但没产出任何清洁文件 → 判失败, 不写检查点"""
    runner = _FakeRunner(success=True)
    qc, cfg = _make_controller(tmp_path, runner)

    assert qc.run_quality_control() is False
    assert not _checkpoint_file(cfg).exists(), "零产出不得写检查点"


def test_stale_checkpoint_with_zero_clean_files_reruns_qc(tmp_path):
    """检查点存在但清洁目录为空(历史毒化检查点) → 作废并重跑, 不得跳过"""
    runner = _FakeRunner(success=True)
    qc, cfg = _make_controller(tmp_path, runner)
    _checkpoint_file(cfg).parent.mkdir(parents=True, exist_ok=True)
    _checkpoint_file(cfg).write_text("2026-09-13 00:00:00")

    result = qc.run_quality_control()

    assert len(runner.calls) == 1, "陈旧检查点必须触发重跑而非直接跳过"
    assert result is False, "重跑后仍零产出 → 判失败"
    assert not _checkpoint_file(cfg).exists(), "作废的检查点应被移除"


def test_healthy_checkpoint_skips_qc(tmp_path):
    """检查点存在且清洁文件在 → 正常跳过(断点续传语义不变)"""
    runner = _FakeRunner(success=True)
    qc, cfg = _make_controller(tmp_path, runner)
    (Path(cfg.clean_fastq_dir) / "a_1_clean.fq.gz").write_bytes(b"x")
    _checkpoint_file(cfg).parent.mkdir(parents=True, exist_ok=True)
    _checkpoint_file(cfg).write_text("2026-09-13 00:00:00")

    assert qc.run_quality_control() is True
    assert len(runner.calls) == 0, "健康检查点应跳过质控"
