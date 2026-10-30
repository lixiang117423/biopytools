"""generate_fof_file 异常样品路径基准测试|Invalid-sample FOF path base tests

根因|Root cause: 异常样品名分支 FOF 写 cleaned_data/... (相对 output_path),
而 kmtricks 子进程按 cwd 解析 FOF 内路径(GATB bank 打开相对 cwd);
cwd != output_path 时报 Unable to open bank (joblog 2026-09-10-18-02)。
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


def _make_processor(input_dir: Path, output_dir: Path):
    mod = _load_pipeline_processor()
    cfg = SimpleNamespace(input_dir=str(input_dir), output_path=output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    return mod.PipelineProcessor(cfg, logging.getLogger("test_kmertools"), None)


def _fof_entries(fof_path: Path):
    """解析FOF条目 (样品名, [文件路径])|Parse FOF entries (sample, [file paths])"""
    entries = []
    for line in fof_path.read_text().splitlines():
        if not line.strip():
            continue
        name, _, files_part = line.partition(":")
        files = [p.strip() for p in files_part.split(";")]
        entries.append((name.strip(), files))
    return entries


def test_invalid_pair_paths_resolve_from_cwd(tmp_path, monkeypatch):
    """异常样品名FOF路径必须相对cwd可解析|Invalid-sample FOF paths must resolve from cwd

    复现真实场景: cwd=输出根目录, output_path=其下 02_kmerdb 子目录。
    """
    monkeypatch.chdir(tmp_path)
    input_dir = tmp_path / "01_data"
    input_dir.mkdir()
    (input_dir / "120.sorted.genome_1.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (input_dir / "120.sorted.genome_2.fq.gz").write_text("@r\nACGT\n+\nII\n")

    output_dir = tmp_path / "02_kmerdb"
    processor = _make_processor(input_dir, output_dir)
    result = processor.generate_fof_file()

    assert result is not None
    entries = _fof_entries(output_dir / "input.fof")
    assert len(entries) == 1
    name, files = entries[0]
    assert name == "120_sorted_genome"
    assert len(files) == 2
    for f in files:
        # kmtricks 子进程按 cwd 解析FOF内路径|kmtricks resolves FOF paths against cwd
        assert Path(f).exists(), f"FOF路径从cwd不可解析|FOF path unresolvable from cwd: {f}"


def test_invalid_single_end_paths_resolve_from_cwd(tmp_path, monkeypatch):
    """单末端异常样品名FOF路径必须相对cwd可解析|Single-end invalid-sample paths resolve from cwd"""
    monkeypatch.chdir(tmp_path)
    input_dir = tmp_path / "01_data"
    input_dir.mkdir()
    (input_dir / "demo.sample.fq.gz").write_text("@r\nACGT\n+\nII\n")

    output_dir = tmp_path / "02_kmerdb"
    processor = _make_processor(input_dir, output_dir)
    result = processor.generate_fof_file()

    assert result is not None
    entries = _fof_entries(output_dir / "input.fof")
    assert len(entries) == 1
    name, files = entries[0]
    assert name == "demo_sample"
    assert len(files) == 1
    assert Path(files[0]).exists(), f"FOF路径从cwd不可解析|FOF path unresolvable from cwd: {files[0]}"


def test_normal_pair_paths_resolve_from_cwd(tmp_path, monkeypatch):
    """正常样品名FOF路径相对cwd可解析(回归)|Normal-sample paths resolve from cwd (regression)"""
    monkeypatch.chdir(tmp_path)
    input_dir = tmp_path / "01_data"
    input_dir.mkdir()
    (input_dir / "sampleA_1.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (input_dir / "sampleA_2.fq.gz").write_text("@r\nACGT\n+\nII\n")

    output_dir = tmp_path / "02_kmerdb"
    processor = _make_processor(input_dir, output_dir)
    result = processor.generate_fof_file()

    assert result is not None
    entries = _fof_entries(output_dir / "input.fof")
    assert len(entries) == 1
    name, files = entries[0]
    assert name == "sampleA"
    for f in files:
        assert Path(f).exists(), f"FOF路径从cwd不可解析|FOF path unresolvable from cwd: {f}"


def test_cleaned_copy_skipped_when_dst_exists(tmp_path, monkeypatch):
    """已存在的清理副本重跑跳过拷贝(断点续传)|Existing cleaned copies skipped on rerun"""
    monkeypatch.chdir(tmp_path)
    input_dir = tmp_path / "01_data"
    input_dir.mkdir()
    (input_dir / "120.sorted.genome_1.fq.gz").write_text("@r\nACGT\n+\nII\n")
    (input_dir / "120.sorted.genome_2.fq.gz").write_text("@r\nACGT\n+\nII\n")

    output_dir = tmp_path / "02_kmerdb"
    mod = _load_pipeline_processor()
    processor = _make_processor(input_dir, output_dir)
    assert processor.generate_fof_file() is not None
    cleaned_r1 = output_dir / "cleaned_data" / "120_sorted_genome_1.fq.gz"
    assert cleaned_r1.exists()

    # 重跑时不再重复拷贝|Rerun must not re-copy existing cleaned files
    calls = []
    monkeypatch.setattr(mod.shutil, "copy2",
                        lambda src, dst: calls.append((str(src), str(dst))))
    assert processor.generate_fof_file() is not None
    assert calls == [], f"重跑仍触发拷贝|Rerun re-copied: {calls}"
