"""输入发现后缀契约测试|Input discovery suffix contract tests

同根因|Same root cause: 命名迁移(95d775d)后 R1/R2_PATTERNS 与 strip/build 后缀
漏掉旧规范 _1.clean.fq.gz/_2.clean.fq.gz(有旧 nohost 却没有旧 clean),
旧命名输入目录将匹配 0 样本。
"""

import importlib.util
import sys
import types
from pathlib import Path

_REPO = Path(__file__).resolve().parents[2]


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


def test_all_r1_pattern_lists_contain_legacy_dot():
    """每个R1模式列表都须含旧命名|Every R1 pattern list must contain legacy naming"""
    mod = _load_data_processing()
    for name in ("R1_PATTERNS_CLEAN", "R1_PATTERNS_SKIP_QC", "R1_PATTERNS_QC"):
        assert any("_1.clean.fq.gz" in p for p in getattr(mod, name)), name


def test_all_r2_pattern_lists_contain_legacy_dot():
    """每个R2模式列表都须含旧命名|Every R2 pattern list must contain legacy naming"""
    mod = _load_data_processing()
    for name in ("R2_PATTERNS_CLEAN", "R2_PATTERNS_SKIP_QC", "R2_PATTERNS_QC"):
        assert any("_2.clean.fq.gz" in p for p in getattr(mod, name)), name


def test_strip_and_build_suffixes_contain_legacy_dot():
    """strip/build 后缀须含旧命名|Strip/build suffixes must contain legacy naming"""
    mod = _load_data_processing()
    assert "_1.clean.fq.gz" in mod.R1_STRIP_SUFFIXES
    assert "_2.clean.fq.gz" in mod.R2_BUILD_SUFFIXES


def test_legacy_dot_glob_matches(tmp_path):
    """旧命名文件实际可被模式命中(软链接同样)|Legacy files actually matched (symlinks too)"""
    import os
    mod = _load_data_processing()
    target = tmp_path / "t.fq.gz"
    target.write_text("@r\nACGT\n+\nII\n")
    os.symlink(target, tmp_path / "sampleN_1.clean.fq.gz")

    hits = [p for pattern in mod.R1_PATTERNS_CLEAN for p in tmp_path.glob(pattern)]
    assert hits


def test_all_r1_pattern_lists_contain_underscore_nohost():
    """每个R1模式列表都须含下划线nohost(mixrace写出命名)|Underscore nohost in every R1 list"""
    mod = _load_data_processing()
    for name in ("R1_PATTERNS_CLEAN", "R1_PATTERNS_SKIP_QC", "R1_PATTERNS_QC"):
        assert any("_1_nohost.fq.gz" in p for p in getattr(mod, name)), name


def test_all_r2_pattern_lists_contain_underscore_nohost():
    """每个R2模式列表都须含下划线nohost|Underscore nohost in every R2 list"""
    mod = _load_data_processing()
    for name in ("R2_PATTERNS_CLEAN", "R2_PATTERNS_SKIP_QC", "R2_PATTERNS_QC"):
        assert any("_2_nohost.fq.gz" in p for p in getattr(mod, name)), name


def test_strip_and_build_suffixes_contain_underscore_nohost():
    """strip/build 后缀须含下划线nohost(否则样本名剥离错)|Underscore nohost in strip/build"""
    mod = _load_data_processing()
    assert "_1_nohost.fq.gz" in mod.R1_STRIP_SUFFIXES
    assert "_2_nohost.fq.gz" in mod.R2_BUILD_SUFFIXES


def test_underscore_nohost_glob_matches(tmp_path):
    """下划线nohost文件实际可被模式命中|Underscore nohost files actually matched"""
    mod = _load_data_processing()
    f = tmp_path / "S1_1_nohost.fq.gz"
    f.write_bytes(b"x")
    hits = [p for pat in mod.R1_PATTERNS_CLEAN for p in tmp_path.glob(pat)]
    assert str(f) in [str(h) for h in hits]
