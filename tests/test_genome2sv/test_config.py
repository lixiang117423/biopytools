"""Genome2SVConfig flank 参数测试|Genome2SVConfig flank parameter tests"""

import pytest

from biopytools.genome2sv.config import Genome2SVConfig


def make_config(tmp_path, **kwargs):
    """构造最小可用配置|Build a minimal valid config"""
    ref_fa = tmp_path / "ref.fa"
    ref_fa.write_text(">c1\nACGT\n")
    q_fa = tmp_path / "q.fa"
    q_fa.write_text(">c1\nACGT\n")
    fof = tmp_path / "samples.fof"
    fof.write_text(f"REF\t{ref_fa}\nQ1\t{q_fa}\n")
    kwargs.setdefault("input_fof", str(fof))
    kwargs.setdefault("ref_sample", "REF")
    kwargs.setdefault("output_dir", str(tmp_path / "out"))
    return Genome2SVConfig(**kwargs)


class TestFlankConfig:
    def test_default_flank_300(self, tmp_path):
        """默认侧翼 300bp|default flank is 300"""
        cfg = make_config(tmp_path)
        assert cfg.flank == 300

    def test_flank_dir_created(self, tmp_path):
        """06_sv_flanks 目录自动创建|06_sv_flanks dir auto-created"""
        cfg = make_config(tmp_path)
        assert cfg.flank_dir.name == "06_sv_flanks"
        assert cfg.flank_dir.exists()

    def test_negative_flank_rejected(self, tmp_path):
        """负侧翼报错|negative flank rejected"""
        cfg = make_config(tmp_path, flank=-1)
        with pytest.raises(ValueError, match="flank"):
            cfg.validate()

    def test_zero_flank_allowed(self, tmp_path):
        """flank=0 合法(仅 SV 本体)|flank=0 allowed (SV interval only)"""
        cfg = make_config(tmp_path, flank=0)
        assert cfg.validate()
