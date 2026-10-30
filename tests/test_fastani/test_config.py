"""FastaniConfig 测试|FastaniConfig tests"""

import pytest

from biopytools.fastani.config import FastaniConfig


def make_config(tmp_path, **kwargs):
    """构造带假fastANI二进制的配置|Build config with a fake fastANI binary"""
    fake_bin = tmp_path / 'bin' / 'fastANI'
    fake_bin.parent.mkdir(parents=True, exist_ok=True)
    fake_bin.write_text('#!/bin/sh\n')
    kwargs.setdefault('output_dir', str(tmp_path / 'out'))
    cfg = FastaniConfig(**kwargs)
    cfg.fastani_path = str(fake_bin)
    return cfg


class TestModeExclusivity:
    """双模式互斥|Dual-mode exclusivity"""

    def test_input_with_query_rejected(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(g), query=str(g), reference=str(g))
        with pytest.raises(ValueError, match='互斥'):
            cfg.validate()

    def test_query_without_reference_rejected(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, query=str(g))
        with pytest.raises(ValueError, match='reference'):
            cfg.validate()

    def test_no_input_rejected(self, tmp_path):
        cfg = make_config(tmp_path)
        with pytest.raises(ValueError, match='至少'):
            cfg.validate()

    def test_all_vs_all_property(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(g))
        assert cfg.all_vs_all is True
        cfg2 = make_config(tmp_path, query=str(g), reference=str(g))
        assert cfg2.all_vs_all is False


class TestPathChecks:
    """路径存在性|Path existence"""

    def test_missing_input_path(self, tmp_path):
        cfg = make_config(tmp_path, input=str(tmp_path / 'nope'))
        with pytest.raises(ValueError, match='不存在'):
            cfg.validate()

    def test_missing_fastani_binary(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(g))
        cfg.fastani_path = str(tmp_path / 'no' / 'fastANI')
        with pytest.raises(ValueError, match='fastANI'):
            cfg.validate()

    def test_tilde_expansion(self, tmp_path, monkeypatch):
        monkeypatch.setenv('HOME', str(tmp_path))
        cfg = FastaniConfig(output_dir='~/fastani_out_test')
        assert '~' not in cfg.output_dir
        assert cfg.output_dir == str(tmp_path / 'fastani_out_test')


class TestNumericBounds:
    """数值边界|Numeric bounds"""

    def test_kmer_over_16_rejected(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(g), kmer=17)
        with pytest.raises(ValueError, match='k'):
            cfg.validate()

    def test_bad_threads_rejected(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(g), threads=0)
        with pytest.raises(ValueError, match='线程'):
            cfg.validate()

    def test_bad_min_fraction_rejected(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(g), min_fraction=1.5)
        with pytest.raises(ValueError, match='min'):
            cfg.validate()

    def test_bad_frag_len_rejected(self, tmp_path):
        g = tmp_path / 'g.fna'; g.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(g), frag_len=0)
        with pytest.raises(ValueError, match='frag'):
            cfg.validate()

    def test_valid_config_passes(self, tmp_path):
        g1 = tmp_path / 'a.fna'; g1.write_text('>c\nACGT\n')
        g2 = tmp_path / 'b.fna'; g2.write_text('>c\nACGT\n')
        cfg = make_config(tmp_path, input=str(tmp_path))
        assert cfg.validate() is True
