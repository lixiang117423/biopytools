"""config 单测|config tests"""

import os

import pytest

from biopytools.anno_curate.config import AnnoCurateConfig


def test_defaults_match_gsaman():
    cfg = AnnoCurateConfig(input='x.gff3', output_dir='/tmp/anno_curate_test_out')
    assert cfg.min_cds_nt == 90
    assert cfg.max_cds_nt == 10500
    assert cfg.utr5_frac_max == 0.15
    assert cfg.utr3_frac_max == 0.30
    assert cfg.relax == 0.5
    assert cfg.check_utr is True
    assert cfg.skip_diagnosis is False


def test_sample_strips_all_suffix_layers():
    cfg = AnnoCurateConfig(input='/tmp/anno_curate_test_in/anno.gff3',
                           output_dir='/tmp/anno_curate_test_out')
    assert cfg.sample == 'anno'


def test_validate_collects_all_errors(tmp_path):
    gff = tmp_path / 'a.gff3'
    gff.write_text('##gff-version 3\n')
    bam = tmp_path / 'b.bam'
    bam.write_bytes(b'BAM\x01')
    cfg = AnnoCurateConfig(
        input=str(tmp_path / 'missing.gff3'),
        output_dir=str(tmp_path / 'out'),
        rna_bams=(str(tmp_path / 'missing.bam'),),
        min_cds_nt=-1, max_cds_nt=50)
    with pytest.raises(ValueError) as ei:
        cfg.validate()
    msg = str(ei.value)
    assert 'missing.gff3' in msg        # 输入不存在
    assert 'missing.bam' in msg         # rna-bam 不存在
    assert 'min-cds-nt' in msg          # 非法阈值
    assert cfg.sample  # validate 不应改变 sample


def test_relax_negative_clamped_to_zero(tmp_path):
    gff = tmp_path / 'a.gff3'
    gff.write_text('##gff-version 3\n')
    cfg = AnnoCurateConfig(input=str(gff), output_dir=str(tmp_path / 'out'),
                           relax=-0.5)
    assert cfg.relax == 0.0   # __post_init__ clamp
