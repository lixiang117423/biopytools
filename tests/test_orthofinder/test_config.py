"""配置新字段测试|Config new-field tests"""

import pytest

from biopytools.orthofinder.config import PangenomeConfig


def make_config(tmp_path, **overrides):
    """构造测试配置|Build a test config"""
    input_dir = tmp_path / "proteins"
    input_dir.mkdir(exist_ok=True)
    (input_dir / "a.fa").write_text(">s\nMK\n")
    defaults = dict(input_dir=str(input_dir), output_dir=str(tmp_path / "out"))
    defaults.update(overrides)
    return PangenomeConfig(**defaults)


def test_new_resume_fields_default_values(tmp_path):
    """续跑三参数默认值|Default values of the three resume parameters"""
    config = make_config(tmp_path)
    assert config.resume_blast_dir is None
    assert config.use_old_version is False
    assert config.orthofinder_extra == ''


def test_resume_blast_dir_expands_home(tmp_path, monkeypatch):
    """resume_blast_dir须在__post_init__展开~|resume_blast_dir must expand ~ in __post_init__"""
    monkeypatch.setenv('HOME', str(tmp_path))
    config = make_config(tmp_path, resume_blast_dir='~/blast_wd')
    assert config.resume_blast_dir == str(tmp_path / 'blast_wd')


def test_validate_rejects_missing_blast_dir(tmp_path):
    """续跑目录不存在应报错|Missing resume directory should fail validation"""
    config = make_config(tmp_path, resume_blast_dir=str(tmp_path / 'nope'))
    with pytest.raises(ValueError, match='BLAST'):
        config.validate()
