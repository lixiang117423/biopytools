"""vcf2deepbsa config 测试|vcf2deepbsa config tests"""

import pytest

from biopytools.vcf2deepbsa.config import Vcf2DeepBsaConfig


class TestConfig:
    """配置类|Config class"""

    def test_output_dir_created(self, tmp_path, make_vcf):
        vcf = make_vcf()
        out = tmp_path / "result_dir"
        cfg = Vcf2DeepBsaConfig(input_vcf=vcf, output_dir=str(out))
        assert out.is_dir()

    def test_validate_missing_input_raises(self, tmp_path):
        cfg = Vcf2DeepBsaConfig(
            input_vcf=str(tmp_path / "no_such.vcf"),
            output_dir=str(tmp_path / "out"))
        with pytest.raises(ValueError, match="输入文件不存在|Input file not found"):
            cfg.validate()

    def test_validate_ok(self, make_vcf, tmp_path):
        cfg = Vcf2DeepBsaConfig(input_vcf=make_vcf(),
                                output_dir=str(tmp_path / "out"))
        assert cfg.validate() is True

    def test_paths_expanded(self, monkeypatch, tmp_path, make_vcf):
        vcf = make_vcf()
        monkeypatch.setenv("HOME", str(tmp_path))
        cfg = Vcf2DeepBsaConfig(input_vcf=vcf, output_dir="~/out_deepbsa")
        assert "~" not in str(cfg.output_dir)
        assert str(cfg.output_dir).startswith(str(tmp_path))
