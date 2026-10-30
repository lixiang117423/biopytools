"""EviAnnConfig 测试|EviAnnConfig tests"""

import pytest

from biopytools.eviann.config import EviAnnConfig


def make_config(tmp_path, **kwargs):
    """构造配置|Build a config"""
    g = tmp_path / "g.fa"
    g.write_text(">c\nACGT\n")
    reads = tmp_path / "reads"
    reads.mkdir()
    (reads / "A_R1.fq").write_text("x")
    (reads / "A_R2.fq").write_text("x")
    kwargs.setdefault("genome", str(g))
    kwargs.setdefault("output_dir", str(tmp_path / "out"))
    kwargs.setdefault("rnaseq_data", str(reads))
    return EviAnnConfig(**kwargs)


class TestUniprotDefault:
    """UniProt 路径默认解析|uniprot default resolution"""

    def test_default_fallback(self, tmp_path, monkeypatch):
        """无配置无环境变量 → 默认 ~/database 路径|no config/env → ~/database default"""
        import biopytools.common.paths as paths
        monkeypatch.setattr(paths, "load_user_config", lambda: {})
        monkeypatch.setenv("HOME", str(tmp_path))
        cfg = make_config(tmp_path)
        assert cfg.uniprot == str(
            tmp_path / "database" / "uniprot" / "uniprot_sprot.fasta")

    def test_config_databases_priority(self, tmp_path, monkeypatch):
        """配置文件 databases 段优先于默认|config databases overrides default"""
        import biopytools.common.paths as paths
        monkeypatch.setattr(
            paths, "load_user_config",
            lambda: {"databases": {"uniprot_sprot": "~/custom/sprot.fasta"}})
        monkeypatch.setenv("HOME", str(tmp_path))
        cfg = make_config(tmp_path)
        assert cfg.uniprot == str(tmp_path / "custom" / "sprot.fasta")

    def test_env_var_priority(self, tmp_path, monkeypatch):
        """环境变量优先|env var takes priority"""
        import biopytools.common.paths as paths
        monkeypatch.setattr(paths, "load_user_config", lambda: {})
        monkeypatch.setenv("UNIPROT_SPROT_PATH", str(tmp_path / "env" / "sprot.fa"))
        cfg = make_config(tmp_path)
        assert cfg.uniprot == str(tmp_path / "env" / "sprot.fa")

    def test_explicit_overrides_default(self, tmp_path, monkeypatch):
        """显式传入覆盖默认|explicit value wins"""
        import biopytools.common.paths as paths
        monkeypatch.setattr(paths, "load_user_config", lambda: {})
        cfg = make_config(tmp_path, uniprot="~/my/sprot.fasta")
        assert cfg.uniprot.endswith("/my/sprot.fasta")


class TestValidate:
    """参数校验|validation"""

    def test_missing_uniprot_rejected(self, tmp_path):
        """默认路径文件不存在 → 明确报错(而非下载失败)|missing file → clear error"""
        cfg = make_config(tmp_path)
        cfg.uniprot = str(tmp_path / "nope" / "sprot.fasta")
        with pytest.raises(ValueError, match="UniProt"):
            cfg.validate()
