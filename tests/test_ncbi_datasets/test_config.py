"""NCBIDatasetsConfig 测试|NCBIDatasetsConfig tests"""

import pytest

from biopytools.ncbi_datasets.config import NCBIDatasetsConfig


def make_config(tmp_path, **kwargs):
    """构造配置|Build a config"""
    taxon = kwargs.pop('taxon', 67593)
    kwargs.setdefault('output_dir', str(tmp_path / 'out'))
    return NCBIDatasetsConfig(taxon=taxon, **kwargs)


class TestTaxon:
    """taxon 参数|taxon parameter"""

    def test_int_coercion(self, tmp_path):
        cfg = make_config(tmp_path, taxon='67593')
        assert cfg.taxon == 67593

    def test_non_numeric_rejected(self, tmp_path):
        with pytest.raises(ValueError, match='正整数'):
            NCBIDatasetsConfig(taxon='abc', output_dir=str(tmp_path / 'out'))

    def test_non_positive_rejected(self, tmp_path):
        cfg = make_config(tmp_path, taxon=0)
        with pytest.raises(ValueError, match='正整数'):
            cfg.validate()


class TestOutputLayout:
    """by-step 输出目录|by-step output layout"""

    def test_step_dirs_created(self, tmp_path):
        cfg = make_config(tmp_path)
        assert cfg.info_dir.is_dir()
        assert cfg.download_dir.is_dir()
        assert cfg.logs_dir.is_dir()

    def test_file_names(self, tmp_path):
        cfg = make_config(tmp_path)
        assert cfg.manifest_file.name == '67593_assemblies.tsv'
        assert cfg.zip_file.name == '67593_genomes.zip'
        assert cfg.dataset_dir.name == '67593.ncbi_dataset'
        assert cfg.log_file.name == '67593_ncbi_datasets.log'
        assert cfg.manifest_file.parent.name == '00_pipeline_info'
        assert cfg.zip_file.parent.name == '01_download'

    def test_tilde_expansion(self, tmp_path, monkeypatch):
        monkeypatch.setenv('HOME', str(tmp_path))
        cfg = NCBIDatasetsConfig(taxon=67593, output_dir='~/out')
        assert str(cfg.output_path).startswith(str(tmp_path))


class TestFilters:
    """筛选参数校验|filter validation"""

    def test_valid_assembly_source(self, tmp_path):
        cfg = make_config(tmp_path, assembly_source='RefSeq')
        cfg.validate()
        assert cfg.assembly_source == 'refseq'

    def test_invalid_assembly_source(self, tmp_path):
        cfg = make_config(tmp_path, assembly_source='genbankx')
        with pytest.raises(ValueError, match='assembly-source'):
            cfg.validate()

    def test_valid_assembly_level(self, tmp_path):
        cfg = make_config(tmp_path, assembly_level='complete,Chromosome')
        cfg.validate()

    def test_invalid_assembly_level(self, tmp_path):
        cfg = make_config(tmp_path, assembly_level='complete,bogus')
        with pytest.raises(ValueError, match='assembly-level'):
            cfg.validate()


class TestToolPath:
    """工具路径解析|tool path resolution"""

    def test_explicit_path_wins(self, tmp_path):
        cfg = make_config(tmp_path, datasets_path='~/mybin/datasets')
        assert cfg.resolved_datasets_path.endswith('/mybin/datasets')

    def test_env_var_used(self, tmp_path, monkeypatch):
        monkeypatch.setenv('DATASETS_PATH', str(tmp_path / 'env' / 'datasets'))
        cfg = make_config(tmp_path)
        assert cfg.resolved_datasets_path == str(tmp_path / 'env' / 'datasets')

    def test_default_fallback(self, tmp_path, monkeypatch):
        monkeypatch.delenv('DATASETS_PATH', raising=False)
        monkeypatch.setenv('HOME', str(tmp_path))
        # 屏蔽真实用户配置(USER_CONFIG_FILE 是模块级常量,monkeypatch HOME 无效)
        # |Neutralize the real user config (USER_CONFIG_FILE is module-level)
        import biopytools.common.paths as paths
        monkeypatch.setattr(paths, 'load_user_config', lambda: {})
        cfg = make_config(tmp_path)
        assert cfg.resolved_datasets_path == str(tmp_path / 'bin' / 'datasets')
