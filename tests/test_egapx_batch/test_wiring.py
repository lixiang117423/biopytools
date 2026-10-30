"""cleanup_work参数三层接线测试|cleanup_work wiring tests

Why: 新参数须同时接通config字段、模块argparse、Click包装三条链路
(Click包装层不透传则CLI用户永远用不到,教训见cli-wrapper-defaults)
|The new option must be wired through config, argparse and the Click
wrapper; tests pin all three passthrough routes.
"""

import pytest
from click.testing import CliRunner

from biopytools.egapx_batch.config import EGAPxBatchConfig
from biopytools.egapx_batch import main as egapx_main
from biopytools.cli.commands import egapx_batch as egapx_command


class FakeGenerator:
    """捕获kwargs的假生成器|Kwargs-capturing fake generator"""

    captured = None

    def __init__(self, **kwargs):
        FakeGenerator.captured = kwargs

    def generate(self):
        return True


@pytest.fixture(autouse=True)
def patch_generator(monkeypatch):
    """替换main模块中的生成器|Replace the generator in main module"""
    FakeGenerator.captured = None
    monkeypatch.setattr(egapx_main, "EGAPxBatchGenerator", FakeGenerator)


@pytest.fixture
def genome_file(tmp_path):
    """最小基因组文件|Minimal genome file"""
    genome = tmp_path / "genome.fa"
    genome.write_text(">chr1\nACGT\n")
    return genome


class TestConfigField:
    """config字段|Config field"""

    def test_cleanup_work_defaults_true(self, genome_file, tmp_path):
        """默认自动清理|Cleanup on by default"""
        config = EGAPxBatchConfig(
            genome=str(genome_file), output_dir=str(tmp_path / "out")
        )
        assert config.cleanup_work is True

    def test_cleanup_work_accepts_false(self, genome_file, tmp_path):
        """可显式关闭|Can be explicitly disabled"""
        config = EGAPxBatchConfig(
            genome=str(genome_file),
            output_dir=str(tmp_path / "out"),
            cleanup_work=False,
        )
        assert config.cleanup_work is False


class TestArgparseWiring:
    """模块argparse接线|Module argparse wiring"""

    def _run_main(self, monkeypatch, argv):
        monkeypatch.setattr("sys.argv", argv)
        with pytest.raises(SystemExit) as excinfo:
            egapx_main.main()
        assert excinfo.value.code == 0
        assert FakeGenerator.captured is not None

    def test_keep_work_flag_maps_false(self, monkeypatch, genome_file, tmp_path):
        """--keep-work映射cleanup_work=False|--keep-work maps to cleanup_work=False"""
        self._run_main(monkeypatch, [
            "egapx_batch.py", "-g", str(genome_file),
            "-o", str(tmp_path / "out"), "--keep-work",
        ])
        assert FakeGenerator.captured["cleanup_work"] is False

    def test_default_maps_true(self, monkeypatch, genome_file, tmp_path):
        """缺省时cleanup_work=True|Default maps to cleanup_work=True"""
        self._run_main(monkeypatch, [
            "egapx_batch.py", "-g", str(genome_file), "-o", str(tmp_path / "out"),
        ])
        assert FakeGenerator.captured["cleanup_work"] is True


class TestCliWiring:
    """Click包装接线|Click wrapper wiring"""

    def _invoke(self, genome_file, tmp_path, extra):
        runner = CliRunner()
        return runner.invoke(
            egapx_command.egapx_batch,
            [
                "-g", str(genome_file),
                "-o", str(tmp_path / "out"),
                "-e", str(tmp_path),        # 绕过默认EGAPx路径存在性校验
                "--local-cache", str(tmp_path),  # 绕过默认local_cache存在性校验|bypass cache check
                "--sif", "",                # 绕过默认SIF存在性校验
                *extra,
            ],
            catch_exceptions=False,
        )

    def test_cli_keep_work_passthrough(self, genome_file, tmp_path):
        """CLI透传--keep-work|CLI passes --keep-work through"""
        result = self._invoke(genome_file, tmp_path, ["--keep-work"])
        assert result.exit_code == 0, result.output
        assert FakeGenerator.captured["cleanup_work"] is False

    def test_cli_default_cleanup(self, genome_file, tmp_path):
        """CLI缺省自动清理|CLI default keeps cleanup on"""
        result = self._invoke(genome_file, tmp_path, [])
        assert result.exit_code == 0, result.output
        assert FakeGenerator.captured["cleanup_work"] is True
