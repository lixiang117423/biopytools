"""primer3 并行参数配置测试|Tests for primer3 parallel parameters in config"""

import pytest


class TestPrimer3ConfigParallel:
    """threads/parallel_threshold 配置与校验|Config fields and validation"""

    @pytest.fixture
    def make_config(self, tmp_path):
        fasta = tmp_path / "input.fa"
        fasta.write_text(">s\nACGT\n")
        core = tmp_path / "primer3_core_bin"
        core.write_text("")

        def _make(**overrides):
            from biopytools.primer3.config import Primer3Config

            params = dict(
                input_fasta=str(fasta),
                output_dir=str(tmp_path / "out"),
                primer3_core_path=str(core),
            )
            params.update(overrides)
            return Primer3Config(**params)

        return _make

    def test_default_threads_and_threshold(self, make_config):
        """默认12线程/500条阈值|Defaults: threads=12, threshold=500"""
        config = make_config()
        assert config.threads == 12
        assert config.parallel_threshold == 500

    def test_validate_rejects_nonpositive_threads(self, make_config):
        """线程数必须为正数|threads must be positive"""
        config = make_config(threads=0)
        with pytest.raises(ValueError, match="线程|thread"):
            config.validate()

    def test_validate_rejects_negative_threshold(self, make_config):
        """阈值不允许负数|threshold must not be negative"""
        config = make_config(parallel_threshold=-1)
        with pytest.raises(ValueError, match="阈值|threshold"):
            config.validate()

    def test_validate_accepts_parallel_params(self, make_config):
        """合法并行参数通过校验|Valid parallel params pass validation"""
        config = make_config(threads=4, parallel_threshold=100)
        assert config.validate() is True

    def test_evaluator_accepts_parallel_params(self, build_evaluator):
        """评估器透传并行参数|Evaluator forwards parallel params"""
        evaluator = build_evaluator(threads=4, parallel_threshold=100)
        assert evaluator.config.threads == 4
        assert evaluator.config.parallel_threshold == 100
