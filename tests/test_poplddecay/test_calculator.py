"""PopLDdecay计算器命令构建测试|Tests for PopLDdecayCalculator command building"""

from unittest.mock import MagicMock

from biopytools.poplddecay.calculator import PopLDdecayCalculator
from biopytools.poplddecay.config import PopLDdecayConfig


def _make_config(tmp_path, **overrides):
    """构造最小可用配置|Build a minimal usable config"""
    input_file = tmp_path / "in.vcf"
    input_file.write_text("##fileformat=VCFv4.2\n")
    params = dict(
        input_file=str(input_file),
        output_prefix=str(tmp_path / "out"),
    )
    params.update(overrides)
    return PopLDdecayConfig(**params)


class TestThreadsWiring:

    def test_build_command_includes_threads(self, tmp_path):
        """命令应包含-T线程数|Command should include -T thread count"""
        config = _make_config(tmp_path, threads=8)
        calculator = PopLDdecayCalculator(config, MagicMock())
        cmd = calculator._build_command()
        cmd_str = " ".join(str(part) for part in cmd)
        assert " -T 8 " in f" {cmd_str} "

    def test_default_threads_is_12(self, tmp_path):
        """线程数默认12(对齐项目规范)|Default threads should be 12"""
        config = _make_config(tmp_path)
        assert config.threads == 12

    def test_validate_rejects_non_positive_threads(self, tmp_path):
        """线程数必须为正|Threads must be positive"""
        config = _make_config(tmp_path, threads=0)
        try:
            config.validate()
            raised = False
        except ValueError as e:
            raised = True
            assert "线程" in str(e) or "thread" in str(e).lower()
        assert raised


class TestPopLDdecay2Paths:

    def test_default_binary_is_poplddecay2(self, tmp_path):
        """默认二进制应为PopLDdecay2|Default binary should be PopLDdecay2"""
        config = _make_config(tmp_path)
        assert "PopLDdecay2" in config.poplddecay_path

    def test_env_var_overrides_binary_path(self, tmp_path, monkeypatch):
        """环境变量可覆盖二进制路径|Env var should override binary path"""
        custom = tmp_path / "custom_binary"
        custom.write_text("")
        monkeypatch.setenv("POPLDDECAY2_PATH", str(custom))
        config = _make_config(tmp_path)
        assert config.poplddecay_path == str(custom)

    def test_default_plot_scripts_from_poplddecay2_dir(self, tmp_path):
        """绘图脚本默认指向PopLDdecay2目录|Plot scripts should default to PopLDdecay2 dir"""
        config = _make_config(tmp_path)
        assert "PopLDdecay2" in config.plot_one_pop_path
        assert "PopLDdecay2" in config.plot_multi_pop_path
