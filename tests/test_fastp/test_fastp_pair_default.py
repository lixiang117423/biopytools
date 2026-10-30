"""fastp pair步骤默认关闭测试|Tests for fastp pair step default-off behavior"""

import sys

import pytest


class TestArgparsePairDefault:
    """argparse 层 pair 开关默认值|argparse-level pair flag defaults"""

    def _parse(self, monkeypatch, extra):
        from biopytools.fastp.main import parse_arguments

        argv = ["fastp.py", "-i", "raw_data", "-o", "clean_data"] + extra
        monkeypatch.setattr(sys, "argv", argv)
        return parse_arguments()

    def test_default_pair_disabled(self, monkeypatch):
        """默认不启用pair修复|Pair repair disabled by default"""
        args = self._parse(monkeypatch, [])
        assert args.enable_pair is False
        assert args.disable_pair is False

    def test_enable_pair_opt_in(self, monkeypatch):
        """--enable-pair显式开启|--enable-pair opts in"""
        args = self._parse(monkeypatch, ["--enable-pair"])
        assert args.enable_pair is True

    def test_disable_pair_still_accepted(self, monkeypatch):
        """--disable-pair保留可用(向后兼容)|--disable-pair kept for backward compat"""
        args = self._parse(monkeypatch, ["--disable-pair"])
        assert args.disable_pair is True
        assert args.enable_pair is False

    def test_disable_pair_wins_over_enable(self, monkeypatch):
        """同时给两个开关时禁用优先|disable wins when both flags given"""
        args = self._parse(monkeypatch, ["--enable-pair", "--disable-pair"])
        assert args.enable_pair is True
        assert args.disable_pair is True


class TestConfigPairDefault:
    """config 层默认值|config-level default"""

    def test_config_default_pair_disabled(self, tmp_path):
        """FastpConfig默认enable_pair=False|FastpConfig defaults enable_pair=False"""
        from biopytools.fastp.config import FastpConfig

        config = FastpConfig(input_dir=str(tmp_path), output_dir=str(tmp_path / "out"))
        assert config.enable_pair is False


class TestClickForwarding:
    """Click 包装层接线|Click wrapper wiring"""

    def _click_defaults(self):
        from biopytools.cli.commands.fastp import fastp

        return {p.name: p.default for p in fastp.params}

    def test_click_enable_pair_default_off(self):
        """Click层--enable-pair默认False|Click --enable-pair defaults False"""
        assert self._click_defaults()["enable_pair"] is False

    def test_click_forwards_enable_pair(self):
        """Click层显式透传--enable-pair|Click forwards --enable-pair explicitly"""
        import inspect

        from biopytools.cli.commands import fastp as fastp_cmd

        source = inspect.getsource(fastp_cmd)
        assert "'--enable-pair'" in source, (
            "Click层必须在enable_pair时向argparse透传--enable-pair|"
            "Click layer must forward --enable-pair to argparse"
        )


class TestPairDisabledBehavior:
    """enable_pair=False 时的流程行为|Pipeline behavior when pair disabled"""

    @pytest.fixture
    def core(self, tmp_path, monkeypatch):
        from biopytools.fastp.config import FastpConfig
        from biopytools.fastp.processing import FastpCore
        from biopytools.fastp.utils import FastpLogger

        class FakeRunner:
            def __init__(self, logger):
                self.calls = []

            def run(self, cmd, description=""):
                self.calls.append(list(cmd))
                return True

            def check_executable(self, path, timeout=10):
                return True

        config = FastpConfig(
            input_dir=str(tmp_path),
            output_dir=str(tmp_path / "out"),
            enable_pair=False,
        )
        logger = FastpLogger(config.output_path).get_logger()
        fake = FakeRunner(logger)
        monkeypatch.setattr(
            "biopytools.fastp.processing.CommandRunner", lambda logger: fake
        )
        return FastpCore(config, logger, fake), fake

    def test_paired_end_runs_fastp_only(self, core, tmp_path):
        """双端模式默认只执行一次fastp,不跑seqkit pair|Paired-end runs fastp once, no seqkit"""
        fastp_core, fake_runner = core
        r1 = tmp_path / "s1_1.fq.gz"
        r2 = tmp_path / "s1_2.fq.gz"
        r1.write_text("")
        r2.write_text("")

        fastp_core.create_output_directories()
        assert fastp_core.process_sample("s1", r1, r2) is True

        # 只调用一次(fastp),且无seqkit字样|Exactly one call (fastp), no seqkit
        assert len(fake_runner.calls) == 1
        assert "seqkit" not in " ".join(fake_runner.calls[0])
        assert "fastp" in fake_runner.calls[0][0]

    def test_no_temp_dirs_when_pair_disabled(self, core):
        """pair关闭时不创建temp_after_fastp/unpaired目录|No temp dirs when pair disabled"""
        fastp_core, _ = core
        fastp_core.create_output_directories()

        assert not (fastp_core.config.output_path / "temp_after_fastp").exists()
        assert not (fastp_core.config.output_path / "unpaired").exists()
