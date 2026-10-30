"""fastp worker线程数钳制测试|Tests for fastp worker thread cap"""

import pytest


class FakeRunner:
    """记录调用的命令执行器替身|Stub runner recording commands"""

    def __init__(self, logger):
        self.logger = logger
        self.calls = []

    def run(self, cmd, description=""):
        self.calls.append(list(cmd))
        return True

    def check_executable(self, path, timeout=10):
        return True


@pytest.fixture
def build_core(tmp_path, monkeypatch):
    from biopytools.fastp.config import FastpConfig
    from biopytools.fastp.processing import FastpCore
    from biopytools.fastp.utils import FastpLogger

    def _build(threads):
        config = FastpConfig(
            input_dir=str(tmp_path),
            output_dir=str(tmp_path / "out"),
            threads=threads,
        )
        logger = FastpLogger(config.output_path).get_logger()
        runner = FakeRunner(logger)
        monkeypatch.setattr(
            "biopytools.fastp.processing.CommandRunner", lambda logger_: runner
        )
        return FastpCore(config, logger, runner), runner

    return _build


def worker_threads_arg(cmd):
    """从命令中提取-w参数值|Extract -w value from command"""
    return cmd[cmd.index("-w") + 1]


class TestFastpThreadCap:
    """fastp -w 钳制|-w clamping

    fastp 1.x 实测 -w>16 会进入死态(0.3核、零输出), 官方建议 worker<=16
    |fastp 1.x stalls with -w>16 (measured); docs cap workers at 16
    """

    def test_threads_88_clamped_to_16(self, build_core):
        """88线程被钳到16|threads=88 clamped to 16"""
        core, runner = build_core(88)
        r1 = core.config.input_path / "s_1.fq.gz"
        r2 = core.config.input_path / "s_2.fq.gz"
        r1.write_text("")
        r2.write_text("")
        core.create_output_directories()

        assert core.process_sample("s", r1, r2) is True
        assert worker_threads_arg(runner.calls[0]) == "16"

    def test_threads_12_unchanged(self, build_core):
        """12线程不改动|threads=12 unchanged"""
        core, runner = build_core(12)
        r1 = core.config.input_path / "s_1.fq.gz"
        r2 = core.config.input_path / "s_2.fq.gz"
        r1.write_text("")
        r2.write_text("")
        core.create_output_directories()

        assert core.process_sample("s", r1, r2) is True
        assert worker_threads_arg(runner.calls[0]) == "12"

    def test_clamp_warning_logged(self, build_core, caplog):
        """钳制时输出WARNING|WARNING logged when clamped"""
        import logging

        with caplog.at_level(logging.WARNING, logger="biopytools.fastp.utils"):
            core, _ = build_core(88)

        assert core.worker_threads == 16
        assert any(
            "16" in record.message and record.levelno == logging.WARNING
            for record in caplog.records
        )

    def test_single_end_also_clamped(self, build_core):
        """单末端模式同样钳制|Single-end mode clamped too"""
        core, runner = build_core(88)
        r1 = core.config.input_path / "s.fq.gz"
        r1.write_text("")
        core.config.single_end = True
        core.create_output_directories()

        assert core.process_sample("s", r1, None) is True
        assert worker_threads_arg(runner.calls[0]) == "16"
