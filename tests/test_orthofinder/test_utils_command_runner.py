"""CommandRunner长命令输出判定测试|CommandRunner long-command output verdict tests

背景|Context: OrthoFinder 3.1.5 的 main() finally 块裸 sys.exit() 会把失败退出码吞成 0,
仅看 returncode 会把停滞自杀的作业误判为成功——必须叠加输出内容判定.
|The bare sys.exit() in OrthoFinder 3.1.5's main() finally block swallows the failure
exit code into 0, so returncode-only verdicts misjudge a stalled self-killed job as
success — an output-content verdict must be layered on top.
"""

import logging
import sys

import pytest

from biopytools.orthofinder.utils import CommandRunner


@pytest.fixture
def runner_and_records(tmp_path):
    """构造CommandRunner与日志捕获|Build CommandRunner with log capture"""
    logger = logging.getLogger("test_orthofinder_command_runner")
    logger.setLevel(logging.DEBUG)
    logger.handlers.clear()
    logger.propagate = False
    records = []

    class CaptureHandler(logging.Handler):
        def emit(self, record):
            records.append(record)

    logger.addHandler(CaptureHandler())
    return CommandRunner(logger, tmp_path), records


def test_clean_exit_code_zero_returns_true(runner_and_records):
    """干净退出码0应返回True|Clean exit code 0 should return True"""
    runner, _ = runner_and_records
    cmd = [sys.executable, '-c', "print('Connected 111 species')"]
    assert runner.run(cmd, 'desc', long_running=True) is True


def test_exit_zero_with_error_marker_returns_false(runner_and_records):
    """退出码0但stdout含ERROR标记应判失败|Exit code 0 with ERROR marker in stdout should fail"""
    runner, _ = runner_and_records
    cmd = [sys.executable, '-c',
           "print('ERROR: Stalled for 200.0s (completed 0/111).')"]
    assert runner.run(cmd, 'desc', long_running=True) is False


def test_exit_zero_with_ansi_wrapped_error_marker_returns_false(runner_and_records):
    """ANSI着色的ERROR标记同样应识别|ANSI-colored ERROR marker should also be detected"""
    runner, _ = runner_and_records
    cmd = [sys.executable, '-c',
           "print('\x1b[1;31mERROR: An error occurred\x1b[0m in gather')"]
    assert runner.run(cmd, 'desc', long_running=True) is False


def test_nonzero_exit_code_returns_false(runner_and_records):
    """非零退出码应返回False(回归守卫)|Non-zero exit code should return False (regression guard)"""
    runner, _ = runner_and_records
    cmd = [sys.executable, '-c', "import sys; print('boom'); sys.exit(3)"]
    assert runner.run(cmd, 'desc', long_running=True) is False


def test_error_marker_line_logged_at_error_level(runner_and_records):
    """ERROR标记行应实时以ERROR级别落日志(进.err可见)|ERROR marker line should be logged at ERROR level in real time"""
    runner, records = runner_and_records
    cmd = [sys.executable, '-c',
           "print('ERROR: Stalled for 200.0s (completed 0/111).')"]
    runner.run(cmd, 'desc', long_running=True)
    error_records = [r for r in records if r.levelno >= logging.ERROR]
    assert any('Stalled' in r.getMessage() for r in error_records)


def test_nonzero_rc_logs_stdout_tail(runner_and_records):
    """非零退出码时应补记stdout尾巴|Non-zero exit code should also log a stdout tail"""
    runner, records = runner_and_records
    cmd = [sys.executable, '-c',
           "print('unique_stdout_marker_line'); import sys; sys.exit(2)"]
    runner.run(cmd, 'desc', long_running=True)
    error_records = [r for r in records if r.levelno >= logging.ERROR]
    assert any('unique_stdout_marker_line' in r.getMessage() for r in error_records)


def test_multiline_output_fully_drained_after_exit(runner_and_records):
    """进程退出后管道残留行应被清空不丢失|Buffered lines left in pipe after exit should be drained"""
    runner, records = runner_and_records
    lines = [f"progress_line_{i}" for i in range(20)]
    code = "\n".join(f"print('{ln}')" for ln in lines)
    cmd = [sys.executable, '-c', code]
    assert runner.run(cmd, 'desc', long_running=True) is True
    # 所有行均应被读到(DEBUG级)|All lines should have been read (DEBUG level)
    debug_text = "\n".join(r.getMessage() for r in records)
    for ln in lines:
        assert ln in debug_text
