"""utils 模块测试|utils module tests"""

import os
import subprocess
from unittest.mock import MagicMock, patch

import pytest

from biopytools.ncbi_datasets.utils import (
    DATASETS_DOWNLOAD_URL,
    ModuleLogger,
    ensure_datasets_tool,
    format_number,
    install_datasets_tool,
    is_tool_ready,
    run_command,
)


def make_logger(tmp_path):
    """构造日志器|Build a logger"""
    return ModuleLogger(str(tmp_path / 'logs' / 'x.ncbi_datasets.log')).get_logger()


class TestFormatNumber:
    """大数字格式化|large number formatting"""

    def test_million_unit(self):
        assert format_number(10_000_000) == '10.00M'

    def test_small_number(self):
        assert format_number(1234) == '1234'


class TestToolInstall:
    """工具安装|tool installation"""

    def test_is_tool_ready(self, tmp_path):
        tool = tmp_path / 'datasets'
        tool.write_text('#!/bin/sh\n')
        assert not is_tool_ready(str(tool))
        os.chmod(tool, 0o755)
        assert is_tool_ready(str(tool))

    @patch('biopytools.ncbi_datasets.utils.subprocess.run')
    def test_install_success(self, mock_run, tmp_path):
        logger = make_logger(tmp_path)
        dest = tmp_path / 'bin' / 'datasets'

        def fake_curl(cmd, **kwargs):
            # 模拟 curl 真实写出文件|Simulate curl actually writing the file
            dest.write_bytes(b'#!/bin/sh\n')
            return subprocess.CompletedProcess(cmd, 0, stdout='', stderr='')

        mock_run.side_effect = fake_curl
        result = install_datasets_tool(str(dest), logger)
        assert result == str(dest)
        assert is_tool_ready(str(dest))
        # 命令必须记录完整 URL|full URL must be in the command
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == 'curl' and DATASETS_DOWNLOAD_URL in cmd

    @patch('biopytools.ncbi_datasets.utils.subprocess.run')
    def test_install_failure(self, mock_run, tmp_path):
        logger = make_logger(tmp_path)
        dest = tmp_path / 'bin' / 'datasets'
        mock_run.return_value = subprocess.CompletedProcess([], 1, stdout='', stderr='boom')
        assert install_datasets_tool(str(dest), logger) is None

    @patch('biopytools.ncbi_datasets.utils.install_datasets_tool')
    def test_ensure_auto_installs(self, mock_install, tmp_path):
        logger = make_logger(tmp_path)
        dest = tmp_path / 'bin' / 'datasets'
        assert not dest.exists()

        def fake_install(path, logger):
            # 模拟安装真实写出可执行文件|Simulate install writing an executable
            dest.parent.mkdir(parents=True)
            dest.write_bytes(b'#!/bin/sh\n')
            os.chmod(dest, 0o755)
            return str(dest)

        mock_install.side_effect = fake_install
        assert ensure_datasets_tool(str(dest), logger) == str(dest)
        mock_install.assert_called_once()

    def test_ensure_existing_skips_install(self, tmp_path):
        logger = make_logger(tmp_path)
        tool = tmp_path / 'datasets'
        tool.write_text('#!/bin/sh\n')
        os.chmod(tool, 0o755)
        with patch('biopytools.ncbi_datasets.utils.install_datasets_tool') as mock_install:
            assert ensure_datasets_tool(str(tool), logger) == str(tool)
            mock_install.assert_not_called()


class TestRunCommand:
    """命令执行|command execution"""

    @patch('biopytools.ncbi_datasets.utils.subprocess.run')
    def test_capture_passthrough(self, mock_run, tmp_path):
        logger = make_logger(tmp_path)
        mock_run.return_value = subprocess.CompletedProcess([], 0, stdout='ok', stderr='')
        result = run_command(['datasets', '--version'], logger, capture=True)
        assert result.stdout == 'ok'
        assert mock_run.call_args.kwargs['capture_output'] is True

    @patch('biopytools.ncbi_datasets.utils.subprocess.run')
    def test_missing_binary_returns_none(self, mock_run, tmp_path):
        logger = make_logger(tmp_path)
        mock_run.side_effect = FileNotFoundError('no such file')
        assert run_command(['datasets'], logger) is None
