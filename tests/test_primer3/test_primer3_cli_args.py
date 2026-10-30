"""primer3 CLI 并行参数接线测试|Tests for primer3 CLI parallel argument wiring"""

import sys

import pytest


class TestArgparseParallelArgs:
    """argparse 层并行参数|argparse-level parallel parameters"""

    def _parse(self, monkeypatch, extra):
        from biopytools.primer3.main import parse_arguments

        argv = ["primer3.py", "-i", "in.fa", "-o", "out_dir"] + extra
        monkeypatch.setattr(sys, "argv", argv)
        return parse_arguments()

    def test_defaults(self, monkeypatch):
        """默认threads=12/阈值500|Defaults: threads=12, threshold=500"""
        args = self._parse(monkeypatch, [])
        assert args.threads == 12
        assert args.parallel_threshold == 500

    def test_accepts_custom_values(self, monkeypatch):
        """可自定义并行参数|Custom values accepted"""
        args = self._parse(monkeypatch, ["-t", "4", "--parallel-threshold", "100"])
        assert args.threads == 4
        assert args.parallel_threshold == 100


class TestClickForwarding:
    """Click 包装层透传(总透传约定)|Click wrapper forwarding convention"""

    def test_click_options_exist_with_matching_defaults(self):
        """Click层option存在且默认值与argparse一致|Click options exist with matching defaults"""
        from biopytools.cli.commands.primer3 import primer3

        defaults = {p.name: p.default for p in primer3.params}
        assert defaults.get("threads") == 12
        assert defaults.get("parallel_threshold") == 500
