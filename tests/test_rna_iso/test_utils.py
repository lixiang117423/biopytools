"""rna_iso 工具函数测试|Tests for rna_iso utils (logger/run_command/versions yml)"""

import logging
import subprocess
from types import SimpleNamespace

import pytest

import biopytools.rna_iso.utils as rna_utils


class TestLogger:
    def test_three_handlers_and_separation(self, tmp_path):
        manager = rna_utils.RnaIsoLogger(tmp_path)
        logger = manager.get_logger()
        levels = []
        for h in logger.handlers:
            levels.append(h.level)
        assert logging.INFO in levels      # stdout
        assert logging.WARNING in levels   # stderr
        assert logging.DEBUG in levels     # file
        assert logger.name == "rna_iso"
        assert logger.propagate is False

    def test_log_file_created(self, tmp_path):
        manager = rna_utils.RnaIsoLogger(tmp_path)
        manager.get_logger().info("测试消息|test message")
        for h in manager.get_logger().handlers:
            h.flush()
        assert (tmp_path / "rna_iso.log").exists()


class TestRunCommand:
    def test_wraps_with_build_conda_and_logs(self, monkeypatch, tmp_path):
        """run_command 走 build_conda_command 并 INFO 记录完整命令|wraps and logs full command"""
        captured = {}

        def fake_build(command, args):
            captured["built"] = (command, args)
            return ["conda", "run", "-n", "isoseq_v.4.0.0", "--no-capture-output", command] + args

        monkeypatch.setattr(rna_utils, "build_conda_command", fake_build)

        def fake_run(cmd, **kwargs):
            captured["run"] = cmd
            captured["kwargs"] = kwargs
            return SimpleNamespace(returncode=0, stdout="ok", stderr="")

        monkeypatch.setattr(subprocess, "run", fake_run)

        records = []

        class FakeLogger:
            def info(self, msg):
                records.append(msg)

            def error(self, msg):
                records.append(msg)

        ok = rna_utils.run_command(
            ["/fake/env/bin/ccs", "in.bam", "out.bam", "--min-passes", "1"],
            FakeLogger(), work_dir=str(tmp_path),
        )
        assert ok is True
        assert captured["built"][0] == "/fake/env/bin/ccs"
        assert captured["run"][0] == "conda"
        assert "--no-capture-output" in captured["run"]
        assert captured["kwargs"]["cwd"] == str(tmp_path)
        # 完整命令被记录(§2.2.1)|full command logged
        joined = " ".join(captured["run"])
        assert any(joined in r for r in records)

    def test_failure_returns_false(self, monkeypatch):
        monkeypatch.setattr(rna_utils, "build_conda_command", lambda c, a: [c] + a)

        def fake_run(cmd, **kwargs):
            raise subprocess.CalledProcessError(1, cmd, stderr="boom")

        monkeypatch.setattr(subprocess, "run", fake_run)

        class FakeLogger:
            def info(self, msg):
                pass

            def error(self, msg):
                pass

        assert rna_utils.run_command(["/fake/bin/tool"], FakeLogger()) is False


class TestSoftwareVersions:
    def test_yml_contains_tools(self, monkeypatch, tmp_path):
        from biopytools.rna_iso.config import RnaIsoConfig

        (tmp_path / "s.fq.gz").write_text(">x\nACGT\n")
        (tmp_path / "ref.fa").write_text(">x\nACGT\n")
        config = RnaIsoConfig(
            reads=[str(tmp_path / "s.fq.gz")], data_type="ont",
            reference=str(tmp_path / "ref.fa"), output_dir=str(tmp_path / "out"),
        )

        def fake_run(cmd, **kwargs):
            return SimpleNamespace(returncode=0, stdout="tool v9.9.9\n", stderr="")

        monkeypatch.setattr(subprocess, "run", fake_run)

        class FakeLogger:
            def info(self, msg):
                pass

        out = tmp_path / "versions.yml"
        rna_utils.generate_software_versions_yml(config, str(out), FakeLogger())
        content = out.read_text()
        for tool in ("ccs", "isoseq3", "isoquant", "samtools"):
            assert tool in content
        assert "v9.9.9" in content
        from biopytools.rna_iso import __version__ as module_version
        assert module_version in content  # 模块版本动态断言|dynamic module version
