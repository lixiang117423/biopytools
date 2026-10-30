"""CLI包装器参数透传测试|CLI wrapper argument passthrough tests"""

import sys

import pytest
from click.testing import CliRunner

import biopytools.cli.commands.rna_iso as rna_cli


@pytest.fixture
def captured_argv(monkeypatch):
    """捕获模块main收到的sys.argv|Capture sys.argv seen by module main"""
    captured = {}

    def fake_main():
        captured["argv"] = list(sys.argv)

    monkeypatch.setattr(rna_cli, "_lazy_import_main", lambda: fake_main)
    return captured


def _invoke(extra_args):
    runner = CliRunner()
    return runner.isolated_filesystem(), runner, extra_args


class TestPassthrough:
    def test_minimal(self, captured_argv):
        _, runner, _ = _invoke([])
        with runner.isolated_filesystem():
            with open("s.fq.gz", "w") as f:
                f.write(">x\nACGT\n")
            result = runner.invoke(rna_cli.rna_iso, ["-i", "s.fq.gz"])
        assert result.exit_code == 0
        argv = captured_argv["argv"]
        assert "--reads" in argv and "s.fq.gz" in argv
        # 默认值不透传|defaults not forwarded
        assert "--data-type" not in argv
        assert "--threads" not in argv
        assert "--engine" not in argv

    def test_full_options_forwarded(self, captured_argv):
        runner = CliRunner()
        with runner.isolated_filesystem():
            for name in ("s.fq.gz", "ref.fa", "anno.gtf", "primers.fa"):
                with open(name, "w") as f:
                    f.write(">x\nACGT\n")
            runner.invoke(rna_cli.rna_iso, [
                "-i", "s.fq.gz", "--data-type", "ont",
                "-g", "ref.fa", "--genedb", "anno.gtf",
                "--engine", "both", "--primers", "primers.fa",
                "--min-passes", "3", "-t", "24", "-p", "S1",
                "-o", "out/",
            ])
        argv = captured_argv["argv"]
        for flag, val in (("--data-type", "ont"), ("--reference", "ref.fa"),
                          ("--genedb", "anno.gtf"), ("--engine", "both"),
                          ("--primers", "primers.fa"), ("--min-passes", "3"),
                          ("--threads", "24"), ("--prefix", "S1"),
                          ("--output-dir", "out/")):
            i = argv.index(flag)
            assert argv[i + 1] == val

    def test_multiple_reads(self, captured_argv):
        runner = CliRunner()
        with runner.isolated_filesystem():
            for name in ("a.fq.gz", "b.fq.gz"):
                with open(name, "w") as f:
                    f.write(">x\nACGT\n")
            runner.invoke(rna_cli.rna_iso, ["-i", "a.fq.gz", "-i", "b.fq.gz"])
        argv = captured_argv["argv"]
        assert "a.fq.gz" in argv and "b.fq.gz" in argv

    def test_explicit_default_thread_semantics(self, captured_argv):
        """显式-t 12(=argparse默认)被哨兵吞掉,但语义等价(最终threads=12)
        |Explicit -t 12 (== argparse default) swallowed by sentinel check,
        semantically equivalent (final threads=12). 与iqtree -t 88不同:那次
        哨兵值与argparse默认不一致才是bug"""
        runner = CliRunner()
        with runner.isolated_filesystem():
            with open("s.fq.gz", "w") as f:
                f.write(">x\nACGT\n")
            result = runner.invoke(rna_cli.rna_iso, ["-i", "s.fq.gz", "-t", "12"])
        assert result.exit_code == 0
        # 被吞但最终值等价|swallowed yet equivalent
        assert "--threads" not in captured_argv["argv"]
