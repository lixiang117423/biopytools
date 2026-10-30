"""CLI包装器genome-size透传测试|CLI wrapper genome-size passthrough tests

事故|Incident: cli/commands/hifi_hic.py 曾用 `if genome_size != "1.45g"` 哨兵判断
是否透传(与 iqtree -t 88 同族bug);hifi_hic_workflow CLI 则完全没有 --genome-size
选项,用户永远只能用模块内写死的 1.45g。
"""

import pytest
from click.testing import CliRunner

import biopytools.cli.commands.hifi_hic as hifi_cli
import biopytools.cli.commands.hifi_hic_workflow as wf_cli


@pytest.fixture
def captured_hifi_args(monkeypatch):
    """捕获hifi_hic模块main收到的sys.argv|Capture sys.argv seen by module main"""
    captured = {}

    def fake_main():
        captured["argv"] = list(__import__("sys").argv)

    monkeypatch.setattr(hifi_cli, "_lazy_import_assembler_main", lambda: fake_main)
    return captured


class TestHifiHicCli:
    def test_genome_size_no_default(self):
        """CLI层默认None(哨兵比较不可靠)|CLI default is None"""
        param = next(p for p in hifi_cli.hifi_hic.params if p.name == "genome_size")
        assert param.default is None

    def test_genome_size_forwarded_when_given(self, captured_hifi_args):
        """显式传入时透传|forwarded when explicitly given"""
        runner = CliRunner()
        with runner.isolated_filesystem():
            with open("hifi.fq", "w") as f:
                f.write(">r\nACGT\n")
            runner.invoke(hifi_cli.hifi_hic, ["-i", "hifi.fq", "-g", "800m"])
        argv = captured_hifi_args["argv"]
        i = argv.index("--genome-size")
        assert argv[i + 1] == "800m"

    def test_genome_size_omitted_when_default(self, captured_hifi_args):
        """未传入时不透传(交给模块auto)|omitted when not given"""
        runner = CliRunner()
        with runner.isolated_filesystem():
            with open("hifi.fq", "w") as f:
                f.write(">r\nACGT\n")
            runner.invoke(hifi_cli.hifi_hic, ["-i", "hifi.fq"])
        assert "--genome-size" not in captured_hifi_args["argv"]


class TestWorkflowCli:
    def test_exposes_genome_size_option(self):
        """workflow CLI 必须暴露 --genome-size|must expose --genome-size"""
        param = next(
            p for p in wf_cli.hifi_hic_workflow.params if p.name == "genome_size"
        )
        assert param.default is None

    def test_genome_size_forwarded_when_given(self, monkeypatch):
        import sys

        import biopytools.hifi_hic_workflow.main as wf_main_mod

        captured = {}

        def fake_main():
            captured["argv"] = list(sys.argv)

        monkeypatch.setattr(wf_main_mod, "main", fake_main)
        runner = CliRunner()
        with runner.isolated_filesystem():
            for name in ("hifi.fq", "r1.fq", "r2.fq", "ref.fa"):
                with open(name, "w") as f:
                    f.write(">r\nACGT\n")
            runner.invoke(
                wf_cli.hifi_hic_workflow,
                ["--hifi", "hifi.fq", "--hic-r1", "r1.fq", "--hic-r2", "r2.fq",
                 "--ref", "ref.fa", "-o", "out", "--genome-size", "900m"],
            )
        assert "--genome-size" in captured["argv"]
        i = captured["argv"].index("--genome-size")
        assert captured["argv"][i + 1] == "900m"
