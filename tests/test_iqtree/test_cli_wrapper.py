"""CLI包装器参数透传测试|CLI wrapper argument passthrough tests

事故|Incident: cli/commands/iqtree.py 曾用 `if threads != 88` 判断是否透传-t,
用户传-t 88时恰好被吞掉,argparse默认12接管(88核算力只跑12线程)。
"""

import pytest
from click.testing import CliRunner

import biopytools.cli.commands.iqtree as iqtree_cli


@pytest.fixture
def captured_argv(monkeypatch):
    """捕获主函数收到的sys.argv|Capture sys.argv seen by main()"""
    seen = {}

    def fake_main():
        seen['argv'] = list(iqtree_cli.sys.argv)
        return 0

    monkeypatch.setattr(iqtree_cli, '_lazy_import_iqtree_main', lambda: fake_main)
    return seen


class TestThreadsPassthrough:
    """-t 透传|-t passthrough"""

    def test_threads_88_forwarded(self, captured_argv, tmp_path):
        # 88曾是被魔法数吞掉的值, 必须回归覆盖
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        result = CliRunner().invoke(
            iqtree_cli.iqtree,
            ['-i', str(aln), '-o', str(tmp_path / 'out'), '-p', 't', '-t', '88'],
        )
        assert result.exit_code == 0, result.output
        argv = captured_argv['argv']
        assert '-t' in argv
        assert argv[argv.index('-t') + 1] == '88'

    def test_threads_24_forwarded(self, captured_argv, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        result = CliRunner().invoke(
            iqtree_cli.iqtree,
            ['-i', str(aln), '-o', str(tmp_path / 'out'), '-p', 't', '-t', '24'],
        )
        assert result.exit_code == 0, result.output
        argv = captured_argv['argv']
        assert '-t' in argv

    def test_threads_default_omitted(self, captured_argv, tmp_path):
        # 默认12与argparse一致, 不传也可(行为等价)
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        result = CliRunner().invoke(
            iqtree_cli.iqtree,
            ['-i', str(aln), '-o', str(tmp_path / 'out'), '-p', 't'],
        )
        assert result.exit_code == 0, result.output
        argv = captured_argv['argv']
        assert '-t' not in argv


class TestSequenceTypePassthrough:
    """--sequence-type 透传|--sequence-type passthrough"""

    def test_sequence_type_dna_forwarded(self, captured_argv, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        result = CliRunner().invoke(
            iqtree_cli.iqtree,
            ['-i', str(aln), '-o', str(tmp_path / 'out'), '-p', 't',
             '--sequence-type', 'DNA'],
        )
        assert result.exit_code == 0, result.output
        argv = captured_argv['argv']
        assert '--sequence-type' in argv
        assert argv[argv.index('--sequence-type') + 1] == 'DNA'

    def test_sequence_type_lowercase_forwarded_upper(self, captured_argv, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        result = CliRunner().invoke(
            iqtree_cli.iqtree,
            ['-i', str(aln), '-o', str(tmp_path / 'out'), '-p', 't',
             '--sequence-type', 'aa'],
        )
        assert result.exit_code == 0, result.output
        argv = captured_argv['argv']
        assert argv[argv.index('--sequence-type') + 1] == 'AA'

    def test_sequence_type_default_omitted(self, captured_argv, tmp_path):
        # 默认None走模块内嗅探, 包装层不传
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        result = CliRunner().invoke(
            iqtree_cli.iqtree,
            ['-i', str(aln), '-o', str(tmp_path / 'out'), '-p', 't'],
        )
        assert result.exit_code == 0, result.output
        assert '--sequence-type' not in captured_argv['argv']
