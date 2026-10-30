"""CLI与argparse双层接线测试|CLI and argparse wiring tests

Why: Click包装层须显式透传新参数(改模块argparse默认值无效的教训),
以测试钉死三条透传链路|The Click wrapper must explicitly pass through new
options; tests pin all three passthrough routes.
"""

import sys
from pathlib import Path

import pytest
from click.testing import CliRunner

from biopytools.cli.commands import orthofinder as orthofinder_command
from biopytools.orthofinder import main as orthofinder_main


@pytest.fixture
def dirs(tmp_path):
    """构造输入与BLAST续跑目录|Build input and BLAST resume directories"""
    input_dir = tmp_path / "proteins"
    input_dir.mkdir()
    (input_dir / "a.fa").write_text(">s\nMK\n")
    blast_dir = tmp_path / "blast_wd"
    blast_dir.mkdir()
    return input_dir, blast_dir


def test_argparse_wires_resume_fields(tmp_path, dirs, monkeypatch):
    """模块argparse须接线三个新参数|Module argparse must wire the three new options"""
    input_dir, blast_dir = dirs
    captured = {}

    class FakeAnalyzer:
        def __init__(self, **kwargs):
            captured.update(kwargs)

        def run_analysis(self):
            captured['ran'] = True

    monkeypatch.setattr(orthofinder_main, 'OrthoFinderPangenomeAnalyzer', FakeAnalyzer)
    monkeypatch.setattr(sys, 'argv', [
        'orthofinder.py', '-i', str(input_dir), '-o', str(tmp_path / 'out'),
        '--resume-from-blast', str(blast_dir),
        '--old-version',
        '--orthofinder-extra', '-M dendroblast',
    ])
    orthofinder_main.main()
    assert captured['resume_blast_dir'] == str(blast_dir)
    assert captured['use_old_version'] is True
    assert captured['orthofinder_extra'] == '-M dendroblast'
    assert captured['ran'] is True


def test_click_passes_resume_options(tmp_path, dirs, monkeypatch):
    """Click包装层须透传三个新参数|Click wrapper must pass through the three new options"""
    input_dir, blast_dir = dirs
    captured = {}

    def fake_main():
        captured['argv'] = list(sys.argv)

    monkeypatch.setattr('biopytools.orthofinder.main.main', fake_main)
    result = CliRunner().invoke(orthofinder_command.orthofinder, [
        '-i', str(input_dir), '-o', str(tmp_path / 'out'),
        '--resume-from-blast', str(blast_dir),
        '--old-version',
        '--orthofinder-extra', '-M dendroblast',
    ])
    assert result.exit_code == 0, result.output
    argv = captured['argv']
    i = argv.index('--resume-from-blast')
    assert argv[i + 1] == str(blast_dir)
    assert '--old-version' in argv
    j = argv.index('--orthofinder-extra')
    assert argv[j + 1] == '-M dendroblast'
