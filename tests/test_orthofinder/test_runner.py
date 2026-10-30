"""OrthoFinderRunner修复行为测试|OrthoFinderRunner fix-behavior tests

覆盖|Covers:
- 复制前产物校验(Orthogroups.tsv+GeneCount.tsv)与现场保留|Pre-copy artifact validation and scene preservation
- basic模式附带--no-fix-files(跳过-og后的GetOrthologues白跑)|--no-fix-files in basic mode
- --resume-from-blast替换-f/绕过已有结果短路/结果目录搜索扩展|-b replaces -f, bypasses short-circuit, extended results search
- 目标已存在时改名保留而非删除|Rename-aside instead of delete when target exists
"""

import logging
from pathlib import Path

import pytest

from biopytools.orthofinder.config import PangenomeConfig
from biopytools.orthofinder.orthofinder_runner import OrthoFinderRunner


@pytest.fixture
def logger():
    lg = logging.getLogger("test_orthofinder_runner")
    lg.setLevel(logging.DEBUG)
    lg.handlers.clear()
    lg.propagate = False
    return lg


class FakeCmdRunner:
    """记录命令并返回预设结果|Record commands and return preset result"""

    def __init__(self, result=True):
        self.result = result
        self.commands = []

    def run(self, cmd, description="", long_running=False):
        self.commands.append((list(cmd), description, long_running))
        return self.result


def make_config(tmp_path, **overrides):
    """构造测试配置|Build a test config"""
    input_dir = tmp_path / "proteins"
    input_dir.mkdir(exist_ok=True)
    (input_dir / "a.fa").write_text(">s\nMK\n")
    defaults = dict(input_dir=str(input_dir), output_dir=str(tmp_path / "out"))
    defaults.update(overrides)
    return PangenomeConfig(**defaults)


def populate_results(of_base: Path, project_name: str, with_key_files=True, tsv=True):
    """在OrthoFinder基目录下直接放置Results目录|Place a Results dir directly under the OrthoFinder base"""
    results = of_base / f"Results_{project_name}_Sep16"
    og_dir = results / "Orthogroups"
    og_dir.mkdir(parents=True, exist_ok=True)
    if with_key_files:
        og_file = og_dir / ("Orthogroups.tsv" if tsv else "Orthogroups.txt")
        og_file.write_text("og\ta\n")
        (og_dir / "Orthogroups.GeneCount.tsv").write_text("og\t1\n")
    return results


@pytest.fixture
def identity_conda(monkeypatch):
    """build_conda_command替换为恒等拼接|Replace build_conda_command with identity concatenation"""
    monkeypatch.setattr(
        'biopytools.orthofinder.orthofinder_runner.build_conda_command',
        lambda path, args: [path] + list(args)
    )


def test_basic_only_command_includes_no_fix_files(tmp_path, identity_conda, logger):
    """basic模式命令须含-og与--no-fix-files|Basic-mode command must include -og and --no-fix-files"""
    # Why: 预置有效结果会触发续跑短路导致命令不执行, 关闭续跑以检视命令
    config = make_config(tmp_path, resume_from_existing=False)
    populate_results(Path(config.input_dir) / "OrthoFinder", config.project_name)
    fake = FakeCmdRunner(result=True)
    runner = OrthoFinderRunner(config, logger, fake)
    assert runner.run_orthofinder() is True
    of_args = fake.commands[0][0][1:]
    assert '-og' in of_args
    assert '--no-fix-files' in of_args


def test_resume_from_blast_replaces_f_flag(tmp_path, identity_conda, logger):
    """-b续跑须替换-f并透传--old-version与额外参数|-b resume must replace -f and pass through --old-version and extra args"""
    blast_wd = tmp_path / "blast" / "OrthoFinder" / "Results_old" / "WorkingDirectory"
    blast_wd.mkdir(parents=True)
    config = make_config(tmp_path, resume_blast_dir=str(blast_wd),
                         use_old_version=True, orthofinder_extra='-M dendroblast')
    populate_results(blast_wd.parent.parent, config.project_name)  # 新Results落在OF基目录|fresh results land in the OF base
    fake = FakeCmdRunner(result=True)
    runner = OrthoFinderRunner(config, logger, fake)
    assert runner.run_orthofinder() is True
    of_args = fake.commands[0][0][1:]
    assert '-b' in of_args
    assert str(blast_wd) in of_args
    assert '-f' not in of_args
    assert '--old-version' in of_args
    assert of_args[-2:] == ['-M', 'dendroblast']


def test_resume_from_blast_bypasses_existing_results_shortcircuit(tmp_path, identity_conda, logger):
    """显式-b续跑不被已有结果短路|Explicit -b resume must not be short-circuited by existing results"""
    blast_wd = tmp_path / "blast" / "OrthoFinder" / "Results_old" / "WorkingDirectory"
    blast_wd.mkdir(parents=True)
    config = make_config(tmp_path, resume_blast_dir=str(blast_wd))
    populate_results(blast_wd.parent.parent, config.project_name)  # 新Results落在OF基目录|fresh results land in the OF base
    # 输出目录预置一份有效旧结果,默认resume逻辑本会直接短路返回|Pre-place valid old results in output
    target = Path(config.output_dir) / "01_orthofinder_analysis" / "orthofinder_results"
    og_dir = target / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text("og\ta\n")
    (og_dir / "Orthogroups.GeneCount.tsv").write_text("og\t1\n")
    fake = FakeCmdRunner(result=True)
    runner = OrthoFinderRunner(config, logger, fake)
    assert runner.run_orthofinder() is True
    assert len(fake.commands) == 1


def test_incomplete_results_fail_without_copying(tmp_path, identity_conda, logger):
    """产物缺失时判失败且保留现场|Incomplete artifacts must fail without copying or deleting"""
    config = make_config(tmp_path)
    results = populate_results(Path(config.input_dir) / "OrthoFinder", config.project_name, with_key_files=False)
    (results / "Log.txt").write_text("ERROR: Stalled for 200.0s (completed 0/111).\n")
    fake = FakeCmdRunner(result=True)  # OF吞码后表面成功|superficially succeeded due to swallowed code
    runner = OrthoFinderRunner(config, logger, fake)
    assert runner.run_orthofinder() is False
    target = Path(config.output_dir) / "01_orthofinder_analysis" / "orthofinder_results"
    assert not target.exists()
    assert results.exists()


def test_orthogroups_txt_fallback_accepted(tmp_path, identity_conda, logger):
    """Orthogroups.txt可作为tsv兜底通过校验|Orthogroups.txt is accepted as tsv fallback"""
    # Why: 关闭续跑短路, 强制走定位+校验路径|Disable short-circuit to force the locate+validate path
    config = make_config(tmp_path, resume_from_existing=False)
    populate_results(Path(config.input_dir) / "OrthoFinder", config.project_name, tsv=False)
    fake = FakeCmdRunner(result=True)
    runner = OrthoFinderRunner(config, logger, fake)
    assert runner.run_orthofinder() is True


def test_existing_target_renamed_instead_of_deleted(tmp_path, identity_conda, logger):
    """目标已存在时默认改名保留|--force未开时旧结果改名保留"""
    config = make_config(tmp_path, resume_from_existing=False)
    populate_results(Path(config.input_dir) / "OrthoFinder", config.project_name)
    target = Path(config.output_dir) / "01_orthofinder_analysis" / "orthofinder_results"
    og_dir = target / "Orthogroups"
    og_dir.mkdir(parents=True)
    (og_dir / "Orthogroups.tsv").write_text("old\n")
    fake = FakeCmdRunner(result=True)
    runner = OrthoFinderRunner(config, logger, fake)
    assert runner.run_orthofinder() is True
    olds = list(target.parent.glob("orthofinder_results.old.*"))
    assert len(olds) == 1
    assert (olds[0] / "Orthogroups" / "Orthogroups.tsv").read_text() == "old\n"
    assert (target / "Orthogroups" / "Orthogroups.tsv").exists()


def test_replace_directory_renames_by_default(tmp_path, logger):
    """_replace_directory默认改名|_replace_directory renames by default"""
    config = make_config(tmp_path)
    runner = OrthoFinderRunner(config, logger, FakeCmdRunner())
    target = tmp_path / "somewhere" / "data"
    target.mkdir(parents=True)
    (target / "m.txt").write_text("old")
    runner._replace_directory(target)
    olds = list(target.parent.glob("data.old.*"))
    assert len(olds) == 1
    assert (olds[0] / "m.txt").read_text() == "old"
    assert not target.exists()


def test_replace_directory_deletes_with_force(tmp_path, logger):
    """--force时_replace_directory直接删除|_replace_directory deletes outright with --force"""
    config = make_config(tmp_path, force_overwrite=True)
    runner = OrthoFinderRunner(config, logger, FakeCmdRunner())
    target = tmp_path / "somewhere" / "data"
    target.mkdir(parents=True)
    (target / "m.txt").write_text("old")
    runner._replace_directory(target)
    assert not target.exists()
    assert list(target.parent.glob("data.old.*")) == []
