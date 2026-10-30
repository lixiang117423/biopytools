"""生成脚本的work目录清理测试|Work-dir cleanup tests for generated scripts

Why: work/是Nextflow中间产物目录(可达数百GB),EGAPx成功后必须自动清理;
失败时保留以支持断点续跑。清理逻辑嵌在生成的shell脚本里,这里同时做
内容断言与真实bash行为验证|The cleanup lives in the generated shell
script; assert both content and real bash behavior.
"""

import subprocess
from pathlib import Path

import pytest

from biopytools.egapx_batch.utils import JobGenerator


@pytest.fixture
def logger():
    """静默logger|Silent logger"""
    import logging

    return logging.getLogger("test_egapx_cleanup")


@pytest.fixture
def script_paths(tmp_path):
    """生成脚本所需的目录参数|Directory args for script generation"""
    work_dir = tmp_path / "chr1" / "work"
    output_dir = tmp_path / "chr1" / "output"
    chr_dir = tmp_path / "chr1"
    yaml_path = tmp_path / "chr1" / "chr1.yaml"
    return {
        "yaml_path": str(yaml_path),
        "work_dir": str(work_dir),
        "output_dir": str(output_dir),
        "report_name": "EGAPx",
        "chr_dir": str(chr_dir),
        "egapx_dir": str(tmp_path),
        "local_cache": str(tmp_path / "local_cache"),
    }


def make_script(logger, script_paths, cleanup_work):
    """生成脚本内容|Generate script content"""
    generator = JobGenerator(logger)
    return generator.generate_script(cleanup_work=cleanup_work, **script_paths)


def run_cleanup_tail(script_content: str, egapx_exit: int, work_dir: Path,
                     tmp_path: Path) -> subprocess.CompletedProcess:
    """提取脚本尾部清理逻辑并真实执行|Extract the cleanup tail and run it for real

    从"egapx_exit=$?"行截取到脚本末尾,去掉首行后以指定退出码执行,
    验证rm -rf守卫与exit透传的真实bash语义|Slice from the exit-capture
    line to the end, drop that line, run with a given exit code to verify
    the rm -rf guard and exit-code passthrough.
    """
    marker = "egapx_exit=$?"
    assert marker in script_content, "脚本缺少退出码捕获行|Missing exit-code capture line"
    tail = script_content.split(marker, 1)[1]
    # 去掉守卫块前的空行与注释,拼上受控退出码|Prepend controlled exit code
    harness = f"egapx_exit={egapx_exit}\n" + tail
    harness_file = tmp_path / "cleanup_harness.sh"
    harness_file.write_text(harness)
    return subprocess.run(
        ["bash", str(harness_file)], capture_output=True, text=True
    )


class TestScriptCleanupContent:
    """清理块内容断言|Cleanup-block content assertions"""

    def test_default_contains_cleanup_block(self, logger, script_paths):
        """默认生成脚本含rm -rf清理块|Default script contains the rm -rf cleanup block"""
        content = make_script(logger, script_paths, cleanup_work=True)
        assert "rm -rf" in content
        assert script_paths["work_dir"] in content

    def test_cleanup_guarded_by_success_and_existence(self, logger, script_paths):
        """rm必须同时受退出码0与目录存在守卫|rm must be guarded by exit 0 and dir existence"""
        content = make_script(logger, script_paths, cleanup_work=True)
        assert "$egapx_exit -eq 0" in content
        assert "-d" in content

    def test_exit_code_passthrough(self, logger, script_paths):
        """脚本必须透传EGAPx退出码|Script must pass through the EGAPx exit code"""
        for cleanup in (True, False):
            content = make_script(logger, script_paths, cleanup_work=cleanup)
            assert "egapx_exit=$?" in content
            assert content.rstrip().endswith("exit $egapx_exit")

    def test_keep_work_has_no_rm(self, logger, script_paths):
        """cleanup_work=False时不得含rm -rf|No rm -rf when cleanup_work=False"""
        content = make_script(logger, script_paths, cleanup_work=False)
        assert "rm -rf" not in content

    @pytest.mark.parametrize("cleanup_work", [True, False])
    def test_script_is_valid_bash(self, logger, script_paths, cleanup_work, tmp_path):
        """生成脚本须通过bash语法检查|Generated script must pass bash -n"""
        content = make_script(logger, script_paths, cleanup_work=cleanup_work)
        script_file = tmp_path / f"egapx_chr1_{cleanup_work}.sh"
        script_file.write_text(content)
        result = subprocess.run(
            ["bash", "-n", str(script_file)], capture_output=True, text=True
        )
        assert result.returncode == 0, result.stderr


class TestScriptCleanupBehavior:
    """清理块真实bash行为|Real bash behavior of the cleanup block"""

    def test_success_removes_work_dir(self, logger, script_paths, tmp_path):
        """成功(退出码0)删除work目录|Exit 0 removes the work dir"""
        work_dir = Path(script_paths["work_dir"])
        work_dir.mkdir(parents=True)
        (work_dir / "intermediate.txt").write_text("x")
        content = make_script(logger, script_paths, cleanup_work=True)

        result = run_cleanup_tail(content, 0, work_dir, tmp_path)

        assert result.returncode == 0
        assert not work_dir.exists()

    def test_failure_keeps_work_dir(self, logger, script_paths, tmp_path):
        """失败(非0)保留work目录供续跑|Non-zero exit keeps the work dir for resume"""
        work_dir = Path(script_paths["work_dir"])
        work_dir.mkdir(parents=True)
        content = make_script(logger, script_paths, cleanup_work=True)

        result = run_cleanup_tail(content, 1, work_dir, tmp_path)

        assert result.returncode == 1
        assert work_dir.exists()

    def test_missing_work_dir_no_error(self, logger, script_paths, tmp_path):
        """work目录不存在时安全跳过|Missing work dir is skipped safely"""
        work_dir = Path(script_paths["work_dir"])  # 不创建|Not created
        content = make_script(logger, script_paths, cleanup_work=True)

        result = run_cleanup_tail(content, 0, work_dir, tmp_path)

        assert result.returncode == 0

    def test_keep_work_keeps_dir_on_success(self, logger, script_paths, tmp_path):
        """keep-work变体成功也不删|Keep-work variant never removes the dir"""
        work_dir = Path(script_paths["work_dir"])
        work_dir.mkdir(parents=True)
        content = make_script(logger, script_paths, cleanup_work=False)

        result = run_cleanup_tail(content, 0, work_dir, tmp_path)

        assert result.returncode == 0
        assert work_dir.exists()
