"""FastK调用参数守卫测试|FastK invocation parameter guard tests

背景|Context: FastK的-t是k-mer表计数下限(非线程),-T才是线程;
历史上threads被同时传给-t和-T,导致ktab只含计数>=threads的高覆盖k-mer,
下游smudgeplot分箱全0(1.71.2修复)。本测试防止该参数错位回归。
FastK's -t is the k-mer table count threshold (NOT threads); -T is threads.
Threads were once passed to both, truncating ktab to count>=threads k-mers
and zeroing downstream smudgeplot bins (fixed in 1.71.2).
"""

import logging
from unittest import mock

import pytest

from biopytools.smudgescope.utils import SmudgeplotRunner


@pytest.fixture
def runner():
    return SmudgeplotRunner(logging.getLogger("test_fastk"))


@pytest.fixture
def captured_fastk_args(monkeypatch, tmp_path):
    """捕获传给FastK的参数|Capture args passed to FastK"""
    # 临时目录指向pytest tmp,避免写系统/tmp|Point temp dir at pytest tmp
    monkeypatch.setenv("TMPDIR", str(tmp_path))
    holder = {}

    def fake_build_conda_command(command, args, preferred_env=None):
        holder["command"] = command
        holder["args"] = list(args)
        return [command] + list(args)

    monkeypatch.setattr(
        "biopytools.smudgescope.utils.build_conda_command", fake_build_conda_command
    )

    def fake_run(cmd, **kwargs):
        holder["full_cmd"] = list(cmd)
        result = mock.MagicMock()
        result.returncode = 0
        result.stdout = ""
        result.stderr = ""
        return result

    monkeypatch.setattr("biopytools.smudgescope.utils.subprocess.run", fake_run)
    return holder


class TestRunFastkParams:
    """run_fastk参数正确性|run_fastk parameter correctness"""

    def test_t_receives_small_count_threshold_not_threads(self, runner, captured_fastk_args, tmp_path):
        """-t必须是小的计数下限,严禁接收threads值|-t must be a small count threshold, never threads"""
        threads = 88
        ok = runner.run_fastk(
            [str(tmp_path / "sample_R1.fastq")],
            str(tmp_path / "fastk_table"),
            kmer_size=21,
            threads=threads,
            memory="16G",
            sample_name="sample",
        )
        assert ok is True
        args = captured_fastk_args["args"]

        # -T承担线程|-T carries threads
        assert f"-T{threads}" in args
        # -t为官方推荐的固定小阈值,与threads无关|-t is the official small fixed threshold, independent of threads
        assert "-t4" in args
        # 防回归:threads值绝不能出现在-t上|Regression guard: threads value must never land on -t
        assert f"-t{threads}" not in args

    def test_threads_value_never_on_t_for_various_thread_counts(self, runner, captured_fastk_args, tmp_path):
        """不同threads值下-t都不跟随|-t stays fixed across different thread counts"""
        for threads in (4, 12, 64, 88):
            # 每轮换新输出路径避免断点续传跳过|Fresh output path per round to dodge resume skip
            table = tmp_path / f"fastk_table_{threads}"
            ok = runner.run_fastk(
                [str(tmp_path / "sample_R1.fastq")],
                str(table),
                kmer_size=31,
                threads=threads,
                memory="32G",
                sample_name="sample",
            )
            assert ok is True
            args = captured_fastk_args["args"]
            assert f"-T{threads}" in args
            # -t恰有一个且恒为4,不随threads变化|Exactly one -t, always 4, independent of threads
            t_flags = [a for a in args if a.startswith("-t")]
            assert t_flags == ["-t4"]

    def test_command_carries_k_m_N_params(self, runner, captured_fastk_args, tmp_path):
        """-k/-M/-N参数完整传递|-k/-M/-N are passed through"""
        ok = runner.run_fastk(
            [str(tmp_path / "sample_R1.fastq")],
            str(tmp_path / "fastk_table"),
            kmer_size=21,
            threads=12,
            memory="16G",
            sample_name="sample",
        )
        assert ok is True
        args = captured_fastk_args["args"]
        assert "-k21" in args
        assert "-M16" in args
        # -N FastK表前缀|-N FastK table prefix
        assert any(a.startswith("-N") for a in args)
