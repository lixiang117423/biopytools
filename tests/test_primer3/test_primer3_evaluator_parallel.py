"""primer3 评估器并行行为测试|Tests for primer3 evaluator parallel behavior"""

from pathlib import Path

import pytest


def merged_marker(evaluator):
    """断点续传 done-marker 路径|Path of the resume done-marker"""
    return Path(evaluator.primer3_output_file)


def ids_from_input(call):
    """获取调用时捕获的SEQUENCE_ID列表|Call-time captured SEQUENCE_IDs"""
    return call["seq_ids"]


def merged_ids(evaluator):
    """从合并产物提取SEQUENCE_ID顺序|Merged marker SEQUENCE_ID order"""
    return [
        line.split("=", 1)[1]
        for line in merged_marker(evaluator).read_text().splitlines()
        if line.startswith("SEQUENCE_ID=")
    ]


class TestParallelDispatch:
    """阈值分流: 单进程 vs 并行|Dispatch: single vs parallel by threshold"""

    def test_below_threshold_runs_single_invocation(
        self, build_evaluator, make_sequences, fake_runner_cls
    ):
        """低于阈值只调一次primer3且直写marker|Below threshold: one call, straight to marker"""
        evaluator = build_evaluator(threads=4, parallel_threshold=10)
        evaluator._run_primer3(make_sequences(6))

        assert len(fake_runner_cls.instances[0].calls) == 1
        call = fake_runner_cls.instances[0].calls[0]
        assert call["output_file"] == str(merged_marker(evaluator))
        assert ids_from_input(call) == [f"seq_{i:04d}" for i in range(6)]

    def test_at_threshold_splits_into_thread_chunks(
        self, build_evaluator, make_sequences, fake_runner_cls
    ):
        """达到阈值按线程数切分并行|At threshold: split by threads"""
        evaluator = build_evaluator(threads=4, parallel_threshold=10)
        evaluator._run_primer3(make_sequences(10))

        assert len(fake_runner_cls.instances[0].calls) == 4
        # 合并产物保持全局顺序|Merged marker keeps global order
        assert merged_ids(evaluator) == [f"seq_{i:04d}" for i in range(10)]

    def test_chunks_capped_by_sequence_count(
        self, build_evaluator, make_sequences, fake_runner_cls
    ):
        """进程数不超过序列数|Process count capped by sequence count"""
        evaluator = build_evaluator(threads=12, parallel_threshold=1)
        evaluator._run_primer3(make_sequences(2))

        assert len(fake_runner_cls.instances[0].calls) == 2

    def test_threads_one_forces_single_process(
        self, build_evaluator, make_sequences, fake_runner_cls
    ):
        """threads=1时即使超阈值也单进程|threads=1 forces single process"""
        evaluator = build_evaluator(threads=1, parallel_threshold=1)
        evaluator._run_primer3(make_sequences(20))

        assert len(fake_runner_cls.instances[0].calls) == 1
        assert merged_ids(evaluator) == [f"seq_{i:04d}" for i in range(20)]

    def test_partition_covers_all_sequences_once(
        self, build_evaluator, make_sequences, fake_runner_cls
    ):
        """切分不重不漏且各chunk内部保序|Partition exact, per-chunk order kept"""
        evaluator = build_evaluator(threads=3, parallel_threshold=5)
        sequences = make_sequences(7)
        evaluator._run_primer3(sequences)

        calls = fake_runner_cls.instances[0].calls
        assert len(calls) == 3

        chunked = [ids_from_input(call) for call in calls]
        all_ids = [seq_id for chunk in chunked for seq_id in chunk]
        expected = [seq_id for seq_id, _ in sequences]
        assert sorted(all_ids) == sorted(expected)
        assert len(set(all_ids)) == len(all_ids)

        # 各chunk内序列相对顺序与全局一致|Within-chunk relative order matches global
        positions = {seq_id: i for i, seq_id in enumerate(expected)}
        for chunk in chunked:
            chunk_pos = [positions[seq_id] for seq_id in chunk]
            assert chunk_pos == sorted(chunk_pos)


class TestParallelFailureAndResume:
    """失败处理与断点续传|Failure handling and resume"""

    def test_chunk_failure_raises_and_no_marker(
        self, build_evaluator, make_sequences, fake_runner_cls, tmp_path
    ):
        """任一chunk失败则抛错且不落marker|Any chunk failure raises, marker absent"""
        evaluator = build_evaluator(threads=3, parallel_threshold=5)
        fake_runner_cls.instances[0].fail_id = "seq_0002"

        with pytest.raises(RuntimeError, match="Primer3运行失败|execution failed"):
            evaluator._run_primer3(make_sequences(6))

        assert not merged_marker(evaluator).exists()

    def test_chunk_failure_cleans_tmp(
        self, build_evaluator, make_sequences, fake_runner_cls
    ):
        """失败后tmp目录无残留|No leftover files in tmp after failure"""
        evaluator = build_evaluator(threads=3, parallel_threshold=5)
        fake_runner_cls.instances[0].fail_id = "seq_0002"

        with pytest.raises(RuntimeError):
            evaluator._run_primer3(make_sequences(6))

        leftovers = [
            p.name
            for p in evaluator.tmp_dir.iterdir()
            if p.is_file()
        ]
        assert leftovers == []

    def test_completed_marker_skips_rerun(
        self, build_evaluator, make_sequences, fake_runner_cls
    ):
        """已有marker时跳过重跑|Existing marker skips rerun"""
        evaluator = build_evaluator(threads=4, parallel_threshold=5)
        merged_marker(evaluator).write_text("previous run\n")

        evaluator._run_primer3(make_sequences(6))

        assert fake_runner_cls.instances[0].calls == []
        assert "previous run" in merged_marker(evaluator).read_text()
