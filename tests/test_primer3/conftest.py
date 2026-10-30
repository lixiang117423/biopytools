"""确保 biopytools 可导入|Ensure biopytools importable from repo root"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[2]))

import pytest


class FakeCommandRunner:
    """CommandRunner替身: 记录调用, 按chunk输入生成可校验的输出|Stub runner

    按 chunk 输入文件中的 SEQUENCE_ID 生成输出, 使合并顺序可断言;
    fail_id 命中时该 chunk 返回失败(模拟 primer3_core 出错)
    |Emits output keyed by the SEQUENCE_IDs of each chunk input so merge
    order is assertable; fail_id makes that chunk fail like a real crash
    """

    instances = None  # 类属性: 收集实例供测试断言|Collect instances for asserts

    def __init__(self, logger):
        self.logger = logger
        self.calls = []
        self.fail_id = None
        if FakeCommandRunner.instances is not None:
            FakeCommandRunner.instances.append(self)

    def run(self, cmd, description="", output_file=None, **kwargs):
        if output_file is None:
            # 版本探测调用|Version probing call
            self.calls.append({"cmd": list(cmd), "output_file": None, "seq_ids": None})
            return True, "primer3 2.6.1", ""

        # 调用时立即捕获输入ID: 生产代码随后会删临时输入文件(§12.4)
        # |Capture IDs at call time: production deletes the temp input (§12.4)
        input_path = Path(cmd[-1])
        seq_ids = [
            line.split("=", 1)[1]
            for line in input_path.read_text().splitlines()
            if line.startswith("SEQUENCE_ID=")
        ]
        self.calls.append({"cmd": list(cmd), "output_file": output_file, "seq_ids": seq_ids})

        if self.fail_id and self.fail_id in seq_ids:
            return False, "", f"mock failure: {self.fail_id}"

        Path(output_file).write_text(
            "".join(f"SEQUENCE_ID={i}\n=\n" for i in seq_ids)
        )
        return True, "", ""


@pytest.fixture
def fake_runner_cls(monkeypatch):
    """替换 evaluator 内的 CommandRunner 并收集实例|Patch CommandRunner, collect instances"""
    FakeCommandRunner.instances = []
    monkeypatch.setattr(
        "biopytools.primer3.evaluator.CommandRunner", FakeCommandRunner
    )
    return FakeCommandRunner


@pytest.fixture
def workspace(tmp_path):
    """临时输入/输出环境|Temp input/output environment"""
    fasta = tmp_path / "input.fa"
    fasta.write_text(">seq_0000\nACGT\n")
    core = tmp_path / "primer3_core_bin"
    core.write_text("")
    return {
        "fasta": str(fasta),
        "core": str(core),
        "out": str(tmp_path / "out"),
    }


@pytest.fixture
def build_evaluator(workspace, fake_runner_cls):
    """评估器工厂(自动替换CommandRunner)|Evaluator factory (CommandRunner patched)"""
    from biopytools.primer3.evaluator import Primer3Evaluator

    def _build(**overrides):
        params = dict(
            input_fasta=workspace["fasta"],
            output_dir=workspace["out"],
            primer3_core_path=workspace["core"],
        )
        params.update(overrides)
        return Primer3Evaluator(**params)

    return _build


@pytest.fixture
def make_sequences():
    """有序测试序列工厂|Ordered test sequence factory"""
    def _make(n_seqs: int):
        return [(f"seq_{i:04d}", "ACGT" * 15) for i in range(n_seqs)]

    return _make
