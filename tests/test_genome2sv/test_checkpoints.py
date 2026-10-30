"""genome2sv 断点续传单测|Checkpoint-resume tests for genome2sv pipeline

merge 断点(merged VCF 存在即跳过 SURVIVOR)与下游输出断点(三产物齐即跳过)。
|Merge checkpoint (skip SURVIVOR when merged VCF exists) and downstream-output
checkpoint (skip when all three artifacts exist).
"""
import os
from types import SimpleNamespace

from biopytools.genome2sv.pipeline import Genome2SVPipeline


def _make_config(tmp_path):
    """构造最小 config 桩(仅下游步骤用到的字段)|Minimal config stub."""
    return SimpleNamespace(
        merged_dir=tmp_path / "03_merged",
        stats_dir=tmp_path / "04_stats",
        sv_seq_dir=tmp_path / "05_sv_sequences",
        reference_fasta=str(tmp_path / "reference" / "ref.fasta"),
        max_dist=1000, min_support=1, survivor_type=1, survivor_strand=1,
        est_dist=1, min_sv_length=50,
        survivor_path="SURVIVOR", bcftools_path="bcftools",
    )


def _make_pipeline(tmp_path):
    """构造 pipeline + 调用记录器|Build pipeline with call recorder."""
    calls = []
    logger = SimpleNamespace(
        info=lambda msg: calls.append(("info", msg)),
        warning=lambda msg: calls.append(("warning", msg)),
        error=lambda msg: calls.append(("error", msg)),
    )
    for d in ("03_merged", "04_stats", "05_sv_sequences"):
        (tmp_path / d).mkdir(parents=True, exist_ok=True)
    return Genome2SVPipeline(_make_config(tmp_path), logger), calls


def _fake_run_recorder(pipeline):
    """替换 _run 为计数桩|Replace _run with counting stub."""
    state = {"count": 0}

    def fake_run(cmd, description=""):
        state["count"] += 1
        return True

    pipeline._run = fake_run
    return state


def test_merge_skipped_when_merged_vcf_exists(tmp_path):
    """merged VCF 已存在时不调 SURVIVOR|No SURVIVOR call when VCF exists."""
    pipeline, calls = _make_pipeline(tmp_path)
    state = _fake_run_recorder(pipeline)
    merged = tmp_path / "03_merged" / "pan_sv_survivor.vcf"
    merged.write_text("##fileformat=VCFv4.2\n")
    ok = pipeline.merge_and_stats({"s1": str(tmp_path / "s1.vcf")})
    assert ok is True
    assert state["count"] == 0
    assert any("跳过" in msg for _, msg in calls)


def test_merge_runs_when_vcf_missing(tmp_path):
    """merged VCF 缺失时执行 SURVIVOR|SURVIVOR runs when VCF missing."""
    pipeline, calls = _make_pipeline(tmp_path)
    state = _fake_run_recorder(pipeline)
    (tmp_path / "s1.vcf").write_text("x\n")
    pipeline.merge_and_stats({"s1": str(tmp_path / "s1.vcf")})
    assert state["count"] == 1


def test_downstream_skipped_when_all_outputs_exist(tmp_path):
    """三产物齐备时跳过下游步骤|Skip downstream when all outputs exist."""
    pipeline, calls = _make_pipeline(tmp_path)
    (tmp_path / "05_sv_sequences" / "pan_sv_sequences.fa").write_text(">x\nACGT\n")
    (tmp_path / "04_stats" / "pav_matrix.tsv").write_text("h\n")
    (tmp_path / "04_stats" / "pav_binary.tsv").write_text("h\n")
    ok = pipeline.generate_downstream_outputs(
        str(tmp_path / "03_merged" / "pan_sv_survivor.vcf"))
    assert ok is True
    assert any("跳过" in msg for _, msg in calls)
    # 已存在的产物不被改写|mtime 不变近似:内容仍为桩
    assert (tmp_path / "04_stats" / "pav_matrix.tsv").read_text() == "h\n"


def test_downstream_regenerates_when_partial(tmp_path):
    """产物不全时重新生成全部|Regenerate all when partial."""
    pipeline, calls = _make_pipeline(tmp_path)
    merged = tmp_path / "03_merged" / "pan_sv_survivor.vcf"
    merged.write_text(
        "##fileformat=VCFv4.2\n"
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\thap1\thap2\n"
        "chrA\t11\tsvim_asm.INS.1\tA\tAGGGTTTAAA\t.\tPASS\t"
        "SVTYPE=INS;END=11;SVLEN=9\tGT\t1/1:.\t./.:.\n")
    (tmp_path / "04_stats" / "pav_matrix.tsv").write_text("stale\n")
    ok = pipeline.generate_downstream_outputs(str(merged))
    assert ok is True
    pav = (tmp_path / "04_stats" / "pav_matrix.tsv").read_text()
    assert "pan_sv.INS.00001" in pav and "\t1\t0" in pav
    fa = (tmp_path / "05_sv_sequences" / "pan_sv_sequences.fa").read_text()
    assert ">pan_sv.INS.00001" in fa and "GGGTTTAAA" in fa


def test_downstream_missing_vcf_warns_not_fails(tmp_path):
    """merged VCF 缺失:WARNING 不中断|Missing merged VCF: warn, not fail."""
    pipeline, calls = _make_pipeline(tmp_path)
    ok = pipeline.generate_downstream_outputs(
        str(tmp_path / "03_merged" / "pan_sv_survivor.vcf"))
    assert ok is True
    assert any("warning" == level for level, _ in calls)
    assert not (tmp_path / "05_sv_sequences" / "pan_sv_sequences.fa").exists()
