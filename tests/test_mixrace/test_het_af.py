"""mixrace AF直接判定分支测试(论文式,不依赖GT)|paper-style AF-based het tests.

契约(Cao et al. 2026 口径): 从联合 VCF 只取 AD/DP(不取 GT);depth 达标位点
alt_frac=altAD/(refAD+altAD),altAD>=min 且 frac 落 [min,1-min] → 杂合;
alt 主导 → 纯合变异;het_rate_af = het/(het+hom_alt),分母不含纯合参考位点
(与 GT-based het_rate 的全部 called 位点分母是两种口径,不得混用)。
|Contract: AD/DP only (no GT); paper denominator het/(het+hom_alt).
"""

import numpy as np
import pytest

from biopytools.mixrace.config import MixraceConfig
from biopytools.mixrace.het_eval import HetData, run_het_eval
from biopytools.mixrace.het_eval_af import (
    build_af_query_cmd,
    compute_af_based_het,
    load_ad_dp,
    merge_af_rows,
    read_af_table,
    run_af_het_eval,
)


def _het_data(sites, samples=("S1",)):
    """构造单样本 HetData;sites=(ref_ad, alt_ad, dp[, gt])|build one-sample HetData."""
    n = len(samples)
    ref = np.array([[s[0] for s in sites]] * n, dtype=np.int32)
    alt = np.array([[s[1] for s in sites]] * n, dtype=np.int32)
    dp = np.array([[s[2] for s in sites]] * n, dtype=np.int32)
    gt = np.array([[s[3] if len(s) > 3 else -1 for s in sites]] * n, dtype=np.int8)
    return HetData(samples=list(samples), gt=gt, ref_ad=ref, alt_ad=alt, dp=dp,
                   chrom=["chr1"] * len(sites), pos=np.arange(1, len(sites) + 1, dtype=np.int32))


class TestComputeAfBasedHet:
    """阈值判定逐 case|per-case threshold calling."""

    def test_moderate_alt_frac_is_het(self):
        # 任务示例:alt_frac=0.3, depth=20, alt_ad=6 → 杂合
        d = _het_data([(14, 6, 20)])
        r = compute_af_based_het(d)[0]
        assert r["n_het_af"] == 1 and r["n_hom_alt_af"] == 0
        assert r["het_rate_af"] == pytest.approx(1.0)

    def test_low_alt_frac_is_hom_ref(self):
        # 任务示例:alt_frac≈0.03, depth=33 → 纯合参考,不进分母
        d = _het_data([(32, 1, 33)])
        r = compute_af_based_het(d)[0]
        assert r["n_het_af"] == 0 and r["n_hom_alt_af"] == 0
        assert r["n_sites_eval"] == 1
        assert np.isnan(r["het_rate_af"])   # 分母 0/0 → nan|zero denominator -> nan

    def test_alt_dominant_is_hom_alt(self):
        # alt_frac>1-min → 纯合变异型,进分母
        d = _het_data([(1, 59, 60)])
        r = compute_af_based_het(d)[0]
        assert r["n_hom_alt_af"] == 1 and r["n_het_af"] == 0
        assert r["het_rate_af"] == pytest.approx(0.0)

    def test_low_depth_site_not_evaluated(self):
        d = _het_data([(5, 3, 8), (14, 6, 20)])
        r = compute_af_based_het(d)[0]
        assert r["n_sites_eval"] == 1
        assert r["n_het_af"] == 1

    def test_frac_in_range_but_alt_ad_below_min_is_hom_ref(self):
        # frac=0.25 在区间内但 altAD=2<3 → 信号不足,按主导等位判纯合参考
        d = _het_data([(6, 2, 20)])
        r = compute_af_based_het(d)[0]
        assert r["n_het_af"] == 0 and r["n_hom_alt_af"] == 0

    def test_frac_boundaries_are_inclusive(self):
        # frac 恰=0.05 与恰=0.95(闭区间)且 altAD 达标 → 均判杂合
        d = _het_data([(57, 3, 60), (3, 57, 60)])
        r = compute_af_based_het(d)[0]
        assert r["n_het_af"] == 2

    def test_paper_denominator_excludes_hom_ref(self):
        # 3杂合+1纯合变异+2纯合参考 → 3/(3+1)=0.75,分母不含纯合参考
        d = _het_data([(14, 6, 20), (32, 1, 33), (1, 59, 60),
                       (19, 1, 20), (28, 12, 40), (14, 6, 20)])
        r = compute_af_based_het(d)[0]
        assert r["n_sites_eval"] == 6
        assert r["n_het_af"] == 3 and r["n_hom_alt_af"] == 1
        assert r["het_rate_af"] == pytest.approx(0.75)

    def test_all_hom_ref_yields_nan_rate(self):
        d = _het_data([(32, 1, 33), (19, 1, 20)])
        r = compute_af_based_het(d)[0]
        assert np.isnan(r["het_rate_af"])

    def test_result_is_gt_independent(self):
        # AF 分支不看 GT:gt 全 -1 与 gt 有值,AD/DP 相同 → 结果一致
        sites = [(14, 6, 20), (1, 59, 60), (32, 1, 33)]
        no_gt = compute_af_based_het(_het_data(sites))[0]
        with_gt = compute_af_based_het(_het_data([(a, b, c, 2) for a, b, c in sites]))[0]
        assert no_gt == with_gt

    def test_custom_thresholds_apply(self):
        # min_alt_ad=2 后原 altAD=2 位点升格为杂合
        d = _het_data([(6, 2, 20)])
        r = compute_af_based_het(d, min_alt_ad=2)[0]
        assert r["n_het_af"] == 1

    def test_per_sample_rows(self):
        d = _het_data([(14, 6, 20)], samples=("S1", "S2"))
        rows = compute_af_based_het(d)
        assert [r["sample"] for r in rows] == ["S1", "S2"]
        assert all(r["n_het_af"] == 1 for r in rows)

    def test_row_fields(self):
        r = compute_af_based_het(_het_data([(14, 6, 20)]))[0]
        assert set(r) == {"sample", "n_sites_eval", "n_het_af", "n_hom_alt_af",
                          "het_rate_af"}


class TestLoadAdDp:
    """AF 长表(位点4列+每样本 AD/DP 两列)解析|AF long-table parsing."""

    def test_headerless_two_cols_per_sample(self):
        # bcftools 原生无表头输出;样本名来自 samples 参数
        f = "/tmp/fake.tsv"
        from pathlib import Path as _P
        import tempfile
        with tempfile.NamedTemporaryFile("w", suffix=".tsv", delete=False) as fh:
            fh.write("chr1\t100\tA\tG\t14,6\t20\t32,1\t33\n"
                     "chr1\t200\tA\tT\t1,59\t60\t14,6\t20\n")
            f = fh.name
        d = load_ad_dp(f, samples=["S1", "S2"])
        assert d.samples == ["S1", "S2"]
        assert len(d.pos) == 2
        assert d.ref_ad[0].tolist() == [14, 1] and d.alt_ad[0].tolist() == [6, 59]
        assert d.dp[1].tolist() == [33, 20]
        assert (d.gt == -1).all()   # AF 长表无 GT|no GT column

    def test_ad_header_line_accepted(self, tmp_path):
        f = tmp_path / "af.tsv"
        f.write_text("CHROM\tPOS\tREF\tALT\tS1_AD\tS1_DP\n"
                     "chr1\t100\tA\tG\t14,6\t20\n", encoding="utf-8")
        d = load_ad_dp(str(f))
        assert d.samples == ["S1"]
        assert d.alt_ad[0].tolist() == [6]

    def test_multiallelic_and_indel_skipped(self, tmp_path):
        f = tmp_path / "af.tsv"
        f.write_text("chr1\t100\tA\tG,T\t14,6\t20\n"     # 多等位|multiallelic
                     "chr1\t200\tA\tAT\t14,6\t20\n"      # indel
                     "chr1\t300\tA\tG\t14,6\t20\n", encoding="utf-8")
        d = load_ad_dp(str(f))
        assert len(d.pos) == 1 and d.pos[0] == 300

    def test_bad_ad_values_tolerated(self, tmp_path):
        f = tmp_path / "af.tsv"
        f.write_text("chr1\t100\tA\tG\t.,.\t20\n", encoding="utf-8")
        d = load_ad_dp(str(f), samples=["S1"])
        assert d.ref_ad[0].tolist() == [0] and d.alt_ad[0].tolist() == [0]


# ---------- 编排层(命令/断点/整合)|orchestration ----------

# 两样本三位点:GT 口径与 AF 口径刻意不同(位点1 S1 GT=0/1 杂合,S2 GT=1/1;
# AF 口径下 S1=3杂合参考边界混合体,S2 位点1 是 0.98 纯合变异)
# |two samples x three sites; GT and AF calibres intentionally diverge.
_GT_TSV = (
    "chr1\t100\tA\tG\t0/1\t14,6\t20\t1/1\t1,59\t60\n"
    "chr1\t200\tA\tT\t0/0\t32,1\t33\t0/1\t28,12\t40\n"
    "chr1\t300\tA\tC\t0/0\t19,1\t20\t0/0\t19,1\t20\n"
)
# S1: frac 0.3→杂合, 0.03→纯合参考, ~0.05/altAD=1→纯合参考 ⇒ het=1, hom_alt=0
# S2: frac 0.983→纯合变异, 0.3→杂合, 0.05/altAD=1→纯合参考 ⇒ het=1, hom_alt=1
_AF_TSV = (
    "chr1\t100\tA\tG\t14,6\t20\t1,59\t60\n"
    "chr1\t200\tA\tT\t32,1\t33\t28,12\t40\n"
    "chr1\t300\tA\tC\t19,1\t20\t19,1\t20\n"
)


class _QuietLogger:
    def info(self, *a, **kw):
        pass

    def warning(self, *a, **kw):
        pass

    def error(self, *a, **kw):
        pass


class _AfRunner:
    """假执行器:query -l 回样本表,query -o 按目标文件名写假长表|fake runner."""

    def __init__(self, samples=("S1", "S2")):
        self.logger = _QuietLogger()
        self.dry_run = False
        self.samples = list(samples)
        self.query_calls = []

    def run(self, cmd, description="", timeout=None):
        args = [str(x) for x in cmd]
        if "-l" in args:
            return True, "\n".join(self.samples) + "\n", ""
        if "-o" in args:
            out = __import__("pathlib").Path(args[args.index("-o") + 1])
            out.parent.mkdir(parents=True, exist_ok=True)
            out.write_text(_AF_TSV if out.name == "af_ad_dp.tsv" else _GT_TSV,
                           encoding="utf-8")
            self.query_calls.append(out.name)
            return True, "", ""
        return True, "", ""

    def run_conda(self, tool_path, args, description="", timeout=None):
        return self.run(["bcftools"] + list(args), description, timeout)


class _AfCkpt:
    def __init__(self, done=()):
        self.done = set(done)

    def exists(self, step):
        return step in self.done

    def create(self, step):
        self.done.add(step)


def _af_config(tmp_path, af=True, min_sites=1):
    # min_sites=1:夹具长表只有3个位点,默认1000小分母守卫会清空全部het_rate_af
    # |tiny fixture: keep the small-denominator guard from wiping rates
    fa = tmp_path / "ref.fa"
    fa.write_text(">chr1\nACGT\n", encoding="utf-8")
    return MixraceConfig(fastq_dir=str(tmp_path), genome=str(fa),
                         output_dir=str(tmp_path / "out"), af_het_eval=af,
                         min_sites=min_sites,
                         run_kraken2=False)  # 污染评估与AF测试无关,跳过DB校验|kraken2 unrelated here


def _verdict_path(cfg):
    from pathlib import Path as _P
    return _P(cfg.output_dir) / "04_het_eval" / "verdict_table.tsv"


class TestBuildAfQueryCmd:
    def test_format_has_ad_dp_no_gt(self):
        cmd = build_af_query_cmd("~/miniforge3/envs/align/bin/bcftools", "in.vcf.gz",
                                 "out.tsv")
        s = " ".join(str(x) for x in cmd)
        assert "%AD" in s and "%DP" in s and "-o" in s
        assert "%GT" not in s   # 刻意不取 GT|never query GT


class TestRunAfHetEval:
    def test_writes_l1_af_table(self, tmp_path):
        cfg = _af_config(tmp_path)
        rows = run_af_het_eval(cfg, _AfRunner(), _AfCkpt(), "in.vcf.gz")
        l1 = tmp_path / "out" / "04_het_eval" / "l1_het_af_based.tsv"
        assert l1.exists()
        by_sample = {r["sample"]: r for r in rows}
        assert by_sample["S1"]["n_het_af"] == 1 and by_sample["S1"]["n_hom_alt_af"] == 0
        assert by_sample["S1"]["het_rate_af"] == pytest.approx(1.0)
        assert by_sample["S2"]["n_het_af"] == 1 and by_sample["S2"]["n_hom_alt_af"] == 1
        assert by_sample["S2"]["het_rate_af"] == pytest.approx(0.5)

    def test_done_checkpoint_skips_query(self, tmp_path):
        cfg = _af_config(tmp_path)
        out = tmp_path / "out" / "04_het_eval"
        out.mkdir(parents=True)
        (out / "l1_het_af_based.tsv").write_text(
            "sample\tn_sites_eval\tn_het_af\tn_hom_alt_af\thet_rate_af\n"
            "S1\t3\t1\t0\t1.0\n", encoding="utf-8")
        runner = _AfRunner()
        rows = run_af_het_eval(cfg, runner, _AfCkpt(done={"af_het_eval"}), "in.vcf.gz")
        assert runner.query_calls == []   # 全部复用,不重跑|nothing re-runs
        assert rows[0]["sample"] == "S1"
        assert rows[0]["het_rate_af"] == pytest.approx(1.0)

    def test_query_checkpoint_reuses_long_table(self, tmp_path):
        cfg = _af_config(tmp_path)
        out = tmp_path / "out" / "04_het_eval"
        out.mkdir(parents=True)
        (out / "af_ad_dp.tsv").write_text(_AF_TSV, encoding="utf-8")
        runner = _AfRunner()
        rows = run_af_het_eval(cfg, runner, _AfCkpt(done={"af_query"}), "in.vcf.gz")
        assert runner.query_calls == []   # 长表复用|table reused
        assert len(rows) == 2 and rows[1]["n_hom_alt_af"] == 1

    def test_query_failure_returns_empty(self, tmp_path):
        cfg = _af_config(tmp_path)

        class _FailRunner(_AfRunner):
            def run(self, cmd, description="", timeout=None):
                args = [str(x) for x in cmd]
                if "-o" in args:
                    return False, "", "boom"
                return True, "S1\n", ""

        assert run_af_het_eval(cfg, _FailRunner(), _AfCkpt(), "in.vcf.gz") == []

    def test_small_denominator_omits_rate(self, tmp_path):
        # 小分母守卫:het+hom_alt < min_sites 时比率是噪声(纯菌分母可低至几十),
        # het_rate_af 置空防误导;计数照常保留——与 judge() 的 min_sites 降级同语义
        # |small-denominator guard: rate omitted as noise; counts kept
        cfg = _af_config(tmp_path, min_sites=1000)
        rows = run_af_het_eval(cfg, _AfRunner(), _AfCkpt(), "in.vcf.gz")
        l1 = (tmp_path / "out" / "04_het_eval" / "l1_het_af_based.tsv").read_text(
            encoding="utf-8").splitlines()
        # S1 分母=1+0、S2 分母=1+1,均 <1000 → TSV 该列为空|both <1000 -> blank
        assert l1[1].endswith("\t") and l1[2].endswith("\t")
        by_sample = {r["sample"]: r for r in rows}
        assert by_sample["S1"]["n_het_af"] == 1            # 计数保留|counts kept
        assert by_sample["S1"]["het_rate_af"] != by_sample["S1"]["het_rate_af"]  # nan


class TestReadAfTable:
    def test_roundtrip_types(self, tmp_path):
        f = tmp_path / "l1_het_af_based.tsv"
        f.write_text("sample\tn_sites_eval\tn_het_af\tn_hom_alt_af\thet_rate_af\n"
                     "S1\t10\t4\t2\t0.6666666666666666\n"
                     "S2\t0\t0\t0\t\n", encoding="utf-8")
        rows = read_af_table(str(f))
        assert rows[0]["n_het_af"] == 4 and rows[0]["n_hom_alt_af"] == 2
        assert rows[0]["het_rate_af"] == pytest.approx(2 / 3)
        assert rows[1]["het_rate_af"] is None   # 空→None|empty -> None

    def test_missing_file_returns_empty(self, tmp_path):
        assert read_af_table(str(tmp_path / "nope.tsv")) == []


class TestMergeAfRows:
    def test_adds_column_by_sample(self):
        verdict_rows = [{"sample": "S1", "verdict": "pure", "het_rate": 0.0002},
                        {"sample": "S2", "verdict": "contaminated", "het_rate": 0.5}]
        af_rows = [{"sample": "S2", "n_sites_eval": 3, "n_het_af": 1,
                    "n_hom_alt_af": 1, "het_rate_af": 0.5}]
        merged = merge_af_rows(verdict_rows, af_rows)
        assert "het_rate_af" not in merged[0]   # AF 无数据的样本不加键|no key added
        assert merged[1]["het_rate_af"] == pytest.approx(0.5)
        assert merged[1]["verdict"] == "contaminated"   # 只加不改|add-only


class TestRunHetEvalIntegration:
    """run_het_eval 整合:默认不加列;--af-het-eval 加列;旧输出断点补跑整合。
    |integration: no column by default; added when enabled; old-output resume.
    """

    def test_with_af_verdict_table_gains_column(self, tmp_path):
        cfg = _af_config(tmp_path, af=True)
        rows = run_het_eval(cfg, _AfRunner(), _AfCkpt(), "in.vcf.gz")
        lines = _verdict_path(cfg).read_text(encoding="utf-8").splitlines()
        head = lines[0].split("\t")
        assert "het_rate_af" in head and "het_rate" in head
        by_sample = {r["sample"]: r for r in rows}
        assert by_sample["S1"]["het_rate_af"] == pytest.approx(1.0)
        assert by_sample["S2"]["het_rate_af"] == pytest.approx(0.5)
        # 判读仍走 GT 口径且不受 AF 影响:与不开 AF 跑出的 verdict 逐一一致
        # |verdicts stay GT-based: identical to a run without the AF branch
        noaf_dir = tmp_path.parent / "noaf_run"
        noaf_dir.mkdir(exist_ok=True)
        cfg_noaf = _af_config(noaf_dir, af=False)
        rows_noaf = run_het_eval(cfg_noaf, _AfRunner(), _AfCkpt(), "in.vcf.gz")
        assert ([r["verdict"] for r in rows] == [r["verdict"] for r in rows_noaf]
                and all(r["verdict"] for r in rows))

    def test_without_af_no_new_column(self, tmp_path):
        cfg = _af_config(tmp_path, af=False)
        rows = run_het_eval(cfg, _AfRunner(), _AfCkpt(), "in.vcf.gz")
        head = _verdict_path(cfg).read_text(encoding="utf-8").splitlines()[0].split("\t")
        assert "het_rate_af" not in head   # 精确列名(防子串误命中het_rate_after_hotspot)
        assert all("het_rate_af" not in r for r in rows)

    def test_resume_old_output_integrates_af_column(self, tmp_path):
        # 旧输出(无 af 列)+ het_eval 断点已建 → 只补跑 AF 并把列整合进已有表
        # |old output + het_eval checkpoint: AF-only rerun merges into existing table
        cfg = _af_config(tmp_path, af=True)
        out = tmp_path / "out" / "04_het_eval"
        out.mkdir(parents=True)
        (out / "verdict_table.tsv").write_text(
            "sample\tverdict\thet_rate\trationale\n"
            "S1\tuncertain\t0.05\told\n"
            "S2\tuncertain\t0.01\told\n", encoding="utf-8")
        runner = _AfRunner()
        rows = run_het_eval(cfg, runner, _AfCkpt(done={"het_eval"}), "in.vcf.gz")
        # GT 长表不重跑,只跑 AF query|GT query not re-run, AF query only
        assert runner.query_calls == ["af_ad_dp.tsv"]
        by_sample = {r["sample"]: r for r in rows}
        assert by_sample["S1"]["het_rate_af"] == pytest.approx(1.0)
        assert by_sample["S2"]["het_rate_af"] == pytest.approx(0.5)
        # 表被重写:既有列保留+新列追加|table rewritten: old cols kept + new appended
        lines = _verdict_path(cfg).read_text(encoding="utf-8").splitlines()
        head = lines[0].split("\t")
        assert head[:3] == ["sample", "verdict", "het_rate"]
        assert "het_rate_af" in head and "rationale" in head
        assert lines[1].split("\t")[0] == "S1"

    def test_full_checkpoint_skip_keeps_af_column(self, tmp_path):
        # 二次重跑:het_eval+af_het_eval 断点齐全 → 读表返回,列已含 het_rate_af
        # |second rerun: all checkpoints -> read table, af column preserved
        cfg = _af_config(tmp_path, af=True)
        rows = run_het_eval(cfg, _AfRunner(), _AfCkpt(), "in.vcf.gz")
        ckpt = _AfCkpt(done={"het_eval", "af_het_eval"})
        runner = _AfRunner()
        rows2 = run_het_eval(cfg, runner, ckpt, "in.vcf.gz")
        assert runner.query_calls == []
        assert {r["sample"]: r["het_rate_af"] for r in rows2} == pytest.approx(
            {"S1": 1.0, "S2": 0.5})


# ---------- config 校验|config validation ----------

class TestAfConfigValidate:
    def _cfg(self, tmp_path, **kw):
        base = dict(fastq_dir=str(tmp_path), genome=str(tmp_path / "ref.fa"),
                    output_dir=str(tmp_path / "out"), run_kraken2=False)  # 跳过kraken2 DB校验|skip kraken2 db check
        base.update(kw)
        return MixraceConfig(**base)

    def test_defaults_on(self, tmp_path):
        (tmp_path / "ref.fa").write_text(">c\nACGT\n", encoding="utf-8")
        cfg = self._cfg(tmp_path)
        assert cfg.af_het_eval is True   # 默认开启|on by default
        assert (cfg.af_het_min_frac, cfg.af_het_min_depth, cfg.af_het_min_alt_ad) \
            == (0.05, 10, 3)

    @pytest.mark.parametrize("bad", [
        {"af_het_min_frac": 0.0}, {"af_het_min_frac": 0.5}, {"af_het_min_frac": 0.9},
        {"af_het_min_depth": 0}, {"af_het_min_alt_ad": 0}])
    def test_out_of_range_raises(self, tmp_path, bad):
        (tmp_path / "ref.fa").write_text(">c\nACGT\n", encoding="utf-8")
        with pytest.raises(ValueError, match="af_het"):
            self._cfg(tmp_path, **bad).validate()

    def test_valid_values_pass(self, tmp_path):
        (tmp_path / "ref.fa").write_text(">c\nACGT\n", encoding="utf-8")
        self._cfg(tmp_path, af_het_eval=True, af_het_min_frac=0.02,
                  af_het_min_depth=20, af_het_min_alt_ad=5).validate()


# ---------- CLI 透传|CLI passthrough ----------

class TestArgvToConfig:
    def _argv(self, tmp_path, extra):
        import sys as _sys
        fa = tmp_path / "ref.fa"
        if not fa.exists():
            fa.write_text(">c\nACGT\n", encoding="utf-8")
        return ["mixrace", "-g", str(fa), "-i", str(tmp_path),
                "-o", str(tmp_path / "out")] + extra

    def test_default_on(self, tmp_path, monkeypatch):
        import sys as _sys
        from biopytools.mixrace.main import _argv_to_config
        monkeypatch.setattr(_sys, "argv", self._argv(tmp_path, []))
        cfg = _argv_to_config()
        assert cfg.af_het_eval is True

    def test_no_af_flag_disables(self, tmp_path, monkeypatch):
        import sys as _sys
        from biopytools.mixrace.main import _argv_to_config
        monkeypatch.setattr(_sys, "argv", self._argv(tmp_path, ["--no-af-het-eval"]))
        cfg = _argv_to_config()
        assert cfg.af_het_eval is False

    def test_thresholds_passthrough(self, tmp_path, monkeypatch):
        import sys as _sys
        from biopytools.mixrace.main import _argv_to_config
        monkeypatch.setattr(_sys, "argv", self._argv(tmp_path, [
            "--af-het-min-frac", "0.02",
            "--af-het-min-depth", "20", "--af-het-min-alt-ad", "5"]))
        cfg = _argv_to_config()
        assert cfg.af_het_eval is True   # 阈值单独给,不影响开关|thresholds don't toggle
        assert cfg.af_het_min_frac == pytest.approx(0.02)
        assert cfg.af_het_min_depth == 20 and cfg.af_het_min_alt_ad == 5


# ---------- 汇总报告对照列|reporter side-by-side column ----------

def _verdict_row(**kw):
    base = {"sample": "S1", "verdict": "contaminated", "subtag": "",
            "advice": "建议: 需再分离纯化|Advice: re-isolate",
            "rationale": "杂合率高|high het",
            "het_rate": 0.05, "robust_rate": 0.04,
            "shared_only_rate": 0.05, "top_partner": "S2",
            "mix_proportion": None, "het_rate_after_hotspot": 0.06,
            "dp_ratio": 1.8, "host_rate": None, "pathogen_map_rate": 0.9,
            "contamination_rate": 0.02, "mean_depth": 30.0, "breadth_1x": 0.97}
    base.update(kw)
    return base


class TestReporterAfColumn:
    def test_tsv_inserts_af_column_after_het_rate(self):
        from biopytools.mixrace.reporter import build_summary_table
        tsv, _ = build_summary_table([_verdict_row(het_rate_af=0.42)])
        head = tsv.splitlines()[0].split("\t")
        assert head[head.index("het_rate") + 1] == "het_rate_af"   # 并排方便对比
        assert "\t42.00%\t" in tsv.splitlines()[1]

    def test_tsv_without_af_keeps_canonical_columns(self):
        from biopytools.mixrace.reporter import build_summary_table
        tsv, _ = build_summary_table([_verdict_row()])
        head = tsv.splitlines()[0].split("\t")
        assert "het_rate_af" not in head
        assert head.index("robust_rate") == head.index("het_rate") + 1

    def test_html_shows_cn_label_and_raw_value(self):
        from biopytools.mixrace.reporter import build_summary_table, build_sample_report
        _, html = build_summary_table([_verdict_row(het_rate_af=0.42)])
        assert "AF口径杂合率" in html and 'data-v="0.42"' in html
        md = build_sample_report("S1", _verdict_row(het_rate_af=0.42), {})
        assert "het_rate_af" in md and "42.00%" in md

    def test_none_af_value_renders_dash(self):
        from biopytools.mixrace.reporter import build_summary_table
        tsv, _ = build_summary_table([_verdict_row(het_rate_af=None)])
        assert "\t—\t" in tsv.splitlines()[1]

    def test_nan_af_value_renders_dash_not_nan_percent(self, tmp_path):
        # 单次全流程:内存行带 nan(小分母守卫/0分母),不得渲染成 "nan%" 或
        # data-v="nan"(JS parseFloat NaN 破坏排序)|in-memory nan -> dash, not "nan%"
        import math
        from biopytools.mixrace.reporter import build_summary_table
        tsv, html = build_summary_table([_verdict_row(het_rate_af=float("nan"))])
        assert "nan" not in tsv and "\t—\t" in tsv.splitlines()[1]
        assert 'data-v="nan"' not in html and "nan%" not in html
