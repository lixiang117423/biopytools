"""mixrace 污染评估测试(step 6 kraken2+bracken)|contamination assessment tests.

契约:逐样本 kraken2(仅 --report,--output 指向 /dev/null 防逐read明细爆stdout)、
bracken(-r 吸附到库内可用 kmer_distrib 读长)、独立汇总表(不动 verdict_summary)。
|Contract: per-sample kraken2 (--report only, --output /dev/null), bracken
(-r snapped to available kmer_distrib), standalone summary tables.
"""

from pathlib import Path

import pytest

from biopytools.mixrace.config import MixraceConfig
from biopytools.mixrace.contamination import (
    KMER_DISTRIB_LENGTHS,
    parse_kraken2_report,
    parse_bracken_report,
    render_summary_tsv,
    render_detail_tsv,
    snap_read_length,
    run_contamination,
)


# ---------- 吸附读长|snap read length ----------

class TestSnapReadLength:
    def test_exact_match_unchanged(self):
        assert snap_read_length(150) == 150

    def test_snaps_to_nearest_below(self):
        # |120-100|=20 < |120-150|=30
        assert snap_read_length(120) == 100

    def test_snaps_to_nearest_above(self):
        assert snap_read_length(180) == 200

    def test_below_min_clamps_to_min(self):
        assert snap_read_length(30) == 50

    def test_above_max_clamps_to_max(self):
        assert snap_read_length(400) == 300


# ---------- kraken2 report 解析|k2 report parsing ----------

_K2_REPORT = (
    "  5.00\t50\t50\tU\t0\tunclassified\n"
    " 95.00\t950\t0\tR\t1\troot\n"
    " 90.00\t900\t10\tK\t2\tBacteria\n"
)

class TestParseKraken2Report:
    def test_totals(self, tmp_path):
        f = tmp_path / "S1.k2.report.txt"
        f.write_text(_K2_REPORT, encoding="utf-8")
        d = parse_kraken2_report(str(f))
        assert d["classified_reads"] == 950
        assert d["unclassified_reads"] == 50
        assert d["total_reads"] == 1000

    def test_missing_unclassified_row_means_zero(self, tmp_path):
        f = tmp_path / "S1.k2.report.txt"
        f.write_text("100.00\t100\t0\tR\t1\troot\n", encoding="utf-8")
        d = parse_kraken2_report(str(f))
        assert d["unclassified_reads"] == 0
        assert d["total_reads"] == 100

    def test_missing_file_returns_empty(self, tmp_path):
        assert parse_kraken2_report(str(tmp_path / "nope.txt")) == {}


# ---------- bracken 解析|bracken parsing ----------

_BRACKEN = (
    "name\ttaxonomy_id\ttaxonomy_lvl\tkraken_assigned_reads\tadded_reads\t"
    "new_est_reads\tfraction_total_reads\n"
    "Plasmodiophora brassicae\t34289\tS\t900\t50\t950\t0.9500\n"
    "Escherichia coli\t562\tS\t30\t0\t30\t0.0300\n"
    "Staphylococcus aureus\t1280\tS\t1\t0\t1\t0.0001\n"
)

class TestParseBrackenReport:
    def test_rows(self, tmp_path):
        f = tmp_path / "S1.bracken.S.txt"
        f.write_text(_BRACKEN, encoding="utf-8")
        rows = parse_bracken_report(str(f))
        assert len(rows) == 3
        top = rows[0]
        assert top["name"] == "Plasmodiophora brassicae"
        assert top["fraction_total_reads"] == pytest.approx(0.95)
        assert top["new_est_reads"] == 950

    def test_missing_file_returns_empty(self, tmp_path):
        assert parse_bracken_report(str(tmp_path / "nope.txt")) == []


# ---------- 汇总表|summary tables ----------

def _per_sample():
    # bracken 原生 fraction_total_reads 分母=classified(950),非全部reads(1000)
    # 汇总口径必须自算 new_est/total,故此处原生值刻意给成分母950的形态
    # |bracken native fractions use classified-only denominator (950);
    # summary must recompute new_est/total (1000)
    return {
        "S1": {
            "k2": {"total_reads": 1000, "classified_reads": 950,
                   "unclassified_reads": 50},
            "bracken": [
                {"name": "Plasmodiophora brassicae", "taxonomy_id": "34289",
                 "taxonomy_lvl": "S", "new_est_reads": 950,
                 "fraction_total_reads": 1.0},
                {"name": "Escherichia coli", "taxonomy_id": "562",
                 "taxonomy_lvl": "S", "new_est_reads": 30,
                 "fraction_total_reads": 0.0316},
            ],
        },
        "S2": {"k2": {}, "bracken": []},   # 失败样本:dry_run 或 bracken 缺失
    }

class TestRenderSummary:
    def test_header_and_top_species(self):
        tsv = render_summary_tsv(_per_sample())
        lines = tsv.strip().splitlines()
        assert lines[0].split("\t") == [
            "sample", "total_reads", "classified_pct", "unclassified_pct",
            "top_species", "top_species_pct", "other_classified_pct",
            "n_species_ge_1pct"]
        row = lines[1].split("\t")
        assert row[0] == "S1"
        assert row[1] == "1000"
        assert row[2] == "95.00"
        assert row[3] == "5.00"
        assert row[4] == "Plasmodiophora brassicae"
        assert row[5] == "95.00"
        assert row[6] == "0.00"          # 95 - 95
        assert row[7] == "2"             # 95% 与 3% 均计入

    def test_failed_sample_row_is_na(self):
        tsv = render_summary_tsv(_per_sample())
        row = tsv.strip().splitlines()[2].split("\t")
        assert row[0] == "S2"
        assert row[2] == "NA"

class TestRenderDetail:
    def test_filters_below_min_fraction(self):
        tsv = render_detail_tsv(_per_sample(), min_fraction=0.001)
        lines = tsv.strip().splitlines()
        assert lines[0].split("\t") == [
            "sample", "name", "taxonomy_id", "taxonomy_lvl",
            "new_est_reads", "fraction_total_reads"]
        # S1 两条(950/1000、30/1000 均留);S2 无
        assert len(lines) == 3
        assert lines[1].split("\t")[1] == "Plasmodiophora brassicae"
        # fraction 列口径=占全部reads(含unclassified),非bracken原生
        # |fraction denominator = all reads incl. unclassified
        assert lines[1].split("\t")[5] == "0.950000"
        assert lines[2].split("\t")[5] == "0.030000"


# ---------- config 校验|config validation ----------

class TestConfigValidate:
    def _cfg(self, tmp_path, **kw):
        fa = tmp_path / "ref.fa"
        fa.write_text(">c1\nACGT\n", encoding="utf-8")
        fq = tmp_path / "fq"
        fq.mkdir(exist_ok=True)
        db = tmp_path / "k2db"
        db.mkdir(exist_ok=True)
        (db / "hash.k2d").write_bytes(b"x")   # validate 只查存在性|existence only
        base = dict(fastq_dir=str(fq), genome=str(fa),
                    output_dir=str(tmp_path / "out"), kraken2_db=str(db))
        base.update(kw)
        return MixraceConfig(**base)

    def test_step6_now_valid(self, tmp_path):
        self._cfg(tmp_path, step=6).validate()

    def test_step7_invalid(self, tmp_path):
        with pytest.raises(ValueError, match="step"):
            self._cfg(tmp_path, step=7).validate()

    def test_kraken2_db_checked_when_enabled(self, tmp_path):
        with pytest.raises(ValueError, match="kraken2"):
            self._cfg(tmp_path, kraken2_db=str(tmp_path / "nope")).validate()

    def test_kraken2_db_skipped_when_disabled(self, tmp_path):
        self._cfg(tmp_path, kraken2_db=str(tmp_path / "nope"),
                  run_kraken2=False).validate()

    def test_bad_bracken_level(self, tmp_path):
        with pytest.raises(ValueError, match="bracken_level"):
            self._cfg(tmp_path, bracken_level="X").validate()


# ---------- 命令构造与断点|commands & checkpoints ----------

class _FakeRunner:
    """假执行器:kraken2/bracken 调用时写出对应假报告|fake runner writing reports."""

    def __init__(self, k2_text=_K2_REPORT, bracken_text=_BRACKEN):
        self.logger = _QuietLogger()
        self.dry_run = False
        self.calls = []
        self.k2_text = k2_text
        self.bracken_text = bracken_text

    def run_conda(self, tool_path, args, description=""):
        self.calls.append((Path(tool_path).name, list(args)))
        argmap = {args[i]: args[i + 1] for i in range(len(args) - 1)
                  if args[i].startswith("-")}
        if Path(tool_path).name == "kraken2":
            Path(argmap["--report"]).write_text(self.k2_text, encoding="utf-8")
        else:
            Path(argmap["-o"]).write_text(self.bracken_text, encoding="utf-8")
        return True, "", ""


class _QuietLogger:
    def info(self, *a, **kw):
        pass

    def warning(self, *a, **kw):
        pass

    def error(self, *a, **kw):
        pass


class _FakeCkpt:
    def __init__(self):
        self.done = set()

    def exists(self, step):
        return step in self.done

    def create(self, step):
        self.done.add(step)


def _contam_config(tmp_path):
    fa = tmp_path / "ref.fa"
    fa.write_text(">c1\nACGT\n", encoding="utf-8")
    return MixraceConfig(fastq_dir=str(tmp_path), genome=str(fa),
                         output_dir=str(tmp_path / "out"),
                         kraken2_db=str(tmp_path))


class TestRunContamination:
    def _samples(self, tmp_path):
        (tmp_path / "S1_1.fq.gz").write_text("x", encoding="utf-8")
        (tmp_path / "S1_2.fq.gz").write_text("x", encoding="utf-8")
        return [{"sample": "S1", "r1": str(tmp_path / "S1_1.fq.gz"),
                 "r2": str(tmp_path / "S1_2.fq.gz")}]

    def test_kraken2_command_shape(self, tmp_path):
        cfg = _contam_config(tmp_path)
        runner = _FakeRunner()
        run_contamination(cfg, runner, _FakeCkpt(), self._samples(tmp_path))
        k2 = [a for n, a in runner.calls if n == "kraken2"][0]
        assert "--paired" in k2
        assert k2[k2.index("--report") + 1].endswith("S1_k2_report.txt")
        # 逐read明细指向 /dev/null,不得写真实文件
        assert k2[k2.index("--output") + 1] == "/dev/null"
        assert "--memory-mapping" not in k2
        assert str(cfg.threads) in k2

    def test_kraken2_memory_mapping_flag(self, tmp_path):
        cfg = _contam_config(tmp_path)
        cfg.kraken_memory_mapping = True
        runner = _FakeRunner()
        run_contamination(cfg, runner, _FakeCkpt(), self._samples(tmp_path))
        k2 = [a for n, a in runner.calls if n == "kraken2"][0]
        assert "--memory-mapping" in k2

    def test_bracken_command_shape_and_snap(self, tmp_path):
        cfg = _contam_config(tmp_path)
        cfg.read_length = 120            # 吸附到 100
        runner = _FakeRunner()
        run_contamination(cfg, runner, _FakeCkpt(), self._samples(tmp_path))
        br = [a for n, a in runner.calls if n == "bracken"][0]
        assert br[br.index("-r") + 1] == "100"
        assert br[br.index("-l") + 1] == "S"
        assert br[br.index("-i") + 1].endswith("S1_k2_report.txt")

    def test_outputs_written(self, tmp_path):
        cfg = _contam_config(tmp_path)
        run_contamination(cfg, _FakeRunner(), _FakeCkpt(), self._samples(tmp_path))
        out = Path(cfg.output_dir) / "08_contamination"
        assert (out / "kraken2" / "S1_k2_report.txt").exists()
        assert (out / "bracken" / "S1_bracken_S.txt").exists()
        summ = (out / "contamination_summary.tsv").read_text(encoding="utf-8")
        assert "S1" in summ and "Plasmodiophora brassicae" in summ
        assert (out / "contamination_detail.tsv").exists()

    def test_checkpoint_skips_kraken2(self, tmp_path):
        cfg = _contam_config(tmp_path)
        out = Path(cfg.output_dir) / "08_contamination" / "kraken2"
        out.mkdir(parents=True)
        (out / "S1_k2_report.txt").write_text(_K2_REPORT, encoding="utf-8")
        ckpt = _FakeCkpt()
        ckpt.done.add("kraken2_S1")
        runner = _FakeRunner()
        run_contamination(cfg, runner, ckpt, self._samples(tmp_path))
        assert not [a for n, a in runner.calls if n == "kraken2"]
        # bracken 未断点,仍要跑
        assert [a for n, a in runner.calls if n == "bracken"]

    def test_empty_samples_no_crash(self, tmp_path):
        cfg = _contam_config(tmp_path)
        run_contamination(cfg, _FakeRunner(), _FakeCkpt(), [])


# ---------- argparse 接线|CLI arg wiring ----------

class TestArgvParse:
    def _argv(self, monkeypatch, *extra):
        import sys
        argv = ["mixrace", "-g", "ref.fa"] + list(extra)
        monkeypatch.setattr(sys, "argv", argv)

    def test_defaults_kraken2_on(self, monkeypatch):
        from biopytools.mixrace.main import _argv_to_config
        self._argv(monkeypatch)
        cfg = _argv_to_config()
        assert cfg.run_kraken2 is True
        assert cfg.kraken_memory_mapping is False
        assert cfg.bracken_level == "S"

    def test_skip_kraken2(self, monkeypatch):
        from biopytools.mixrace.main import _argv_to_config
        self._argv(monkeypatch, "--skip-kraken2")
        assert _argv_to_config().run_kraken2 is False

    def test_kraken2_db_and_memory_mapping(self, monkeypatch):
        from biopytools.mixrace.main import _argv_to_config
        self._argv(monkeypatch, "--kraken2-db", "/data/k2",
                   "--kraken-memory-mapping")
        cfg = _argv_to_config()
        assert cfg.kraken2_db == "/data/k2"
        assert cfg.kraken_memory_mapping is True

    def test_bracken_level_and_step6(self, monkeypatch):
        from biopytools.mixrace.main import _argv_to_config
        self._argv(monkeypatch, "--bracken-level", "G", "--step", "6")
        cfg = _argv_to_config()
        assert cfg.bracken_level == "G"
        assert cfg.step == 6


# ---------- step6 调度|dispatch ----------

class TestStep6Dispatch:
    def _setup(self, tmp_path):
        fa = tmp_path / "ref.fa"
        fa.write_text(">c1\nACGT\n", encoding="utf-8")
        fq = tmp_path / "fq"
        fq.mkdir(exist_ok=True)
        db = tmp_path / "k2db"
        db.mkdir(exist_ok=True)
        (db / "hash.k2d").write_bytes(b"x")
        qc = tmp_path / "out" / "01_qc"
        qc.mkdir(parents=True)
        (qc / "S1_1_clean.fq.gz").write_text("x", encoding="utf-8")
        (qc / "S1_2_clean.fq.gz").write_text("x", encoding="utf-8")
        return MixraceConfig(fastq_dir=str(fq), genome=str(fa),
                             output_dir=str(tmp_path / "out"), step=6,
                             kraken2_db=str(db))

    def test_step6_runs_contamination(self, tmp_path, monkeypatch):
        import biopytools.mixrace.main as m
        monkeypatch.setattr(m, "write_software_versions", lambda *a, **kw: None)
        cfg = self._setup(tmp_path)
        runner = _FakeRunner()
        ckpt = _FakeCkpt()
        m.run_pipeline(cfg, runner, ckpt, _QuietLogger())
        summ = Path(cfg.output_dir) / "08_contamination" / "contamination_summary.tsv"
        assert summ.exists()
        assert "S1" in summ.read_text(encoding="utf-8")

    def test_disabled_kraken2_skips_step6(self, tmp_path, monkeypatch):
        import biopytools.mixrace.main as m
        monkeypatch.setattr(m, "write_software_versions", lambda *a, **kw: None)
        cfg = self._setup(tmp_path)
        cfg.run_kraken2 = False
        m.run_pipeline(cfg, _FakeRunner(), _FakeCkpt(), _QuietLogger())
        assert not (Path(cfg.output_dir) / "08_contamination").exists()


# ---------- 版本信息|software versions ----------

class TestSoftwareVersions:
    def test_kraken2_bracken_probed(self, tmp_path, monkeypatch):
        import biopytools.mixrace.utils as u
        from types import SimpleNamespace

        def fake_run(cmd, **kw):
            name = cmd[-2] if isinstance(cmd, list) else str(cmd)
            text = "Kraken version 2.17.1" if "kraken" in str(cmd) else "Bracken v3.0.1"
            return SimpleNamespace(stdout=text, stderr="")
        monkeypatch.setattr(u.subprocess, "run", fake_run)
        cfg = _contam_config(tmp_path)
        out = tmp_path / "versions.yml"
        u.write_software_versions(cfg, _QuietLogger(), str(out))
        import yaml
        data = yaml.safe_load(out.read_text(encoding="utf-8"))
        assert "kraken2" in data["tools"]
        assert "bracken" in data["tools"]
        assert data["parameters"]["run_kraken2"] is True
        assert "kraken2_db" in data["parameters"]
