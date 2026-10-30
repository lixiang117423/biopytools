"""RnaIsoConfig 嗅探与校验测试|Tests for input sniffing and validation

isoseq3 26.2 CLI 适配后语义|Semantics after isoseq3 26.2 adaptation:
- isoseq3_supported: "full"(subreads/ccs.bam) / "none"(reads fasta、ONT)
- 流程: ccs → refine(去引物) → isoquant/cluster2
"""

import os

import pytest

from biopytools.rna_iso.config import (
    CLONTECH_PRIMERS_FA,
    RnaIsoConfig,
    detect_input_type,
)


class TestDetectInputType:
    def test_subreads(self):
        assert detect_input_type(["m1.subreads.bam", "m2.subreads.bam"]) == "pacbio_subreads"

    def test_ccs_bam(self):
        assert detect_input_type(["sample.ccs.bam"]) == "pacbio_ccs"

    def test_hifi_bam(self):
        assert detect_input_type(["sample.hifi_reads.bam"]) == "pacbio_ccs"

    def test_reads_fasta(self):
        assert detect_input_type(["s1.fa.gz", "s2.fasta"]) == "reads"

    def test_reads_fastq(self):
        assert detect_input_type(["s1.fq", "s2.fastq.gz"]) == "reads"

    def test_mixed_raises(self):
        with pytest.raises(ValueError):
            detect_input_type(["m1.subreads.bam", "s1.fq.gz"])

    def test_unknown_raises(self):
        with pytest.raises(ValueError):
            detect_input_type(["s1.txt"])


def _cfg(tmp_path, **kw):
    """构造最小RnaIsoConfig并校验(ONT reads场景)|Build minimal RnaIsoConfig and validate"""
    args = dict(
        reads=[str(tmp_path / "s.fq.gz")],
        data_type="ont",
        reference=str(tmp_path / "ref.fa"),
        output_dir=str(tmp_path / "out"),
    )
    args.update(kw)
    for name in ("s.fq.gz", "ref.fa"):
        (tmp_path / name).write_text(">x\nACGT\n")
    config = RnaIsoConfig(**args)
    config.validate()
    return config


class TestDerived:
    def test_ont_reads_kind(self, tmp_path):
        cfg = _cfg(tmp_path)
        assert cfg.input_kind == "reads"
        assert cfg.needs_ccs is False
        assert cfg.needs_refine is False
        assert cfg.isoquant_data_type == "nanopore"
        assert cfg.isoseq3_supported == "none"

    def test_subreads_full_chain(self, tmp_path):
        (tmp_path / "m.subreads.bam").write_text("x")
        cfg = _cfg(tmp_path, reads=[str(tmp_path / "m.subreads.bam")], data_type=None)
        assert cfg.input_kind == "pacbio_subreads"
        assert cfg.needs_ccs is True
        assert cfg.needs_refine is True
        assert cfg.isoseq3_supported == "full"
        assert cfg.isoquant_data_type == "pacbio_ccs"

    def test_ccs_bam_skips_ccs_only_refine(self, tmp_path):
        (tmp_path / "m.ccs.bam").write_text("x")
        cfg = _cfg(tmp_path, reads=[str(tmp_path / "m.ccs.bam")], data_type=None)
        assert cfg.input_kind == "pacbio_ccs"
        assert cfg.needs_ccs is False
        assert cfg.needs_refine is True
        assert cfg.isoseq3_supported == "full"

    def test_pacbio_reads_isoseq3_none(self, tmp_path):
        (tmp_path / "s.fa").write_text(">x\nACGT\n")
        cfg = _cfg(tmp_path, reads=[str(tmp_path / "s.fa")], data_type="pacbio")
        assert cfg.input_kind == "reads"
        assert cfg.isoseq3_supported == "none"
        assert cfg.isoquant_data_type == "pacbio_ccs"

    def test_reads_fasta_needs_datatype(self, tmp_path):
        (tmp_path / "s.fa").write_text(">x\nACGT\n")
        with pytest.raises(ValueError, match="data-type"):
            _cfg(tmp_path, reads=[str(tmp_path / "s.fa")], data_type=None)


class TestValidate:
    def test_ok(self, tmp_path):
        assert _cfg(tmp_path).validate()

    def test_isoquant_needs_reference(self, tmp_path):
        with pytest.raises(ValueError, match="reference"):
            _cfg(tmp_path, reference=None)

    def test_isoseq3_only_no_reference_required(self, tmp_path):
        (tmp_path / "m.subreads.bam").write_text("x")
        cfg = _cfg(
            tmp_path,
            reads=[str(tmp_path / "m.subreads.bam")],
            data_type=None,
            reference=None,
            engine="isoseq3",
        )
        assert cfg.validate()

    def test_bad_engine(self, tmp_path):
        with pytest.raises(ValueError, match="engine"):
            _cfg(tmp_path, engine="bogus")

    def test_reads_missing(self, tmp_path):
        with pytest.raises(ValueError, match="reads"):
            _cfg(tmp_path, reads=[str(tmp_path / "nope.fq.gz")])

    def test_isoseq3_unsupported_on_reads(self, tmp_path):
        with pytest.raises(ValueError, match="isoseq3"):
            _cfg(tmp_path, engine="isoseq3")

    def test_isoseq3_unsupported_on_ont(self, tmp_path):
        with pytest.raises(ValueError, match="isoseq3"):
            _cfg(tmp_path, engine="both", data_type="ont")

    def test_primers_need_two_records(self, tmp_path):
        bad = tmp_path / "bad.fa"
        bad.write_text(">a\nACGT\n")
        with pytest.raises(ValueError, match="primers"):
            _cfg(tmp_path, primers=str(bad))

    def test_primers_file_missing(self, tmp_path):
        with pytest.raises(ValueError, match="primers"):
            _cfg(tmp_path, primers=str(tmp_path / "nope.fa"))

    def test_genedb_missing(self, tmp_path):
        with pytest.raises(ValueError, match="genedb"):
            _cfg(tmp_path, genedb=str(tmp_path / "nope.gtf"))

    def test_bad_threads(self, tmp_path):
        with pytest.raises(ValueError, match="线程|threads"):
            _cfg(tmp_path, threads=0)


class TestDirectoryInput:
    """-i 目录支持:自动展开reads文件|Directory input: auto-expand reads files"""

    def _make_dir(self, tmp_path):
        import os
        d = tmp_path / "clean"
        d.mkdir()
        for name in ("S2_subreads_clean.fq.gz", "S1_subreads_clean.fq.gz"):
            (d / name).write_text(">x\nACGT\n")
        (d / "fastp_reports").mkdir()
        (d / "fastp_reports" / "report.html").write_text("x")
        return d

    def test_dir_expands_sorted_reads(self, tmp_path):
        d = self._make_dir(tmp_path)
        cfg = _cfg(tmp_path, reads=[str(d)])
        basenames = [os.path.basename(r) for r in cfg.reads]
        # 排序+只取reads文件(忽略子目录/html)|sorted, reads files only (skip subdirs)
        assert basenames == ["S1_subreads_clean.fq.gz", "S2_subreads_clean.fq.gz"]

    def test_dir_no_reads_raises(self, tmp_path):
        empty = tmp_path / "empty"
        empty.mkdir()
        with pytest.raises(ValueError, match="目录|directory"):
            _cfg(tmp_path, reads=[str(empty)])

    def test_mixed_dir_and_file(self, tmp_path):
        d = self._make_dir(tmp_path)
        cfg = _cfg(tmp_path, reads=[str(d), str(tmp_path / "s.fq.gz")])
        assert len(cfg.reads) == 3

    def test_dir_in_input_kind_reads(self, tmp_path):
        d = self._make_dir(tmp_path)
        cfg = _cfg(tmp_path, reads=[str(d)])
        assert cfg.input_kind == "reads"
        assert cfg.validate()


class TestPrimersFasta:
    def test_builtin_written_to_tmp(self, tmp_path):
        cfg = _cfg(tmp_path)
        target = cfg.write_primers_fasta()
        assert target == cfg.tmp_dir + "/primers.fasta"
        content = open(target).read()
        assert content == CLONTECH_PRIMERS_FA
        assert "primer_5p" in content

    def test_user_primers_copied(self, tmp_path):
        user = tmp_path / "user_primers.fa"
        user.write_text(">p5\nACGT\n>p3\nTTTT\n")
        cfg = _cfg(tmp_path, primers=str(user))
        target = cfg.write_primers_fasta()
        assert open(target).read() == ">p5\nACGT\n>p3\nTTTT\n"

    def test_clontech_constant(self):
        assert "primer_5p" in CLONTECH_PRIMERS_FA
        assert "primer_3p" in CLONTECH_PRIMERS_FA
        assert "AAGCAGTGGTATCAACGCAGAGTACATGGG" in CLONTECH_PRIMERS_FA
