"""build_command 传 -st 的测试|Tests for -st handling in TreeBuilder.build_command

复现事故|Reproduced incident: psoja_365.min4.phy(非ACGT≈15%)在IQ-TREE 3.1.3下
无-st即报 "Unknown sequence type" 退出码2(2026-09-02, joblogs/...iqtree-ps.err)。
"""

import logging

import pytest

from biopytools.iqtree.config import TreeConfig
from biopytools.iqtree.analysis import TreeBuilder


def _make_builder(tmp_path, input_file, **kwargs):
    """构造最小TreeBuilder(不执行命令)|Build minimal TreeBuilder without running"""
    cfg_kwargs = dict(
        input_file=str(input_file),
        output_dir=str(tmp_path / "out"),
        iqtree_path="/nonexistent/bin/iqtree",
    )
    cfg_kwargs.update(kwargs)
    config = TreeConfig(**cfg_kwargs)
    logger = logging.getLogger("test_iqtree_st")
    return TreeBuilder(config, logger, cmd_runner=None)


class TestBuildCommandSequenceType:
    """build_command 的 -st 接线|-st wiring in build_command"""

    def test_explicit_sequence_type_dna(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        builder = _make_builder(tmp_path, aln, sequence_type='DNA')
        cmd = builder.build_command()
        idx = cmd.index('-st')
        assert cmd[idx + 1] == 'DNA'

    def test_explicit_sequence_type_aa(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nMKVL\n")
        builder = _make_builder(tmp_path, aln, sequence_type='AA')
        cmd = builder.build_command()
        idx = cmd.index('-st')
        assert cmd[idx + 1] == 'AA'

    def test_sniffed_dna_alignment_gets_st_dna(self, tmp_path):
        # 事故场景: 简并码富集DNA比对, 未显式指定→嗅探必须补上-st DNA
        aln = tmp_path / "aln.min4.phy"
        aln.write_text(
            "2 40\n"
            "s1            ACGTRYSWKMBDHVNACGTACGTRYSWKMBDHVNACGT\n"
            "s2            ACGTNNNNACGTACGTRYSWACGTNNNNACGTACGT\n"
        )
        builder = _make_builder(tmp_path, aln)
        cmd = builder.build_command()
        idx = cmd.index('-st')
        assert cmd[idx + 1] == 'DNA'

    def test_sniffed_protein_gets_st_aa(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">p1\nMKVLPEQILNF\n>p2\nMKVLPEQILNY\n")
        builder = _make_builder(tmp_path, aln)
        cmd = builder.build_command()
        idx = cmd.index('-st')
        assert cmd[idx + 1] == 'AA'

    def test_unknown_format_omits_st(self, tmp_path):
        # NEXUS不嗅探→不传-st, 保持IQ-TREE自动检测(现状行为)
        aln = tmp_path / "aln.nex"
        aln.write_text("#NEXUS\nBEGIN DATA;\nMATRIX\ns1 ACGT\n;\nEND;\n")
        builder = _make_builder(tmp_path, aln)
        cmd = builder.build_command()
        assert '-st' not in cmd


class TestTreeConfigSequenceType:
    """TreeConfig.sequence_type 校验|Validation"""

    def test_invalid_sequence_type_rejected(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        builder = _make_builder(tmp_path, aln, sequence_type='RNA')
        with pytest.raises(ValueError, match='序列类型'):
            builder.config.validate()

    def test_valid_sequence_types_pass(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        for st in ('DNA', 'AA', 'BIN', 'MORPH'):
            builder = _make_builder(tmp_path, aln, sequence_type=st)
            assert builder.config.validate() is True

    def test_lowercase_normalized(self, tmp_path):
        # 小写输入归一为大写后传给IQ-TREE
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT\n")
        builder = _make_builder(tmp_path, aln, sequence_type='dna')
        assert builder.config.sequence_type == 'DNA'
