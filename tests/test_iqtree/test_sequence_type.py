"""sniff_sequence_type 嗅探函数测试|Tests for sniff_sequence_type

背景|Context: IQ-TREE 3.x 自动检测对简并码/N富集(非ACGT≥10%)的DNA比对
误报 "Unknown sequence type"(3.1.3 实测)，模块需显式传 -st。
"""

import pytest

from biopytools.iqtree.utils import sniff_sequence_type


class TestSniffFasta:
    """FASTA 嗅探|FASTA sniffing"""

    def test_dna_with_ambiguity_codes(self, tmp_path):
        # 含 R/Y/S/W/K/M 简并码的DNA比对必须判为DNA
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGTRYSWKMBDHVN\n>s2\nTTTACGTAAAACGTN\n")
        assert sniff_sequence_type(str(aln)) == 'DNA'

    def test_plain_dna(self, tmp_path):
        aln = tmp_path / "aln.fa"
        aln.write_text(">s1\nACGTACGT\n>s2\nACGTACGA\n")
        assert sniff_sequence_type(str(aln)) == 'DNA'

    def test_lowercase_dna(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nacgtnnn\n>s2\nacgtacg\n")
        assert sniff_sequence_type(str(aln)) == 'DNA'

    def test_protein(self, tmp_path):
        aln = tmp_path / "protein.fasta"
        aln.write_text(">p1\nMKVLPEQILNF\n>p2\nMKVLPEQILNY\n")
        assert sniff_sequence_type(str(aln)) == 'AA'

    def test_dna_with_gaps_and_question(self, tmp_path):
        aln = tmp_path / "aln.fasta"
        aln.write_text(">s1\nACGT---??N\n>s2\nACGTACGTTTN\n")
        assert sniff_sequence_type(str(aln)) == 'DNA'


class TestSniffPhylip:
    """PHYLIP 嗅探|PHYLIP sniffing"""

    def test_sequential_dna_with_ambiguity(self, tmp_path):
        # 复刻实际事故文件形态：头行+名字+填充空格+长序列(含简并码)
        aln = tmp_path / "aln.phy"
        aln.write_text(
            "3 20\n"
            "41            ACGTRYSWKMBDHVNACGTA\n"
            "SRR10321576   ACGTACGTAAAACGTNNNNN\n"
            "T30-4         ACGTACGTAAAACGTNNNNA\n"
        )
        assert sniff_sequence_type(str(aln)) == 'DNA'

    def test_name_with_protein_chars_still_dna(self, tmp_path):
        # 名字含蛋白特异字母(E/F/I/L/P/Q/Z)不得污染判定——名字区必须被剥离
        aln = tmp_path / "aln.phy"
        aln.write_text(
            "2 12\n"
            "EFPLQZ_1      ACGTACGTAAAC\n"
            "EFPLQZ_2      ACGTACGTAAAT\n"
        )
        assert sniff_sequence_type(str(aln)) == 'DNA'

    def test_protein_phylip(self, tmp_path):
        aln = tmp_path / "protein.phy"
        aln.write_text(
            "2 10\n"
            "sp1           MKVLPEQILN\n"
            "sp2           MKVLPEQILD\n"
        )
        assert sniff_sequence_type(str(aln)) == 'AA'

    def test_interleaved_second_block_without_names(self, tmp_path):
        # 交错格式第二块行首无名字, 整行是序列, 不得误剥前10列
        aln = tmp_path / "aln.phy"
        aln.write_text(
            "2 16\n"
            "sp1           ACGTACGTACGTACGT\n"
            "sp2           ACGTACGTACGTACGG\n"
            "ACGTRYSW\n"
            "ACGTNNNN\n"
        )
        assert sniff_sequence_type(str(aln)) == 'DNA'


class TestSniffFallback:
    """未知格式回退|Fallback for unknown formats"""

    def test_nexus_returns_none(self, tmp_path):
        # NEXUS 复杂格式不嗅探, 回退IQ-TREE自动检测(现状行为)
        aln = tmp_path / "aln.nex"
        aln.write_text("#NEXUS\nBEGIN DATA;\nMATRIX\nsp1 ACGT\n;\nEND;\n")
        assert sniff_sequence_type(str(aln)) is None

    def test_unknown_suffix_returns_none(self, tmp_path):
        aln = tmp_path / "aln.txt"
        aln.write_text("ACGTACGT\n")
        assert sniff_sequence_type(str(aln)) is None

    def test_nonexistent_file_returns_none(self, tmp_path):
        # 文件不存在时不抛异常(后续validate会拦截), 返回None回退自动检测
        assert sniff_sequence_type(str(tmp_path / "nope.fasta")) is None
