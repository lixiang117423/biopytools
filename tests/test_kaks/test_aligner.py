"""密码子比对器测试|Codon aligner tests

背景|Context: kaks 模块原 prepare_input_file 对不等长序列直接 [:min_len] 截断,
含 indel 的配对喂给 KaKs_Calculator 后全部错位(项目28抽样36.1%配对受害)。
本测试定义「翻译→逐对蛋白比对→回译密码子比对」的标准流程行为。
"""

import pytest

from biopytools.kaks.aligner import (
    CodonAligner,
    PairAlignmentResult,
    back_translate,
    translate_cds,
)


class TestTranslateCds:
    """CDS翻译|CDS translation"""

    def test_valid_cds_with_terminal_stop(self):
        """末端终止子被剥离,蛋白与截短CDS一一对应|Terminal stop stripped, protein matches stopless CDS"""
        protein, cds_stripped, reason = translate_cds("ATGGCTAAATTGTAA")
        assert reason == ""
        assert protein == "MAKL"
        assert cds_stripped == "ATGGCTAAATTG"

    def test_valid_cds_without_terminal_stop(self):
        """无末端终止子的CDS正常翻译|CDS without terminal stop translates fine"""
        protein, cds_stripped, reason = translate_cds("ATGGCTAAATTG")
        assert reason == ""
        assert protein == "MAKL"
        assert cds_stripped == "ATGGCTAAATTG"

    def test_lowercase_input_normalized(self):
        """小写输入归一为大写|Lowercase input normalized to upper"""
        protein, _, reason = translate_cds("atggctaaattgtaa")
        assert reason == ""
        assert protein == "MAKL"

    def test_internal_stop_reports_reason(self):
        """内部终止子返回原因,不产蛋白|Internal stop returns reason, no protein"""
        protein, cds_stripped, reason = translate_cds("ATGTAAAAGGCTTAA")
        assert reason == "internal_stop"
        assert protein is None
        assert cds_stripped is None

    def test_length_not_multiple_of_three(self):
        """长度非3倍数返回原因|Length not multiple of 3 returns reason"""
        _, _, reason = translate_cds("ATGGCTAAATT")
        assert reason == "length_not_multiple_of_3"

    def test_empty_sequence(self):
        """空序列返回原因|Empty sequence returns reason"""
        _, _, reason = translate_cds("")
        assert reason == "empty_sequence"

    def test_ambiguous_codon_becomes_x_and_preserved(self):
        """NNN密码子译为X且CDS原样保留(无损)|NNN codon translates to X, CDS preserved losslessly"""
        protein, cds_stripped, reason = translate_cds("ATGNNNTTGTAA")
        assert reason == ""
        assert protein == "MXL"
        assert cds_stripped == "ATGNNNTTG"


class TestBackTranslate:
    """蛋白比对回译|Protein alignment back-translation"""

    def test_gap_becomes_triplet(self):
        """蛋白gap扩展为---三连|Protein gap expands to --- triplet"""
        codons = back_translate("MK-L", "ATGAAATTG")
        assert codons == "ATGAAA---TTG"

    def test_roundtrip_recovers_cds(self):
        """去掉---应还原原始CDS|Removing --- recovers original CDS"""
        aligned_prot = "MK-L"
        cds = "ATGAAATTG"
        codons = back_translate(aligned_prot, cds)
        assert codons.replace("---", "") == cds

    def test_length_mismatch_raises(self):
        """蛋白列数与cds密码子数不符时报错|Protein columns vs codon count mismatch raises"""
        with pytest.raises(ValueError):
            back_translate("MK", "ATGAAATTGTTA")


class TestCodonAlignerAlignPair:
    """逐对端到端比对|Pairwise end-to-end alignment"""

    def setup_method(self):
        self.aligner = CodonAligner()

    def test_indel_pair_equal_length_with_gap(self):
        """含indel配对产出等长密码子比对且带---|Indel pair yields equal-length codon alignment with ---"""
        cds1 = "ATGGCTAAATTGTAA"   # MA KL
        cds2 = "ATGGCTTTGTAA"      # MA-L (缺AAA, 3bp缺失)
        result = self.aligner.align_pair("s1", "s2", "p1", cds1, cds2)
        assert isinstance(result, PairAlignmentResult)
        assert result.skip_reason == ""
        assert len(result.codon_aln1) == len(result.codon_aln2)
        assert "---" in result.codon_aln1 or "---" in result.codon_aln2
        assert result.gap_codons >= 1
        assert result.aln_codons == len(result.codon_aln1) // 3

    def test_gap_side_not_counted_as_substitution(self):
        """回译去掉gap后应还原各自CDS(去终止子)|Back-translation minus gaps recovers each stopless CDS"""
        cds1 = "ATGGCTAAATTGTAA"
        cds2 = "ATGGCTTTGTAA"
        result = self.aligner.align_pair("s1", "s2", "p1", cds1, cds2)
        assert result.codon_aln1.replace("---", "") == "ATGGCTAAATTG"
        assert result.codon_aln2.replace("---", "") == "ATGGCTTTG"

    def test_identical_pair_no_gaps_full_identity(self):
        """相同序列无gap,identity=1|Identical sequences yield no gaps, identity=1"""
        cds = "ATGGCTAAATTGTAA"
        result = self.aligner.align_pair("s1", "s2", "p1", cds, cds)
        assert result.skip_reason == ""
        assert result.gap_codons == 0
        assert result.protein_identity == 1.0

    def test_substitution_lowers_identity(self):
        """氨基酸替换降低identity|Amino acid substitution lowers identity"""
        cds1 = "ATGGCTAAATTGTAA"   # MAKL
        cds2 = "ATGGCTAAAGAATAA"   # MAKE (K保留, L->E)
        result = self.aligner.align_pair("s1", "s2", "p1", cds1, cds2)
        assert 0.0 < result.protein_identity < 1.0
        assert result.gap_codons == 0

    def test_internal_stop_skips_pair(self):
        """内部终止子配对被跳过并记录原因|Pair with internal stop is skipped with reason"""
        cds1 = "ATGGCTAAATTGTAA"
        cds2 = "ATGTAAAAGGCTTAA"   # 内部TAA
        result = self.aligner.align_pair("s1", "s2", "p1", cds1, cds2)
        assert result.skip_reason == "internal_stop"
        assert result.codon_aln1 == ""
        assert result.codon_aln2 == ""

    def test_result_carries_pair_identifiers(self):
        """结果携带配对标识与原始长度|Result carries pair identifiers and raw lengths"""
        cds1 = "ATGGCTAAATTGTAA"
        cds2 = "ATGGCTTTGTAA"
        result = self.aligner.align_pair("id1", "id2", "pairname", cds1, cds2)
        assert result.seq1_id == "id1"
        assert result.seq2_id == "id2"
        assert result.pair_name == "pairname"
        assert result.cds1_len == len(cds1)
        assert result.cds2_len == len(cds2)

    def test_alignment_is_deterministic(self):
        """同输入两次比对结果一致|Same input aligns identically twice"""
        cds1 = "ATGGCTAAATTGCGTTAA"
        cds2 = "ATGGCTCGTTAA"
        r1 = self.aligner.align_pair("s1", "s2", "p1", cds1, cds2)
        r2 = self.aligner.align_pair("s1", "s2", "p1", cds1, cds2)
        assert r1.codon_aln1 == r2.codon_aln1
        assert r1.codon_aln2 == r2.codon_aln2
