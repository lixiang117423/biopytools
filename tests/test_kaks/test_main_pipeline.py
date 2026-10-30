"""主流程分块与合并测试|Main pipeline chunking and merge tests

并行设计: 配对按线程数分块,每块独立(比对→AXT→KaKs_Calculator),结果合并后统一后处理。
"""

import pandas as pd
import pytest

from biopytools.kaks.aligner import PairAlignmentResult
from biopytools.kaks.main import KaKsAnalyzer, chunk_list


class TestChunkList:
    """配对分块|Pair chunking"""

    def test_even_split(self):
        chunks = chunk_list(list(range(10)), 2)
        assert [len(c) for c in chunks] == [5, 5]
        assert [x for c in chunks for x in c] == list(range(10))

    def test_uneven_split_keeps_all(self):
        chunks = chunk_list(list(range(10)), 3)
        assert sum(len(c) for c in chunks) == 10
        assert [x for c in chunks for x in c] == list(range(10))

    def test_more_chunks_than_items(self):
        chunks = chunk_list([1, 2], 8)
        assert [x for c in chunks for x in c] == [1, 2]
        assert all(len(c) > 0 for c in chunks)

    def test_empty_input(self):
        assert chunk_list([], 4) == []


class TestMergeChunkResults:
    """分块结果合并|Chunk result merging"""

    def _make_analyzer(self, tmp_path):
        """构造分析器但不运行|Build analyzer without running"""
        return KaKsAnalyzer(
            fasta1="f1.fa", fasta2="f2.fa", pairs="p.txt",
            output_dir=str(tmp_path), kaks_path="KaKs_Calculator",
            check_installation=False,  # 不依赖本机装 KaKs_Calculator|no local binary needed
        )

    def test_alignment_stats_joined_by_pair(self, tmp_path):
        """比对统计列按pair_name并入结果|Alignment stats joined by pair_name"""
        analyzer = self._make_analyzer(tmp_path)
        chunk_df = pd.DataFrame({
            "Sequence": ["p1"], "Method": ["GMYN"], "Ka": [0.1], "Ks": [0.2],
            "Ka/Ks": [0.5], "Pair_ID": ["p1"],
        })
        records = [
            PairAlignmentResult(
                seq1_id="a", seq2_id="b", pair_name="p1",
                codon_aln1="ATG", codon_aln2="ATG",
                protein_identity=0.95, aln_codons=1, gap_codons=0,
                cds1_len=3, cds2_len=3,
            )
        ]
        merged = analyzer._merge_chunk_results([chunk_df], records)
        row = merged[merged["Pair_ID"] == "p1"].iloc[0]
        assert row["protein_identity"] == 0.95
        assert row["seq1_id"] == "a"
        assert row["seq2_id"] == "b"

    def test_skipped_pair_kept_with_reason(self, tmp_path):
        """被跳过的配对保留行并标记原因|Skipped pairs kept with reason"""
        analyzer = self._make_analyzer(tmp_path)
        chunk_df = pd.DataFrame({
            "Sequence": ["p_ok"], "Method": ["GMYN"], "Ka": [0.1], "Ks": [0.2],
            "Ka/Ks": [0.5], "Pair_ID": ["p_ok"],
        })
        records = [
            PairAlignmentResult(seq1_id="a", seq2_id="b", pair_name="p_ok"),
            PairAlignmentResult(seq1_id="a", seq2_id="c", pair_name="p_bad",
                                skip_reason="internal_stop"),
        ]
        merged = analyzer._merge_chunk_results([chunk_df], records)
        assert len(merged) == 2
        bad = merged[merged["Pair_ID"] == "p_bad"].iloc[0]
        assert bad["skip_reason"] == "internal_stop"
        assert pd.isna(bad["Ka/Ks"])

    def test_pair_missing_from_results_gets_nan_row(self, tmp_path):
        """KaKs漏算的配对也保留行|Pairs missing from KaKs output still get a row"""
        analyzer = self._make_analyzer(tmp_path)
        chunk_df = pd.DataFrame({
            "Sequence": ["p1"], "Method": ["GMYN"], "Ka": [0.1], "Ks": [0.2],
            "Ka/Ks": [0.5], "Pair_ID": ["p1"],
        })
        records = [
            PairAlignmentResult(seq1_id="a", seq2_id="b", pair_name="p1"),
            PairAlignmentResult(seq1_id="a", seq2_id="d", pair_name="p_kaks_dropped"),
        ]
        merged = analyzer._merge_chunk_results([chunk_df], records)
        assert len(merged) == 2
        dropped = merged[merged["Pair_ID"] == "p_kaks_dropped"].iloc[0]
        assert pd.isna(dropped["Ka"])

    def test_strict_mode_records_without_stats(self, tmp_path):
        """严格模式记录无identity但有长度|Strict-mode records lack identity but carry lengths"""
        analyzer = self._make_analyzer(tmp_path)
        chunk_df = pd.DataFrame({
            "Sequence": ["p1"], "Method": ["GMYN"], "Ka": [0.1], "Ks": [0.2],
            "Ka/Ks": [0.5], "Pair_ID": ["p1"],
        })
        records = [
            PairAlignmentResult(seq1_id="a", seq2_id="b", pair_name="p1",
                                codon_aln1="ATGGCT", codon_aln2="ATGGCG",
                                aln_codons=2, cds1_len=6, cds2_len=6),
        ]
        merged = analyzer._merge_chunk_results([chunk_df], records)
        row = merged.iloc[0]
        assert row["aln_codons"] == 2
        assert row["skip_reason"] == ""
