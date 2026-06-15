"""候选基因Junction Reads分析模块|Candidate Gene Junction Reads Analysis Module"""

import os
import logging
from dataclasses import dataclass
from typing import Dict, List, Tuple, Optional

from .config import GeneRnaseqCheckConfig
from .parse_gff import EffectorRecord


@dataclass
class JunctionInfo:
    """Junction信息|Junction information"""
    donor: int       # 0-based donor position (inclusive)
    acceptor: int    # 0-based acceptor position (exclusive)
    read_count: int
    matches_annotation: bool = False
    canonical_splice: bool = False


@dataclass
class GeneJunctions:
    """基因Junction分析结果|Gene junction analysis result"""
    gene_id: str
    is_multi_exon: bool
    matched_junction_count: int
    novel_junctions: List[JunctionInfo]
    novel_junction_count: int
    total_junction_reads: int


def _get_junctions_from_read(read) -> List[Tuple[int, int]]:
    """从read的CIGAR提取junction坐标|Extract junction coordinates from read CIGAR

    Returns:
        [(donor_0based, acceptor_0based_exclusive), ...]
        donor: 剪接供体位点（0-based, inclusive）
        acceptor: 剪接受体位点（0-based, exclusive, 即第一个比对碱基的位置）
    """
    junctions = []
    ref_pos = read.reference_start  # 0-based

    if not read.cigartuples:
        return junctions

    for length, op in read.cigartuples:
        if op == 0 or op == 2 or op == 3 or op == 7 or op == 8:
            # M/D/N/=/X 消耗参考序列|Consume reference
            ref_pos += length
        elif op == 3:
            # N (intron skip) - ref_pos 已经跳过了 intron
            # donor = 跳跃前的最后一个位置 (0-based: ref_pos_before_N)
            # acceptor = 跳跃后的第一个位置 (0-based: ref_pos_after_N)
            donor = ref_pos - length  # N之前的位置
            acceptor = ref_pos          # N之后的位置
            junctions.append((donor, acceptor))
        # op 1 (I) 和 op 4 (S) 不消耗参考序列|I and S don't consume reference

    return junctions


def _junction_matches_annotation(
    junction: Tuple[int, int],
    introns: List[Tuple[int, int]],
    tolerance: int = 5,
) -> bool:
    """检查junction是否与注释的intron匹配|Check if junction matches annotated intron

    Args:
        junction: (donor_0based, acceptor_0based_exclusive)
        introns: [(donor_1based, acceptor_1based), ...] GFF3 1-based
        tolerance: 容差碱基数|Tolerance in bp

    Returns:
        True if any annotated intron matches within tolerance
    """
    j_donor, j_acceptor = junction
    # 转换到1-based: donor_1based = donor_0based + 1, acceptor_1based = acceptor_0based
    j_donor_1 = j_donor + 1
    j_acceptor_1 = j_acceptor  # 0-based exclusive = 1-based inclusive

    for intron_donor, intron_acceptor in introns:
        if (abs(j_donor_1 - intron_donor) <= tolerance and
                abs(j_acceptor_1 - intron_acceptor) <= tolerance):
            return True
    return False


def _check_canonical_splice(fasta, chrom: str, donor_0: int, acceptor_0: int) -> bool:
    """检查剪接位点是否为经典GT-AG|Check if splice site is canonical GT-AG"""
    try:
        donor_seq = fasta.fetch(chrom, donor_0, donor_0 + 2).upper()
        acceptor_seq = fasta.fetch(chrom, acceptor_0 - 2, acceptor_0).upper()

        canonical_pairs = [
            (donor_seq == 'GT' and acceptor_seq == 'AG'),
            (donor_seq == 'GC' and acceptor_seq == 'AG'),
            (donor_seq == 'AT' and acceptor_seq == 'AC'),
        ]
        return any(canonical_pairs)
    except Exception:
        return False


class JunctionAnalyzer:
    """Junction Reads分析器|Junction Reads Analyzer"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def analyze_gene_junctions(
        self,
        bam_files: List[str],
        records: Dict[str, EffectorRecord],
        genome_fa: str = None,
    ) -> Dict[str, GeneJunctions]:
        """分析所有目标基因的junction reads|Analyze junction reads for all target genes

        Args:
            bam_files: BAM文件列表|BAM file list
            records: 基因注释记录|Gene annotation records
            genome_fa: 基因组FASTA（用于检查经典剪接位点）|Genome FASTA for canonical splice check

        Returns:
            基因ID到GeneJunctions的映射|Gene ID to GeneJunctions mapping
        """
        import pysam

        # 打开基因组FASTA（可选）|Open genome FASTA (optional)
        fasta = None
        if genome_fa and os.path.exists(genome_fa):
            try:
                fasta = pysam.FastaFile(genome_fa)
            except Exception:
                self.logger.warning("无法打开基因组FASTA，跳过经典剪接检查|Cannot open genome FASTA, skipping canonical check")

        results = {}
        total_genes = len(records)
        processed = 0

        for gene_id, rec in records.items():
            processed += 1
            if processed % 50 == 0 or processed == total_genes:
                self.logger.info(
                    f"Junction分析进度|Junction analysis progress: {processed}/{total_genes}"
                )

            if not rec.is_multi_exon:
                results[gene_id] = GeneJunctions(
                    gene_id=gene_id,
                    is_multi_exon=False,
                    matched_junction_count=0,
                    novel_junctions=[],
                    novel_junction_count=0,
                    total_junction_reads=0,
                )
                continue

            # 收集所有BAM的junction reads|Collect junction reads from all BAMs
            junction_counts = {}  # (donor, acceptor) -> read_count

            for bam_path in bam_files:
                try:
                    bam = pysam.AlignmentFile(bam_path, "rb")
                except Exception:
                    continue

                for read in bam.fetch(rec.chrom, rec.gene_start - 1, rec.gene_end):
                    if read.is_unmapped or read.is_secondary or read.is_supplementary:
                        continue
                    for junction in _get_junctions_from_read(read):
                        junction_counts[junction] = junction_counts.get(junction, 0) + 1
                bam.close()

            # 分类junctions|Classify junctions
            supporting_count = 0
            novel_junctions = []
            total_reads = 0

            for (donor, acceptor), count in junction_counts.items():
                matches = _junction_matches_annotation(
                    (donor, acceptor), rec.introns, self.config.junction_tolerance
                )
                canonical = False
                if fasta and not matches:
                    canonical = _check_canonical_splice(fasta, rec.chrom, donor, acceptor)

                info = JunctionInfo(
                    donor=donor, acceptor=acceptor, read_count=count,
                    matches_annotation=matches, canonical_splice=canonical,
                )

                if matches:
                    supporting_count += count
                else:
                    novel_junctions.append(info)
                total_reads += count

            results[gene_id] = GeneJunctions(
                gene_id=gene_id,
                is_multi_exon=True,
                matched_junction_count=supporting_count,
                novel_junctions=novel_junctions,
                novel_junction_count=len(novel_junctions),
                total_junction_reads=total_reads,
            )

        if fasta:
            fasta.close()

        # 保存结果|Save results
        out_dir = os.path.join(self.config.output_dir, '03_junctions')
        os.makedirs(out_dir, exist_ok=True)
        out_file = os.path.join(out_dir, 'junction_analysis.tsv')
        with open(out_file, 'w') as f:
            f.write("gene_id\tis_multi_exon\tmatched_junction_reads\t"
                    "novel_junction_count\ttotal_junction_reads\t"
                    "novel_junction_details\n")
            for gid, junc in sorted(results.items()):
                if junc.is_multi_exon and junc.novel_junctions:
                    details = ';'.join(
                        f"({j.donor},{j.acceptor}):{j.read_count}"
                        f"{'_canon' if j.canonical_splice else ''}"
                        for j in junc.novel_junctions
                    )
                else:
                    details = 'NA'
                f.write(f"{gid}\t{junc.is_multi_exon}\t{junc.matched_junction_count}\t"
                        f"{junc.novel_junction_count}\t{junc.total_junction_reads}\t"
                        f"{details}\n")

        self.logger.info(f"Junction分析完成|Junction analysis completed: {out_file}")
        return results
