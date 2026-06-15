"""候选基因覆盖度分析模块|Candidate Gene Coverage Analysis Module"""

import os
import logging
from dataclasses import dataclass, field
from typing import Dict, List, Tuple

from .config import GeneRnaseqCheckConfig
from .utils import CommandRunner, build_conda_command
from .parse_gff import EffectorRecord, ExonInfo, build_exon_bed, build_gene_bed, build_flanking_bed


@dataclass
class ExonCoverage:
    """单个exon的覆盖度信息|Single exon coverage information"""
    gene_id: str
    exon_index: int
    start: int       # GFF3 1-based
    end: int
    breadth: float   # 0-1, 被覆盖碱基比例
    mean_depth: float


@dataclass
class GeneCoverage:
    """基因整体覆盖度|Gene overall coverage"""
    gene_id: str
    overall_breadth: float    # 0-1
    overall_mean_depth: float
    exons: List[ExonCoverage]
    upstream_mean_depth: float = 0.0
    downstream_mean_depth: float = 0.0
    min_exon_breadth: float = 1.0
    min_exon_depth: float = float('inf')


def _parse_bedtools_coverage_line(line: str) -> Tuple[str, float, float]:
    """解析bedtools coverage输出行|Parse bedtools coverage output line

    bedtools coverage -mean 输出格式:
    chrom  start  end  name  score  strand  num_reads  length  breadth  mean_depth

    Returns:
        (name, breadth, mean_depth)
    """
    cols = line.strip().split('\t')
    if len(cols) < 10:
        return cols[3], 0.0, 0.0
    name = cols[3]
    try:
        breadth = float(cols[9]) / 100.0 if '%' in cols[9] else float(cols[9])
        if cols[9].endswith('.00'):
            breadth = float(cols[9]) / 100.0
        else:
            breadth = float(cols[9])
        mean_depth = float(cols[10]) if len(cols) > 10 else 0.0
    except (ValueError, IndexError):
        breadth = 0.0
        mean_depth = 0.0

    # bedtools coverage -mean 第10列是 breadth_pct，第11列是 meandepth_per_base
    # 但某些版本列顺序不同，做安全解析
    try:
        # 尝试倒数第二列作为breadth, 最后一列作为depth
        if len(cols) >= 12:
            breadth = float(cols[-2])
            mean_depth = float(cols[-1])
        elif len(cols) >= 10:
            # 列顺序: ... numBasesCovered  numBases  fraction  meandepth
            if '/' in cols[9]:
                parts = cols[9].split('/')
                breadth = float(parts[0]) / float(parts[1]) if float(parts[1]) > 0 else 0.0
                mean_depth = float(cols[10]) if len(cols) > 10 else 0.0
            else:
                breadth = float(cols[9]) if float(cols[9]) <= 1.0 else float(cols[9]) / 100.0
                mean_depth = float(cols[10]) if len(cols) > 10 else 0.0
    except (ValueError, ZeroDivisionError):
        pass

    return name, breadth, mean_depth


class CoverageAnalyzer:
    """覆盖度分析器|Coverage Analyzer"""

    def __init__(self, config: GeneRnaseqCheckConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def compute_coverage(
        self,
        records: Dict[str, EffectorRecord],
        bam_files: List[str],
        chrom_lengths: Dict[str, int],
    ) -> Dict[str, GeneCoverage]:
        """计算所有效应子的覆盖度|Compute coverage for all target genes

        Args:
            records: 基因注释记录|Gene annotation records
            bam_files: BAM文件列表|BAM file list
            chrom_lengths: 染色体长度|Chromosome lengths

        Returns:
            基因ID到GeneCoverage的映射|Gene ID to GeneCoverage mapping
        """
        cov_dir = os.path.join(self.config.output_dir, '02_coverage')
        os.makedirs(cov_dir, exist_ok=True)

        # 生成BED文件|Generate BED files
        exon_bed = os.path.join(cov_dir, 'effector_exons.bed')
        gene_bed = os.path.join(cov_dir, 'effector_genes.bed')
        up_bed = os.path.join(cov_dir, 'effector_upstream.bed')
        dn_bed = os.path.join(cov_dir, 'effector_downstream.bed')

        build_exon_bed(records, exon_bed)
        build_gene_bed(records, gene_bed)
        build_flanking_bed(records, up_bed, self.config.flanking_window,
                           chrom_lengths, 'upstream')
        build_flanking_bed(records, dn_bed, self.config.flanking_window,
                           chrom_lengths, 'downstream')

        self.logger.info(f"BED文件已生成|BED files generated in {cov_dir}")

        strand_flag = ['-s'] if self.config.strandness != 'unstranded' else []

        # 对每个BAM运行bedtools coverage|Run bedtools coverage for each BAM
        exon_cov_data = {}   # gene_id -> list of (breadth, depth) per exon
        gene_cov_data = {}   # gene_id -> (breadth, depth)
        up_cov_data = {}     # gene_id -> depth
        dn_cov_data = {}     # gene_id -> depth

        for bam_file in bam_files:
            sample_name = os.path.splitext(os.path.basename(bam_file))[0]
            self.logger.info(f"计算覆盖度|Computing coverage: {sample_name}")

            # Exon级别覆盖|Exon-level coverage
            exon_sample = self._run_bedtools_coverage(
                exon_bed, bam_file, strand_flag, f"Exon覆盖度|Exon coverage: {sample_name}"
            )
            for gene_id, breadth, depth in exon_sample:
                exon_cov_data.setdefault(gene_id, []).append((breadth, depth))

            # Gene级别覆盖（合并所有exon）|Gene-level merged coverage
            gene_sample = self._run_bedtools_coverage(
                gene_bed, bam_file, strand_flag, f"Gene覆盖度|Gene coverage: {sample_name}"
            )
            for gene_id, breadth, depth in gene_sample:
                # 多样本取最大值|Take max across samples
                if gene_id not in gene_cov_data or breadth > gene_cov_data[gene_id][0]:
                    gene_cov_data[gene_id] = (breadth, depth)

            # 上游覆盖|Upstream coverage
            up_sample = self._run_bedtools_coverage(
                up_bed, bam_file, strand_flag, f"上游覆盖|Upstream coverage: {sample_name}",
                mean_only=True,
            )
            for gene_id, breadth, depth in up_sample:
                if gene_id not in up_cov_data or depth > up_cov_data[gene_id]:
                    up_cov_data[gene_id] = depth

            # 下游覆盖|Downstream coverage
            dn_sample = self._run_bedtools_coverage(
                dn_bed, bam_file, strand_flag, f"下游覆盖|Downstream coverage: {sample_name}",
                mean_only=True,
            )
            for gene_id, breadth, depth in dn_sample:
                if gene_id not in dn_cov_data or depth > dn_cov_data[gene_id]:
                    dn_cov_data[gene_id] = depth

        # 组装结果|Assemble results
        results = {}
        for gene_id, rec in records.items():
            exon_covs = []
            for idx, exon in enumerate(rec.exons):
                exon_key = gene_id
                data_points = exon_cov_data.get(exon_key, [])
                if data_points:
                    # 取最大breadth和最大depth across samples|Take max across samples
                    max_breadth = max(b for b, d in data_points)
                    max_depth = max(d for b, d in data_points)
                else:
                    max_breadth = 0.0
                    max_depth = 0.0
                exon_covs.append(ExonCoverage(
                    gene_id=gene_id,
                    exon_index=idx,
                    start=exon.start,
                    end=exon.end,
                    breadth=max_breadth,
                    mean_depth=max_depth,
                ))

            g_breadth, g_depth = gene_cov_data.get(gene_id, (0.0, 0.0))
            up_depth = up_cov_data.get(gene_id, 0.0)
            dn_depth = dn_cov_data.get(gene_id, 0.0)

            min_exon_breadth = min((e.breadth for e in exon_covs), default=0.0)
            min_exon_depth = min((e.mean_depth for e in exon_covs), default=0.0)

            results[gene_id] = GeneCoverage(
                gene_id=gene_id,
                overall_breadth=g_breadth,
                overall_mean_depth=g_depth,
                exons=exon_covs,
                upstream_mean_depth=up_depth,
                downstream_mean_depth=dn_depth,
                min_exon_breadth=min_exon_breadth,
                min_exon_depth=min_exon_depth,
            )

            # 零覆盖快速标记|Quick flag zero coverage
            if g_breadth < 0.001 and up_depth < 0.5 and dn_depth < 0.5:
                self.logger.info(
                    f"基因 {gene_id} 覆盖度接近零|Gene {gene_id} has near-zero coverage"
                )

        # 保存中间结果|Save intermediate results
        self._save_summary(results, os.path.join(cov_dir, 'coverage_summary.tsv'))
        self._save_per_exon(results, os.path.join(cov_dir, 'coverage_per_exon.tsv'))

        self.logger.info(f"覆盖度分析完成，共 {len(results)} 个基因|Coverage analysis completed for {len(results)} genes")
        return results

    def _run_bedtools_coverage(
        self,
        bed_file: str,
        bam_file: str,
        strand_flags: list,
        description: str,
        mean_only: bool = False,
    ) -> List[Tuple[str, float, float]]:
        """运行bedtools coverage|Run bedtools coverage

        Returns:
            [(name, breadth, mean_depth), ...]
        """
        if not os.path.exists(bed_file):
            return []

        args = ['coverage', '-a', bed_file, '-b', bam_file, '-mean']
        args.extend(strand_flags)

        cmd = build_conda_command(self.config.bedtools_path, args)
        success, stdout, stderr = self.cmd_runner.run(cmd, description)

        if not success:
            self.logger.warning(f"bedtools coverage失败|bedtools coverage failed: {description}")
            return []

        results = []
        for line in stdout.strip().split('\n'):
            if not line:
                continue
            name, breadth, depth = _parse_bedtools_coverage_line(line)
            results.append((name, breadth, depth))
        return results

    def _save_summary(self, coverage_data: Dict[str, GeneCoverage], output_path: str):
        """保存覆盖度汇总|Save coverage summary"""
        with open(output_path, 'w') as f:
            f.write("gene_id\toverall_breadth\toverall_mean_depth\t"
                    "min_exon_breadth\tupstream_mean_depth\tdownstream_mean_depth\n")
            for gid, cov in sorted(coverage_data.items()):
                f.write(f"{gid}\t{cov.overall_breadth:.4f}\t{cov.overall_mean_depth:.2f}\t"
                        f"{cov.min_exon_breadth:.4f}\t{cov.upstream_mean_depth:.2f}\t"
                        f"{cov.downstream_mean_depth:.2f}\n")

    def _save_per_exon(self, coverage_data: Dict[str, GeneCoverage], output_path: str):
        """保存per-exon覆盖度|Save per-exon coverage"""
        with open(output_path, 'w') as f:
            f.write("gene_id\texon_index\tstart\tend\tbreadth\tmean_depth\n")
            for gid, cov in sorted(coverage_data.items()):
                for exon in cov.exons:
                    f.write(f"{gid}\t{exon.exon_index}\t{exon.start}\t{exon.end}\t"
                            f"{exon.breadth:.4f}\t{exon.mean_depth:.2f}\n")
