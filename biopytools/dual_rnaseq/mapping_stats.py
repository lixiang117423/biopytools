"""
互作转录组比对统计模块|Dual RNA-seq Alignment Statistics Module

统计每个样品中reads精确比对到两个基因组上的数量和百分比
Count and percentage of reads precisely aligned to two genomes for each sample
"""

import os
import pysam
import pandas as pd
from typing import Dict, List, Tuple


class MappingStatistics:
    """比对统计器|Alignment Statistics Calculator"""

    def __init__(self, config, logger):
        """
        初始化比对统计器|Initialize alignment statistics calculator

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger

    def calculate_sample_stats(self, species1_bam: str, species2_bam: str,
                               sample_name: str) -> Dict:
        """
        计算单个样本的比对统计|Calculate alignment statistics for a single sample

        Args:
            species1_bam: 物种1的BAM文件|Species 1 BAM file
            species2_bam: 物种2的BAM文件|Species 2 BAM file
            sample_name: 样本名称|Sample name

        Returns:
            统计结果字典|Statistics result dictionary
        """
        self.logger.info(f"计算样本比对统计|Calculating sample alignment statistics: {sample_name}")

        # 读取BAM文件|Read BAM files
        bam1 = pysam.AlignmentFile(species1_bam, "rb")
        bam2 = pysam.AlignmentFile(species2_bam, "rb")

        # 统计变量|Statistics variables
        total_reads = 0
        species1_only = 0  # 仅比对到物种1|Only aligned to species 1
        species2_only = 0  # 仅比对到物种2|Only aligned to species 2
        both_mapped = 0     # 同时比对到两个物种|Aligned to both species
        neither_mapped = 0  # 两个都没比对上|Aligned to neither

        # 存储reads信息|Store reads information
        reads1 = {}  # read_name -> (is_mapped, mapq)
        reads2 = {}

        # 读取物种1的reads|Read species 1 reads
        for read in bam1.fetch():
            read_name = read.query_name
            is_mapped = not read.is_unmapped and not read.is_secondary and not read.is_supplementary
            mapq = read.mapping_quality if is_mapped else 0

            if read_name not in reads1:
                reads1[read_name] = (is_mapped, mapq)

        # 读取物种2的reads|Read species 2 reads
        for read in bam2.fetch():
            read_name = read.query_name
            is_mapped = not read.is_unmapped and not read.is_secondary and not read.is_supplementary
            mapq = read.mapping_quality if is_mapped else 0

            if read_name not in reads2:
                reads2[read_name] = (is_mapped, mapq)

        # 获取所有唯一的reads|Get all unique reads
        all_reads = set(reads1.keys()) | set(reads2.keys())
        total_reads = len(all_reads)

        # 统计比对情况|Count alignment status
        for read_name in all_reads:
            in_species1 = read_name in reads1
            in_species2 = read_name in reads2

            if not in_species1 and not in_species2:
                # 这个分支理论上不会执行，因为reads至少来自一个文件
                # This branch should theoretically not execute as reads come from at least one file
                neither_mapped += 1
            elif in_species1 and not in_species2:
                # 仅在物种1|Only in species 1
                is_mapped, mapq = reads1[read_name]
                if is_mapped and self._is_high_quality(mapq):
                    species1_only += 1
            elif in_species2 and not in_species1:
                # 仅在物种2|Only in species 2
                is_mapped, mapq = reads2[read_name]
                if is_mapped and self._is_high_quality(mapq):
                    species2_only += 1
            else:
                # 两个都在|In both
                is_mapped1, mapq1 = reads1[read_name]
                is_mapped2, mapq2 = reads2[read_name]

                species1_high_quality = is_mapped1 and self._is_high_quality(mapq1)
                species2_high_quality = is_mapped2 and self._is_high_quality(mapq2)

                if species1_high_quality and species2_high_quality:
                    both_mapped += 1
                elif species1_high_quality:
                    species1_only += 1
                elif species2_high_quality:
                    species2_only += 1

        # 关闭文件|Close files
        bam1.close()
        bam2.close()

        # 计算百分比|Calculate percentages
        total_mapped = species1_only + species2_only + both_mapped
        neither_mapped = total_reads - total_mapped

        stats = {
            'sample_name': sample_name,
            'total_reads': total_reads,
            'species1_only': species1_only,
            'species2_only': species2_only,
            'both_mapped': both_mapped,
            'neither_mapped': neither_mapped,
            'total_mapped': total_mapped,
            'species1_only_pct': (species1_only / total_reads * 100) if total_reads > 0 else 0,
            'species2_only_pct': (species2_only / total_reads * 100) if total_reads > 0 else 0,
            'both_mapped_pct': (both_mapped / total_reads * 100) if total_reads > 0 else 0,
            'neither_mapped_pct': (neither_mapped / total_reads * 100) if total_reads > 0 else 0,
            'total_mapped_pct': (total_mapped / total_reads * 100) if total_reads > 0 else 0,
        }

        # 输出统计信息|Output statistics
        self.logger.info(f"样本|Sample: {sample_name}")
        self.logger.info(f"  总reads数|Total reads: {stats['total_reads']:,}")
        self.logger.info(f"  仅比对到{self.config.species1_name}|Only aligned to {self.config.species1_name}: {stats['species1_only']:,} ({stats['species1_only_pct']:.2f}%)")
        self.logger.info(f"  仅比对到{self.config.species2_name}|Only aligned to {self.config.species2_name}: {stats['species2_only']:,} ({stats['species2_only_pct']:.2f}%)")
        self.logger.info(f"  同时比对到两者|Aligned to both: {stats['both_mapped']:,} ({stats['both_mapped_pct']:.2f}%)")
        self.logger.info(f"  未比对上|Not aligned: {stats['neither_mapped']:,} ({stats['neither_mapped_pct']:.2f}%)")
        self.logger.info(f"  总比对数|Total aligned: {stats['total_mapped']:,} ({stats['total_mapped_pct']:.2f}%)")

        return stats

    def _is_high_quality(self, mapq: int) -> bool:
        """
        判断MAPQ是否高质量|Check if MAPQ is high quality

        Args:
            mapq: Mapping quality score

        Returns:
            是否高质量|Whether high quality
        """
        if self.config.unique_only:
            return mapq >= self.config.min_mapq
        else:
            return mapq >= 0

    def generate_mapping_summary(self, sample_stats: List[Dict], output_dir: str) -> bool:
        """
        生成比对统计汇总报告|Generate alignment statistics summary report

        Args:
            sample_stats: 所有样本的统计结果列表|List of statistics for all samples
            output_dir: 输出目录|Output directory

        Returns:
            是否成功|Whether successful
        """
        try:
            # 创建DataFrame|Create DataFrame
            df = pd.DataFrame(sample_stats)

            # 重新排列列顺序|Reorder columns
            columns = [
                'sample_name',
                'total_reads',
                'species1_only',
                'species1_only_pct',
                'species2_only',
                'species2_only_pct',
                'both_mapped',
                'both_mapped_pct',
                'neither_mapped',
                'neither_mapped_pct',
                'total_mapped',
                'total_mapped_pct'
            ]

            df = df[columns]

            # 重命名列|Rename columns
            df.columns = [
                'Sample Name|样本名称',
                f'Total Reads|总reads数',
                f'Only {self.config.species1_name}|仅{self.config.species1_name}',
                f'% Only {self.config.species1_name}|%仅{self.config.species1_name}',
                f'Only {self.config.species2_name}|仅{self.config.species2_name}',
                f'% Only {self.config.species2_name}|%仅{self.config.species2_name}',
                'Both Mapped|同时比对',
                '% Both Mapped|%同时比对',
                'Not Mapped|未比对',
                '% Not Mapped|%未比对',
                'Total Mapped|总比对数',
                '% Total Mapped|%总比对'
            ]

            # 保存为TSV|Save as TSV
            output_file = os.path.join(output_dir, 'mapping_statistics.tsv')
            df.to_csv(output_file, sep='\t', index=False, float_format='%.2f')

            self.logger.info(f"比对统计报告已保存|Mapping statistics report saved: {output_file}")

            # 同时生成纯英文版本的报告（方便后续分析）|Also generate English-only report
            output_file_en = os.path.join(output_dir, 'mapping_statistics_en.tsv')
            df_en = df.copy()
            df_en.columns = [
                'Sample',
                'Total_Reads',
                f'Only_{self.config.species1_name}',
                f'Pct_Only_{self.config.species1_name}',
                f'Only_{self.config.species2_name}',
                f'Pct_Only_{self.config.species2_name}',
                'Both_Mapped',
                'Pct_Both_Mapped',
                'Not_Mapped',
                'Pct_Not_Mapped',
                'Total_Mapped',
                'Pct_Total_Mapped'
            ]
            df_en.to_csv(output_file_en, sep='\t', index=False, float_format='%.2f')

            self.logger.info(f"英文版比对统计报告已保存|English mapping statistics report saved: {output_file_en}")

            return True

        except Exception as e:
            self.logger.error(f"生成比对统计报告时出错|Error generating mapping statistics report: {e}")
            return False

    def process_all_samples(self, bam_files: List[Tuple[str, str, str]], output_dir: str) -> bool:
        """
        处理所有样本的比对统计|Process mapping statistics for all samples

        Args:
            bam_files: BAM文件列表，每个元素为(sample_name, species1_bam, species2_bam)|List of BAM files, each element is (sample_name, species1_bam, species2_bam)
            output_dir: 输出目录|Output directory

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始计算所有样本的比对统计|Starting to calculate mapping statistics for all samples")
        self.logger.info("=" * 60)

        all_stats = []

        for sample_name, species1_bam, species2_bam in bam_files:
            if not os.path.exists(species1_bam):
                self.logger.warning(f"物种1 BAM文件不存在|Species 1 BAM file does not exist: {species1_bam}")
                continue

            if not os.path.exists(species2_bam):
                self.logger.warning(f"物种2 BAM文件不存在|Species 2 BAM file does not exist: {species2_bam}")
                continue

            stats = self.calculate_sample_stats(species1_bam, species2_bam, sample_name)
            all_stats.append(stats)

        if not all_stats:
            self.logger.error("没有可用的统计结果|No valid statistics results")
            return False

        # 生成汇总报告|Generate summary report
        success = self.generate_mapping_summary(all_stats, output_dir)

        if success:
            self.logger.info("=" * 60)
            self.logger.info("所有样本的比对统计计算完成|Mapping statistics calculation completed for all samples")
            self.logger.info("=" * 60)

        return success
