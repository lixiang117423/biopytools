"""
BAM文件批量统计处理器|BAM File Batch Statistics Processor

输出两个文件:
  1. {prefix}.summary.tsv        - 全局统计长表（样本为行，指标为列）
  2. {prefix}.per_chromosome.tsv  - 染色体级别统计
"""

import json
import subprocess
from pathlib import Path
from multiprocessing import Pool, cpu_count

import pandas as pd
from tqdm import tqdm

from .sample_stats import SampleStatsGenerator
from .utils import (
    CommandRunner,
    get_sample_name,
)


class BAMBatchProcessor:
    """BAM文件批量处理器|BAM File Batch Processor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def analyze_single_bam(self, bam_file: str) -> dict:
        """分析单个BAM文件|Analyze a single BAM file"""
        sample_name = get_sample_name(bam_file)
        try:
            self._ensure_bam_index(bam_file)

            cmd_runner = CommandRunner(
                self.logger, threads=self.config.threads
            )
            generator = SampleStatsGenerator(
                self.config, self.logger, cmd_runner
            )
            return generator.analyze_sample(bam_file)

        except Exception as e:
            self.logger.error(
                f"处理失败|Processing failed for {sample_name}: {e}"
            )
            return {
                'sample_name': sample_name,
                'bam_file': bam_file,
                'error': str(e),
            }

    def process_bam_files(
        self, bam_files: list, max_workers: int = None
    ) -> list:
        """批量处理BAM文件|Batch process BAM files"""
        if max_workers is None:
            max_workers = self.config.max_workers

        self.logger.info(
            f"找到{len(bam_files)}个BAM文件|"
            f"Found {len(bam_files)} BAM files"
        )
        self.logger.info(
            f"使用{max_workers}个进程并行处理|"
            f"Using {max_workers} processes for parallel processing"
        )

        all_results = []
        try:
            with Pool(processes=max_workers) as pool:
                results_iter = pool.imap_unordered(
                    self.analyze_single_bam, bam_files
                )
                all_results = list(tqdm(
                    results_iter,
                    total=len(bam_files),
                    desc="处理BAM文件|Processing BAM files",
                ))
        except Exception as e:
            self.logger.error(
                f"并行处理出错，降级到串行|"
                f"Parallel processing error, falling back to serial: {e}"
            )
            all_results = []
            for bam_file in tqdm(
                bam_files, desc="串行处理BAM文件|Serial processing"
            ):
                all_results.append(self.analyze_single_bam(bam_file))

        successful = [r for r in all_results if 'error' not in r]
        failed = [r for r in all_results if 'error' in r]

        if failed:
            self.logger.error(
                f"以下样本处理失败|The following samples failed:"
            )
            for r in failed:
                self.logger.error(
                    f"  - {r['sample_name']}: {r['error']}"
                )

        if not successful:
            self.logger.error("未能成功处理任何文件|Failed to process any files")
            return []

        self.logger.info(
            f"成功处理{len(successful)}个，失败{len(failed)}个|"
            f"Successfully processed {len(successful)}, "
            f"failed {len(failed)}"
        )
        return successful

    def generate_report(
        self, results: list, output_file: str, output_dir: str, prefix: str
    ) -> bool:
        """生成统计报告|Generate statistics reports"""
        try:
            summary_file = Path(output_file)
            chrom_file = Path(output_dir) / f"{prefix}.per_chromosome.tsv"

            summary_df = self._build_summary_df(results)
            summary_df.to_csv(summary_file, sep='\t', index=False)
            self.logger.info(
                f"全局统计报告已保存|Summary report saved: {summary_file}"
            )

            chrom_df = self._build_chromosome_df(results)
            if not chrom_df.empty:
                chrom_df.to_csv(chrom_file, sep='\t', index=False)
                self.logger.info(
                    f"染色体统计报告已保存|"
                    f"Chromosome report saved: {chrom_file}"
                )

            genome_json = Path(output_dir) / f"{prefix}.genome_stats.json"
            self._save_genome_stats(results, genome_json)

            return True

        except Exception as e:
            self.logger.error(f"生成报告失败|Failed to generate report: {e}")
            return False

    def _build_summary_df(self, results: list) -> pd.DataFrame:
        """构建全局统计长表|Build summary tidy DataFrame"""
        rows = []

        for sample in results:
            row = {'Sample': sample['sample_name']}

            alignment = sample.get('alignment_stats', {})
            if alignment:
                row['Total_Reads'] = alignment.get('total_reads')
                row['Mapped_Reads'] = alignment.get('mapped_reads')
                row['Unmapped_Reads'] = alignment.get('unmapped_reads')
                row['Map_Rate(%)'] = alignment.get('mapping_rate')
                row['Proper_Pair(%)'] = alignment.get('proper_pair_rate')
                row['Secondary_Reads'] = alignment.get('secondary_reads')
                row['Supplementary_Reads'] = alignment.get(
                    'supplementary_reads'
                )
                row['Singletons'] = alignment.get('singletons')
                row['Unique_Mapped'] = alignment.get('unique_mapped')
                row['Multi_Mapped'] = alignment.get('multi_mapped')

            coverage = sample.get('coverage_stats', {})
            if coverage:
                row['Avg_Depth'] = coverage.get('mean_coverage')
                row['Median_Depth'] = coverage.get('median_coverage')
                row['Depth_Std'] = coverage.get('coverage_std')
                row['Depth_CV'] = coverage.get('coverage_cv')
                row['Cov_Bases_Rate(%)'] = round(
                    coverage.get('coverage_rate', 0) * 100, 2
                ) if coverage.get('coverage_rate') is not None else None

            dup = sample.get('duplicate_stats', {})
            if dup:
                row['Duplicates_Rate(%)'] = dup.get('duplicate_rate')
                row['Duplicates_Marked'] = dup.get('duplicates_marked')

            ins = sample.get('insert_stats', {})
            if ins and ins.get('insert_stats_available', False):
                row['Mean_Insert_Size'] = ins.get('mean_insert_size')
                row['Median_Insert_Size'] = ins.get('median_insert_size')
                row['Insert_Size_Std'] = ins.get('insert_size_std')

            seq = sample.get('sequence_stats', {})
            if seq:
                row['Avg_Read_Length'] = seq.get('average_read_length')
                row['Avg_GC(%)'] = seq.get('average_gc_content')

            var = sample.get('variation_stats', {})
            if var:
                row['Mismatch_Rate(%)'] = var.get('mismatch_rate')
                row['SoftClip_Rate(%)'] = var.get('softclip_rate')
                row['SoftClip_Reads_Rate(%)'] = var.get(
                    'reads_with_softclip_rate'
                )
                row['GC_Bias(%)'] = var.get('gc_bias')
                row['GC_Deviation'] = var.get('gc_deviation')

            # 污染预警|Contamination warning
            row['Contam_Warning'] = self._check_contamination(sample)

            rows.append(row)

        df = pd.DataFrame(rows)
        if not df.empty:
            df = df.sort_values('Sample').reset_index(drop=True)
        return df

    def _check_contamination(self, sample: dict) -> str:
        """检查潜在污染|Check potential contamination"""
        warnings = []

        alignment = sample.get('alignment_stats', {})
        if alignment:
            map_rate = alignment.get('mapping_rate', 100)
            if map_rate < 70:
                warnings.append(f"LowMapRate({map_rate:.1f}%)")

        coverage = sample.get('coverage_stats', {})
        if coverage:
            chrom_cov = coverage.get('chromosome_coverage', {})
            total_mapped = alignment.get('mapped_reads', 1) if alignment else 1

            # 检查非主要染色体的异常高比例
            if chrom_cov and total_mapped > 0:
                chrom_read_counts = []
                for chrom, stats in chrom_cov.items():
                    if 'mean_coverage' in stats:
                        chrom_read_counts.append(stats)

                chrom_read_counts.sort(
                    key=lambda x: x.get('mean_coverage', 0), reverse=True
                )

                # 如果染色体数量 > 8，检查top8之外的染色体
                if len(chrom_read_counts) > 8:
                    top8_mean = sum(
                        c['mean_coverage'] for c in chrom_read_counts[:8]
                    ) / 8
                    for chrom, stats in chrom_cov.items():
                        if stats not in chrom_read_counts[:8]:
                            chrom_mean = stats.get('mean_coverage', 0)
                            if chrom_mean > top8_mean * 0.5:
                                warnings.append(
                                    f"{chrom}({chrom_mean:.1f}x)"
                                )

        if warnings:
            return ';'.join(warnings)
        return 'Normal'

    def _build_chromosome_df(self, results: list) -> pd.DataFrame:
        """构建染色体级别统计表|Build per-chromosome DataFrame"""
        rows = []

        for sample in results:
            sample_name = sample['sample_name']

            coverage = sample.get('coverage_stats', {})
            if not coverage:
                continue

            chrom_cov = coverage.get('chromosome_coverage', {})
            alignment = sample.get('alignment_stats', {})
            total_mapped = alignment.get('mapped_reads', 1) if alignment else 1

            for chrom, stats in chrom_cov.items():
                row = {
                    'Sample': sample_name,
                    'Chromosome': chrom,
                    'Coverage_Depth': stats.get('mean_coverage'),
                    'Coverage_Rate(%)': round(
                        stats.get('coverage_rate', 0) * 100, 2
                    ) if stats.get('coverage_rate') is not None else None,
                    'Windows': stats.get('windows'),
                    'Covered_Windows': stats.get('covered_windows'),
                }
                rows.append(row)

        if not rows:
            return pd.DataFrame()

        df = pd.DataFrame(rows)
        return df.sort_values(['Sample', 'Chromosome']).reset_index(drop=True)

    def _save_genome_stats(self, results: list, output_file: Path):
        """保存基因组级别汇总|Save genome-level summary"""
        from .genome_stats import GenomeStatsGenerator

        generator = GenomeStatsGenerator(self.config, self.logger)
        genome_stats = generator.generate_genome_stats(results)

        with open(output_file, 'w', encoding='utf-8') as f:
            json.dump(genome_stats, f, indent=2, ensure_ascii=False)

        self.logger.info(
            f"基因组统计已保存|Genome stats saved: {output_file}"
        )

    def _ensure_bam_index(self, bam_file: str):
        """确保BAM索引文件存在|Ensure BAM index file exists"""
        bam_path = Path(bam_file)
        bai_file = bam_path.with_suffix('.bam.bai')

        if bai_file.exists():
            return True

        self.logger.info(
            f"正在创建BAM索引|Creating BAM index: {bai_file.name}"
        )
        cmd = f"{self.config.samtools_path} index {bam_file}"
        self.logger.info(f"命令|Command: {cmd}")

        try:
            subprocess.run(
                cmd, shell=True, capture_output=True, text=True, check=True,
            )
            self.logger.info(
                f"BAM索引创建成功|BAM index created: {bai_file.name}"
            )
            return True
        except (subprocess.CalledProcessError, FileNotFoundError) as e:
            self.logger.error(
                f"创建BAM索引失败|Failed to create BAM index: {e}"
            )
            return False
