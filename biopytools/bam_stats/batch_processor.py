"""
BAM文件批量统计处理器|BAM File Batch Statistics Processor
"""

import re
import subprocess
import pandas as pd
from pathlib import Path
from tqdm import tqdm
from typing import Dict, List
from multiprocessing import Pool, cpu_count


class BAMBatchProcessor:
    """BAM文件批量处理器|BAM File Batch Processor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def parse_flagstat(self, flagstat_output: str) -> Dict:
        """从samtools flagstat的输出中解析关键指标|Parse key metrics from samtools flagstat output"""
        stats = {}

        total_reads_match = re.search(r'(\d+) \+ \d+ in total', flagstat_output)
        mapped_reads_match = re.search(r'(\d+) \+ \d+ mapped', flagstat_output)
        properly_paired_match = re.search(r'(\d+) \+ \d+ properly paired', flagstat_output)

        total_reads = int(total_reads_match.group(1)) if total_reads_match else 0
        mapped_reads = int(mapped_reads_match.group(1)) if mapped_reads_match else 0
        properly_paired = int(properly_paired_match.group(1)) if properly_paired_match else 0

        stats['Total_Reads'] = total_reads
        stats['Mapped_Reads'] = mapped_reads
        stats['Properly_Paired_Reads'] = properly_paired

        if total_reads > 0:
            stats['Overall_Mapping_Rate(%)'] = round((mapped_reads / total_reads) * 100, 2)
            stats['Properly_Paired_Rate(%)'] = round((properly_paired / total_reads) * 100, 2)
        else:
            stats['Overall_Mapping_Rate(%)'] = 0.0
            stats['Properly_Paired_Rate(%)'] = 0.0

        return stats

    def analyze_single_bam(self, bam_file: Path) -> Dict:
        """对单个BAM文件执行分析|Analyze a single BAM file"""
        try:
            sample_name = bam_file.name.replace('.sorted.bam', '').replace('.bam', '')
            results = {'Sample': sample_name}

            self.logger.info(f"分析文件|Analyzing file: {bam_file.name}")

            # 检查并创建索引|Check and create index
            if not self._ensure_bam_index(bam_file):
                return {'Sample': sample_name, 'Error': 'Failed to create BAM index'}

            # 1. 运行samtools flagstat|Run samtools flagstat
            self.logger.debug(f"运行samtools flagstat|Running samtools flagstat: {bam_file.name}")
            flagstat_proc = subprocess.run(['samtools', 'flagstat', str(bam_file)],
                                        capture_output=True, text=True, check=True)
            results.update(self.parse_flagstat(flagstat_proc.stdout))

            # 2. 运行samtools coverage|Run samtools coverage
            self.logger.debug(f"运行samtools coverage|Running samtools coverage: {bam_file.name}")
            coverage_proc = subprocess.run(['samtools', 'coverage', str(bam_file)],
                                        capture_output=True, text=True, check=True)

            # 3. 运行samtools idxstats|Run samtools idxstats
            self.logger.debug(f"运行samtools idxstats|Running samtools idxstats: {bam_file.name}")
            idxstats_proc = subprocess.run(['samtools', 'idxstats', str(bam_file)],
                                        capture_output=True, text=True, check=True)

            # 解析覆盖度信息|Parse coverage information
            coverage_lines = coverage_proc.stdout.strip().splitlines()
            # 寻找最后一行非注释行作为整体统计|Find last non-comment line as overall statistics
            summary_line = None
            for line in reversed(coverage_lines):
                if not line.startswith('#'):
                    summary_line = line
                    break

            if summary_line:
                overall_cov_parts = summary_line.strip().split('\t')
                if len(overall_cov_parts) > 6:
                    results['Overall_Covered_Bases_Rate(%)'] = float(overall_cov_parts[5])
                    results['Overall_Coverage_Depth'] = float(overall_cov_parts[6])

            # 解析每条染色体的统计|Parse per-chromosome statistics
            total_mapped = results.get('Mapped_Reads', 1)
            for line in idxstats_proc.stdout.strip().splitlines():
                parts = line.split('\t')
                if len(parts) < 3:
                    continue
                chrom, _, mapped = parts[0], parts[1], parts[2]
                if chrom == '*':
                    continue
                results[f'{chrom}_Mapped_Reads'] = int(mapped)
                if total_mapped > 0:
                    results[f'{chrom}_Mapped_Reads_Percentage(%)'] = round((int(mapped) / total_mapped) * 100, 2)

            # 解析每条染色体的覆盖度|Parse per-chromosome coverage
            for line in coverage_lines:
                if line.startswith('#'):
                    continue
                parts = line.strip().split('\t')
                chrom = parts[0]
                if len(parts) > 6 and chrom != "genome" and chrom != summary_line.split('\t')[0]:  # 排除整体行
                    results[f'{chrom}_Coverage_Depth'] = float(parts[6])
                    results[f'{chrom}_Covered_Bases_Rate(%)'] = float(parts[5])

            self.logger.info(f"完成分析|Analysis completed: {bam_file.name}")
            return results

        except subprocess.CalledProcessError as e:
            error_msg = e.stderr.strip() if e.stderr else str(e)
            self.logger.error(f"处理失败|Processing failed for {bam_file.name}: {error_msg}")
            return {'Sample': bam_file.name, 'Error': error_msg}
        except Exception as e:
            self.logger.error(f"处理失败|Processing failed for {bam_file.name}: {str(e)}")
            return {'Sample': bam_file.name, 'Error': str(e)}

    def process_bam_files(self, bam_files: List[Path], processes: int = None) -> List[Dict]:
        """批量处理BAM文件|Batch process BAM files"""
        if processes is None:
            processes = cpu_count()

        self.logger.info(f"找到{len(bam_files)}个BAM文件|Found {len(bam_files)} BAM files")
        self.logger.info(f"使用{processes}个进程并行处理|Using {processes} processes for parallel processing")

        all_results = []

        try:
            with Pool(processes=processes) as pool:
                results_iterator = pool.imap_unordered(self.analyze_single_bam, bam_files)
                all_results = list(tqdm(results_iterator, total=len(bam_files),
                                       desc="处理BAM文件|Processing BAM files"))

        except Exception as e:
            self.logger.error(f"并行处理出错|Parallel processing error: {str(e)}")
            # 降级到串行处理|Fallback to serial processing
            self.logger.info("降级到串行处理|Fallback to serial processing")
            all_results = []
            for bam_file in tqdm(bam_files, desc="串行处理BAM文件|Serial processing BAM files"):
                result = self.analyze_single_bam(bam_file)
                all_results.append(result)

        successful_results = [res for res in all_results if 'Error' not in res]
        failed_samples = [res for res in all_results if 'Error' in res]

        if failed_samples:
            self.logger.error(f"\n以下样本处理失败|The following samples failed to process:")
            for res in failed_samples:
                self.logger.error(f"  - {res['Sample']}: {res['Error']}")

        if not successful_results:
            self.logger.error("\n未能成功处理任何文件，程序退出|Failed to process any files, exiting")
            return []

        self.logger.info(f"成功处理{len(successful_results)}个文件，失败{len(failed_samples)}个|Successfully processed {len(successful_results)} files, failed {len(failed_samples)}")

        return successful_results

    def generate_report(self, results: List[Dict], output_file: str) -> bool:
        """生成统计报告|Generate statistics report"""
        try:
            self.logger.info("正在汇总所有成功的结果并生成报告|Summarizing all successful results and generating report")

            df = pd.DataFrame(results)
            df = df.set_index('Sample').T

            # 重新排列统计指标的顺序|Reorder statistics indices
            core_indices = [
                'Total_Reads', 'Mapped_Reads', 'Overall_Mapping_Rate(%)',
                'Properly_Paired_Reads', 'Properly_Paired_Rate(%)',
                'Overall_Coverage_Depth', 'Overall_Covered_Bases_Rate(%)'
            ]
            chrom_indices = sorted([idx for idx in df.index if idx not in core_indices])
            final_indices = core_indices + chrom_indices
            df = df.reindex(index=final_indices)

            # 根据文件扩展名选择输出格式|Choose output format based on file extension
            if Path(output_file).suffix.lower() == '.xlsx':
                df.to_excel(output_file, index=True)
                self.logger.info(f"Excel报告已保存|Excel report saved: {output_file}")
            else:
                df.to_csv(output_file, index=True)
                self.logger.info(f"CSV报告已保存|CSV report saved: {output_file}")

            return True

        except Exception as e:
            self.logger.error(f"生成报告失败|Failed to generate report: {str(e)}")
            return False

    def _ensure_bam_index(self, bam_file: Path):
        """确保BAM索引文件存在|Ensure BAM index file exists"""
        bai_file = bam_file.with_suffix('.bam.bai')

        if bai_file.exists():
            self.logger.debug(f"BAM索引文件已存在|BAM index file already exists: {bai_file}")
            return True

        try:
            self.logger.info(f"正在创建BAM索引文件|Creating BAM index file: {bai_file}")
            result = subprocess.run(['samtools', 'index', str(bam_file)],
                                  capture_output=True, text=True, check=True)
            self.logger.info(f"BAM索引文件创建成功|BAM index file created successfully: {bai_file}")
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"创建BAM索引失败|Failed to create BAM index: {e.stderr}")
            return False
        except FileNotFoundError:
            self.logger.error("samtools命令未找到|samtools command not found")
            return False
