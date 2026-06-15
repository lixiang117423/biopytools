"""
Kraken2核心分析模块|Kraken2 Core Analysis Module
"""

import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .config import Kraken2Config
from .utils import (
    CommandRunner,
    PairFinder,
    build_conda_command,
    format_number,
    parse_kraken2_report,
)


class Kraken2Analyzer:
    """Kraken2分析器|Kraken2 Analyzer"""

    def __init__(self, config: Kraken2Config, logger):
        self.config = config
        self.logger = logger
        self.cmd_runner = CommandRunner(logger)

        # 输出子目录|Output subdirectories
        self.kraken2_dir = os.path.join(self.config.output_dir, '01_kraken2')
        self.bracken_dir = os.path.join(self.config.output_dir, '02_bracken')
        self.summary_dir = os.path.join(self.config.output_dir, '03_summary')

        for d in (self.kraken2_dir, self.bracken_dir, self.summary_dir):
            os.makedirs(d, exist_ok=True)

        self._pairs: Optional[Dict[str, Tuple[str, str]]] = None

    def _find_pairs(self) -> Dict[str, Tuple[str, str]]:
        """检测配对FASTQ文件|Detect paired FASTQ files"""
        self.logger.info(
            f"扫描输入目录|Scanning input directory: {self.config.input_dir}"
        )

        finder = PairFinder(
            self.config.input_dir,
            self.config.r1_suffix,
            self.config.r2_suffix
        )
        pairs = finder.find_pairs()

        if not pairs:
            self.logger.error(
                "未找到配对的FASTQ文件|No paired FASTQ files found"
            )
            self.logger.error(
                f"请确认文件使用后缀|Please confirm files use suffixes: "
                f"{self.config.r1_suffix} / {self.config.r2_suffix}"
            )
            sys.exit(1)

        self.logger.info(
            f"找到|Found {len(pairs)} 个配对样本|paired samples"
        )
        for name in sorted(pairs.keys()):
            r1, r2 = pairs[name]
            self.logger.info(f"  {name}: {os.path.basename(r1)} + {os.path.basename(r2)}")

        return pairs

    def _is_step_completed(self, output_file: str) -> bool:
        """检查步骤是否已完成|Check if step is completed"""
        return os.path.exists(output_file)

    def _run_kraken2(self, sample_name: str, r1: str, r2: str) -> bool:
        """运行Kraken2分类|Run Kraken2 classification"""
        kraken_output = os.path.join(
            self.kraken2_dir, f"{sample_name}.kraken"
        )
        report_output = os.path.join(
            self.kraken2_dir, f"{sample_name}.kraken_report.txt"
        )

        if (self._is_step_completed(kraken_output)
                and self._is_step_completed(report_output)):
            self.logger.info(
                f"跳过已完成样本|Skipping completed sample: {sample_name}"
            )
            return True

        self.logger.info(
            f"开始Kraken2分类|Starting Kraken2 classification: {sample_name}"
        )

        args = [
            '--db', self.config.db_path,
            '--threads', str(self.config.threads),
            '--paired', r1, r2,
            '--gzip-compressed',
            '--output', kraken_output,
            '--report', report_output,
            '--use-names',
        ]

        if self.config.confidence > 0:
            args.extend(['--confidence', str(self.config.confidence)])

        cmd = build_conda_command('kraken2', args)

        success, stdout, stderr = self.cmd_runner.run(
            cmd,
            f"Kraken2分类|Kraken2 classification: {sample_name}"
        )

        if not success:
            self.logger.error(
                f"Kraken2分类失败|Kraken2 classification failed: {sample_name}"
            )
            return False

        if stderr:
            self.logger.info(f"Kraken2输出|Kraken2 output: {stderr.strip()}")

        return True

    def _run_bracken(self, sample_name: str) -> bool:
        """运行Bracken丰度估算|Run Bracken abundance estimation"""
        if not self.config.run_bracken:
            return True

        bracken_report = os.path.join(
            self.kraken2_dir, f"{sample_name}.kraken_report.txt"
        )
        bracken_output = os.path.join(
            self.bracken_dir, f"{sample_name}.bracken.txt"
        )
        bracken_new_report = os.path.join(
            self.bracken_dir, f"{sample_name}.bracken_report.txt"
        )

        if (self._is_step_completed(bracken_output)
                and self._is_step_completed(bracken_new_report)):
            self.logger.info(
                f"跳过已完成Bracken|Skipping completed Bracken: {sample_name}"
            )
            return True

        if not self._is_step_completed(bracken_report):
            self.logger.warning(
                f"跳过Bracken: Kraken2报告不存在|"
                f"Skipping Bracken: Kraken2 report not found for {sample_name}"
            )
            return False

        self.logger.info(
            f"开始Bracken丰度估算|Starting Bracken estimation: {sample_name}"
        )

        args = [
            '-d', self.config.db_path,
            '-i', bracken_report,
            '-o', bracken_output,
            '-w', bracken_new_report,
            '-r', str(self.config.read_len),
            '-l', self.config.bracken_level,
            '-t', str(self.config.bracken_threshold),
        ]

        cmd = build_conda_command('bracken', args)

        success, stdout, stderr = self.cmd_runner.run(
            cmd,
            f"Bracken丰度估算|Bracken estimation: {sample_name}"
        )

        if not success:
            self.logger.error(
                f"Bracken丰度估算失败|Bracken estimation failed: {sample_name}"
            )
            return False

        return True

    def _generate_summary(self) -> bool:
        """生成汇总报告|Generate summary report"""
        self.logger.info("开始生成汇总报告|Starting summary report generation")

        summary_file = os.path.join(
            self.summary_dir, 'kraken2_summary.tsv'
        )
        species_file = os.path.join(
            self.summary_dir, 'bracken_species.tsv'
        )

        summary_lines = [
            'Sample\tTotal_Reads\tClassified\tClassified(%)\t'
            'Unclassified\tUnclassified(%)\tTop_Species'
        ]

        all_species = {}

        for sample_name in sorted(self._pairs.keys()):
            report = os.path.join(
                self.kraken2_dir, f"{sample_name}.kraken_report.txt"
            )

            if not os.path.exists(report):
                self.logger.warning(
                    f"跳过汇总: 报告文件不存在|"
                    f"Skipping summary: report not found for {sample_name}"
                )
                continue

            parsed = parse_kraken2_report(report)

            top_species = 'N/A'
            if parsed['species']:
                top_species = parsed['species'][0]['name']

            summary_lines.append(
                f"{sample_name}\t"
                f"{format_number(parsed['total_reads'])}\t"
                f"{format_number(parsed['classified_reads'])}\t"
                f"{parsed['classified_pct']:.2f}%\t"
                f"{format_number(parsed['unclassified_reads'])}\t"
                f"{parsed['unclassified_pct']:.2f}%\t"
                f"{top_species}"
            )

            # 收集Bracken结果|Collect Bracken results
            bracken_file = os.path.join(
                self.bracken_dir, f"{sample_name}.bracken.txt"
            )
            if os.path.exists(bracken_file):
                all_species[sample_name] = self._parse_bracken(bracken_file)

        with open(summary_file, 'w') as f:
            f.write('\n'.join(summary_lines) + '\n')

        self.logger.info(f"汇总报告已保存|Summary report saved: {summary_file}")

        # 生成Bracken物种汇总|Generate Bracken species summary
        if all_species:
            self._write_species_summary(all_species, species_file)

        return True

    def _parse_bracken(self, bracken_file: str) -> List[dict]:
        """解析Bracken输出文件|Parse Bracken output file"""
        species_list = []
        with open(bracken_file, 'r') as f:
            header = f.readline()
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 6:
                    species_list.append({
                        'name': fields[0],
                        'taxid': fields[1],
                        'kraken_assigned': int(fields[2]),
                        'bracken_assigned': int(fields[3]),
                        'new_est_reads': int(fields[4]),
                        'fraction_total_reads': float(fields[5]),
                    })
        return species_list

    def _write_species_summary(
        self, all_species: Dict[str, List[dict]], output_file: str
    ):
        """写入物种丰度汇总表|Write species abundance summary table"""
        all_names = set()
        for entries in all_species.values():
            for entry in entries:
                all_names.add(entry['name'])

        sorted_names = sorted(all_names)

        with open(output_file, 'w') as f:
            header = 'Species\t' + '\t'.join(sorted(
                all_species.keys()
            ))
            f.write(header + '\n')

            for name in sorted_names:
                row = [name]
                for sample_name in sorted(all_species.keys()):
                    val = 0
                    for entry in all_species[sample_name]:
                        if entry['name'] == name:
                            val = entry['fraction_total_reads']
                            break
                    row.append(f"{val:.6f}")
                f.write('\t'.join(row) + '\n')

        self.logger.info(
            f"物种丰度汇总已保存|Species summary saved: {output_file}"
        )

    def run_pipeline(self) -> bool:
        """运行完整分析流程|Run complete analysis pipeline"""
        self.logger.info("开始Kraken2分析流程|Starting Kraken2 analysis pipeline")

        # 步骤1: 检测配对FASTQ|Step 1: Detect paired FASTQ files
        self._pairs = self._find_pairs()

        failed_samples = []

        # 步骤2: 逐样本运行Kraken2|Step 2: Run Kraken2 per sample
        for i, (sample_name, (r1, r2)) in enumerate(
            sorted(self._pairs.items()), 1
        ):
            self.logger.info(
                f"处理样本 [{i}/{len(self._pairs)}]|"
                f"Processing sample [{i}/{len(self._pairs)}]: {sample_name}"
            )

            if not self._run_kraken2(sample_name, r1, r2):
                failed_samples.append(sample_name)
                continue

            if not self._run_bracken(sample_name):
                self.logger.warning(
                    f"Bracken失败，Kraken2结果仍可用|"
                    f"Bracken failed, Kraken2 results still available: "
                    f"{sample_name}"
                )

        # 步骤3: 生成汇总报告|Step 3: Generate summary
        self._generate_summary()

        # 汇总结果|Summary results
        total = len(self._pairs)
        success = total - len(failed_samples)

        self.logger.info(
            f"Kraken2分析完成|Kraken2 analysis completed: "
            f"{success}/{total} 成功|succeeded"
        )

        if failed_samples:
            self.logger.warning(
                f"失败样本|Failed samples: {', '.join(failed_samples)}"
            )

        return len(failed_samples) == 0
