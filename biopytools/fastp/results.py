"""
FASTP结果处理模块|FASTP Results Processing Module
"""

import os
from pathlib import Path


class SummaryGenerator:
    """总结生成器|Summary Generator"""

    def __init__(self, config, logger):
        """
        初始化总结生成器|Initialize summary generator

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger

    def generate_summary_report(self, successful_count: int, failed_count: int, total_samples: int):
        """
        生成总结报告|Generate summary report

        Args:
            successful_count: 成功处理的样本数|Number of successfully processed samples
            failed_count: 失败的样本数|Number of failed samples
            total_samples: 总样本数|Total number of samples

        Returns:
            生成是否成功|Whether generation succeeded
        """
        report_file = self.config.output_path / "fastp_processing_summary.txt"

        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("FASTP质控处理总结报告|FASTP Quality Control Processing Summary Report\n")
                f.write("=" * 70 + "\n\n")

                # 处理统计|Processing statistics
                f.write("处理统计|Processing Statistics:\n")
                f.write(f"  - 总样本数|Total samples: {total_samples:,}\n")
                f.write(f"  - 成功处理|Successfully processed: {successful_count:,}\n")
                f.write(f"  - 失败样本|Failed samples: {failed_count:,}\n")
                f.write(f"  - 成功率|Success rate: {(successful_count/total_samples)*100:.1f}%\n\n")

                # 输入输出信息|Input/Output information
                f.write("文件路径|File Paths:\n")
                f.write(f"  - 输入目录|Input directory: {self.config.input_dir}\n")
                f.write(f"  - 输出目录|Output directory: {self.config.output_dir}\n")
                f.write(f"  - 报告目录|Report directory: {self.config.report_path}\n\n")

                # 处理参数|Processing parameters
                f.write("处理参数|Processing Parameters:\n")
                f.write(f"  - Fastp路径|Fastp path: {self.config.fastp_path}\n")
                f.write(f"  - 线程数|Thread count: {self.config.threads:,}\n")
                f.write(f"  - 质量阈值|Quality threshold: {self.config.quality_threshold}\n")
                f.write(f"  - 最小长度|Minimum length: {self.config.min_length}\n")
                f.write(f"  - 不合格百分比|Unqualified percentage: {self.config.unqualified_percent}\n")
                f.write(f"  - N碱基限制|N base limit: {self.config.n_base_limit}\n")
                f.write(f"  - Read1后缀|Read1 suffix: {self.config.read1_suffix}\n")
                f.write(f"  - Read2后缀|Read2 suffix: {self.config.read2_suffix}\n\n")

                # 输出文件说明|Output file description
                f.write("输出文件说明|Output File Description:\n")
                f.write("  - *_1.clean.fq.gz: 质控后的Read1文件|Quality-controlled Read1 files\n")
                f.write("  - *_2.clean.fq.gz: 质控后的Read2文件|Quality-controlled Read2 files\n")
                f.write("  - fastp_reports/*.html: 质控报告HTML文件|Quality control HTML reports\n")
                f.write("  - fastp_reports/*.json: 质控报告JSON文件|Quality control JSON reports\n")

            self.logger.info(f"总结报告已生成|Summary report generated: {report_file}")
            return True

        except Exception as e:
            self.logger.error(f"生成总结报告时出错|Error generating summary report: {e}")
            return False
