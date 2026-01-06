"""
ANNOVAR结果处理模块|ANNOVAR Results Processing Module
"""

import os
from pathlib import Path

class SummaryGenerator:
    """总结生成器|Summary Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_summary_report(self):
        """生成总结报告|Generate summary report"""
        report_file = os.path.join(self.config.output_dir, "annotation_summary.txt")

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("ANNOVAR注释总结报告|ANNOVAR Annotation Summary Report\n")
            f.write("=" * 60 + "\n\n")

            # 输入文件信息|Input file information
            f.write("输入文件|Input Files:\n")
            f.write(f"  GFF3文件|GFF3 file: {self.config.gff3_file}\n")
            f.write(f"  基因组文件|Genome file: {self.config.genome_file}\n")
            f.write(f"  VCF文件|VCF file: {self.config.vcf_file}\n\n")

            # 输出文件信息|Output file information
            f.write("输出文件|Output Files:\n")
            if hasattr(self.config, 'output_files'):
                for file in self.config.output_files:
                    f.write(f"  {file}\n")
            else:
                f.write("  无输出文件记录|No output files recorded\n")
            f.write("\n")

            # 配置参数|Configuration parameters
            f.write("配置参数|Configuration Parameters:\n")
            f.write(f"  质量阈值|Quality threshold: {self.config.qual_threshold}\n")
            f.write(f"  构建版本|Build version: {self.config.build_ver}\n")
            f.write(f"  输出目录|Output directory: {self.config.output_dir}\n")
            f.write(f"  ANNOVAR路径|ANNOVAR path: {self.config.annovar_path}\n")
            f.write(f"  数据库路径|Database path: {self.config.database_path}\n")
            f.write(f"  跳过VCF过滤|Skip VCF filter: {'是|Yes' if self.config.skip_vcf_filter else '否|No'}\n")
            f.write(f"  跳过GFF修复|Skip GFF fix: {'是|Yes' if self.config.skip_gff_fix else '否|No'}\n")

            # 执行步骤|Executed steps
            if self.config.step:
                f.write(f"  执行步骤|Executed step: {self.config.step}\n")
            else:
                f.write(f"  执行步骤|Executed steps: 完整流程|Full pipeline\n")

            f.write(f"\n注释流程完成时间|Annotation completed at: {self._get_current_time()}\n")

        self.logger.info(f"总结报告已生成|Summary report generated: {report_file}")

    def _get_current_time(self):
        """获取当前时间|Get current time"""
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
