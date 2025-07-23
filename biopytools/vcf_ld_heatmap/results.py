"""
VCF LD热图结果处理模块 | VCF LD Heatmap Results Processing Module
"""

import os
from pathlib import Path

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self):
        """生成总结报告 | Generate summary report"""
        report_file = os.path.join(self.config.output_dir, "ld_analysis_summary.txt")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("VCF LD热图分析总结报告 | VCF LD Heatmap Analysis Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            # 输入文件信息 | Input file information
            f.write("输入文件 | Input Files:\n")
            f.write(f"  - VCF文件 | VCF file: {self.config.vcf_file}\n")
            if self.config.region:
                f.write(f"  - 分析区域 | Analysis region: {self.config.region}\n")
            if self.config.samples:
                f.write(f"  - 指定样本数 | Specified samples: {len(self.config.samples)}\n")
            f.write("\n")
            
            # 分析参数 | Analysis parameters
            f.write("分析参数 | Analysis Parameters:\n")
            f.write(f"  - MAF阈值 | MAF threshold: {self.config.maf}\n")
            f.write(f"  - 最大SNP数 | Maximum SNPs: {self.config.max_snps}\n")
            f.write(f"  - LD阈值 | LD threshold: {self.config.ld_threshold}\n")
            f.write(f"  - 图形尺寸 | Figure size: {self.config.figsize}\n")
            f.write(f"  - DPI: {self.config.dpi}\n")
            f.write(f"  - 颜色映射 | Color map: {self.config.colormap}\n")
            f.write(f"  - 仅显示上三角 | Triangle only: {'是 | Yes' if self.config.triangle_only else '否 | No'}\n")
            f.write("\n")
            
            # 输出文件 | Output files
            f.write("输出文件 | Output Files:\n")
            f.write(f"  - 热图文件 | Heatmap file: {self.config.output_file}\n")
            if self.config.save_matrix:
                f.write(f"  - LD矩阵文件 | LD matrix file: {self.config.save_matrix}\n")
            f.write(f"  - 日志文件 | Log file: ld_analysis.log\n")
            f.write(f"  - 总结报告 | Summary report: {report_file}\n")
        
        self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
