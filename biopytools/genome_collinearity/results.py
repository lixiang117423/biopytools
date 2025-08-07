"""
结果汇总模块 | Results Summary Module
"""

import os
from pathlib import Path
from datetime import datetime

class SummaryGenerator:
    """结果汇总生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, alignment_files: dict, syri_files: dict) -> None:
        """生成总结报告 | Generate summary report"""
        
        report_file = Path(self.config.output_dir) / "analysis_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            # 基本信息 | Basic information
            f.write("基因组共线性分析总结报告 | Genome Collinearity Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            f.write(f"分析时间 | Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"输出目录 | Output directory: {self.config.output_dir}\n\n")
            
            # 分析参数 | Analysis parameters
            f.write("分析参数 | Analysis Parameters:\n")
            f.write("-" * 40 + "\n")
            f.write(f"样本数量 | Number of samples: {len(self.config.sample_list)}\n")
            f.write(f"样本顺序 | Sample order: {' -> '.join(self.config.sample_list)}\n")
            f.write(f"线程数 | Threads: {self.config.threads}\n")
            f.write(f"Minimap2预设 | Minimap2 preset: {self.config.minimap2_preset}\n")
            if self.config.chromosome:
                f.write(f"指定染色体 | Specified chromosome: {self.config.chromosome}\n")
            else:
                f.write(f"分析范围 | Analysis scope: 全基因组 | Whole genome\n")
            f.write(f"图像格式 | Image format: {self.config.plotsr_format}\n")
            f.write(f"图像尺寸 | Image size: {self.config.figure_width}x{self.config.figure_height}\n\n")
            
            # 比对结果 | Alignment results
            f.write("基因组比对结果 | Genome Alignment Results:\n")
            f.write("-" * 40 + "\n")
            for pair_name, bam_file in alignment_files.items():
                ref_sample, query_sample = pair_name.split('_')
                f.write(f"  {ref_sample} vs {query_sample}: {os.path.basename(bam_file)}\n")
            f.write(f"总比对数 | Total alignments: {len(alignment_files)}\n\n")
            
            # SyRI分析结果 | SyRI analysis results
            f.write("SyRI结构变异分析结果 | SyRI Structural Variation Analysis Results:\n")
            f.write("-" * 40 + "\n")
            for pair_name, syri_file in syri_files.items():
                ref_sample, query_sample = pair_name.split('_')
                f.write(f"  {ref_sample} vs {query_sample}: {os.path.basename(syri_file)}\n")
            f.write(f"总SyRI分析数 | Total SyRI analyses: {len(syri_files)}\n\n")
            
            # 输出文件 | Output files
            f.write("主要输出文件 | Main Output Files:\n")
            f.write("-" * 40 + "\n")
            f.write(f"  - alignments/: BAM比对文件 | BAM alignment files\n")
            f.write(f"  - syri_results/: SyRI分析结果 | SyRI analysis results\n")
            f.write(f"  - plots/: 可视化图表 | Visualization plots\n")
            f.write(f"  - genomes.txt: 基因组配置文件 | Genome configuration file\n")
            
            if self.config.chromosome:
                f.write(f"  - collinearity_{self.config.chromosome}.{self.config.plotsr_format}: 染色体共线性图 | Chromosome collinearity plot\n")
            else:
                f.write(f"  - collinearity_all.{self.config.plotsr_format}: 全基因组共线性图 | Whole genome collinearity plot\n")
            
            f.write(f"\n分析日志 | Analysis log: collinearity_analysis.log\n")
        
        self.logger.info(f"📋 总结报告已生成 | Summary report generated: {report_file}")
