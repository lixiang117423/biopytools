# """
# 结果处理模块 | Results Processing Module
# """

# import os
# import pandas as pd
# from pathlib import Path
# from typing import List

# class ResultsProcessor:
#     """结果处理器 | Results Processor"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def process_diversity_results(self, diversity_files: List[str]):
#         """处理多样性分析结果 | Process diversity analysis results"""
#         self.logger.info("处理多样性分析结果 | Processing diversity analysis results")
        
#         for file_prefix in diversity_files:
#             # 处理π值结果
#             pi_file = f"{file_prefix}.windowed.pi"
#             if os.path.exists(pi_file):
#                 self._format_output_file(pi_file, "pi")
            
#             # 处理Tajima's D结果
#             tajima_file = f"{file_prefix}.Tajima.D"
#             if os.path.exists(tajima_file):
#                 self._format_output_file(tajima_file, "tajima_d")
    
#     def process_fst_results(self, fst_files: List[str]):
#         """处理Fst分析结果 | Process Fst analysis results"""
#         self.logger.info("处理Fst分析结果 | Processing Fst analysis results")
        
#         for file_prefix in fst_files:
#             fst_file = f"{file_prefix}.windowed.weir.fst"
#             if os.path.exists(fst_file):
#                 self._format_output_file(fst_file, "fst")
    
#     def _format_output_file(self, input_file: str, analysis_type: str):
#         """格式化输出文件 | Format output file"""
#         try:
#             df = pd.read_csv(input_file, sep='\t')
            
#             # 根据输出格式保存 | Save according to output format
#             base_name = Path(input_file).stem
#             output_file = self.config.output_path / f"{base_name}_{analysis_type}.{self.config.output_format}"
            
#             if self.config.output_format == 'csv':
#                 df.to_csv(output_file, index=False)
#             elif self.config.output_format == 'tsv':
#                 df.to_csv(output_file, sep='\t', index=False)
#             elif self.config.output_format == 'json':
#                 df.to_json(output_file, orient='records', indent=2)
#             else:  # txt
#                 df.to_csv(output_file, sep='\t', index=False)
            
#             self.logger.info(f"格式化输出文件 | Formatted output file: {output_file}")
            
#         except Exception as e:
#             self.logger.error(f"处理文件失败 | Failed to process file {input_file}: {e}")

# class SummaryGenerator:
#     """总结报告生成器 | Summary Report Generator"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_summary(self):
#         """生成分析总结报告 | Generate analysis summary report"""
#         self.logger.info("生成分析总结报告 | Generating analysis summary report")
        
#         summary_file = self.config.output_path / "analysis_summary.txt"
        
#         with open(summary_file, 'w', encoding='utf-8') as f:
#             f.write("群体遗传分析总结报告 | Population Genetics Analysis Summary Report\n")
#             f.write("="*80 + "\n")
#             f.write(f"分析时间 | Analysis time: {pd.Timestamp.now()}\n")
#             f.write(f"输入VCF文件 | Input VCF file: {self.config.vcf_file}\n")
#             f.write(f"输出目录 | Output directory: {self.config.output_dir}\n\n")
            
#             # 分析参数 | Analysis parameters
#             f.write("分析参数 | Analysis Parameters:\n")
#             f.write("-" * 40 + "\n")
#             f.write(f"  窗口大小 | Window sizes: {self.config.window_sizes} bp\n")
#             f.write(f"  窗口重叠率 | Window overlap: {self.config.window_overlap * 100}%\n")
#             f.write(f"  MAF阈值 | MAF threshold: {self.config.maf}\n")
#             f.write(f"  缺失率阈值 | Missing rate threshold: {self.config.missing_rate}\n")
#             f.write(f"  HWE p值阈值 | HWE p-value threshold: {self.config.hwe_pvalue}\n")
#             f.write(f"  输出格式 | Output format: {self.config.output_format}\n\n")
            
#             # 执行的分析 | Performed analyses
#             f.write("执行的分析 | Performed Analyses:\n")
#             f.write("-" * 40 + "\n")
#             if self.config.calculate_pi:
#                 f.write("  ✓ π (nucleotide diversity) - 核苷酸多样性\n")
#             if self.config.calculate_theta_w:
#                 f.write("  ✓ θw (Watterson's estimator) - Watterson估计量\n")
#             if self.config.calculate_tajima_d:
#                 f.write("  ✓ Tajima's D - Tajima's D中性检验\n")
#             if self.config.calculate_fst:
#                 f.write("  ✓ Fst (population differentiation) - 群体分化指数\n")
#             if self.config.calculate_ibd:
#                 f.write("  ✓ IBD (identity-by-descent) - 同源性分析\n")
#             if self.config.calculate_ld:
#                 f.write("  ✓ LD (linkage disequilibrium) - 连锁不平衡\n")
#             if self.config.calculate_ne:
#                 f.write("  ✓ Ne (effective population size) - 有效群体大小\n")
        
#         self.logger.info(f"总结报告已生成 | Summary report generated: {summary_file}")

"""
结果处理模块 | Results Processing Module
"""

import os
import pandas as pd
from pathlib import Path
from typing import List

class ResultsProcessor:
    """结果处理器 | Results Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def process_all_results(self):
        """处理所有分析结果 | Process all analysis results"""
        self.logger.info("处理所有分析结果 | Processing all analysis results")
        
        # 查找并处理所有结果文件 | Find and process all result files
        output_path = self.config.output_path
        
        # 处理多样性结果文件 | Process diversity result files
        diversity_files = list(output_path.glob("diversity_window_*"))
        if diversity_files:
            self.logger.info(f"找到 {len(diversity_files)} 个多样性结果文件 | Found {len(diversity_files)} diversity result files")
            diversity_prefixes = [str(f) for f in diversity_files]
            self.process_diversity_results(diversity_prefixes)
        
        # 处理Fst结果文件 | Process Fst result files
        fst_files = list(output_path.glob("fst_*"))
        if fst_files:
            self.logger.info(f"找到 {len(fst_files)} 个Fst结果文件 | Found {len(fst_files)} Fst result files")
            fst_prefixes = [str(f) for f in fst_files]
            self.process_fst_results(fst_prefixes)
        
        # 处理其他结果文件 | Process other result files
        self._process_other_results()
    
    def process_diversity_results(self, diversity_files: List[str]):
        """处理多样性分析结果 | Process diversity analysis results"""
        self.logger.info("处理多样性分析结果 | Processing diversity analysis results")
        
        for file_prefix in diversity_files:
            # 处理π值结果 | Process π value results
            pi_file = f"{file_prefix}.windowed.pi"
            if os.path.exists(pi_file):
                self._format_output_file(pi_file, "pi")
            
            # 处理Tajima's D结果 | Process Tajima's D results
            tajima_file = f"{file_prefix}.Tajima.D"
            if os.path.exists(tajima_file):
                self._format_output_file(tajima_file, "tajima_d")
    
    def process_fst_results(self, fst_files: List[str]):
        """处理Fst分析结果 | Process Fst analysis results"""
        self.logger.info("处理Fst分析结果 | Processing Fst analysis results")
        
        for file_prefix in fst_files:
            fst_file = f"{file_prefix}.windowed.weir.fst"
            if os.path.exists(fst_file):
                self._format_output_file(fst_file, "fst")
    
    def _process_other_results(self):
        """处理其他分析结果 | Process other analysis results"""
        self.logger.info("检查其他分析结果文件 | Checking other analysis result files")
        
        output_path = self.config.output_path
        
        # 检查IBD结果 | Check IBD results
        ibd_files = list(output_path.glob("*ibd*"))
        if ibd_files:
            self.logger.info(f"找到IBD结果文件: {[f.name for f in ibd_files]}")
        
        # 检查LD结果 | Check LD results
        ld_files = list(output_path.glob("*ld*"))
        if ld_files:
            self.logger.info(f"找到LD结果文件: {[f.name for f in ld_files]}")
        
        # 检查Ne结果 | Check Ne results
        ne_files = list(output_path.glob("*ne*"))
        if ne_files:
            self.logger.info(f"找到Ne结果文件: {[f.name for f in ne_files]}")
    
    def _format_output_file(self, input_file: str, analysis_type: str):
        """格式化输出文件 | Format output file"""
        try:
            if not os.path.exists(input_file):
                self.logger.warning(f"文件不存在: {input_file} | File does not exist: {input_file}")
                return
            
            df = pd.read_csv(input_file, sep='\t')
            
            # 根据输出格式保存 | Save according to output format
            base_name = Path(input_file).stem
            output_file = self.config.output_path / f"{base_name}_{analysis_type}.{self.config.output_format}"
            
            if self.config.output_format == 'csv':
                df.to_csv(output_file, index=False)
            elif self.config.output_format == 'tsv':
                df.to_csv(output_file, sep='\t', index=False)
            elif self.config.output_format == 'json':
                df.to_json(output_file, orient='records', indent=2)
            else:  # txt
                df.to_csv(output_file, sep='\t', index=False)
            
            self.logger.info(f"格式化输出文件 | Formatted output file: {output_file}")
            
        except Exception as e:
            self.logger.error(f"处理文件失败 | Failed to process file {input_file}: {e}")

class SummaryGenerator:
    """总结报告生成器 | Summary Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self):
        """生成分析总结报告 | Generate analysis summary report"""
        self.logger.info("生成分析总结报告 | Generating analysis summary report")
        
        summary_file = self.config.output_path / "analysis_summary.txt"
        
        try:
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write("群体遗传分析总结报告 | Population Genetics Analysis Summary Report\n")
                f.write("="*80 + "\n")
                f.write(f"分析时间 | Analysis time: {pd.Timestamp.now()}\n")
                f.write(f"输入VCF文件 | Input VCF file: {self.config.vcf_file}\n")
                f.write(f"输出目录 | Output directory: {self.config.output_dir}\n\n")
                
                # 分析参数 | Analysis parameters
                f.write("分析参数 | Analysis Parameters:\n")
                f.write("-" * 40 + "\n")
                f.write(f"  窗口大小 | Window sizes: {self.config.window_sizes} bp\n")
                f.write(f"  窗口重叠率 | Window overlap: {self.config.window_overlap * 100}%\n")
                f.write(f"  MAF阈值 | MAF threshold: {self.config.maf}\n")
                f.write(f"  缺失率阈值 | Missing rate threshold: {self.config.missing_rate}\n")
                f.write(f"  HWE p值阈值 | HWE p-value threshold: {self.config.hwe_pvalue}\n")
                f.write(f"  输出格式 | Output format: {self.config.output_format}\n")
                f.write(f"  跳过质控 | Skip QC: {getattr(self.config, 'skip_qc', False)}\n\n")
                
                # 执行的分析 | Performed analyses
                f.write("执行的分析 | Performed Analyses:\n")
                f.write("-" * 40 + "\n")
                if self.config.calculate_pi:
                    f.write("  ✓ π (nucleotide diversity) - 核苷酸多样性\n")
                if self.config.calculate_theta_w:
                    f.write("  ✓ θw (Watterson's estimator) - Watterson估计量\n")
                if self.config.calculate_tajima_d:
                    f.write("  ✓ Tajima's D - Tajima's D中性检验\n")
                if self.config.calculate_fst:
                    f.write("  ✓ Fst (population differentiation) - 群体分化指数\n")
                if self.config.calculate_ibd:
                    f.write("  ✓ IBD (identity-by-descent) - 同源性分析\n")
                if self.config.calculate_ld:
                    f.write("  ✓ LD (linkage disequilibrium) - 连锁不平衡\n")
                if self.config.calculate_ne:
                    f.write("  ✓ Ne (effective population size) - 有效群体大小\n")
                
                # 输出文件统计 | Output file statistics
                f.write("\n输出文件统计 | Output File Statistics:\n")
                f.write("-" * 40 + "\n")
                self._write_file_statistics(f)
            
            self.logger.info(f"总结报告已生成 | Summary report generated: {summary_file}")
            return str(summary_file)
            
        except Exception as e:
            self.logger.error(f"生成总结报告失败 | Failed to generate summary report: {e}")
            return None
    
    def _write_file_statistics(self, f):
        """写入文件统计信息 | Write file statistics"""
        try:
            output_path = self.config.output_path
            all_files = list(output_path.glob("*"))
            
            file_types = {
                'pi': len(list(output_path.glob("*.windowed.pi"))),
                'tajima_d': len(list(output_path.glob("*.Tajima.D"))),
                'fst': len(list(output_path.glob("*.weir.fst"))),
                'ld': len(list(output_path.glob("*ld*"))),
                'ibd': len(list(output_path.glob("*ibd*"))),
                'log': len(list(output_path.glob("*.log"))),
                'other': len([f for f in all_files if f.is_file()])
            }
            
            for file_type, count in file_types.items():
                if count > 0:
                    f.write(f"  {file_type}: {count} 个文件 | files\n")
            
            total_size = sum(f.stat().st_size for f in all_files if f.is_file())
            f.write(f"\n  总文件大小 | Total file size: {total_size / (1024*1024):.2f} MB\n")
            
        except Exception as e:
            f.write(f"  统计信息获取失败 | Failed to get statistics: {e}\n")