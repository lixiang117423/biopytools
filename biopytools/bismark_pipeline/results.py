"""
Bismark流程结果处理模块 | Bismark Pipeline Results Module
"""
from datetime import datetime
from pathlib import Path
from . import __version__

class ResultsManager:
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    # def split_cx_report(self, cx_report_path: Path):
    #     """
    #     将CX_report.txt文件按上下文 (CG, CHG, CHH) 拆分，并使用友好命名。
    #     Splits the CX_report.txt file by context (CG, CHG, CHH) with user-friendly names.
    #     """
    #     # *** FIX START: Correctly handle 'CG' context from Bismark output ***
    #     base_name = cx_report_path.name.rsplit('.', 2)[0] if '.CX_report.' in cx_report_path.name else cx_report_path.stem
        
    #     # Map file context codes to user-friendly names for output files
    #     context_map = {
    #         'CG': 'CpG',
    #         'CHG': 'CHG',
    #         'CHH': 'CHH'
    #     }
        
    #     output_files = {
    #         file_context: open(cx_report_path.with_name(f"{base_name}.{friendly_name}.txt"), 'w')
    #         for file_context, friendly_name in context_map.items()
    #     }
        
    #     counts = {context: 0 for context in context_map.keys()}
    #     counts['other'] = 0

    #     try:
    #         with open(cx_report_path, 'r') as infile:
    #             for line in infile:
    #                 parts = line.strip().split('\t')
    #                 if len(parts) >= 6:
    #                     context_from_file = parts[5]  # This will be 'CG', 'CHG', or 'CHH'
    #                     if context_from_file in output_files:
    #                         output_files[context_from_file].write(line)
    #                         counts[context_from_file] += 1
    #                     else:
    #                         counts['other'] += 1
            
    #         self.logger.info("    拆分统计 | Split statistics:")
    #         for file_context, friendly_name in context_map.items():
    #             count = counts[file_context]
    #             if count > 0:
    #                 self.logger.info(f"      - {friendly_name} (context {file_context}): {count} 行 | lines")
    #         if counts['other'] > 0:
    #             self.logger.warning(f"      - 其他 | Other contexts: {counts['other']} 行 | lines")
        
    #     except Exception as e:
    #         self.logger.error(f"  ❌ 拆分CX报告时出错 | Error splitting CX report: {e}")
        
    #     finally:
    #         for f in output_files.values():
    #             f.close()
        # *** FIX END ***
    
    # 这是新的、经过严格调试的正确代码
    def split_cx_report(self, cx_report_path: Path):
        """
        将CX_report.txt文件按上下文 (CG, CHG, CHH) 拆分，并使用友好命名。
        Splits the CX_report.txt file by context (CG, CHG, CHH) with user-friendly names.
        """
        # *** FIX START: Robustly extract the base name ***
        # Example filename: FZY4205_1_clean_bismark_bt2_pe.CX_report.txt
        # We want to get "FZY4205_1_clean_bismark_bt2_pe"
        original_filename = cx_report_path.name
        if '.CX_report.txt' in original_filename:
            base_name = original_filename.replace('.CX_report.txt', '')
        else:
            # Fallback for unexpected naming, less robust
            base_name = cx_report_path.stem 
        # Now base_name is correctly "FZY4205_1_clean_bismark_bt2_pe"
        # *** FIX END ***

        # 定义 "报告中的上下文" -> "输出文件名中的标签" 的映射
        # Map from context in report -> label in output filename
        context_map = {
            'CG': 'CpG',
            'CHG': 'CHG',
            'CHH': 'CHH'
        }
        
        # 根据报告中的上下文 ('CG', 'CHG', 'CHH') 创建文件句柄
        # The keys of this dictionary MUST match the context strings from the report file
        output_files = {
            report_context: open(cx_report_path.with_name(f"{base_name}.{friendly_name}.txt"), 'w')
            for report_context, friendly_name in context_map.items()
        }
        
        # 初始化计数器
        counts = {context: 0 for context in context_map.keys()}
        counts['other'] = 0

        try:
            with open(cx_report_path, 'r') as infile:
                for line in infile:
                    parts = line.strip().split('\t')
                    if len(parts) >= 6:
                        context_from_file = parts[5]  # 这将是 'CG', 'CHG', 或 'CHH'
                        
                        # 直接使用文件中的上下文作为key来查找对应的文件句柄
                        if context_from_file in output_files:
                            output_files[context_from_file].write(line)
                            counts[context_from_file] += 1
                        else:
                            counts['other'] += 1
            
            self.logger.info("    拆分统计 | Split statistics:")
            for report_context, friendly_name in context_map.items():
                count = counts[report_context]
                # 即使数量为0也打印，方便调试
                self.logger.info(f"      - {friendly_name} (context {report_context}): {count} 行 | lines")
            if counts['other'] > 0:
                self.logger.warning(f"      - 其他 | Other contexts: {counts['other']} 行 | lines")
        
        except Exception as e:
            self.logger.error(f"  ❌ 拆分CX报告时出错 | Error splitting CX report: {e}")
        
        finally:
            # 确保所有文件都被关闭
            for f in output_files.values():
                f.close()

    def generate_summary_report(self):
        report_file = Path(self.config.output_dir) / "bismark_pipeline_summary.txt"
        self.logger.info(f"生成分析总结报告 | Generating analysis summary report: {report_file}")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*60 + "\n")
            f.write("Bismark 甲基化分析总结报告 | Bismark Methylation Analysis Summary Report\n")
            f.write("="*60 + "\n\n")
            f.write(f"分析时间 | Analysis Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"分析版本 | Pipeline Version: {__version__}\n\n")
            f.write("--- 输入与输出 | Inputs & Outputs ---\n")
            f.write(f"基因组文件 | Genome File: {self.config.genome_fa}\n")
            f.write(f"原始数据目录 | Raw Data Directory: {self.config.raw_dir}\n")
            f.write(f"主输出目录 | Main Output Directory: {self.config.output_dir}\n\n")
            f.write("--- 关键参数 | Key Parameters ---\n")
            f.write(f"文件匹配模式 (R1) | File Pattern (R1): *{self.config.pattern}\n")
            f.write(f"线程数 | Threads: {self.config.threads}\n")
            f.write(f"排序缓存 | Sort Buffer: {self.config.sort_buffer}\n")
            f.write(f"忽略重叠Reads | Ignore Overlapping Reads: {'是 | Yes' if self.config.no_overlap else '否 | No'}\n\n")
            f.write("--- 结果文件位置 | Result File Locations ---\n")
            f.write(f"基因组索引目录 | Genome Index Directory: {self.config.genome_dir}\n")
            f.write(f"比对BAM文件 | Alignment BAM Files: {self.config.mapping_dir}\n")
            f.write(f"甲基化结果 | Methylation Results: {self.config.result_dir}\n")
            f.write("  (包括原始CX报告和按CpG/CHG/CHH拆分的文件 | Includes original CX reports and files split by CpG/CHG/CHH)\n")
            f.write(f"中间/报告文件 | Intermediate/Report Files: {self.config.tmp_dir}\n")
            f.write(f"运行日志 | Run Logs: {Path(self.config.output_dir) / 'logs'}\n")
        
        self.logger.info("✅ 总结报告已生成 | Summary report has been generated.")
