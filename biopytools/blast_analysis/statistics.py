# """
# 📈 BLAST统计分析模块
# """

# import pandas as pd
# from pathlib import Path
# from datetime import datetime

# class StatisticsGenerator:
#     """📈 BLAST统计生成器"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def generate_statistics_report(self, summary_file: str) -> str:
#         """生成统计报告"""
#         self.logger.info("生成BLAST统计报告")
        
#         stats_file = self.config.output_path / "blast_statistics.txt"
        
#         try:
#             # 读取汇总数据
#             if Path(summary_file).exists() and Path(summary_file).stat().st_size > 0:
#                 df = pd.read_csv(summary_file, sep='\t', encoding='utf-8')
#             else:
#                 df = pd.DataFrame()
            
#             with open(stats_file, 'w', encoding='utf-8') as f:
#                 self._write_basic_info(f)
#                 self._write_summary_stats(f, df)
#                 self._write_file_stats(f, df)
#                 self._write_evalue_distribution(f, df)
#                 self._write_identity_distribution(f, df)
#                 self._write_coverage_distribution(f, df)
            
#             self.logger.info(f"统计报告已生成: {stats_file}")
#             return str(stats_file)
            
#         except Exception as e:
#             self.logger.error(f"生成统计报告时出错: {e}")
#             return ""
    
#     def _write_basic_info(self, f):
#         """写入基本信息"""
#         f.write("BLAST比对统计报告\n")
#         f.write("=" * 80 + "\n")
#         f.write(f"分析日期: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
#         f.write(f"BLAST程序: {self.config.blast_type.upper()}\n")
        
#         if self.config.input_path:
#             f.write(f"输入路径: {self.config.input_path}\n")
#         elif self.config.sample_map_file:
#             f.write(f"样品映射文件: {self.config.sample_map_file}\n")
        
#         f.write(f"目标基因文件: {self.config.target_file}\n")
#         f.write(f"使用线程数: {self.config.threads}\n")
#         f.write(f"E-value阈值: {self.config.evalue}\n")
#         f.write(f"最小相似度: {self.config.min_identity}%%\n")
#         f.write(f"最小覆盖度: {self.config.min_coverage}%%\n")
#         f.write("\n")
    
#     def _write_summary_stats(self, f, df: pd.DataFrame):
#         """写入汇总统计"""
#         f.write("汇总统计\n")
#         f.write("-" * 40 + "\n")
        
#         if df.empty:
#             f.write("无比对结果\n\n")
#             return
        
#         total_hits = len(df)
#         unique_files = df['输入文件'].nunique() if '输入文件' in df.columns else 0
#         unique_samples = df['样品名称'].nunique() if '样品名称' in df.columns else 0
        
#         f.write(f"总比对命中数: {total_hits}\n")
#         f.write(f"涉及输入文件数: {unique_files}\n")
#         f.write(f"涉及样品数: {unique_samples}\n")
        
#         if total_hits > 0:
#             avg_identity = df['序列相似度(%)'].mean() if '序列相似度(%)' in df.columns else 0
#             avg_coverage = df['目标序列覆盖度(%)'].mean() if '目标序列覆盖度(%)' in df.columns else 0
#             f.write(f"平均序列相似度: {avg_identity:.2f}%%\n")
#             f.write(f"平均目标覆盖度: {avg_coverage:.2f}%%\n")
        
#         f.write("\n")
    
#     def _write_file_stats(self, f, df: pd.DataFrame):
#         """写入文件统计"""
#         f.write("各输入文件命中统计\n")
#         f.write("-" * 40 + "\n")
        
#         if df.empty or '输入文件' not in df.columns:
#             f.write("无数据可统计\n\n")
#             return
        
#         file_stats = df['输入文件'].value_counts().head(20)
        
#         f.write("输入文件\t命中数\n")
#         for file_name, count in file_stats.items():
#             f.write(f"{file_name}\t{count}\n")
        
#         f.write("\n")
    
#     def _write_evalue_distribution(self, f, df: pd.DataFrame):
#         """写入E-value分布"""
#         f.write("E-value分布统计\n")
#         f.write("-" * 40 + "\n")
        
#         if df.empty or 'E-value' not in df.columns:
#             f.write("无E-value数据\n\n")
#             return
        
#         f.write("E-value范围\t命中数\n")
        
#         evalue_ranges = [
#             ("≤1e-50", lambda x: x <= 1e-50),
#             ("1e-50 ~ 1e-30", lambda x: 1e-50 < x <= 1e-30),
#             ("1e-30 ~ 1e-20", lambda x: 1e-30 < x <= 1e-20),
#             ("1e-20 ~ 1e-10", lambda x: 1e-20 < x <= 1e-10),
#             ("1e-10 ~ 1e-5", lambda x: 1e-10 < x <= 1e-5),
#             (">1e-5", lambda x: x > 1e-5)
#         ]
        
#         for range_name, condition in evalue_ranges:
#             count = df[df['E-value'].apply(condition)].shape[0]
#             f.write(f"{range_name}\t{count}\n")
        
#         f.write("\n")
    
#     def _write_identity_distribution(self, f, df: pd.DataFrame):
#         """写入相似度分布"""
#         f.write("序列相似度分布\n")
#         f.write("-" * 40 + "\n")
        
#         if df.empty or '序列相似度(%)' not in df.columns:
#             f.write("无相似度数据\n\n")
#             return
        
#         f.write("相似度范围\t命中数\n")
        
#         identity_ranges = [
#             ("≥95%%", lambda x: x >= 95),
#             ("90-95%%", lambda x: 90 <= x < 95),
#             ("80-90%%", lambda x: 80 <= x < 90),
#             ("70-80%%", lambda x: 70 <= x < 80),
#             ("<70%%", lambda x: x < 70)
#         ]
        
#         for range_name, condition in identity_ranges:
#             count = df[df['序列相似度(%)'].apply(condition)].shape[0]
#             f.write(f"{range_name}\t{count}\n")
        
#         f.write("\n")
    
#     def _write_coverage_distribution(self, f, df: pd.DataFrame):
#         """写入覆盖度分布"""
#         f.write("目标序列覆盖度分布\n")
#         f.write("-" * 40 + "\n")
        
#         if df.empty or '目标序列覆盖度(%)' not in df.columns:
#             f.write("无覆盖度数据\n\n")
#             return
        
#         f.write("覆盖度范围\t命中数\n")
        
#         coverage_ranges = [
#             ("≥90%%", lambda x: x >= 90),
#             ("80-90%%", lambda x: 80 <= x < 90),
#             ("70-80%%", lambda x: 70 <= x < 80),
#             ("50-70%%", lambda x: 50 <= x < 70),
#             ("<50%%", lambda x: x < 50)
#         ]
        
#         for range_name, condition in coverage_ranges:
#             count = df[df['目标序列覆盖度(%)'].apply(condition)].shape[0]
#             f.write(f"{range_name}\t{count}\n")
        
#         f.write("\n")

"""
📈 BLAST统计分析模块 | BLAST Statistics Analysis Module
"""

import pandas as pd
from pathlib import Path
from datetime import datetime

class StatisticsGenerator:
    """📈 BLAST统计生成器 | BLAST Statistics Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_statistics_report(self, summary_file: str) -> str:
        """生成统计报告 | Generate Statistics Report"""
        self.logger.info("生成BLAST统计报告 | Generating BLAST statistics report")
        
        stats_file = self.config.output_path / "blast_statistics.txt"
        
        try:
            # 读取汇总数据 | Read summary data
            if Path(summary_file).exists() and Path(summary_file).stat().st_size > 0:
                # 使用双语列名读取 | Read with bilingual column names
                df = pd.read_csv(summary_file, sep='\t', encoding='utf-8')
            else:
                df = pd.DataFrame()
            
            with open(stats_file, 'w', encoding='utf-8') as f:
                self._write_basic_info(f)
                self._write_summary_stats(f, df)
                self._write_file_stats(f, df)
                self._write_evalue_distribution(f, df)
                self._write_identity_distribution(f, df)
                self._write_coverage_distribution(f, df)
            
            self.logger.info(f"统计报告已生成 | Statistics report generated: {stats_file}")
            return str(stats_file)
            
        except Exception as e:
            self.logger.error(f"生成统计报告时出错 | Error generating statistics report: {e}")
            return ""
    
    def _write_basic_info(self, f):
        """写入基本信息 | Write Basic Information"""
        f.write("BLAST比对统计报告 | BLAST Alignment Statistics Report\n")
        f.write("=" * 80 + "\n")
        f.write(f"分析日期 | Analysis Date: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
        f.write(f"BLAST程序 | BLAST Program: {self.config.blast_type.upper()}\n")
        
        if self.config.input_path:
            f.write(f"输入路径 | Input Path: {self.config.input_path}\n")
        elif self.config.sample_map_file:
            f.write(f"样品映射文件 | Sample Map File: {self.config.sample_map_file}\n")
        
        f.write(f"目标基因文件 | Target Gene File: {self.config.target_file}\n")
        f.write(f"使用线程数 | Threads Used: {self.config.threads}\n")
        f.write(f"E-value阈值 | E-value Threshold: {self.config.evalue}\n")
        f.write(f"最小相似度 | Minimum Identity: {self.config.min_identity}%%\n")
        f.write(f"最小覆盖度 | Minimum Coverage: {self.config.min_coverage}%%\n")
        f.write("\n")
    
    def _write_summary_stats(self, f, df: pd.DataFrame):
        """写入汇总统计 | Write Summary Statistics"""
        f.write("汇总统计 | Summary Statistics\n")
        f.write("-" * 40 + "\n")
        
        if df.empty:
            f.write("无比对结果 | No alignment results\n\n")
            return
        
        # 使用双语列名 | Use bilingual column names
        input_file_col = '输入文件 | Input File'
        sample_name_col = '样品名称 | Sample Name'
        identity_col = '序列相似度(%) | Sequence Identity (%)'
        coverage_col = '目标序列覆盖度(%) | Target Sequence Coverage (%)'

        total_hits = len(df)
        unique_files = df[input_file_col].nunique() if input_file_col in df.columns else 0
        unique_samples = df[sample_name_col].nunique() if sample_name_col in df.columns else 0
        
        f.write(f"总比对命中数 | Total Alignment Hits: {total_hits}\n")
        f.write(f"涉及输入文件数 | Number of Input Files Involved: {unique_files}\n")
        f.write(f"涉及样品数 | Number of Samples Involved: {unique_samples}\n")
        
        if total_hits > 0:
            avg_identity = df[identity_col].mean() if identity_col in df.columns else 0
            avg_coverage = df[coverage_col].mean() if coverage_col in df.columns else 0
            f.write(f"平均序列相似度 | Average Sequence Identity: {avg_identity:.2f}%%\n")
            f.write(f"平均目标覆盖度 | Average Target Coverage: {avg_coverage:.2f}%%\n")
        
        f.write("\n")
    
    def _write_file_stats(self, f, df: pd.DataFrame):
        """写入文件统计 | Write File Statistics"""
        f.write("各输入文件命中统计 | Hit Statistics per Input File\n")
        f.write("-" * 40 + "\n")
        
        input_file_col = '输入文件 | Input File'
        if df.empty or input_file_col not in df.columns:
            f.write("无数据可统计 | No data to count\n\n")
            return
        
        file_stats = df[input_file_col].value_counts().head(20)
        
        f.write("输入文件 | Input File\t命中数 | Hit Count\n")
        for file_name, count in file_stats.items():
            f.write(f"{file_name}\t{count}\n")
        
        f.write("\n")
    
    def _write_evalue_distribution(self, f, df: pd.DataFrame):
        """写入E-value分布 | Write E-value Distribution"""
        f.write("E-value分布统计 | E-value Distribution Statistics\n")
        f.write("-" * 40 + "\n")
        
        if df.empty or 'E-value' not in df.columns:
            f.write("无E-value数据 | No E-value data\n\n")
            return
        
        f.write("E-value范围 | E-value Range\t命中数 | Hit Count\n")
        
        evalue_ranges = [
            ("≤1e-50", lambda x: x <= 1e-50),
            ("1e-50 ~ 1e-30", lambda x: 1e-50 < x <= 1e-30),
            ("1e-30 ~ 1e-20", lambda x: 1e-30 < x <= 1e-20),
            ("1e-20 ~ 1e-10", lambda x: 1e-20 < x <= 1e-10),
            ("1e-10 ~ 1e-5", lambda x: 1e-10 < x <= 1e-5),
            (">1e-5", lambda x: x > 1e-5)
        ]
        
        for range_name, condition in evalue_ranges:
            count = df[df['E-value'].apply(condition)].shape[0]
            f.write(f"{range_name}\t{count}\n")
        
        f.write("\n")
    
    def _write_identity_distribution(self, f, df: pd.DataFrame):
        """写入相似度分布 | Write Identity Distribution"""
        f.write("序列相似度分布 | Sequence Identity Distribution\n")
        f.write("-" * 40 + "\n")
        
        identity_col = '序列相似度(%) | Sequence Identity (%)'
        if df.empty or identity_col not in df.columns:
            f.write("无相似度数据 | No identity data\n\n")
            return
        
        f.write("相似度范围 | Identity Range\t命中数 | Hit Count\n")
        
        identity_ranges = [
            ("≥95%%", lambda x: x >= 95),
            ("90-95%%", lambda x: 90 <= x < 95),
            ("80-90%%", lambda x: 80 <= x < 90),
            ("70-80%%", lambda x: 70 <= x < 80),
            ("<70%%", lambda x: x < 70)
        ]
        
        for range_name, condition in identity_ranges:
            count = df[df[identity_col].apply(condition)].shape[0]
            f.write(f"{range_name}\t{count}\n")
        
        f.write("\n")
    
    def _write_coverage_distribution(self, f, df: pd.DataFrame):
        """写入覆盖度分布 | Write Coverage Distribution"""
        f.write("目标序列覆盖度分布 | Target Sequence Coverage Distribution\n")
        f.write("-" * 40 + "\n")
        
        coverage_col = '目标序列覆盖度(%) | Target Sequence Coverage (%)'
        if df.empty or coverage_col not in df.columns:
            f.write("无覆盖度数据 | No coverage data\n\n")
            return
        
        f.write("覆盖度范围 | Coverage Range\t命中数 | Hit Count\n")
        
        coverage_ranges = [
            ("≥90%%", lambda x: x >= 90),
            ("80-90%%", lambda x: 80 <= x < 90),
            ("70-80%%", lambda x: 70 <= x < 80),
            ("50-70%%", lambda x: 50 <= x < 70),
            ("<50%%", lambda x: x < 50)
        ]
        
        for range_name, condition in coverage_ranges:
            count = df[df[coverage_col].apply(condition)].shape[0]
            f.write(f"{range_name}\t{count}\n")
        
        f.write("\n")