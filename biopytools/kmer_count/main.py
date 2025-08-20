"""
🚀 主程序模块 | Main Program Module
"""

import os
import sys
import argparse
import shutil
import pandas as pd
from typing import List
from pathlib import Path
from .utils import setup_logger, check_dependencies, CommandRunner
from .file_processor import FileProcessor
from .jellyfish_processor import JellyfishProcessor
from .data_processor import DataProcessor
from .window_analyzer import WindowAnalyzer
from .config import KmerCountConfig

class KmerCountAnalyzer:
    """🧬 K-mer计数分析器 | K-mer Count Analyzer"""
    
    def __init__(self, config):
        self.config = config
        self.logger = setup_logger(config.output_dir, config.verbose)
        self.cmd_runner = CommandRunner(self.logger)
        
        # 初始化处理器 | Initialize processors
        self.file_processor = FileProcessor(config, self.logger, self.cmd_runner)
        self.jellyfish_processor = JellyfishProcessor(config, self.logger, self.cmd_runner)
        self.data_processor = DataProcessor(config, self.logger)
        self.window_analyzer = WindowAnalyzer(config, self.logger)

    # def _merge_sample_files(self, sample_files: List[Path]) -> pd.DataFrame:
    #     """🔗 合并样品文件 | Merge sample files"""
    #     self.logger.info("🔗 读取并合并样品文件 | Reading and merging sample files")
        
    #     # 读取第一个文件作为基础
    #     base_df = pd.read_csv(sample_files[0], sep='\t')
    #     self.logger.info(f"📊 基础文件: {len(base_df)} 行")
        
    #     # 逐个合并其他文件
    #     for sample_file in sample_files[1:]:
    #         sample_df = pd.read_csv(sample_file, sep='\t')
    #         sample_name = sample_file.stem.replace('_kmer_abundance', '')
            
    #         # 只保留样品丰度列
    #         sample_cols = ['kmer'] + [col for col in sample_df.columns if col not in ['chr', 'start', 'end', 'kmer', 'kmer_id', 'score', 'strand']]
    #         merge_df = sample_df[sample_cols]
            
    #         base_df = pd.merge(base_df, merge_df, on='kmer', how='left')
    #         base_df = base_df.fillna(0)
    #         self.logger.info(f"📊 合并后: {len(base_df)} 行")
        
    #     return base_df

    def _merge_sample_files(self, sample_files: List[Path]) -> pd.DataFrame:
        """🔗 合并样品文件 | Merge sample files"""
        self.logger.info("🔗 基于BED文件重新构建最终矩阵 | Rebuilding final matrix based on BED file")
        
        # 以BED文件为基础构建最终矩阵
        if self.config.bed_file:
            final_df = self.data_processor.bed_info.copy()
            # 创建unique key：chr+start+end+kmer
            final_df['unique_key'] = final_df['chr'] + '_' + final_df['start'].astype(str) + '_' + final_df['end'].astype(str) + '_' + final_df['kmer']
        else:
            # 🔥 修改：如果没有BED文件，使用完整的k-mer ID-序列对
            final_df = pd.DataFrame(self.data_processor.kmer_pairs)
            self.logger.info(f"📊 基础矩阵(基于k-mer对): {len(final_df)} 行")
        
        self.logger.info(f"📊 基础矩阵: {len(final_df)} 行")
        
        # 为每个样本添加丰度列
        for sample_file in sample_files:
            sample_df = pd.read_csv(sample_file, sep='\t')
            sample_name = sample_file.stem.replace('_kmer_abundance', '')
            
            # 🔥 统一使用unique_key进行映射
            unique_abundance = dict(zip(sample_df['unique_key'], sample_df[sample_name]))
            final_df[sample_name] = final_df['unique_key'].map(unique_abundance).fillna(0).astype(int)
            
            self.logger.info(f"📊 添加样本{sample_name}: {len(final_df)} 行")
        
        # 🔥 关键修复：使用unique_key去除重复行
        original_len = len(final_df)
        final_df = final_df.drop_duplicates(subset=['unique_key'], keep='first')
        deduplicated_len = len(final_df)
        
        if original_len != deduplicated_len:
            self.logger.warning(f"🗑️ 检测到重复行，已去重: {original_len} → {deduplicated_len}")
        
        # 删除辅助列unique_key
        final_df = final_df.drop('unique_key', axis=1)
        
        self.logger.info(f"✅ 最终矩阵: {len(final_df)} 行")
        return final_df
    
    # def run_analysis(self):
    #     """🚀 运行完整分析 | Run complete analysis"""
    #     try:
    #         self.logger.info("="*60)
    #         self.logger.info("🚀 开始K-mer丰度分析 | Starting K-mer abundance analysis")
    #         self.logger.info("="*60)
            
    #         # 1. 检查依赖 | Check dependencies
    #         self.logger.info("🔍 步骤1: 检查依赖软件 | Step 1: Checking dependencies")
    #         if not check_dependencies(self.config.jellyfish_path, self.logger):
    #             raise RuntimeError("❌ 依赖检查失败 | Dependency check failed")
            
    #         # 2. 设置临时目录 | Setup temporary directory
    #         self.logger.info("📁 步骤2: 设置临时目录 | Step 2: Setting up temporary directory")
    #         self.config.setup_temp_dir()
    #         self.logger.info(f"📂 临时目录: {self.config.temp_dir}")
            
    #         # 3. 查找FASTQ文件 | Find FASTQ files
    #         self.logger.info("🔍 步骤3: 查找FASTQ文件 | Step 3: Finding FASTQ files")
    #         samples = self.file_processor.find_fastq_files()
            
    #         # 4. 准备k-mer库 | Prepare k-mer library
    #         self.logger.info("🧬 步骤4: 准备k-mer库 | Step 4: Preparing k-mer library")
    #         kmer_lib = self.file_processor.prepare_kmer_library()
            
    #         # 5. 解析注释信息 | Parse annotation information
    #         self.logger.info("📋 步骤5: 解析注释信息 | Step 5: Parsing annotation information")
    #         self.data_processor.parse_kmer_library()
    #         self.data_processor.parse_bed_file()
            
    #         # # 6. 处理每个样本 | Process each sample
    #         # self.logger.info("🔄 步骤6: 处理样本 | Step 6: Processing samples")
    #         # sample_results = []
            
    #         # for sample_name, r1_file, r2_file in samples:
    #         #     self.logger.info(f"👤 处理样本 | Processing sample: {sample_name}")
                
    #         #     # 6.1 解压缩文件 | Decompress files
    #         #     fastq_files = self.file_processor.decompress_files([r1_file, r2_file])
                
    #         #     # 6.2 k-mer计数 | K-mer counting
    #         #     jf_file = self.jellyfish_processor.count_kmers(sample_name, fastq_files)
                
    #         #     # 6.3 查询k-mer丰度 | Query k-mer abundance
    #         #     count_file = self.jellyfish_processor.query_kmers(sample_name, jf_file, kmer_lib)
                
    #         #     # 6.4 解析结果 | Parse results
    #         #     sample_df = self.data_processor.parse_jellyfish_output(count_file, sample_name)
    #         #     sample_results.append(sample_df)
                
    #         #     # 清理样本临时文件 | Clean sample temporary files
    #         #     if not self.config.keep_temp:
    #         #         for f in fastq_files:
    #         #             if f and os.path.exists(f):
    #         #                 os.remove(f)
    #         #         if not self.config.keep_binary and jf_file.exists():
    #         #             os.remove(jf_file)
            
    #         # # 7. 合并样本结果 | Merge sample results
    #         # self.logger.info("🔗 步骤7: 合并样本结果 | Step 7: Merging sample results")
    #         # merged_df = self.data_processor.merge_sample_results(sample_results)

    #         # 6. 处理每个样本 | Process each sample
    #         self.logger.info("🔄 步骤6: 处理样本 | Step 6: Processing samples")
    #         sample_files = []  # 存储每个样品的结果文件路径

    #         # for sample_name, r1_file, r2_file in samples:
    #         #     self.logger.info(f"👤 处理样本 | Processing sample: {sample_name}")
                
    #         #     # 6.1 解压缩文件 | Decompress files
    #         #     fastq_files = self.file_processor.decompress_files([r1_file, r2_file])
                
    #         #     # 6.2 k-mer计数 | K-mer counting
    #         #     jf_file = self.jellyfish_processor.count_kmers(sample_name, fastq_files)
    #         for sample_info in samples:
    #             if len(sample_info) == 4:  # 新格式：(sample_name, file1, file2_or_none, file_type)
    #                 sample_name, file1, file2, file_type = sample_info
    #                 input_files = [file1, file2] if file2 else [file1]
    #             else:  # 兼容旧格式
    #                 sample_name, file1, file2 = sample_info
    #                 input_files = [file1, file2] if file2 else [file1]
    #                 file_type = 'fastq'  # 默认为fastq
                
    #             self.logger.info(f"👤 处理样本 | Processing sample: {sample_name} ({file_type})")
                
    #             # 6.1 准备文件
    #             prepared_files = self.file_processor.decompress_files(input_files)
                
    #             # 6.2 k-mer计数
    #             jf_file = self.jellyfish_processor.count_kmers(sample_name, prepared_files, file_type)
                
    #             # 6.3 查询k-mer丰度 | Query k-mer abundance
    #             count_file = self.jellyfish_processor.query_kmers(sample_name, jf_file, kmer_lib)
                
    #             # 6.4 解析结果 | Parse results
    #             sample_df = self.data_processor.parse_jellyfish_output(count_file, sample_name)
                
    #             # 6.5 立即与BED合并并保存 | Immediately merge with BED and save
    #             # merged_sample_df = self.data_processor.merge_single_sample_with_bed(sample_df, sample_name)
    #             # sample_file = self.config.output_dir / f"{sample_name}_kmer_abundance.tsv"
    #             # 6.5 立即与BED合并并保存 | Immediately merge with BED and save
    #             merged_sample_df = self.data_processor.merge_single_sample_with_bed(sample_df, sample_name)

    #             # 创建each_sample子目录
    #             each_sample_dir = self.config.output_dir / "each_sample"
    #             each_sample_dir.mkdir(exist_ok=True)

    #             sample_file = each_sample_dir / f"{sample_name}_kmer_abundance.tsv"

    #             merged_sample_df.to_csv(sample_file, sep='\t', index=False)
    #             sample_files.append(sample_file)
    #             self.logger.info(f"✅ 样品{sample_name}结果已保存: {sample_file}")
                
    #             # 清理解压文件，但保留jf和count文件
    #             for f in fastq_files:
    #                 if f and os.path.exists(f):
    #                     os.remove(f)

    #         # 7. 合并所有样本结果 | Merge all sample results
    #         self.logger.info("🔗 步骤7: 合并所有样本结果 | Step 7: Merging all sample results")
    #         final_df = self._merge_sample_files(sample_files)
            
    #         # 8. 添加注释信息 | Add annotation information
    #         # self.logger.info("📝 步骤8: 添加注释信息 | Step 8: Adding annotation information")
    #         # final_df = self.data_processor.add_annotation_info(merged_df)
    #         # 8. 保存最终结果 | Save final results
    #         self.logger.info("💾 步骤8: 保存最终结果 | Step 8: Saving final results")
    #         output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
    #         final_df.to_csv(output_file, sep='\t', index=False)
    #         self.logger.info(f"✅ 最终结果已保存 | Final results saved: {output_file}")
                        
    #         # 9. 保存主要结果 | Save main results
    #         self.logger.info("💾 步骤9: 保存结果 | Step 9: Saving results")
    #         output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
    #         final_df.to_csv(output_file, sep='\t', index=False)
    #         self.logger.info(f"✅ 主要结果已保存 | Main results saved: {output_file}")
            
    #         # 10. 滑动窗口分析 | Sliding window analysis
    #         if self.config.bed_file and self.config.window_size:
    #             self.logger.info("🪟 步骤10: 滑动窗口分析 | Step 10: Sliding window analysis")
    #             window_df = self.window_analyzer.sliding_window_analysis(final_df)
    #             if window_df is not None:
    #                 window_file = self.config.output_dir / "sliding_window_analysis.tsv"
    #                 window_df.to_csv(window_file, sep='\t', index=False)
    #                 self.logger.info(f"✅ 滑动窗口结果已保存 | Sliding window results saved: {window_file}")
            
    #         # 11. 创建0/1矩阵 | Create 0/1 matrix
    #         if self.config.keep_binary:
    #             self.logger.info("🔢 步骤11: 创建0/1矩阵 | Step 11: Creating 0/1 matrix")
    #             binary_df = self.window_analyzer.create_binary_matrix(final_df)
    #             binary_file = self.config.output_dir / "kmer_binary_matrix.tsv"
    #             binary_df.to_csv(binary_file, sep='\t', index=False)
    #             self.logger.info(f"✅ 0/1矩阵已保存 | 0/1 matrix saved: {binary_file}")
            
    #         # 12. 生成统计报告 | Generate statistics report
    #         self._generate_summary_report(final_df)
            
    #         self.logger.info("="*60)
    #         self.logger.info("🎉 K-mer丰度分析完成! | K-mer abundance analysis completed!")
    #         self.logger.info("="*60)
            
    #     except Exception as e:
    #         self.logger.error(f"💥 分析过程中出现错误 | Error during analysis: {e}")
    #         raise
    #     finally:
    #         # 清理临时目录 | Clean up temporary directory
    #         if self.config.temp_dir and not self.config.keep_temp:
    #             try:
    #                 shutil.rmtree(self.config.temp_dir)
    #                 self.logger.info("🧹 临时目录已清理 | Temporary directory cleaned")
    #             except Exception as e:
    #                 self.logger.warning(f"⚠️ 清理临时目录失败 | Failed to clean temporary directory: {e}")

    def run_analysis(self):
        """🚀 运行完整分析 | Run complete analysis"""
        try:
            self.logger.info("="*60)
            self.logger.info("🚀 开始K-mer丰度分析 | Starting K-mer abundance analysis")
            self.logger.info("="*60)
            
            # 1. 检查依赖 | Check dependencies
            self.logger.info("🔍 步骤1: 检查依赖软件 | Step 1: Checking dependencies")
            if not check_dependencies(self.config.jellyfish_path, self.logger):
                raise RuntimeError("❌ 依赖检查失败 | Dependency check failed")
            
            # 2. 设置临时目录 | Setup temporary directory
            self.logger.info("📁 步骤2: 设置临时目录 | Step 2: Setting up temporary directory")
            self.config.setup_temp_dir()
            self.logger.info(f"📂 临时目录: {self.config.temp_dir}")
            
            # 3. 🔥 修改：查找输入文件（支持FASTQ和FASTA）| Find input files (support FASTQ and FASTA)
            self.logger.info("🔍 步骤3: 查找输入文件 | Step 3: Finding input files")
            samples = self.file_processor.find_input_files()  # 🔥 改名
            
            # 4. 准备k-mer库 | Prepare k-mer library
            self.logger.info("🧬 步骤4: 准备k-mer库 | Step 4: Preparing k-mer library")
            kmer_lib = self.file_processor.prepare_kmer_library()
            
            # 5. 解析注释信息 | Parse annotation information
            self.logger.info("📋 步骤5: 解析注释信息 | Step 5: Parsing annotation information")
            self.data_processor.parse_kmer_library()
            self.data_processor.parse_bed_file()
            
            # 6. 🔥 修改：处理每个样本（支持文件类型信息）| Process each sample (support file type info)
            self.logger.info("🔄 步骤6: 处理样本 | Step 6: Processing samples")
            sample_files = []
            
            for sample_info in samples:
                # 🔥 新格式：(sample_name, file1, file2_or_none, file_type)
                if len(sample_info) == 4:
                    sample_name, file1, file2, file_type = sample_info
                    input_files = [file1] if file2 is None else [file1, file2]
                else:
                    # 兼容旧格式：(sample_name, file1, file2_or_none)
                    sample_name, file1, file2 = sample_info
                    input_files = [file1] if file2 is None else [file1, file2]
                    file_type = 'fastq'  # 默认为FASTQ
                
                self.logger.info(f"👤 处理样本 | Processing sample: {sample_name} ({file_type})")
                
                # 6.1 准备文件
                prepared_files = self.file_processor.decompress_files(input_files)
                
                # 6.2 🔥 修改：k-mer计数（传递文件类型）| K-mer counting (pass file type)
                jf_file = self.jellyfish_processor.count_kmers(sample_name, prepared_files, file_type)
                
                # 6.3 查询k-mer丰度 | Query k-mer abundance
                count_file = self.jellyfish_processor.query_kmers(sample_name, jf_file, kmer_lib)
                
                # 6.4 解析结果 | Parse results
                sample_df = self.data_processor.parse_jellyfish_output(count_file, sample_name)
                
                # 6.5 立即与BED合并并保存 | Immediately merge with BED and save
                merged_sample_df = self.data_processor.merge_single_sample_with_bed(sample_df, sample_name)
                
                # 创建each_sample子目录
                each_sample_dir = self.config.output_dir / "each_sample"
                each_sample_dir.mkdir(exist_ok=True)
                
                sample_file = each_sample_dir / f"{sample_name}_kmer_abundance.tsv"
                merged_sample_df.to_csv(sample_file, sep='\t', index=False)
                sample_files.append(sample_file)
                self.logger.info(f"✅ 样品{sample_name}结果已保存: {sample_file}")
                
                # 清理解压文件，但保留jf和count文件
                for f in prepared_files:
                    if f and os.path.exists(f) and f != file1 and f != file2:  # 只删除临时解压文件
                        os.remove(f)
            
            # 7. 合并所有样本结果 | Merge all sample results
            self.logger.info("🔗 步骤7: 合并所有样本结果 | Step 7: Merging all sample results")
            final_df = self._merge_sample_files(sample_files)
            
            # 8. 保存最终结果 | Save final results
            self.logger.info("💾 步骤8: 保存最终结果 | Step 8: Saving final results")
            output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
            final_df.to_csv(output_file, sep='\t', index=False)
            self.logger.info(f"✅ 最终结果已保存 | Final results saved: {output_file}")
                        
            # 其余步骤保持不变...
            # 9. 滑动窗口分析等后续处理
            
        except Exception as e:
            self.logger.error(f"💥 分析过程中出现错误 | Error during analysis: {e}")
            raise
    
    def _generate_summary_report(self, final_df: pd.DataFrame):
        """📊 生成统计摘要报告 | Generate summary statistics report"""
        self.logger.info("📊 生成统计摘要报告 | Generating summary statistics report")
        
        # 获取样本列 | Get sample columns
        sample_columns = [col for col in final_df.columns 
                         if col not in ['chr', 'start', 'end', 'kmer', 'kmer_id'] 
                         and not col.startswith('tmp_')]
        
        report_lines = []
        report_lines.append("🧬 K-mer丰度分析统计报告 | K-mer Abundance Analysis Statistics Report")
        report_lines.append("="*80)
        report_lines.append("")
        
        # 基本统计 | Basic statistics
        report_lines.append("📈 基本统计信息 | Basic Statistics:")
        report_lines.append(f"  🧬 总k-mer数量 | Total k-mers: {len(final_df):,}")
        report_lines.append(f"  👥 样本数量 | Number of samples: {len(sample_columns)}")
        report_lines.append(f"  📏 K-mer长度 | K-mer size: {self.config.kmer_size}")
        report_lines.append("")
        
        # 样本统计 | Sample statistics
        report_lines.append("👤 样本统计信息 | Sample Statistics:")
        for sample in sample_columns:
            total_abundance = final_df[sample].sum()
            present_kmers = (final_df[sample] > 0).sum()
            avg_abundance = final_df[sample].mean()
            report_lines.append(f"  📊 {sample}:")
            report_lines.append(f"    💯 总丰度 | Total abundance: {total_abundance:,}")
            report_lines.append(f"    ✅ 存在k-mer数 | Present k-mers: {present_kmers:,}")
            report_lines.append(f"    📊 平均丰度 | Average abundance: {avg_abundance:.2f}")
            report_lines.append(f"    📈 存在比例 | Presence ratio: {present_kmers/len(final_df):.2%}")
        report_lines.append("")
        
        # 保存报告 | Save report
        report_file = self.config.output_dir / "analysis_summary.txt"
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(report_lines))
        
        self.logger.info(f"✅ 统计报告已保存 | Summary report saved: {report_file}")

def create_parser():
    """🔧 创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='🧬 K-mer丰度分析工具 | K-mer abundance analysis tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
🌟 示例用法 | Example usage:
  %(prog)s -i /data/fastq -p "*_1.fq.gz" -k kmers.fasta -o results/
  %(prog)s -i ./samples -p "*_R1.fastq" -k kmers.fasta -b kmers.bed -w 500000 -o results/
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('📋 必需参数 | Required arguments')
    # required.add_argument('-i', '--input', required=True,
    #                     help='📁 FASTQ文件输入目录 | FASTQ files input directory')
    # required.add_argument('-p', '--pattern', required=True,
    #                     help='📁 FASTQ文件模式，如*_1.fq.gz | FASTQ file pattern, e.g. *_1.fq.gz')
    required.add_argument('-i', '--input', required=True,
                    help='📁 输入文件目录 | Input files directory')
    required.add_argument('-p', '--pattern', required=True,
                        help='📁 文件模式，支持FASTQ和FASTA格式，如*_1.fq.gz或*.fasta | File pattern, support FASTQ and FASTA formats, e.g. *_1.fq.gz or *.fasta')
    required.add_argument('-k', '--kmer-lib', required=True,
                        help='🧬 K-mer库文件(FASTA格式) | K-mer library file (FASTA format)')
    required.add_argument('-o', '--output', required=True,
                        help='📂 输出目录 | Output directory')
    
    # 可选参数 | Optional arguments
    optional = parser.add_argument_group('⚙️ 可选参数 | Optional arguments')
    optional.add_argument('-b', '--bed-file',
                        help='📋 BED文件路径 | BED file path')
    optional.add_argument('-m', '--kmer-size', type=int, default=51,
                        help='📏 K-mer长度 (默认: %(default)s) | K-mer size (default: %(default)s)')
    optional.add_argument('-s', '--hash-size', default='1000M',
                        help='🗂️ 哈希表大小 (默认: %(default)s) | Hash table size (default: %(default)s)')
    optional.add_argument('-t', '--threads', type=int, default=8,
                        help='🧵 线程数 (默认: %(default)s) | Number of threads (default: %(default)s)')
    optional.add_argument('-w', '--window-size', type=int, default=500000,
                        help='🪟 滑动窗口大小bp (默认: %(default)s) | Sliding window size in bp (default: %(default)s)')
    optional.add_argument('--step-size', type=int,
                        help='👣 滑动窗口步长bp (默认: window-size/5) | Sliding window step size in bp (default: window-size/5)')
    optional.add_argument('-C', '--canonical', action='store_true',
                        help='🔄 统计正向和反向互补链 | Count both forward and reverse complement')
    optional.add_argument('--keep-temp', action='store_true',
                        help='💾 保留临时文件 | Keep temporary files')
    optional.add_argument('--keep-binary', action='store_true',
                        help='🔢 保留0/1存在缺失矩阵 | Keep 0/1 presence/absence matrix')
    optional.add_argument('--jellyfish-path', default='jellyfish',
                        help='🐙 Jellyfish程序路径 (默认: %(default)s) | Jellyfish program path (default: %(default)s)')
    optional.add_argument('-v', '--verbose', action='store_true',
                        help='📝 详细输出 | Verbose output')
    
    return parser


def main():
    """🚀 主函数 - 这是run_kmer_count脚本的入口点 | Main function - Entry point for run_kmer_count script"""
    parser = create_parser()
    args = parser.parse_args()
    
    try:
        # 创建配置 | Create configuration
        config = KmerCountConfig(args)
        
        # 创建分析器并运行 | Create analyzer and run
        analyzer = KmerCountAnalyzer(config)
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        print("\n⏹️ 分析被用户中断 | Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 错误 | Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
