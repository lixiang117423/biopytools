"""
🚀 主程序模块 | Main Program Module
"""

import os
import sys
import argparse
import shutil
import pandas as pd
import subprocess
from typing import List
from pathlib import Path
from .utils import setup_logger, check_dependencies, CommandRunner
from .file_processor import FileProcessor
from .jellyfish_processor import JellyfishProcessor
from .data_processor import DataProcessor
from .window_analyzer import WindowAnalyzer
from .config import KmerCountConfig

# 自定义异常类
class FileIntegrityError(Exception):
    """💥 文件完整性错误"""
    pass

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
    #     """合并样品文件（统一strand为'.'）"""
    #     self.logger.info("合并所有样本结果")
        
    #     if self.config.bed_file:
    #         # 只使用正向序列作为基础矩阵
    #         base_bed = self.data_processor.bed_info[self.data_processor.bed_info['strand'] == '+'].copy()
    #         final_df = base_bed.copy()
    #         final_df['unique_key'] = final_df['chr'] + '_' + final_df['start'].astype(str) + '_' + final_df['end'].astype(str) + '_' + final_df['kmer']
    #     else:
    #         final_df = pd.DataFrame(self.data_processor.kmer_pairs)
        
    #     # 为每个样本添加丰度列
    #     for sample_file in sample_files:
    #         sample_df = pd.read_csv(sample_file, sep='\t')
    #         sample_name = sample_file.stem.replace('_kmer_abundance', '')
            
    #         unique_abundance = dict(zip(sample_df['unique_key'], sample_df[sample_name]))
    #         final_df[sample_name] = final_df['unique_key'].map(unique_abundance).fillna(0).astype(int)
        
    #     # 统一将strand设为'.'
    #     final_df['strand'] = '.'
        
    #     # 删除unique_key列
    #     final_df = final_df.drop('unique_key', axis=1)
        
    #     self.logger.info(f"最终矩阵: {len(final_df)} 行")
    #     return final_df

    def _merge_sample_files(self, sample_files: List[Path]) -> pd.DataFrame:
        """合并样品文件（处理空文件）"""
        self.logger.info("合并所有样本结果")
        
        if self.config.bed_file:
            # 只使用正向序列作为基础矩阵
            base_bed = self.data_processor.bed_info[self.data_processor.bed_info['strand'] == '+'].copy()
            final_df = base_bed.copy()
            final_df['unique_key'] = final_df['chr'] + '_' + final_df['start'].astype(str) + '_' + final_df['end'].astype(str) + '_' + final_df['kmer']
        else:
            final_df = pd.DataFrame(self.data_processor.kmer_pairs)
        
        # 为每个样本添加丰度列
        valid_samples = 0
        for sample_file in sample_files:
            sample_name = sample_file.stem.replace('_kmer_abundance', '')
            
            try:
                # 检查文件是否为空或太小
                if sample_file.stat().st_size < 10:  # 小于10字节认为是空文件
                    self.logger.warning(f"样本文件为空或过小，跳过: {sample_name}")
                    final_df[sample_name] = 0  # 全部设为0
                    continue
                
                # 尝试读取文件
                sample_df = pd.read_csv(sample_file, sep='\t')
                
                # 检查是否有数据行
                if len(sample_df) == 0:
                    self.logger.warning(f"样本文件无数据行: {sample_name}")
                    final_df[sample_name] = 0
                    continue
                
                # 检查是否有required列
                if 'unique_key' not in sample_df.columns:
                    self.logger.warning(f"样本文件缺少unique_key列: {sample_name}")
                    final_df[sample_name] = 0
                    continue
                
                # 正常处理
                unique_abundance = dict(zip(sample_df['unique_key'], sample_df[sample_name]))
                final_df[sample_name] = final_df['unique_key'].map(unique_abundance).fillna(0).astype(int)
                valid_samples += 1
                self.logger.info(f"成功合并样本: {sample_name}")
                
            except Exception as e:
                self.logger.error(f"读取样本文件失败: {sample_name}, 错误: {e}")
                final_df[sample_name] = 0  # 出错时设为0
        
        # 统一将strand设为'.'
        final_df['strand'] = '.'
        
        # 删除unique_key列
        final_df = final_df.drop('unique_key', axis=1)
        
        self.logger.info(f"最终矩阵: {len(final_df)} 行, 有效样本: {valid_samples}")
        return final_df
    
    
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
            
            # 3. 查找输入文件 | Find input files
            self.logger.info("🔍 步骤3: 查找输入文件 | Step 3: Finding input files")
            samples = self.file_processor.find_fastq_files()
            
            # 4. 准备k-mer库 | Prepare k-mer library
            self.logger.info("🧬 步骤4: 准备k-mer库 | Step 4: Preparing k-mer library")
            kmer_lib = self.file_processor.prepare_kmer_library()
            
            # 5. 解析注释信息 | Parse annotation information
            self.logger.info("📋 步骤5: 解析注释信息 | Step 5: Parsing annotation information")
            self.data_processor.parse_kmer_library()
            self.data_processor.parse_bed_file()

            # 6. 处理每个样本 | Process each sample
            self.logger.info("🔄 步骤6: 处理样本 | Step 6: Processing samples")
            sample_files = []  # 存储每个样品的结果文件路径

            # 将file_processor引用传给jellyfish_processor，用于文件完整性检查
            self.jellyfish_processor.file_processor = self.file_processor

            total_samples = len(samples)
            processed_samples = 0
            skipped_samples = 0

            for sample_name, r1_file, r2_file in samples:
                try:
                    self.logger.info(f"👤 处理样本 | Processing sample: {sample_name}")
                    
                    # 6.1 准备文件 | Prepare files
                    fastq_files = self.file_processor.decompress_files([r1_file, r2_file])
                    
                    # 6.2 k-mer计数 | K-mer counting
                    jf_file = self.jellyfish_processor.count_kmers(sample_name, fastq_files)
                    
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
                    processed_samples += 1
                    self.logger.info(f"✅ 样品{sample_name}结果已保存: {sample_file}")
                    
                    # 清理临时文件，但保留jf和count文件
                    for f in fastq_files:
                        if f and os.path.exists(f):
                            # 只删除临时解压文件，不删除原始文件
                            if str(self.config.temp_dir) in f:
                                os.remove(f)

                except Exception as e:
                    skipped_samples += 1
                    error_msg = str(e)
                    
                    # 检查是否是SIGPIPE相关错误
                    if ("killed by signal 13" in error_msg or 
                        "Some generator commands failed" in error_msg or
                        "SIGPIPE" in error_msg.upper()):
                        self.logger.error(f"💥 样本 {sample_name} 遇到SIGPIPE错误，可能是文件损坏")
                        self.logger.warning(f"💔 错误详情: killed by signal 13 (SIGPIPE)")
                    elif "FileIntegrityError" in error_msg or "文件完整性" in error_msg:
                        self.logger.error(f"💥 样本 {sample_name} 文件完整性检查失败")
                        self.logger.warning(f"💔 错误详情: {error_msg}")
                    else:
                        self.logger.error(f"💥 样本 {sample_name} 处理失败: {error_msg}")
                    
                    self.logger.warning(f"⏭️ 跳过该样本，继续处理其他样本")
                    continue

            # 输出处理总结
            self.logger.info("="*50)
            self.logger.info(f"📊 样本处理总结 | Sample Processing Summary:")
            self.logger.info(f"  📈 总样本数: {total_samples}")
            self.logger.info(f"  ✅ 成功处理: {processed_samples}")
            self.logger.info(f"  ⏭️ 跳过样本: {skipped_samples}")
            if skipped_samples > 0:
                self.logger.warning(f"  ⚠️ 跳过比例: {skipped_samples/total_samples:.1%}")
            self.logger.info("="*50)

            # 检查是否有成功处理的样本
            if not sample_files:
                raise RuntimeError("❌ 没有样本成功处理，分析无法继续 | No samples processed successfully, analysis cannot continue")

            # 7. 合并所有样本结果 | Merge all sample results
            self.logger.info("🔗 步骤7: 合并所有样本结果 | Step 7: Merging all sample results")
            final_df = self._merge_sample_files(sample_files)
            
            # 8. 保存最终结果 | Save final results
            self.logger.info("💾 步骤8: 保存最终结果 | Step 8: Saving final results")
            output_file = self.config.output_dir / "kmer_abundance_matrix.tsv"
            final_df.to_csv(output_file, sep='\t', index=False)
            self.logger.info(f"✅ 最终结果已保存 | Final results saved: {output_file}")
            
            # 9. 滑动窗口分析 | Sliding window analysis
            if self.config.bed_file and self.config.window_size:
                self.logger.info("🪟 步骤9: 滑动窗口分析 | Step 9: Sliding window analysis")
                window_df = self.window_analyzer.sliding_window_analysis(final_df)
                if window_df is not None:
                    window_file = self.config.output_dir / "sliding_window_analysis.tsv"
                    window_df.to_csv(window_file, sep='\t', index=False)
                    self.logger.info(f"✅ 滑动窗口结果已保存 | Sliding window results saved: {window_file}")
            
            # 10. 创建0/1矩阵 | Create 0/1 matrix
            if self.config.keep_binary:
                self.logger.info("🔢 步骤10: 创建0/1矩阵 | Step 10: Creating 0/1 matrix")
                binary_df = self.window_analyzer.create_binary_matrix(final_df)
                binary_file = self.config.output_dir / "kmer_binary_matrix.tsv"
                binary_df.to_csv(binary_file, sep='\t', index=False)
                self.logger.info(f"✅ 0/1矩阵已保存 | 0/1 matrix saved: {binary_file}")
            
            # 11. 生成统计报告 | Generate statistics report
            self._generate_summary_report(final_df, processed_samples, skipped_samples)
            
            self.logger.info("="*60)
            if skipped_samples > 0:
                self.logger.warning(f"⚠️ K-mer丰度分析完成，但有{skipped_samples}个样本被跳过 | K-mer abundance analysis completed, but {skipped_samples} samples were skipped")
            else:
                self.logger.info("🎉 K-mer丰度分析完成! | K-mer abundance analysis completed!")
            self.logger.info("="*60)
            
        except Exception as e:
            self.logger.error(f"💥 分析过程中出现错误 | Error during analysis: {e}")
            raise
        finally:
            # 清理临时目录 | Clean up temporary directory
            if self.config.temp_dir and not self.config.keep_temp:
                try:
                    shutil.rmtree(self.config.temp_dir)
                    self.logger.info("🧹 临时目录已清理 | Temporary directory cleaned")
                except Exception as e:
                    self.logger.warning(f"⚠️ 清理临时目录失败 | Failed to clean temporary directory: {e}")

    def _generate_summary_report(self, final_df: pd.DataFrame, processed_samples: int = 0, skipped_samples: int = 0):
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
        
        # 处理统计 | Processing statistics
        if processed_samples > 0 or skipped_samples > 0:
            total_samples = processed_samples + skipped_samples
            report_lines.append("📈 处理统计信息 | Processing Statistics:")
            report_lines.append(f"  📊 总样本数 | Total samples: {total_samples}")
            report_lines.append(f"  ✅ 成功处理 | Successfully processed: {processed_samples}")
            report_lines.append(f"  ⏭️ 跳过样本 | Skipped samples: {skipped_samples}")
            if total_samples > 0:
                report_lines.append(f"  📈 成功率 | Success rate: {processed_samples/total_samples:.1%}")
            report_lines.append("")
        
        # 基本统计 | Basic statistics
        report_lines.append("📈 基本统计信息 | Basic Statistics:")
        report_lines.append(f"  🧬 总k-mer数量 | Total k-mers: {len(final_df):,}")
        report_lines.append(f"  👥 最终样本数量 | Final sample count: {len(sample_columns)}")
        report_lines.append(f"  📏 K-mer长度 | K-mer size: {self.config.kmer_size}")
        report_lines.append("")

        # 样本统计 | Sample statistics
        if sample_columns:
            report_lines.append("👤 样本统计信息 | Sample Statistics:")
            for sample in sample_columns:
                try:
                    # 🔥 添加这一行：确保列是数值类型
                    sample_data = pd.to_numeric(final_df[sample], errors='coerce').fillna(0)
                    
                    total_abundance = sample_data.sum()
                    present_kmers = (sample_data > 0).sum()
                    avg_abundance = sample_data.mean()
                    
                    report_lines.append(f"  📊 {sample}:")
                    report_lines.append(f"    💯 总丰度 | Total abundance: {total_abundance:,}")
                    report_lines.append(f"    ✅ 存在k-mer数 | Present k-mers: {present_kmers:,}")
                    report_lines.append(f"    📊 平均丰度 | Average abundance: {avg_abundance:.2f}")
                    report_lines.append(f"    📈 存在比例 | Presence ratio: {present_kmers/len(final_df):.2%}")
                    
                except Exception as e:
                    self.logger.warning(f"⚠️ 样本 {sample} 统计计算失败: {e}")
                    report_lines.append(f"  📊 {sample}: 统计计算失败")
            
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
  %(prog)s -i /data/fasta -p "*.fa" -k kmers.fasta -o results/
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('📋 必需参数 | Required arguments')
    required.add_argument('-i', '--input', required=True,
                        help='📁 输入文件目录 | Input files directory')
    required.add_argument('-p', '--pattern', required=True,
                        help='📁 文件模式，支持FASTQ和FASTA格式，如*_1.fq.gz、*.fa | File pattern, support FASTQ and FASTA formats, e.g. *_1.fq.gz, *.fa')
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