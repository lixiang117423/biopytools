"""
BAM统计分析主程序模块 🚀 | BAM Statistics Analysis Main Module
"""

import argparse
import json
import sys
from pathlib import Path
from typing import List, Dict, Any
from concurrent.futures import ThreadPoolExecutor, as_completed
import threading

from .config import BAMStatsConfig
from .utils import BAMStatsLogger, CommandRunner, check_dependencies
from .sample_stats import SampleStatsGenerator
from .genome_stats import GenomeStatsGenerator

class BAMStatsAnalyzer:
    """BAM统计分析主类 🔬 | Main BAM Statistics Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = BAMStatsConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = BAMStatsLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path, self.config.threads)
        
        # 初始化分析器 | Initialize analyzers
        self.sample_analyzer = SampleStatsGenerator(self.config, self.logger, self.cmd_runner)
        self.genome_analyzer = GenomeStatsGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的BAM统计分析流程 🚀⚡ | Run complete BAM statistics analysis pipeline"""
        try:
            self.logger.info("🚀 开始BAM统计分析流程 | Starting BAM statistics analysis pipeline")
            
            # 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # 输出分析配置 | Output analysis configuration
            self._log_analysis_config()
            
            # 分析所有样品 | Analyze all samples
            all_sample_stats = []
            
            if self.config.parallel_samples and len(self.config.bam_files) > 1:
                # 多线程并行分析多个文件 | Multi-threaded parallel analysis for multiple files
                self.logger.info(f"🔥 启用多线程并行分析 {len(self.config.bam_files)} 个文件，使用 {self.config.max_workers} 个工作线程 | Enabling multi-threaded parallel analysis for {len(self.config.bam_files)} files with {self.config.max_workers} worker threads")
                all_sample_stats = self._run_parallel_analysis()
            else:
                # 单文件多线程分析或顺序分析 | Single-file multi-threaded or sequential analysis
                if len(self.config.bam_files) == 1:
                    self.logger.info("🚀 单文件多线程分析模式 | Single-file multi-threaded analysis mode")
                else:
                    self.logger.info("📝 使用单线程顺序分析 | Using single-threaded sequential analysis")
                all_sample_stats = self._run_sequential_analysis()
            
            if not all_sample_stats:
                raise RuntimeError("❌ 没有成功分析的样品 | No samples were successfully analyzed")
            
            # 生成基因组级别统计 | Generate genome-level statistics
            self.logger.info("🌍 生成基因组级别统计 | Generating genome-level statistics")
            genome_stats = self.genome_analyzer.generate_genome_stats(all_sample_stats)
            
            # 输出结果 | Output results
            self._output_results(all_sample_stats, genome_stats)
            
            # 生成总结报告 | Generate summary report
            self._generate_summary_report(all_sample_stats, genome_stats)
            
            self.logger.info("✅ BAM统计分析完成 | BAM statistics analysis completed successfully")
            
        except Exception as e:
            self.logger.error(f"分析过程中发生错误 | Error during analysis: {e}")
            raise
    
    def _run_parallel_analysis(self) -> List[Dict[str, Any]]:
        """运行并行分析 🔥 | Run parallel analysis"""
        all_sample_stats = []
        completed_count = 0
        total_files = len(self.config.bam_files)
        
        # 创建线程锁用于安全地更新进度 | Create thread lock for safe progress updates
        progress_lock = threading.Lock()
        
        def analyze_sample_wrapper(bam_file: str) -> Dict[str, Any]:
            """线程安全的样品分析包装器 | Thread-safe sample analysis wrapper"""
            nonlocal completed_count
            
            try:
                # 创建每个线程独立的分析器 | Create independent analyzer for each thread
                thread_analyzer = SampleStatsGenerator(self.config, self.logger, self.cmd_runner)
                sample_stats = thread_analyzer.analyze_sample(bam_file)
                
                with progress_lock:
                    completed_count += 1
                    self.logger.info(f"✅ 完成进度 {completed_count}/{total_files} | Progress: {completed_count}/{total_files} - {Path(bam_file).name}")
                
                return sample_stats
                
            except Exception as e:
                with progress_lock:
                    completed_count += 1
                    self.logger.error(f"❌ 样品分析失败 {completed_count}/{total_files} | Sample analysis failed {completed_count}/{total_files} - {bam_file}: {e}")
                return None
        
        # 使用线程池执行并行分析 | Use thread pool for parallel analysis
        with ThreadPoolExecutor(max_workers=self.config.max_workers) as executor:
            # 提交所有任务 | Submit all tasks
            future_to_bam = {
                executor.submit(analyze_sample_wrapper, bam_file): bam_file 
                for bam_file in self.config.bam_files
            }
            
            # 收集结果 | Collect results
            for future in as_completed(future_to_bam):
                bam_file = future_to_bam[future]
                try:
                    sample_stats = future.result()
                    if sample_stats is not None:
                        all_sample_stats.append(sample_stats)
                except Exception as e:
                    self.logger.error(f"❌ 线程执行异常 | Thread execution error for {bam_file}: {e}")
        
        return all_sample_stats
    
    def _run_sequential_analysis(self) -> List[Dict[str, Any]]:
        """运行顺序分析 📝 | Run sequential analysis"""
        all_sample_stats = []
        
        for i, bam_file in enumerate(self.config.bam_files, 1):
            self.logger.info(f"🔄 处理文件 {i}/{len(self.config.bam_files)} | Processing file {i}/{len(self.config.bam_files)}: {Path(bam_file).name}")
            
            try:
                sample_stats = self.sample_analyzer.analyze_sample(bam_file)
                all_sample_stats.append(sample_stats)
            except Exception as e:
                self.logger.error(f"❌ 样品分析失败 | Sample analysis failed for {bam_file}: {e}")
                continue
        
        return all_sample_stats
    
    def _log_analysis_config(self):
        """记录分析配置 📊 | Log analysis configuration"""
        self.logger.info("🔧 分析配置 | Analysis Configuration:")
        self.logger.info(f"  📁 输入路径 | Input path: {self.config.input_path}")
        self.logger.info(f"  📄 找到BAM文件数量 | BAM files found: {len(self.config.bam_files)}")
        self.logger.info(f"  📂 输出目录 | Output directory: {self.config.output_dir}")
        self.logger.info(f"  🎯 最小MAPQ | Min MAPQ: {self.config.min_mapq}")
        self.logger.info(f"  💎 最小碱基质量 | Min base quality: {self.config.min_base_quality}")
        self.logger.info(f"  🪟 滑窗大小 | Window size: {self.config.window_size//1000}kb")
        self.logger.info(f"  👣 滑窗步长 | Step size: {self.config.step_size//1000}kb")
        self.logger.info(f"  🔥 线程数 | Threads: {self.config.threads}")
        self.logger.info(f"  ⚡ 最大并发样品数 | Max concurrent samples: {self.config.max_workers}")
        self.logger.info(f"  🚀 并行分析 | Parallel analysis: {'启用' if self.config.parallel_samples else '禁用'} | {'Enabled' if self.config.parallel_samples else 'Disabled'}")
        
        # 记录分析模块状态 | Log analysis module status
        modules_status = []
        if not self.config.skip_alignment_stats:
            modules_status.append("🎯 比对统计 | Alignment")
        if not self.config.skip_coverage_stats:
            modules_status.append("📈 覆盖度统计 | Coverage")
        if not self.config.skip_sequence_stats:
            modules_status.append("🧬 序列特征 | Sequence")
        if not self.config.skip_insert_stats:
            modules_status.append("📏 插入片段 | Insert size")
        if not self.config.skip_duplicate_stats:
            modules_status.append("🔄 重复序列 | Duplicates")
        
        self.logger.info(f"  🧩 启用的分析模块 | Enabled modules: {', '.join(modules_status)}")
        
        # 显示预期性能信息 | Show expected performance info
        if len(self.config.bam_files) > 1 and self.config.parallel_samples:
            estimated_batches = (len(self.config.bam_files) + self.config.max_workers - 1) // self.config.max_workers
            self.logger.info(f"  ⏱️ 预计批次数 | Expected batches: {estimated_batches}")
            self.logger.info(f"  🔧 性能优化: 使用 {self.config.max_workers} 个并发线程处理 {len(self.config.bam_files)} 个文件")
        elif len(self.config.bam_files) == 1:
            self.logger.info(f"  🚀 单文件多线程优化: 使用 {self.config.threads} 个线程加速分析工具")
    
    def _output_results(self, sample_stats: List[Dict[str, Any]], genome_stats: Dict[str, Any]):
        """输出结果文件 | Output result files"""
        output_path = Path(self.config.output_dir)
        
        # 输出样品级别统计 | Output sample-level statistics
        sample_stats_file = output_path / f"{self.config.prefix}_sample_statistics.json"
        with open(sample_stats_file, 'w', encoding='utf-8') as f:
            json.dump(sample_stats, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"💾 样品统计已保存 | Sample statistics saved: {sample_stats_file}")
        
        # 输出基因组级别统计 | Output genome-level statistics
        genome_stats_file = output_path / f"{self.config.prefix}_genome_statistics.json"
        with open(genome_stats_file, 'w', encoding='utf-8') as f:
            json.dump(genome_stats, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"💾 基因组统计已保存 | Genome statistics saved: {genome_stats_file}")
        
        # 输出人类可读的总结 | Output human-readable summary
        summary_file = output_path / f"{self.config.prefix}_analysis_summary.txt"
        self._write_readable_summary(summary_file, sample_stats, genome_stats)
        
        self.logger.info(f"📄 分析总结已保存 | Analysis summary saved: {summary_file}")
    
    def _write_readable_summary(self, summary_file: Path, sample_stats: List[Dict[str, Any]], genome_stats: Dict[str, Any]):
        """写入人类可读的总结文件 | Write human-readable summary file"""
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("BAM统计分析总结报告 | BAM Statistics Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 基本信息 | Basic information
            f.write("基本信息 | Basic Information:\n")
            f.write("-" * 40 + "\n")
            f.write(f"分析的样品数量 | Samples analyzed: {len(sample_stats)}\n")
            f.write(f"输出目录 | Output directory: {self.config.output_dir}\n")
            f.write(f"分析时间 | Analysis time: {genome_stats.get('analysis_info', {}).get('analysis_timestamp', 'N/A')}\n\n")
            
            # 基因组级别统计 | Genome-level statistics
            f.write("基因组级别统计汇总 | Genome-level Statistics Summary:\n")
            f.write("-" * 50 + "\n")
            
            alignment_summary = genome_stats.get('alignment_summary', {})
            if alignment_summary:
                f.write("比对统计 | Alignment Statistics:\n")
                f.write(f"  总读段数 | Total reads: {alignment_summary.get('total_reads_across_samples', 'N/A'):,}\n")
                f.write(f"  平均比对率 | Average mapping rate: {alignment_summary.get('average_mapping_rate', 'N/A'):.2f}%\n")
                f.write(f"  比对率范围 | Mapping rate range: {alignment_summary.get('min_mapping_rate', 'N/A'):.2f}% - {alignment_summary.get('max_mapping_rate', 'N/A'):.2f}%\n\n")
            
            coverage_summary = genome_stats.get('coverage_summary', {})
            if coverage_summary:
                f.write("覆盖度统计 | Coverage Statistics:\n")
                f.write(f"  平均覆盖深度 | Average coverage: {coverage_summary.get('average_coverage_across_samples', 'N/A'):.2f}X\n")
                f.write(f"  覆盖深度范围 | Coverage range: {coverage_summary.get('min_coverage', 'N/A'):.2f}X - {coverage_summary.get('max_coverage', 'N/A'):.2f}X\n\n")
            
            # 样品级别统计 | Sample-level statistics
            f.write("样品级别统计详情 | Sample-level Statistics Details:\n")
            f.write("-" * 50 + "\n")
            
            for i, stats in enumerate(sample_stats, 1):
                f.write(f"{i}. 样品 | Sample: {stats.get('sample_name', 'Unknown')}\n")
                
                alignment = stats.get('alignment_stats', {})
                if alignment:
                    f.write(f"   总读段数 | Total reads: {alignment.get('total_reads', 'N/A'):,}\n")
                    f.write(f"   比对率 | Mapping rate: {alignment.get('mapping_rate', 'N/A'):.2f}%\n")
                
                coverage = stats.get('coverage_stats', {})
                if coverage:
                    f.write(f"   平均覆盖深度 | Average coverage: {coverage.get('mean_coverage', 'N/A'):.2f}X\n")
                
                insert_stats = stats.get('insert_stats', {})
                if insert_stats and insert_stats.get('insert_stats_available', False):
                    f.write(f"   平均插入大小 | Average insert size: {insert_stats.get('mean_insert_size', 'N/A'):.0f}bp\n")
                
                f.write("\n")
            
            # 分析说明 | Analysis notes
            f.write("分析说明 | Analysis Notes:\n")
            f.write("-" * 30 + "\n")
            f.write("• 详细的JSON格式统计数据请查看对应的JSON文件 | Detailed statistics in JSON format are available in corresponding JSON files\n")
            f.write("• 基因组统计文件包含跨样品的汇总信息 | Genome statistics file contains cross-sample summary information\n")
            f.write("• 样品统计文件包含每个样品的详细分析结果 | Sample statistics file contains detailed analysis results for each sample\n")
    
    def _generate_summary_report(self, sample_stats: List[Dict[str, Any]], genome_stats: Dict[str, Any]):
        """生成总结报告 | Generate summary report"""
        output_path = Path(self.config.output_dir)
        report_file = output_path / f"{self.config.prefix}_analysis_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("BAM统计分析完成报告 | BAM Statistics Analysis Completion Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("输出文件列表 | Output Files:\n")
            f.write("-" * 30 + "\n")
            
            # 主要输出文件 | Main output files
            f.write("主要统计文件 | Main Statistics Files:\n")
            f.write(f"  - {self.config.prefix}_sample_statistics.json: 样品级别统计 | Sample-level statistics\n")
            f.write(f"  - {self.config.prefix}_genome_statistics.json: 基因组级别统计 | Genome-level statistics\n")
            f.write(f"  - {self.config.prefix}_analysis_summary.txt: 人类可读总结 | Human-readable summary\n")
            
            # 日志文件 | Log files
            f.write(f"\n日志文件 | Log Files:\n")
            f.write(f"  - bam_stats.log: 分析过程日志 | Analysis process log\n")
            
            f.write(f"\n输出目录 | Output directory: {self.config.output_dir}\n")
        
        self.logger.info(f"📋 分析报告已生成 | Analysis report generated: {report_file}")

def create_argument_parser():
    """创建参数解析器 | Create argument parser"""
    parser = argparse.ArgumentParser(
        description="BAM文件统计分析工具 🔬📊 | BAM File Statistics Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法 | Example Usage:
  # 分析单个BAM文件（多线程模式）| Analyze single BAM file (multi-threaded mode)
  python -m bam_stats -i sample.bam -o results --threads 88
  
  # 分析文件夹中的所有BAM文件 | Analyze all BAM files in directory  
  python -m bam_stats --input /path/to/bam/files/ --output results
  
  # 自定义参数分析（高性能模式）| Custom parameter analysis (high performance mode)
  python -m bam_stats -i data/ -o results --min-mapq 30 --threads 88
  
  # 自定义滑窗参数 | Custom sliding window parameters
  python -m bam_stats -i sample.bam -o results --window-size 2000000 --step-size 200000
  
  # 限制并发样品数以控制资源使用 | Limit concurrent samples to control resource usage
  python -m bam_stats -i data/ -o results --threads 88 --max-workers 8
  
  # 禁用并行处理（内存受限环境）| Disable parallel processing (memory-constrained environment)
  python -m bam_stats -i data/ -o results --no-parallel
  
  # 跳过某些分析模块 | Skip certain analysis modules
  python -m bam_stats -i sample.bam -o results --skip-coverage --skip-duplicates
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', dest='input_path', required=True,
                        help='📁 输入BAM文件或包含BAM文件的目录 | Input BAM file or directory containing BAM files')
    
    parser.add_argument('-o', '--output', dest='output_dir', required=True,
                        help='📂 输出目录 | Output directory')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-p', '--prefix', default='sample',
                        help='输出文件前缀 | Output file prefix (default: sample)')
    
    parser.add_argument('-r', '--reference', dest='reference_file',
                        help='参考基因组文件 | Reference genome file')
    
    parser.add_argument('-b', '--bed', dest='bed_file',
                        help='BED文件，用于限制分析区域 | BED file for restricting analysis regions')
    
    # 质量参数 | Quality parameters
    parser.add_argument('--min-mapq', type=int, default=20,
                        help='最小MAPQ阈值 | Minimum MAPQ threshold (default: 20)')
    
    parser.add_argument('--min-base-quality', type=int, default=20,
                        help='最小碱基质量阈值 | Minimum base quality threshold (default: 20)')
    
    parser.add_argument('--max-insert-size', type=int, default=1000,
                        help='📏 最大插入片段大小 | Maximum insert size (default: 1000)')
    
    # 覆盖度计算参数 📊 | Coverage calculation parameters
    parser.add_argument('--window-size', type=int, default=1000000,
                        help='🪟 滑窗大小（bp）| Sliding window size in bp (default: 1000000 = 1Mb)')
    
    parser.add_argument('--step-size', type=int, default=100000,
                        help='👣 滑窗步长（bp）| Sliding window step size in bp (default: 100000 = 100kb)')
    
    # 分析模块开关 | Analysis module switches
    parser.add_argument('--skip-alignment', action='store_true', dest='skip_alignment_stats',
                        help='跳过比对统计分析 | Skip alignment statistics analysis')
    
    parser.add_argument('--skip-coverage', action='store_true', dest='skip_coverage_stats',
                        help='跳过覆盖度统计分析 | Skip coverage statistics analysis')
    
    parser.add_argument('--skip-sequence', action='store_true', dest='skip_sequence_stats',
                        help='跳过序列特征分析 | Skip sequence feature analysis')
    
    parser.add_argument('--skip-insert', action='store_true', dest='skip_insert_stats',
                        help='跳过插入片段分析 | Skip insert size analysis')
    
    parser.add_argument('--skip-duplicates', action='store_true', dest='skip_duplicate_stats',
                        help='跳过重复序列分析 | Skip duplicate analysis')
    
    parser.add_argument('--skip-chromosomes', action='store_true', dest='skip_chromosome_stats',
                        help='跳过染色体统计分析 | Skip chromosome statistics analysis')
    
    # 可视化选项 | Visualization options
    parser.add_argument('--plots', action='store_true', dest='generate_plots',
                        help='生成统计图表 | Generate statistical plots')
    
    parser.add_argument('--plot-format', choices=['png', 'pdf', 'svg'], default='png',
                        help='图表格式 | Plot format (default: png)')
    
    # 工具路径 | Tool paths
    parser.add_argument('--samtools-path', default='samtools',
                        help='samtools路径 | samtools path (default: samtools)')
    
    parser.add_argument('--bedtools-path', default='bedtools',
                        help='bedtools路径 | bedtools path (default: bedtools)')
    
    # 性能参数 🚀 | Performance parameters
    parser.add_argument('-t', '--threads', type=int, default=88,
                        help='🔥 线程数 | Number of threads (default: 88)')
    
    parser.add_argument('--max-workers', type=int, default=16,
                        help='⚡ 最大并发样品处理数 | Maximum concurrent sample processing (default: 16)')
    
    parser.add_argument('--no-parallel', action='store_false', dest='parallel_samples',
                        help='🚫 禁用并行样品分析 | Disable parallel sample analysis')
    
    parser.add_argument('--memory', dest='memory_limit', default='4G',
                        help='💾 内存限制 | Memory limit (default: 4G)')
    
    return parser

def main():
    """主函数 | Main function"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    # 创建分析器实例 | Create analyzer instance
    try:
        analyzer = BAMStatsAnalyzer(**vars(args))
        analyzer.run_analysis()
        
        print("\n🎉 分析完成！| Analysis completed!")
        print(f"📊 结果已保存到 | Results saved to: {args.output_dir}")
        
    except KeyboardInterrupt:
        print("\n⚠️ 用户中断分析 | Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n❌ 分析失败 | Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
