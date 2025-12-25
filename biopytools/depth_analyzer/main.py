"""
🧬 覆盖度分析主程序模块 | Depth Analysis Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import DepthConfig
from .utils import DepthLogger, CommandRunner, check_dependencies
from .depth_processor import DepthProcessor
from .statistics import DepthStatistics

class DepthAnalyzer:
    """🧬 覆盖度分析主类 | Main Depth Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = DepthConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        output_dir = Path(self.config.output_file).parent
        self.logger_manager = DepthLogger(output_dir)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, output_dir)
        
        # 初始化各个处理器 | Initialize processors
        self.depth_processor = DepthProcessor(self.config, self.logger, self.cmd_runner)
        self.statistics = DepthStatistics(self.config, self.logger)
    
    def check_dependencies(self):
        """🔍 检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """🚀 运行完整的覆盖度分析流程 | Run complete depth analysis pipeline"""
        try:
            self.logger.info("🧬 开始覆盖度分析流程 | Starting depth analysis pipeline")
            self.logger.info(f"📋 配置信息 | Configuration:")
            self.logger.info(f"  📁 输入文件数量 | Input files count: {len(self.config.input_files)}")
            self.logger.info(f"  🧬 目标染色体 | Target chromosome: {self.config.chromosome}")
            self.logger.info(f"  📍 分析区间 | Analysis region: {self.config.region}")
            self.logger.info(f"  🔢 线程数 | Threads: {self.config.threads}")
            self.logger.info(f"  📄 输出文件 | Output file: {self.config.output_file}")
            if self.config.enable_window_analysis:
                self.logger.info(f"  📊 窗口分析 | Window analysis: 启用 | Enabled")
                self.logger.info(f"    📏 窗口大小 | Window size: {self.config.window_size}bp")
                self.logger.info(f"    👣 步长 | Step size: {self.config.window_step}bp")
                window_output = self.depth_processor._get_window_output_file()
                self.logger.info(f"    📄 窗口结果文件 | Window results file: {window_output}")
            
            # Step 1: 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # Step 2: 处理覆盖度分析 | Process depth analysis
            if not self.depth_processor.process_all_files():
                raise RuntimeError("覆盖度分析失败 | Depth analysis failed")
            
            # Step 3: 生成统计报告 | Generate statistics report
            self.statistics.generate_statistics()
            
            # 分析完成 | Analysis completed
            self.logger.info("🎉 覆盖度分析流程已完成 | Depth analysis pipeline completed successfully")
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_file}")
            
            if self.config.enable_window_analysis:
                window_output = self.depth_processor._get_window_output_file()
                self.logger.info(f"窗口分析结果 | Window analysis results: {window_output}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """🚀 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 BAM/SAM文件覆盖度分析脚本 (模块化版本，自动处理索引) | BAM/SAM Depth Analysis Script (Modular Version, Auto-handle Index)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
📋 示例 | Examples:
  %(prog)s -i sample.bam -o depth_results.txt
  %(prog)s -i sample1.bam sample2.bam -o results.txt -c chr12 -r 136491092:138554123
  %(prog)s -i /path/to/bam_files/ -o ./results/ --chromosome all --threads 16
  %(prog)s -i data.bam -o ./ -c chr1 -r 1000:5000 -q 20 -Q 30
  %(prog)s -i /data/bam_dir/ /data/extra.bam -o /output/dir/
  %(prog)s -i sample.bam -o ./ --enable-windows --window-size 1000 --window-step 500
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', nargs='+', required=True, 
                       help='🎯 输入BAM/SAM文件路径或文件夹 (支持多个文件/文件夹，文件夹时自动识别所有BAM/SAM文件) | Input BAM/SAM file paths or directories (supports multiple files/directories, auto-detect all BAM/SAM files in directories)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./depth_results.txt', 
                       help='📄 输出文件路径或目录 (目录时自动生成文件名) | Output file path or directory (auto-generate filename when directory)')
    
    # 分析参数 | Analysis parameters
    parser.add_argument('-c', '--chromosome', default='all', 
                       help='🧬 目标染色体名称 (默认all表示所有染色体) | Target chromosome name (default all means all chromosomes)')
    parser.add_argument('-r', '--region', default='all', 
                       help='📍 染色体区间，格式: start:end (如 100:1235, 1-based坐标) | Chromosome region, format: start:end (e.g., 100:1235, 1-based coordinates)')
    parser.add_argument('-t', '--threads', type=int, default=88, 
                       help='🔢 线程数 | Number of threads')
    
    # samtools参数 | samtools parameters
    parser.add_argument('-q', '--quality', type=int, default=0, 
                       help='⭐ 最小碱基质量阈值 | Minimum base quality threshold')
    parser.add_argument('-Q', '--mapping-quality', type=int, default=0, 
                       help='🎯 最小比对质量阈值 | Minimum mapping quality threshold')
    
    # 工具路径 | Tool paths
    parser.add_argument('--samtools-path', default='samtools', 
                       help='🔧 samtools软件路径 | samtools software path')
    
    # 输出格式 | Output format
    parser.add_argument('--output-format', default='txt', choices=['txt'], 
                       help='📄 输出格式 | Output format')
    parser.add_argument('--compress', action='store_true', 
                       help='🗜️ 压缩输出文件 | Compress output file')
    
    # 滑窗分析参数 | Sliding window analysis parameters
    parser.add_argument('--enable-windows', action='store_true',
                       help='📊 启用滑窗分析 | Enable sliding window analysis')
    parser.add_argument('--window-size', type=int, default=1000,
                       help='📏 窗口大小(bp) | Window size (bp)')
    parser.add_argument('--window-step', type=int, default=0,
                       help='👣 窗口步长(bp)，0表示无重叠，小于窗口大小则重叠 | Window step size (bp), 0 means no overlap, less than window size means overlap')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = DepthAnalyzer(
        input_files=args.input,
        output_file=args.output,
        chromosome=args.chromosome,
        region=args.region,
        threads=args.threads,
        samtools_path=args.samtools_path,
        quality_threshold=args.quality,
        mapping_quality=args.mapping_quality,
        output_format=args.output_format,
        compress_output=args.compress,
        enable_window_analysis=args.enable_windows,
        window_size=args.window_size,
        window_step=args.window_step
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
