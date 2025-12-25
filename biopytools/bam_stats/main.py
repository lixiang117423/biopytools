"""
BAM文件批量统计报告生成器主程序 | BAM File Batch Statistics Report Generator Main Program
"""

import argparse
import sys
from pathlib import Path

# 确保可以导入模块 | Ensure module can be imported
sys.path.insert(0, str(Path(__file__).parent.parent))

from .core import BAMAnalyzer


def parse_arguments():
    """解析命令行参数 | Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="BAM 文件批量统计报告生成器 📊 (并行处理版 + Excel/CSV输出) | BAM File Batch Statistics Report Generator 📊 (Parallel Processing + Excel/CSV Output)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:

  # 基本用法 | Basic usage
  python -m biopytools.bam_stats -i ./bam_files -o report.xlsx

  # 使用指定进程数 | Use specified number of processes
  python -m biopytools.bam_stats -i ./bam_files -o report.csv -p 4

  # 指定日志目录 | Specify log directory
  python -m biopytools.bam_stats -i ./bam_files -o report.xlsx --log-dir ./logs

  # 使用相对路径 | Use relative paths
  python -m biopytools.bam_stats -i bam_directory -o results/report.xlsx
        """
    )

    parser.add_argument('-i', '--input-dir',
                        required=True,
                        help='📁 包含 BAM 文件的输入文件夹路径 | Input directory path containing BAM files')

    parser.add_argument('-o', '--output-file',
                        required=True,
                        help='📄 输出的报告文件名 (支持 .xlsx 和 .csv 格式) | Output report file name (supports .xlsx and .csv formats)')

    parser.add_argument('-p', '--processes',
                        type=int,
                        default=None,
                        help='🔢 使用的并行进程数 (默认为全部可用核心) | Number of parallel processes to use (default: all available cores)')

    parser.add_argument('--log-dir',
                        default='.',
                        help='📁 日志文件输出目录 | Log file output directory (default: current directory)')

    return parser.parse_args()


def main():
    """主函数 | Main function"""
    try:
        # 解析参数 | Parse arguments
        args = parse_arguments()

        # 配置分析器 | Configure analyzer
        analyzer_config = {
            'input_dir': args.input_dir,
            'output_file': args.output_file,
            'processes': args.processes,
            'log_dir': args.log_dir
        }

        # 创建并运行分析器 | Create and run analyzer
        analyzer = BAMAnalyzer(**analyzer_config)
        success = analyzer.run_analysis()

        if success:
            print(f"\n🎉 分析完成！报告已保存至 | Analysis completed! Report saved to: {args.output_file}")
            sys.exit(0)
        else:
            print("\n❌ 分析失败 | Analysis failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n🛑 用户中断操作 | User interrupted operation")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 程序错误 | Program error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()