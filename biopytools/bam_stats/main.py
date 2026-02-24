"""
BAM文件批量统计报告生成器主程序|BAM File Batch Statistics Report Generator Main Program
"""

import argparse
import sys
from pathlib import Path

# 确保可以导入模块|Ensure module can be imported
sys.path.insert(0, str(Path(__file__).parent.parent))

from .core import BAMAnalyzer


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="BAM文件批量统计报告生成器|BAM File Batch Statistics Report Generator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i ./bam_files -o report.xlsx
  %(prog)s -i sample.bam -o report.xlsx
        '''
    )

    parser.add_argument('-i', '--input-dir',
                        required=True,
                        help='BAM文件或包含BAM文件的目录|BAM file or directory containing BAM files')

    parser.add_argument('-o', '--output-file',
                        required=True,
                        help='输出的报告文件名|Output report file name (supports .xlsx and .csv formats)')

    parser.add_argument('-p', '--processes',
                        type=int,
                        default=None,
                        help='使用的并行进程数|Number of parallel processes to use')

    parser.add_argument('--log-dir',
                        default='.',
                        help='日志文件输出目录|Log file output directory')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    try:
        # 解析参数|Parse arguments
        args = parse_arguments()

        # 配置分析器|Configure analyzer
        analyzer_config = {
            'input_dir': args.input_dir,
            'output_file': args.output_file,
            'processes': args.processes,
            'log_dir': args.log_dir
        }

        # 创建并运行分析器|Create and run analyzer
        analyzer = BAMAnalyzer(**analyzer_config)
        success = analyzer.run_analysis()

        if success:
            print(f"分析完成|Analysis completed! 报告已保存至|Report saved to: {args.output_file}")
            sys.exit(0)
        else:
            print("分析失败|Analysis failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("用户中断操作|User interrupted operation")
        sys.exit(1)
    except Exception as e:
        print(f"程序错误|Program error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
