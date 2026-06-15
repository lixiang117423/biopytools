"""
BAM文件批量统计分析主程序|BAM File Batch Statistics Analysis Main Program
"""

import argparse
import sys
from pathlib import Path

from .config import BAMStatsConfig
from .core import BAMAnalyzer


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="BAM文件批量统计分析工具|BAM File Batch Statistics Analysis Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i ./bam_files -o result.summary.tsv
  %(prog)s -i sample.bam -o result.summary.tsv -g reference.fa
        '''
    )

    parser.add_argument(
        '-i', '--input', required=True,
        help='BAM文件或目录|BAM file or directory containing BAM files'
    )
    parser.add_argument(
        '-o', '--output', default='bam_stats.summary.tsv',
        help='输出文件|Output file (default: bam_stats.summary.tsv)'
    )
    parser.add_argument(
        '-t', '--threads', type=int, default=12,
        help='samtools线程数|Samtools threads (default: 12)'
    )
    parser.add_argument(
        '-p', '--processes', type=int, default=16,
        help='并行处理的样本数|Max parallel samples (default: 16)'
    )
    parser.add_argument(
        '-g', '--reference',
        help='参考基因组FASTA|Reference genome FASTA (for GC bias)'
    )
    parser.add_argument(
        '--bed-file',
        help='目标区域BED文件|Target regions BED file'
    )
    parser.add_argument(
        '--min-mapq', type=int, default=20,
        help='最小MAPQ阈值|Minimum MAPQ threshold (default: 20)'
    )
    parser.add_argument(
        '--max-insert', type=int, default=1000,
        help='最大插入片段|Maximum insert size (default: 1000)'
    )
    parser.add_argument(
        '--skip-alignment', action='store_true',
        help='跳过比对统计|Skip alignment statistics'
    )
    parser.add_argument(
        '--skip-coverage', action='store_true',
        help='跳过覆盖度统计|Skip coverage statistics'
    )
    parser.add_argument(
        '--skip-sequence', action='store_true',
        help='跳过序列特征统计|Skip sequence feature statistics'
    )
    parser.add_argument(
        '--skip-insert', action='store_true',
        help='跳过插入片段统计|Skip insert size statistics'
    )
    parser.add_argument(
        '--skip-duplicate', action='store_true',
        help='跳过重复统计|Skip duplicate statistics'
    )
    parser.add_argument(
        '--skip-variation', action='store_true',
        help='跳过变异统计|Skip variation statistics'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    try:
        args = parse_arguments()

        config = BAMStatsConfig(
            input_path=args.input,
            output_file=args.output,
            reference_file=args.reference,
            bed_file=args.bed_file,
            min_mapq=args.min_mapq,
            max_insert_size=args.max_insert,
            threads=args.threads,
            max_workers=args.processes,
            skip_alignment=args.skip_alignment,
            skip_coverage=args.skip_coverage,
            skip_sequence=args.skip_sequence,
            skip_insert=args.skip_insert,
            skip_duplicate=args.skip_duplicate,
            skip_variation=args.skip_variation,
        )

        analyzer = BAMAnalyzer(config)
        success = analyzer.run_analysis()

        if success:
            print(
                f"\n分析完成|Analysis completed! "
                f"结果保存至|Results saved to: {config.output_file}"
            )
            sys.exit(0)
        else:
            print("\n分析失败|Analysis failed")
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n用户中断操作|User interrupted operation")
        sys.exit(1)
    except Exception as e:
        print(f"程序错误|Program error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
