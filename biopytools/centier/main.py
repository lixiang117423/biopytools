"""
CentIER着丝粒鉴定命令行入口|CentIER Centromere Identification Command Line Entry
"""

import argparse
import sys
import json
from pathlib import Path
from typing import Optional

from .config import CentIERConfig
from .centier_analyzer import CentIERAnalyzer
from .utils import CentIERLogger, check_dependencies


class CentIERRunner:
    """CentIER着丝粒鉴定运行器|CentIER Centromere Identification Runner"""

    def __init__(self, **kwargs):
        """
        初始化运行器|Initialize runner

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = CentIERConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        log_file = self.config.output_path / '99_logs' / 'centier.log'
        self.logger_manager = CentIERLogger(log_file, log_level="INFO")
        self.logger = self.logger_manager.get_logger()

        # 初始化分析器|Initialize analyzer
        self.analyzer = CentIERAnalyzer(self.config, self.logger)

    def run(self):
        """运行着丝粒鉴定分析|Run centromere identification analysis"""
        return self.analyzer.run()

    def get_summary(self):
        """获取分析结果摘要|Get analysis result summary"""
        return self.analyzer.get_summary()


def parse_arguments():
    """
    解析命令行参数|Parse command line arguments

    Returns:
        argparse.Namespace: 解析后的参数|Parsed arguments
    """
    parser = argparse.ArgumentParser(
        description="CentIER着丝粒鉴定工具|CentIER Centromere Identification Tool\n"
                    "用于T2T基因组组装的着丝粒识别和注释|"
                    "Identify and annotate centromeric regions in T2T-assembled genomes",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '-i', '--genome',
        required=True,
        help='基因组FASTA文件|Genome FASTA file path'
    )

    # 路径配置|Path configuration
    parser.add_argument(
        '-o', '--output-dir',
        default='./centier_output',
        help='输出目录|Output directory (default: ./centier_output)'
    )

    parser.add_argument(
        '--centier-path',
        default='~/software/CentIER/CentIER-main',
        help='CentIER软件路径|CentIER software path (default: ~/software/CentIER/CentIER-main)'
    )

    # 可选文件|Optional files
    parser.add_argument(
        '--gff',
        help='GFF/GTF注释文件|GFF/GTF annotation file path (optional)'
    )

    parser.add_argument(
        '--matrix1',
        help='Hi-C矩阵文件(100000分辨率)|Hi-C matrix file at 100000 resolution (optional)'
    )

    parser.add_argument(
        '--matrix2',
        help='Hi-C矩阵文件(200000分辨率)|Hi-C matrix file at 200000 resolution (optional)'
    )

    parser.add_argument(
        '--bed1',
        help='Hi-C BED文件(对应matrix1)|Hi-C BED file for matrix1 (optional)'
    )

    parser.add_argument(
        '--bed2',
        help='Hi-C BED文件(对应matrix2)|Hi-C BED file for matrix2 (optional)'
    )

    # 分析参数|Analysis parameters
    parser.add_argument(
        '-t', '--threads',
        type=int,
        default=12,
        help='线程数|Number of threads (default: 12)'
    )

    parser.add_argument(
        '-k', '--kmer-size',
        type=int,
        default=21,
        help='K-mer大小|K-mer size (default: 21)'
    )

    parser.add_argument(
        '-c', '--center-tolerance',
        type=int,
        default=15,
        help='中心容差|Center tolerance (default: 15)'
    )

    parser.add_argument(
        '--step-len',
        type=int,
        default=10000,
        help='步长|Step length (default: 10000)'
    )

    parser.add_argument(
        '--mul-cents',
        action='store_true',
        help='保留所有潜在的着丝粒区域|Retain all potential centromeric regions'
    )

    parser.add_argument(
        '--mingap',
        type=int,
        default=2,
        help='最小Gap值|Minimum gap value n*100000 (default: 2)'
    )

    parser.add_argument(
        '--signal-threshold',
        type=float,
        default=0.7,
        help='信号阈值|Signal threshold (default: 0.7)'
    )

    # 步骤控制|Step control
    parser.add_argument(
        '-s', '--step',
        type=int,
        choices=[1, 2, 3, 4, 5, 6],
        help='运行指定步骤|Run only specified step (1-6)'
    )

    # 其他选项|Other options
    parser.add_argument(
        '--skip-dependency-check',
        action='store_true',
        help='跳过依赖检查|Skip dependency check'
    )

    parser.add_argument(
        '--summary',
        action='store_true',
        help='输出分析结果摘要|Output analysis result summary'
    )

    return parser.parse_args()


def main():
    """
    主函数|Main function
    """
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = CentIERConfig(
            genome_fasta=args.genome,
            centier_path=args.centier_path,
            output_dir=args.output_dir,
            gff_annotation=args.gff,
            matrix1=args.matrix1,
            matrix2=args.matrix2,
            bed1=args.bed1,
            bed2=args.bed2,
            threads=args.threads,
            kmer_size=args.kmer_size,
            center_tolerance=args.center_tolerance,
            step_len=args.step_len,
            mul_cents=args.mul_cents,
            mingap=args.mingap,
            signal_threshold=args.signal_threshold,
            step=args.step
        )

        # 验证配置|Validate configuration
        config.validate()

        # 创建日志管理器|Create logger manager
        log_file = config.output_path / '99_logs' / 'centier.log'
        logger_manager = CentIERLogger(log_file, log_level="INFO")
        logger = logger_manager.get_logger()

        # 检查依赖|Check dependencies
        if not args.skip_dependency_check:
            if not check_dependencies(config, logger):
                logger.error("依赖检查失败|Dependency check failed")
                sys.exit(1)

        # 运行分析|Run analysis
        analyzer = CentIERAnalyzer(config, logger)
        result = analyzer.run()

        # 输出摘要|Output summary
        if args.summary and result.get('success'):
            summary = analyzer.get_summary()
            summary_file = config.output_path / 'centier_summary.json'
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            logger.info(f"结果摘要已保存|Result summary saved to: {summary_file}")

        # 退出|Exit
        if result.get('success'):
            logger.info("分析完成|Analysis completed successfully")
            sys.exit(0)
        else:
            logger.error("分析失败|Analysis failed")
            sys.exit(1)

    except ValueError as e:
        print(f"配置错误|Configuration Error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
