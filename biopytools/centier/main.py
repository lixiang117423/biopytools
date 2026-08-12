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

    # Hi-C FASTQ 自动模式参数|Hi-C FASTQ auto mode parameters
    hic_group = parser.add_argument_group('Hi-C FASTQ 自动模式|Hi-C FASTQ auto mode')
    hic_group.add_argument(
        '-1', '--fastq-r1',
        help='Hi-C R1 FASTQ(提供即启用自动模式)|Hi-C R1 FASTQ (enables auto mode)'
    )
    hic_group.add_argument(
        '-2', '--fastq-r2',
        help='Hi-C R2 FASTQ|Hi-C R2 FASTQ'
    )
    hic_group.add_argument(
        '-g', '--genome-id',
        help='基因组ID(bowtie2索引命名,默认从基因组文件名推导)|'
             'Genome ID for bowtie2 index naming (default: derived from genome filename)'
    )
    hic_group.add_argument(
        '--restriction-enzyme',
        default='MboI',
        help='限制性内切酶|Restriction enzyme (default: MboI)'
    )
    hic_group.add_argument(
        '--bowtie2-idx',
        help='Bowtie2索引路径(默认自动建)|Bowtie2 index path (auto-built if not given)'
    )
    hic_group.add_argument(
        '--bin-sizes',
        default='100000 20000',
        help='HiC-Pro bin大小(空格分隔)|HiC-Pro bin sizes (default: 100000 20000)'
    )
    hic_group.add_argument(
        '--max-memory',
        type=int,
        default=200,
        help='HiC-Pro最大内存(GB)|HiC-Pro max memory in GB (default: 200)'
    )
    hic_group.add_argument(
        '--force-hicpro',
        action='store_true',
        help='强制重跑HiC-Pro|Force rerun HiC-Pro'
    )
    hic_group.add_argument(
        '--hic-matrix-type',
        choices=['raw', 'iced'],
        default='raw',
        help='Hi-C矩阵类型|Hi-C matrix type (default: raw)'
    )
    hic_group.add_argument(
        '--strict-chrname',
        action='store_true',
        help='染色体命名不符ChrN时中止|Abort if chr naming not ChrN'
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
        # 构建配置参数字典|Build configuration parameter dict
        kwargs = {
            'genome_fasta': args.genome,
            'centier_path': args.centier_path,
            'output_dir': args.output_dir,
            'gff_annotation': args.gff,
            'matrix1': args.matrix1,
            'matrix2': args.matrix2,
            'bed1': args.bed1,
            'bed2': args.bed2,
            'threads': args.threads,
            'kmer_size': args.kmer_size,
            'center_tolerance': args.center_tolerance,
            'step_len': args.step_len,
            'mul_cents': args.mul_cents,
            'mingap': args.mingap,
            'signal_threshold': args.signal_threshold,
            'fastq_r1': args.fastq_r1,
            'fastq_r2': args.fastq_r2,
            'genome_id': args.genome_id,
            'restriction_enzyme': args.restriction_enzyme,
            'bowtie2_idx': args.bowtie2_idx,
            'bin_sizes': args.bin_sizes,
            'max_memory_gb': args.max_memory,
            'force_hicpro': args.force_hicpro,
            'hic_matrix_type': args.hic_matrix_type,
            'strict_chrname': args.strict_chrname,
            'step': args.step
        }

        # 使用CentIERRunner统一管理配置、日志、分析和依赖检查|
        # Use CentIERRunner to manage config, logging, analysis, and dependency check
        runner = CentIERRunner(**kwargs)

        # 检查依赖|Check dependencies
        if not args.skip_dependency_check:
            if not check_dependencies(runner.config, runner.logger):
                runner.logger.error("依赖检查失败|Dependency check failed")
                sys.exit(1)

        # 运行分析|Run analysis
        result = runner.run()

        # 输出摘要|Output summary
        if args.summary and result.get('success'):
            summary = runner.get_summary()
            summary_file = runner.config.output_path / 'centier_summary.json'
            with open(summary_file, 'w') as f:
                json.dump(summary, f, indent=2)
            runner.logger.info(f"结果摘要已保存|Result summary saved to: {summary_file}")

        # 退出|Exit
        if result.get('success'):
            runner.logger.info("分析完成|Analysis completed successfully")
            sys.exit(0)
        else:
            runner.logger.error("分析失败|Analysis failed")
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
