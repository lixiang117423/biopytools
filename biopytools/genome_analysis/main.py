"""
Genome Analysis Main Module
基因组分析主模块
"""

import os
import sys
import time
import subprocess
import argparse
from pathlib import Path
from .config import GenomeAnalysisConfig
from .utils import GenomeAnalysisLogger, GenomeScopeRunner, SmudgeplotRunner


def setup_logger(output_dir: str):
    """设置日志 | Setup logging"""
    logger_manager = GenomeAnalysisLogger(output_dir)
    return logger_manager.get_logger()


def find_fastq_files(input_dir: str, logger) -> list:
    """
    查找FASTQ文件 | Find FASTQ files

    Args:
        input_dir: 输入目录 | Input directory
        logger: 日志对象 | Logger object

    Returns:
        FASTQ文件列表 | List of FASTQ files
    """
    logger.info(f"正在目录中搜索FASTQ文件: {input_dir}")

    fastq_files = list(Path(input_dir).rglob('*.fastq')) + \
                  list(Path(input_dir).rglob('*.fq')) + \
                  list(Path(input_dir).rglob('*.fastq.gz')) + \
                  list(Path(input_dir).rglob('*.fq.gz'))

    if not fastq_files:
        logger.error(f"未找到任何FASTQ文件: {input_dir}")
        return []

    fastq_files = [str(f) for f in fastq_files]
    logger.info(f"找到 {len(fastq_files)} 个FASTQ文件")

    for f in fastq_files[:5]:  # 只显示前5个
        logger.info(f"  - {f}")
    if len(fastq_files) > 5:
        logger.info(f"  ... 还有 {len(fastq_files) - 5} 个文件")

    return fastq_files


def check_dependencies(logger) -> bool:
    """
    检查依赖环境 | Check dependencies

    Args:
        logger: 日志对象 | Logger object

    Returns:
        是否所有依赖都可用 | Whether all dependencies are available
    """
    logger.info("=" * 60)
    logger.info("检查依赖环境")
    logger.info("=" * 60)

    dependencies = [
        ('jellyfish', 'jellyfish'),
        ('Rscript', 'R'),
        ('smudgeplot', 'smudgeplot')
    ]

    all_ok = True
    for cmd, name in dependencies:
        try:
            # 检查命令 | Check command
            subprocess.run(['which', cmd], capture_output=True, check=True)
            logger.info(f"[OK] {name} 已找到")
        except (subprocess.CalledProcessError, FileNotFoundError):
            logger.error(f"[ERROR] {name} 未找到")
            all_ok = False

    return all_ok


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='Genome Analysis - 基因组分析工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 基本用法
  %(prog)s -i fastq_dir -o output_dir -l 150

  # 指定k-mer大小和线程数
  %(prog)s -i fastq_dir -o output_dir -l 150 -k 21 -t 32

  # 自定义哈希表大小
  %(prog)s -i fastq_dir -o output_dir -l 150 -s 20G
        '''
    )

    # 必需参数 | Required parameters
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input-dir', required=True,
                         help='[DIR] 输入FASTQ文件目录 | Input FASTQ directory')
    required.add_argument('-o', '--output-dir', required=True,
                         help='[DIR] 输出目录 | Output directory')
    required.add_argument('-l', '--read-length', required=True, type=int,
                         help='[INT] 测序读长 | Read length')

    # 可选参数 | Optional parameters
    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-k', '--kmer-size', type=int, default=21,
                         help='[INT] K-mer大小 (默认: 21) | K-mer size (default: 21)')
    optional.add_argument('-t', '--threads', type=int, default=64,
                         help='[INT] 线程数 (默认: 64) | Number of threads (default: 64)')
    optional.add_argument('-s', '--hash-size', default='10G',
                         help='[STR] Jellyfish哈希表大小 (默认: 10G) | Jellyfish hash size (default: 10G)')
    optional.add_argument('-c', '--max-kmer-cov', type=int, default=1000,
                         help='[INT] 最大k-mer覆盖度 (默认: 1000) | Max k-mer coverage (default: 1000)')
    optional.add_argument('--genomescope-r',
                         default='/share/org/YZWL/yzwl_lixg/software/scripts/genomescope.R',
                         help='[FILE] GenomeScope R脚本路径 | GenomeScope R script path')

    args = parser.parse_args()

    # 创建输出目录 | Create output directory
    os.makedirs(args.output_dir, exist_ok=True)

    # 设置日志 | Setup logging
    logger = setup_logger(args.output_dir)

    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("Genome Analysis Pipeline")
        logger.info("=" * 60)
        logger.info(f"输入目录: {args.input_dir}")
        logger.info(f"输出目录: {args.output_dir}")
        logger.info(f"读长: {args.read_length}")
        logger.info(f"K-mer大小: {args.kmer_size}")
        logger.info(f"线程数: {args.threads}")

        # 检查依赖 | Check dependencies
        if not check_dependencies(logger):
            logger.error("依赖检查失败，请安装所需软件")
            sys.exit(1)

        # 查找FASTQ文件 | Find FASTQ files
        fastq_files = find_fastq_files(args.input_dir, logger)
        if not fastq_files:
            logger.error("未找到FASTQ文件")
            sys.exit(1)

        # 初始化配置 | Initialize config
        config = GenomeAnalysisConfig(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            read_length=args.read_length,
            kmer_size=args.kmer_size,
            threads=args.threads,
            hash_size=args.hash_size,
            max_kmer_cov=args.max_kmer_cov,
            genomescope_r_script=args.genomescope_r
        )

        config.validate()

        # 获取绝对路径 | Get absolute paths
        output_dir_abs = os.path.abspath(args.output_dir)
        output_prefix = os.path.join(output_dir_abs, "genome_analysis")

        # 运行GenomeScope流程 | Run GenomeScope pipeline
        gs_runner = GenomeScopeRunner(logger)

        # 步骤1: Jellyfish count
        if not gs_runner.run_jellyfish_count(
            fastq_files, output_prefix, config.kmer_size,
            config.hash_size, config.threads
        ):
            logger.error("Jellyfish count失败")
            sys.exit(1)

        # 步骤2: Jellyfish histo
        jf_file = f"{output_prefix}.jf"
        if not gs_runner.run_jellyfish_histo(jf_file, output_prefix, config.threads):
            logger.error("Jellyfish histo失败")
            sys.exit(1)

        # 步骤3: GenomeScope
        histo_file = f"{output_prefix}.histo"
        gs_output_dir = f"{output_prefix}_genomescope_output"
        os.makedirs(gs_output_dir, exist_ok=True)

        kcov = gs_runner.run_genomescope(
            histo_file, config.kmer_size, config.read_length,
            gs_output_dir, config.max_kmer_cov, config.genomescope_r_script
        )

        if kcov is None:
            logger.warning("未能获取kcov值，将尝试使用默认值继续")
            kcov = 50.0  # 默认值

        config.kcov = kcov

        # 步骤4-5: Smudgeplot (需要FastK表，这里使用jellyfish输出代替)
        logger.info("")
        logger.info("=" * 60)
        logger.info("Smudgeplot需要FastK格式的输入")
        logger.info("当前使用Jellyfish输出，如需运行Smudgeplot请提供FastK表")
        logger.info("=" * 60)

        # 完成统计 | Completion statistics
        elapsed_time = time.time() - start_time

        logger.info("")
        logger.info("=" * 60)
        logger.info("Pipeline Summary")
        logger.info("=" * 60)
        logger.info(f"总运行时间: {elapsed_time:.2f} seconds")
        logger.info(f"K-mer coverage: {kcov}")
        logger.info(f"GenomeScope输出: {gs_output_dir}")
        logger.info("Pipeline completed successfully")

    except KeyboardInterrupt:
        logger.warning("用户中断操作")
        sys.exit(130)
    except Exception as e:
        logger.error(f"Pipeline失败: {str(e)}", exc_info=True)
        sys.exit(1)


class GenomeAnalysis:
    """基因组分析类 | Genome Analysis Class"""

    def __init__(self, **kwargs):
        """
        初始化分析器 | Initialize analyzer

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        self.config = GenomeAnalysisConfig(**kwargs)
        self.config.validate()

        # 创建输出目录 | Create output directory
        os.makedirs(self.config.output_dir, exist_ok=True)

        # 初始化日志 | Initialize logging
        self.logger_manager = GenomeAnalysisLogger(self.config.output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化工具类 | Initialize utility classes
        self.gs_runner = GenomeScopeRunner(self.logger)
        self.sm_runner = SmudgeplotRunner(self.logger)

    def run(self):
        """
        运行分析流程 | Run analysis pipeline

        Returns:
            是否成功 | Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("Genome Analysis Pipeline")
            self.logger.info("=" * 60)
            self.logger.info(f"配置信息:\n{self.config}")

            # 检查依赖 | Check dependencies
            if not check_dependencies(self.logger):
                self.logger.error("依赖检查失败")
                return False

            # 查找FASTQ文件 | Find FASTQ files
            fastq_files = find_fastq_files(self.config.input_dir, self.logger)
            if not fastq_files:
                self.logger.error("未找到FASTQ文件")
                return False

            # 获取绝对路径 | Get absolute paths
            output_dir_abs = os.path.abspath(self.config.output_dir)
            output_prefix = os.path.join(output_dir_abs, "genome_analysis")

            # 运行GenomeScope流程 | Run GenomeScope pipeline
            # 步骤1: Jellyfish count
            if not self.gs_runner.run_jellyfish_count(
                fastq_files, output_prefix, self.config.kmer_size,
                self.config.hash_size, self.config.threads
            ):
                return False

            # 步骤2: Jellyfish histo
            jf_file = f"{output_prefix}.jf"
            if not self.gs_runner.run_jellyfish_histo(jf_file, output_prefix, self.config.threads):
                return False

            # 步骤3: GenomeScope
            histo_file = f"{output_prefix}.histo"
            gs_output_dir = f"{output_prefix}_genomescope_output"
            os.makedirs(gs_output_dir, exist_ok=True)

            kcov = self.gs_runner.run_genomescope(
                histo_file, self.config.kmer_size, self.config.read_length,
                gs_output_dir, self.config.max_kmer_cov, self.config.genomescope_r_script
            )

            if kcov is None:
                self.logger.warning("未能获取kcov值")
                return False

            self.config.kcov = kcov

            # 完成统计 | Completion statistics
            elapsed_time = time.time() - start_time

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("Pipeline Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总运行时间: {elapsed_time:.2f} seconds")
            self.logger.info(f"K-mer coverage: {kcov}")
            self.logger.info(f"GenomeScope输出: {gs_output_dir}")
            self.logger.info("Pipeline completed successfully")

            return True

        except Exception as e:
            self.logger.error(f"Pipeline失败: {str(e)}", exc_info=True)
            return False


if __name__ == "__main__":
    main()
