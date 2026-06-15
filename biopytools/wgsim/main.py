"""
Wgsim主模块|Wgsim Main Module
"""

import os
import sys
import time
import argparse
from .config import WgsimConfig
from .utils import (
    WgsimLogger,
    get_input_files,
    format_number,
    run_wgsim,
    fix_quality_lines,
    compress_file,
)


def setup_logger(output_dir: str):
    logger_manager = WgsimLogger(output_dir)
    return logger_manager.get_logger()


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Wgsim基因组测序数据模拟|Wgsim genome sequencing simulation',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i genome.fna -o output_dir
  %(prog)s -i genome_dir/ -o output_dir -N 10000000 -1 150
        '''
    )

    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', required=True,
                         help='[FILE/DIR] 输入基因组文件或目录|Input genome file or directory')
    required.add_argument('-o', '--output-dir', required=True,
                         help='[DIR] 输出目录|Output directory')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('-N', '--num-reads', type=int, default=50000000,
                         help='[INT] 模拟reads数量，默认50000000|Number of reads to simulate, default 50000000')
    optional.add_argument('-1', '--read-length', type=int, default=150,
                         help='[INT] reads长度，默认150|Read length, default 150')
    optional.add_argument('-s', '--seed', type=int, default=0,
                         help='[INT] 随机种子，默认0|Random seed, default 0')
    optional.add_argument('-e', '--error-rate', type=float, default=0.020,
                         help='[FLOAT] 测序错误率，默认0.020|Sequencing error rate, default 0.020')
    optional.add_argument('-r', '--mutation-rate', type=float, default=0.001,
                         help='[FLOAT] 突变率，默认0.001|Mutation rate, default 0.001')
    optional.add_argument('-d', '--outer-distance', type=int, default=500,
                         help='[INT] 外部距离，默认500|Outer distance, default 500')
    optional.add_argument('-D', '--inner-distance', type=int, default=0,
                         help='[INT] 内部距离，默认0|Inner distance, default 0')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    logger = setup_logger(args.output_dir)
    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("Wgsim基因组测序数据模拟|Wgsim Genome Sequencing Simulation")
        logger.info("=" * 60)

        config = WgsimConfig(
            input_path=args.input,
            output_dir=args.output_dir,
            num_reads=args.num_reads,
            read_length=args.read_length,
            seed=args.seed,
            error_rate=args.error_rate,
            mutation_rate=args.mutation_rate,
            outer_distance=args.outer_distance,
            inner_distance=args.inner_distance,
        )
        config.validate()

        runner = WgsimRunner(config, logger)
        success = runner.run()

        if success:
            elapsed = time.time() - start_time
            logger.info("=" * 60)
            logger.info(f"全部完成|All finished, 耗时|elapsed: {elapsed:.2f}s ({elapsed/60:.2f}min)")
            logger.info("=" * 60)
            return 0
        else:
            logger.error("部分样本模拟失败|Some samples simulation failed")
            return 1

    except KeyboardInterrupt:
        logger.warning("用户中断操作|Operation interrupted by user")
        return 130
    except ValueError as e:
        logger.error(f"配置错误|Configuration error: {e}")
        return 1
    except Exception as e:
        logger.error(f"运行失败|Run failed: {e}", exc_info=True)
        return 1


class WgsimRunner:
    """Wgsim运行器|Wgsim Runner"""

    def __init__(self, config: WgsimConfig, logger=None):
        self.config = config
        if logger is None:
            logger_manager = WgsimLogger(config.output_dir)
            self.logger = logger_manager.get_logger()
        else:
            self.logger = logger

    def run(self) -> bool:
        """运行模拟|Run simulation"""
        input_files = get_input_files(self.config.input_path, self.config.extensions)

        if not input_files:
            self.logger.error(f"未找到输入文件(扩展名: {self.config.extensions})|No input files found (extensions: {self.config.extensions})")
            return False

        self.logger.info(f"找到|Found {len(input_files)} 个基因组文件|genome file(s)")
        self.logger.info(f"模拟参数|Simulation parameters:")
        self.logger.info(f"  reads数量|num reads: {format_number(self.config.num_reads)}")
        self.logger.info(f"  reads长度|read length: {self.config.read_length}bp")
        self.logger.info(f"  错误率|error rate: {self.config.error_rate}")
        self.logger.info(f"  突变率|mutation rate: {self.config.mutation_rate}")
        self.logger.info(f"  外部距离|outer distance: {self.config.outer_distance}")
        self.logger.info(f"  内部距离|inner distance: {self.config.inner_distance}")
        self.logger.info(f"  随机种子|seed: {self.config.seed}")

        total = len(input_files)
        success_count = 0
        skip_count = 0

        for i, genome_file in enumerate(input_files, 1):
            base = os.path.splitext(os.path.basename(genome_file))[0]
            out1_gz = os.path.join(self.config.output_dir, f"{base}_1.fq.gz")
            out2_gz = os.path.join(self.config.output_dir, f"{base}_2.fq.gz")

            self.logger.info(f"[{i}/{total}] 处理|Processing: {base}")

            if os.path.exists(out1_gz) and os.path.exists(out2_gz):
                self.logger.info(f"跳过(已存在)|Skip (exists): {base}")
                skip_count += 1
                continue

            tmp1 = os.path.join(self.config.output_dir, f"{base}_1.fq")
            tmp2 = os.path.join(self.config.output_dir, f"{base}_2.fq")

            ok = run_wgsim(
                self.logger,
                self.config.wgsim_path,
                genome_file,
                tmp1,
                tmp2,
                self.config.num_reads,
                self.config.read_length,
                seed=self.config.seed,
                error_rate=self.config.error_rate,
                mutation_rate=self.config.mutation_rate,
                outer_distance=self.config.outer_distance,
                inner_distance=self.config.inner_distance,
            )

            if not ok:
                self.logger.error(f"模拟失败|Simulation failed: {base}")
                for f in [tmp1, tmp2]:
                    if os.path.exists(f):
                        os.remove(f)
                continue

            if fix_quality_lines(self.logger, tmp1) and fix_quality_lines(self.logger, tmp2):
                self.logger.info(f"质量行已修复|Quality lines fixed: {base}")
            else:
                self.logger.error(f"质量行修复失败|Quality line fix failed: {base}")
                for f in [tmp1, tmp2]:
                    if os.path.exists(f):
                        os.remove(f)
                continue

            if compress_file(self.logger, tmp1) and compress_file(self.logger, tmp2):
                self.logger.info(f"完成|Done: {base}")
                success_count += 1
            else:
                self.logger.error(f"压缩失败|Compression failed: {base}")
                for f in [tmp1, tmp2, tmp1 + '.gz', tmp2 + '.gz']:
                    if os.path.exists(f):
                        os.remove(f)

        self.logger.info(f"结果统计|Summary: 成功|success={success_count}, 跳过|skip={skip_count}, 失败|fail={total - success_count - skip_count}/{total}")

        return (total - success_count - skip_count) == 0


if __name__ == "__main__":
    sys.exit(main())
