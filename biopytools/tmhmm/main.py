"""
TMHMM主模块|TMHMM Main Module
"""

import os
import sys
import time
import argparse
from .config import TmhmmConfig
from .utils import TmhmmLogger, run_tmhmm, parse_tmhmm_output, write_clean_tsv


def setup_logger(output_dir: str):
    logger_manager = TmhmmLogger(output_dir)
    return logger_manager.get_logger()


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='TMHMM跨膜螺旋预测|TMHMM transmembrane helix prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i proteins.fa -o output_dir
        '''
    )

    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', required=True,
                         help='[FILE] 输入蛋白质FASTA文件|Input protein FASTA file')
    required.add_argument('-o', '--output-dir', required=True,
                         help='[DIR] 输出目录|Output directory')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--noplot', action='store_true', default=True,
                         help='[FLAG] 不生成图形，默认开启|No plot generation, default on')
    optional.add_argument('--plot', action='store_true', default=False,
                         help='[FLAG] 生成图形|Generate plots')
    optional.add_argument('--prefix',
                         default=None,
                         help='[STR] 输出文件前缀(默认使用输入文件名)|Output file prefix')

    args = parser.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    logger = setup_logger(args.output_dir)
    start_time = time.time()

    try:
        logger.info("=" * 60)
        logger.info("TMHMM跨膜螺旋预测|TMHMM Transmembrane Helix Prediction")
        logger.info("=" * 60)

        noplot = not args.plot

        config = TmhmmConfig(
            input_file=args.input,
            output_dir=args.output_dir,
            noplot=noplot,
            output_prefix=args.prefix,
        )
        config.validate()

        runner = TmhmmRunner(config, logger)
        success = runner.run()

        if success:
            elapsed = time.time() - start_time
            logger.info("=" * 60)
            logger.info(f"全部完成|All finished, 耗时|elapsed: {elapsed:.2f}s")
            logger.info("=" * 60)
            return 0
        else:
            logger.error("预测失败|Prediction failed")
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


class TmhmmRunner:
    """TMHMM运行器|TMHMM Runner"""

    def __init__(self, config: TmhmmConfig, logger=None):
        self.config = config
        if logger is None:
            logger_manager = TmhmmLogger(config.output_dir)
            self.logger = logger_manager.get_logger()
        else:
            self.logger = logger

    def run(self) -> bool:
        """运行预测|Run prediction"""
        raw_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_tmhmm_raw.txt"
        )
        clean_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_tmhmm.tsv"
        )

        self.logger.info(f"输入文件|Input file: {self.config.input_file}")

        if os.path.exists(raw_file) and os.path.exists(clean_file):
            self.logger.info(f"跳过(已存在)|Skip (exists): {clean_file}")
            return True

        if not os.path.exists(raw_file):
            ok = run_tmhmm(
                self.logger,
                self.config.tmhmm_path,
                self.config.input_file,
                raw_file,
                noplot=self.config.noplot,
            )
            if not ok:
                return False
            self.logger.info(f"原始输出|Raw output: {raw_file}")

        records = parse_tmhmm_output(raw_file)
        write_clean_tsv(records, clean_file)

        tmh0 = sum(1 for r in records if r['pred_tmhs'] == 0)
        tmh1 = sum(1 for r in records if r['pred_tmhs'] == 1)
        self.logger.info(f"整理完成|Clean TSV: {clean_file}")
        self.logger.info(f"蛋白总数|Total proteins: {len(records)}")
        self.logger.info(f"0个TMH: {tmh0}, 1个TMH: {tmh1}, >=2个TMH: {len(records) - tmh0 - tmh1}")

        return True


if __name__ == "__main__":
    sys.exit(main())
