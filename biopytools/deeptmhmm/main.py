"""
DeepTMHMM 1.0跨膜螺旋/信号肽预测主模块|DeepTMHMM 1.0 Main Module
"""

import argparse
import os
import sys
import time

from .config import DeeptmhmmConfig
from .utils import (
    DeeptmhmmLogger,
    run_deeptmhmm,
    parse_deeptmhmm_output,
    write_clean_tsv,
)


class DeeptmhmmPredictor:
    """DeepTMHMM预测主类|DeepTMHMM Prediction Main Class"""

    def __init__(self, **kwargs):
        # 初始化并验证配置|Init and validate config
        self.config = DeeptmhmmConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Init logging
        self.logger_manager = DeeptmhmmLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

    def _output_paths(self):
        """计算各输出文件路径|Compute output file paths"""
        prefix = self.config.output_prefix
        base = self.config.output_dir
        return {
            'tsv': os.path.join(base, f'{prefix}_deeptmhmm_summary.tsv'),
            'three_line': os.path.join(base, f'{prefix}_deeptmhmm_topologies.3line'),
            'gff3': os.path.join(base, f'{prefix}_deeptmhmm_tmr.gff3'),
        }

    def run(self) -> bool:
        """运行预测|Run prediction"""
        cfg = self.config
        log = self.logger
        paths = self._output_paths()
        start = time.time()

        log.info("=" * 60)
        log.info("DeepTMHMM跨膜螺旋/信号肽预测|DeepTMHMM TM helix & signal peptide prediction")
        log.info("=" * 60)
        log.info(f"输入文件|Input file: {cfg.input_file}")
        log.info(f"输出目录|Output dir: {cfg.output_dir}")
        log.info(f"conda环境|conda env: {cfg.conda_env}")
        log.info(f"工具目录|Tool dir: {cfg.deeptmhmm_dir}")

        # 断点续传: 主输出和原始输出齐全则跳过|Checkpoint resume
        if (os.path.exists(paths['tsv'])
                and os.path.exists(paths['three_line'])
                and os.path.exists(paths['gff3'])):
            log.info(f"跳过(已完成)|Skip (already done): {os.path.basename(paths['tsv'])}")
            return True

        # 运行predict.py(含临时目录搬运)|Run predict.py (handles temp dir move)
        ok = run_deeptmhmm(log, cfg)
        if not ok:
            return False

        # 检查必要原始输出|Check required raw outputs
        if not (os.path.exists(paths['three_line']) and os.path.exists(paths['gff3'])):
            log.error(
                f"必要输出缺失|Missing required outputs: "
                f"{paths['three_line']} / {paths['gff3']}"
            )
            return False

        # 解析并写干净TSV|Parse and write clean TSV
        records = parse_deeptmhmm_output(
            paths['three_line'], paths['gff3'], log
        )
        write_clean_tsv(records, paths['tsv'])

        # 统计|Stats
        tmh0 = sum(1 for r in records if r['pred_tmhs'] == 0)
        tmh_ge1 = len(records) - tmh0
        sp_count = sum(1 for r in records if r['signal_peptide'] != 'no')
        log.info(f"整理完成|Clean TSV: {paths['tsv']}")
        log.info(f"蛋白总数|Total proteins: {len(records)}")
        log.info(
            f"0个TMH: {tmh0}, >=1个TMH: {tmh_ge1}, "
            f"含信号肽|with signal peptide: {sp_count}"
        )

        log.info("=" * 60)
        log.info(f"全部完成|All finished, 耗时|elapsed: {time.time() - start:.2f}s")
        log.info("=" * 60)
        return True


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='DeepTMHMM跨膜螺旋/信号肽预测|DeepTMHMM TM helix & signal peptide prediction',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='示例|Examples: biopytools deeptmhmm -i proteins.fa -o output_dir'
    )

    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', required=True,
                          help='[FILE] 输入蛋白质FASTA文件|Input protein FASTA file')
    required.add_argument('-o', '--output-dir', required=True,
                          help='[DIR] 输出目录|Output directory')

    optional = parser.add_argument_group('optional arguments')
    optional.add_argument('--prefix', default=None,
                          help='[STR] 输出文件前缀(默认输入文件名)|Output file prefix')
    optional.add_argument('--conda-env', default='deeptmhmm_v.1.0',
                          help='[STR] conda环境名|conda env name')
    optional.add_argument(
        '--deeptmhmm-dir',
        default='~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0',
        help='[STR] DeepTMHMM安装目录|DeepTMHMM install directory'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        predictor = DeeptmhmmPredictor(
            input_file=args.input,
            output_dir=args.output_dir,
            output_prefix=args.prefix,
            conda_env=args.conda_env,
            deeptmhmm_dir=args.deeptmhmm_dir,
        )
        success = predictor.run()
        sys.exit(0 if success else 1)

    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)
    except KeyboardInterrupt:
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
