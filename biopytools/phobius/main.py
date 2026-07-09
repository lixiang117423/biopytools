"""
Phobius跨膜拓扑+信号肽预测主模块|Phobius TM topology & signal peptide prediction main module
"""

import argparse
import os
import sys
import time

from .config import PhobiusConfig
from .utils import (
    PhobiusLogger,
    run_phobius,
    parse_short_output,
    parse_long_output,
    build_tsv_records,
    write_clean_tsv,
)


class PhobiusPredictor:
    """Phobius预测主类|Phobius predictor main class"""

    def __init__(self, **kwargs):
        # 初始化配置(含validate)|Init config (with validate)
        self.config = PhobiusConfig(**kwargs)
        self.config.validate()
        # 初始化日志|Init logging
        self.logger_manager = PhobiusLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

    def run(self) -> bool:
        """运行预测(带断点续传)|Run prediction (with checkpoint resume)"""
        cfg = self.config
        prefix = cfg.output_prefix
        short_file = os.path.join(cfg.output_dir, f"{prefix}.phobius.short.txt")
        long_file = os.path.join(cfg.output_dir, f"{prefix}.phobius.long.txt")
        tsv_file = os.path.join(cfg.output_dir, f"{prefix}.phobius.tsv")

        self.logger.info("=" * 60)
        self.logger.info("Phobius跨膜拓扑+信号肽预测|Phobius TM topology & signal peptide prediction")
        self.logger.info("=" * 60)
        self.logger.info(f"输入文件|Input file: {cfg.input_file}")
        self.logger.info(f"输出目录|Output dir: {cfg.output_dir}")
        self.logger.info(f"输出前缀|Prefix: {prefix}")

        # 步骤1: short(断点续传)|Step 1: short (checkpoint)
        if os.path.exists(short_file):
            self.logger.info(f"跳过(已存在)|Skip (exists): {short_file}")
        else:
            if not run_phobius(self.logger, cfg.phobius_path, "-short", cfg.input_file, short_file):
                return False

        # 步骤2: long(断点续传)|Step 2: long (checkpoint)
        if os.path.exists(long_file):
            self.logger.info(f"跳过(已存在)|Skip (exists): {long_file}")
        else:
            if not run_phobius(self.logger, cfg.phobius_path, "-long", cfg.input_file, long_file):
                return False

        # 步骤3: 解析合并写TSV|Step 3: parse, merge, write TSV
        short_map = parse_short_output(short_file)
        long_map = parse_long_output(long_file)
        records = build_tsv_records(short_map, long_map)
        write_clean_tsv(records, tsv_file)

        # 统计|Stats
        n_total = len(records)
        n_sp = sum(1 for r in records if r["sp"] == "Y")
        n_tm = sum(1 for r in records if r["tm"] > 0)
        self.logger.info(f"整理完成|Clean TSV: {tsv_file}")
        self.logger.info(f"蛋白总数|Total proteins: {n_total}")
        self.logger.info(f"有信号肽|With signal peptide: {n_sp}")
        self.logger.info(f"有跨膜区(TM>0)|With TM helices: {n_tm}")
        self.logger.info("=" * 60)
        self.logger.info("Phobius预测完成|Phobius prediction completed")
        self.logger.info("=" * 60)
        return True


def main():
    """命令行入口|CLI entry"""
    parser = argparse.ArgumentParser(
        description="Phobius跨膜拓扑+信号肽预测|Phobius TM topology & signal peptide prediction",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Examples: biopytools phobius -i proteins.fa -o output_dir/",
    )
    parser.add_argument("-i", "--input", required=True,
                        help="[FILE] 输入蛋白质FASTA|Input protein FASTA")
    parser.add_argument("-o", "--output-dir", required=True,
                        help="[DIR] 输出目录|Output directory")
    parser.add_argument("--prefix", default=None,
                        help="[STR] 输出前缀(默认输入文件名)|Output prefix (default: input filename)")
    parser.add_argument("--phobius-path",
                        default="~/miniforge3/envs/phobius_v.1.0.1/bin/phobius.pl",
                        help="phobius.pl路径|phobius.pl path")
    args = parser.parse_args()

    try:
        predictor = PhobiusPredictor(
            input_file=args.input,
            output_dir=args.output_dir,
            output_prefix=args.prefix,
            phobius_path=args.phobius_path,
        )
        success = predictor.run()
        sys.exit(0 if success else 1)
    except KeyboardInterrupt:
        sys.exit(130)
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"运行失败|Run failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
