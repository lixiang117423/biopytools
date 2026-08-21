"""xpclr 命令行入口|xpclr CLI entry.

解析参数 → 构造 config → validate → XpclrCalculator.run。
|argparse → config → validate → run.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path
from typing import List, Optional

from .calculator import XpclrCalculator
from .config import XpclrConfig
from .utils import ModuleLogger


def build_parser() -> argparse.ArgumentParser:
    psr = argparse.ArgumentParser(
        prog="xpclr",
        description="XP-CLR 跨群体选择信号扫描(逐染色体自动循环+全基因组合并)"
                    "|XP-CLR cross-population selection scan")
    psr.add_argument("-i", "--input", dest="input_vcf", required=True,
                     help="bgzip VCF(须带 .tbi/.csi 索引)|bgzipped VCF with tabix index")
    psr.add_argument("-a", "--samples-a", required=True,
                     help="群体A样本列表文件(每行一个ID)|Pop A sample list (one ID per line)")
    psr.add_argument("-b", "--samples-b", required=True,
                     help="群体B样本列表文件|Pop B sample list")
    psr.add_argument("-o", "--output-dir", default="xpclr_out",
                     help="输出目录|Output directory (default: xpclr_out)")
    psr.add_argument("--label", default="popA_vs_popB",
                     help="结果文件前缀|Result file prefix (default: popA_vs_popB)")
    psr.add_argument("--chroms", default=None,
                     help="逗号分隔染色体列表(默认VCF全部contig)|Comma-separated "
                          "chromosomes (default: all contigs in VCF)")
    psr.add_argument("--size", type=int, default=20000,
                     help="窗口大小bp|Window size bp (default: 20000)")
    psr.add_argument("--step", type=int, default=20000,
                     help="滑窗步长bp|Step size bp (default: 20000)")
    psr.add_argument("--maxsnps", type=int, default=200,
                     help="窗口最大SNP数|Max SNPs per window (default: 200)")
    psr.add_argument("--minsnps", type=int, default=10,
                     help="窗口最小SNP数|Min SNPs per window (default: 10)")
    psr.add_argument("--ld", type=float, default=0.95,
                     help="LD加权截断|LD cutoff (default: 0.95)")
    psr.add_argument("--phased", action="store_true", default=False,
                     help="数据已phased(更精确r2)|Data phased (more precise r2)")
    psr.add_argument("--rrate", type=float, default=1e-8,
                     help="每碱基重组率|Recombination rate per base (default: 1e-8)")
    psr.add_argument("--top-n", type=int, default=50,
                     help="Top候选窗口数|Top candidate windows (default: 50)")
    psr.add_argument("--xpclr-path",
                     default="~/miniforge3/envs/selective_sweep/bin/xpclr",
                     help="xpclr可执行路径|xpclr executable path")
    psr.add_argument("--log-level", default="INFO",
                     help="日志级别|Log level (default: INFO)")
    return psr


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    config = XpclrConfig(
        input_vcf=args.input_vcf,
        samples_a=args.samples_a,
        samples_b=args.samples_b,
        output_dir=args.output_dir,
        label=args.label,
        chroms=[c.strip() for c in args.chroms.split(",")] if args.chroms else None,
        size=args.size, step=args.step, maxsnps=args.maxsnps, minsnps=args.minsnps,
        ld=args.ld, phased=args.phased, rrate=args.rrate, top_n=args.top_n,
        xpclr_path=args.xpclr_path,
    )
    try:
        config.validate()
    except ValueError as e:
        print(f"[xpclr] 参数校验失败|Validation failed:\n{e}", file=sys.stderr)
        sys.exit(2)
    log_file = str(Path(config.output_dir) / "99_logs" / "xpclr.log")
    logger = ModuleLogger(log_file=log_file, log_level=args.log_level).get_logger()
    logger.info("xpclr 启动|xpclr start")
    ok = XpclrCalculator(config, logger).run()
    if not ok:
        logger.error("流程失败|Pipeline failed")
        sys.exit(1)
    logger.info("xpclr 完成|xpclr done")
    return 0


if __name__ == "__main__":
    main()
