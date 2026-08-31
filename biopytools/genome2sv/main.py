"""genome2sv 命令行入口|genome2sv CLI entry"""
import argparse
import sys

from .config import Genome2SVConfig
from .pipeline import Genome2SVPipeline
from .utils import ModuleLogger, check_dependencies


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="genome2sv",
        description="assembly-to-assembly SV calling|"
                    "Assembly-to-assembly SV calling (minimap2 + svim-asm + SURVIVOR)")
    parser.add_argument("-i", "--input", required=True,
                        help="样本清单 fof(name<TAB>path)|Sample manifest fof")
    parser.add_argument("-r", "--ref", required=True,
                        help="参考样本名(fof 第一列)|Reference sample name")
    parser.add_argument("-o", "--output-dir", default="./genome2sv_output",
                        help="输出目录|Output directory")
    parser.add_argument("-t", "--threads", type=int, default=12,
                        help="线程数(每样本)|Threads per sample (default 12)")
    parser.add_argument("--preset", default="asm5",
                        choices=["asm5", "asm10", "asm20"],
                        help="minimap2 预设|minimap2 preset (default asm5)")
    parser.add_argument("--svim-mode", default="haploid",
                        choices=["haploid", "diploid"],
                        help="svim-asm 模式|svim-asm mode (default haploid)")
    parser.add_argument("--max-dist", type=int, default=1000,
                        help="SURVIVOR 断点最大距离(bp)|SURVIVOR max breakpoint dist (default 1000)")
    parser.add_argument("--min-sv-length", type=int, default=50,
                        help="SURVIVOR 最小 SV 长度(bp)|SURVIVOR min SV length (default 50)")
    parser.add_argument("--survivor-type", type=int, default=1, choices=[0, 1],
                        help="SV 类型一致(1)/任意(0)|Require same type (default 1)")
    parser.add_argument("--survivor-strand", type=int, default=1, choices=[0, 1],
                        help="链方向一致(1)/任意(0)|Require same strand (default 1)")
    parser.add_argument("--min-support", type=int, default=1,
                        help="SURVIVOR 最小支持调用数|SURVIVOR min supporting callers (default 1)")
    parser.add_argument("--flank", type=int, default=300,
                        help="SV 上下游侧翼长度 bp(默认300)|SV flank length bp (default 300)")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别|Log level (default INFO)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    try:
        config = Genome2SVConfig(
            input_fof=args.input,
            ref_sample=args.ref,
            output_dir=args.output_dir,
            threads=args.threads,
            preset=args.preset,
            svim_mode=args.svim_mode,
            max_dist=args.max_dist,
            min_sv_length=args.min_sv_length,
            survivor_type=args.survivor_type,
            survivor_strand=args.survivor_strand,
            min_support=args.min_support,
            flank=args.flank,
            log_level=args.log_level,
        )
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    log_file = str(config.logs_dir / "genome2sv.log")
    logger = ModuleLogger(log_file=log_file, log_level=config.log_level).get_logger()
    logger.info("=" * 60)
    logger.info("genome2sv 启动|genome2sv started")
    logger.info(f"样本清单|Manifest: {config.input_fof}")
    logger.info(f"参考样本|Reference: {config.ref_sample}")
    logger.info(f"输出目录|Output: {config.output_dir}")

    if not check_dependencies(config, logger):
        logger.error("依赖检查未通过,终止|Dependency check failed; abort")
        return 1

    pipeline = Genome2SVPipeline(config, logger)
    return pipeline.run()


if __name__ == "__main__":
    sys.exit(main())
