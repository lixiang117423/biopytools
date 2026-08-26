"""reads2tree 命令行入口|reads2tree CLI entry"""
import argparse
import sys

from ..common.paths import expand_path
from .calculator import Reads2TreeCalculator
from .config import Reads2TreeConfig
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="reads2tree",
        description="fastq 目录→物种进化树(WASTER from raw reads)|FASTQ dir to "
                    "species tree (WASTER from raw reads)")
    parser.add_argument("-i", "--input", required=True,
                        help="fastq 目录(自动识别双端 _R1/_R2、_1/_2 等,可.gz)|"
                             "FASTQ dir (auto-detect paired-end _R1/_R2, _1/_2, "
                             ".gz ok)")
    parser.add_argument("-o", "--output-dir", default="./reads2tree_output",
                        help="输出目录(默认./reads2tree_output)|Output directory")
    parser.add_argument("-t", "--threads", type=int, default=12,
                        help="线程数(默认12)|Threads (default 12)")
    parser.add_argument("--root", default="",
                        help="外群物种名(出有根树)|Outgroup species name")
    parser.add_argument("--branch-length", action="store_true",
                        help="追加 waster_branchlength 枝长计算|Also compute "
                             "branch lengths on the fixed topology")
    parser.add_argument("--samples-map", default="",
                        help="个体→物种映射文件(个体名<TAB>物种名)|"
                             "individual-to-species map")
    parser.add_argument("--merge", action="store_true",
                        help="重叠双端 reads 用 BBMerge 合并(默认 cat 拼接)|"
                             "BBMerge overlapping paired reads (default: cat)")
    parser.add_argument("--waster-path", default=None,
                        help="waster 路径(默认~/software/ASTER/bin/waster)|waster path")
    parser.add_argument("--bbmerge-path", default=None,
                        help="bbmerge.sh 路径(--merge 时用)|bbmerge.sh path (for --merge)")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别(默认INFO)|Log level (default INFO)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    kwargs = dict(
        input_dir=args.input,
        output_dir=args.output_dir,
        threads=args.threads,
        root=args.root,
        branch_length=args.branch_length,
        samples_map=args.samples_map,
        merge=args.merge,
        log_level=args.log_level,
    )
    if args.waster_path:
        kwargs["waster_path"] = expand_path(args.waster_path)
    if args.bbmerge_path:
        kwargs["bbmerge_path"] = expand_path(args.bbmerge_path)
    try:
        config = Reads2TreeConfig(**kwargs)
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    logger = ModuleLogger(log_file=str(config.logs_dir / "reads2tree.log"),
                          log_level=config.log_level).get_logger()
    logger.info("=" * 60)
    logger.info("reads2tree 启动|reads2tree started")
    logger.info(f"输入目录|Input dir: {config.input_dir}")
    logger.info(f"输出目录|Output: {config.output_dir}")
    logger.info(f"双端样本|Paired samples: {len(config.paired)}, "
                f"单端|single-end: {len(config.singles)}")
    return Reads2TreeCalculator(config, logger).run()


if __name__ == "__main__":
    sys.exit(main())
