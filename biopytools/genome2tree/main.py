"""genome2tree 命令行入口|genome2tree CLI entry"""
import argparse
import sys

from ..common.paths import expand_path
from .calculator import Genome2TreeCalculator
from .config import Genome2TreeConfig
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="genome2tree",
        description="基因组目录→物种进化树(waster)|Genome dir to species tree (waster)")
    parser.add_argument("-i", "--input", required=True,
                        help="基因组目录(每样本一个 fasta,可.gz)|Genome dir "
                             "(one fasta per sample, .gz ok)")
    parser.add_argument("-o", "--output-dir", default="./genome2tree_output",
                        help="输出目录(默认./genome2tree_output)|Output directory")
    parser.add_argument("-t", "--threads", type=int, default=12,
                        help="线程数(默认12)|Threads (default 12)")
    parser.add_argument("--root", default="",
                        help="外群物种名(出有根树)|Outgroup species name")
    parser.add_argument("--branch-length", action="store_true",
                        help="追加 waster_branchlength 枝长计算|Also compute "
                             "branch lengths on the fixed topology")
    parser.add_argument("--samples-map", default="",
                        help="个体→物种映射文件(个体stem<TAB>物种名)|"
                             "individual-to-species map")
    parser.add_argument("--waster-path", default=None,
                        help="waster 路径(默认~/software/ASTER/bin/waster)|waster path")
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
        log_level=args.log_level,
    )
    if args.waster_path:
        kwargs["waster_path"] = expand_path(args.waster_path)
    try:
        config = Genome2TreeConfig(**kwargs)
        config.validate()
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    logger = ModuleLogger(log_file=str(config.logs_dir / "genome2tree.log"),
                          log_level=config.log_level).get_logger()
    logger.info("=" * 60)
    logger.info("genome2tree 启动|genome2tree started")
    logger.info(f"输入目录|Input dir: {config.input_dir}")
    logger.info(f"输出目录|Output: {config.output_dir}")
    return Genome2TreeCalculator(config, logger).run()


if __name__ == "__main__":
    sys.exit(main())
