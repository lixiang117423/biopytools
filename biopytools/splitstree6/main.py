"""splitstree6 命令行入口|SplitsTree6 CLI entry"""
import argparse
import sys

from ..common.paths import expand_path
from .calculator import Splitstree6Calculator
from .config import Splitstree6Config, EXPORT_FORMATS
from .utils import ModuleLogger


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="splitstree6",
        description="SplitsTree6 免比对建网/建树(VCF 默认输入,多格式导出)"
                    "|SplitsTree6 network/tree (VCF input by default, "
                    "multi-format export)")
    parser.add_argument("-i", "--input", required=True,
                        help="输入数据:.vcf/.vcf.gz(默认自动转距离矩阵)或其他 "
                             "SplitsTree6 支持格式(fasta/nexus/phylip/newick/…)"
                             "|Input: VCF (auto-converted) or any SplitsTree6-"
                             "readable format")
    parser.add_argument("-o", "--output-dir", default="./splitstree6_output",
                        help="输出目录(默认./splitstree6_output)|Output directory")
    parser.add_argument("-e", "--export-formats", default=None,
                        help=f"输出格式,逗号分隔(默认 Newick,Nexus,GML)|Export "
                             f"formats comma-separated (default: Newick,Nexus,"
                             f"GML)。可选|valid: {', '.join(EXPORT_FORMATS)}")
    parser.add_argument("--input-format", default="",
                        help="指定输入格式(默认自动识别)|Input format (default "
                             "auto-detect)")
    parser.add_argument("-w", "--workflow", default="",
                        help="自定义 .stree6 工作流(默认内置 NeighborNet 模板)|"
                             "Custom .stree6 workflow (default: built-in "
                             "NeighborNet template)")
    parser.add_argument("-n", "--node-name", default="Splits",
                        help="导出节点名(默认 Splits 网络节点)|Node to export "
                             "(default: Splits)")
    parser.add_argument("-t", "--threads", type=int, default=12,
                        help="线程数(默认12)|Threads (default 12)")
    parser.add_argument("--tools-dir", default=None,
                        help="splitstree6-tools 目录(jars 所在)|splitstree6-tools dir")
    parser.add_argument("--xvfb-path", default=None,
                        help="Xvfb 路径(JavaFX 需要虚拟显示)|Xvfb path")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别(默认INFO)|Log level (default INFO)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    kwargs = dict(
        input=args.input,
        output_dir=args.output_dir,
        export_formats=args.export_formats,
        threads=args.threads,
        log_level=args.log_level,
    )
    if args.input_format:
        kwargs["input_format"] = args.input_format
    if args.workflow:
        kwargs["workflow"] = expand_path(args.workflow)
    if args.node_name:
        kwargs["node_name"] = args.node_name
    if args.tools_dir:
        kwargs["splitstree_tools_dir"] = expand_path(args.tools_dir)
    if args.xvfb_path:
        kwargs["xvfb_path"] = expand_path(args.xvfb_path)

    try:
        config = Splitstree6Config(**kwargs)
    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        return 1

    logger = ModuleLogger(log_file=str(config.log_file),
                          log_level=config.log_level).get_logger()
    logger.info("=" * 60)
    logger.info("splitstree6 启动|splitstree6 started")
    logger.info(f"输入|Input: {config.input}")
    logger.info(f"输出目录|Output: {config.output_dir}")
    try:
        return Splitstree6Calculator(config, logger).run()
    except ValueError as e:
        logger.error(str(e))
        return 1


if __name__ == "__main__":
    sys.exit(main())
