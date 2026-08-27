"""splitstree6 CLI 包装器|SplitsTree6 Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...splitstree6.main import main as splitstree6_main
        return splitstree6_main
    except ImportError as e:
        click.echo(f"导入错误|Import error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """是否帮助请求|Is help request"""
    return any(a in {"-h", "--help"} for a in sys.argv)


def _validate_path_exists(path):
    """校验路径存在|Validate path exists"""
    if not _is_help_request() and path and not os.path.exists(os.path.expanduser(path)):
        raise click.BadParameter(f"路径不存在|Path does not exist: {path}")
    return path


@click.command(short_help="SplitsTree6 建网建树(VCF 输入,多格式导出)|"
                          "SplitsTree6 network/tree (VCF input, multi-format export)")
@click.option("-i", "--input", "input_path", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="输入数据 .vcf/.vcf.gz(自动转距离矩阵)或其他 SplitsTree6 格式"
                   "|Input: VCF (auto-converted) or any SplitsTree6-readable file")
@click.option("-o", "--output-dir", default="./splitstree6_output",
              help="输出目录(默认./splitstree6_output)|Output directory")
@click.option("-e", "--export-formats", default=None,
              help="输出格式,逗号分隔(默认 Newick,Nexus,GML)|Export formats "
                   "comma-separated (default Newick,Nexus,GML)。可选|valid: "
                   "Newick,Nexus,GML,PlainText,Phylip,FastA,Clustal")
@click.option("-w", "--workflow", default=None,
              help="自定义 .stree6 工作流(默认内置 NeighborNet 模板)|Custom .stree6 "
                   "workflow (default built-in NeighborNet template)")
@click.option("-n", "--node-name", default="Splits",
              help="导出节点名(默认 Splits 网络节点)|Node to export (default Splits)")
@click.option("--input-format", default="",
              help="指定输入格式(默认自动识别)|Input format (default auto-detect)")
@click.option("-t", "--threads", default=12, type=int,
              help="线程数(默认12)|Threads (default 12)")
@click.option("--tools-dir", default=None,
              help="splitstree6-tools 目录(jars 所在)|splitstree6-tools dir")
@click.option("--xvfb-path", default=None,
              help="Xvfb 路径|Xvfb path (JavaFX requires a display)")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def splitstree6(input_path, output_dir, export_formats, workflow, node_name,
                input_format, threads, tools_dir, xvfb_path, log_level):
    """SplitsTree6 免比对构建进化网络/树(VCF 默认输入,多格式导出)
    |Alignment-free phylogenetic networks/trees via SplitsTree6
    (VCF-default input, multi-format export).

    示例|Examples: biopytools splitstree6 -i variants.vcf -o splitstree_out/
    """
    splitstree6_main = _lazy_import_main()
    args = ["splitstree6.py",
            "-i", input_path,
            "-o", output_dir,
            "-t", str(threads),
            "-n", node_name,
            "--log-level", log_level]
    if export_formats:
        args.extend(["-e", export_formats])
    if workflow:
        args.extend(["-w", workflow])
    if input_format:
        args.extend(["--input-format", input_format])
    if tools_dir:
        args.extend(["--tools-dir", tools_dir])
    if xvfb_path:
        args.extend(["--xvfb-path", xvfb_path])
    original = sys.argv
    sys.argv = args
    try:
        rc = splitstree6_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
