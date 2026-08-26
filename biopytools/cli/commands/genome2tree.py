"""genome2tree CLI 包装器|genome2tree Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...genome2tree.main import main as genome2tree_main
        return genome2tree_main
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


@click.command(short_help="基因组目录免比对构建物种树(waster)|"
                          "Alignment-free species tree from genome dir (waster)")
@click.option("-i", "--input", "input_dir", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="基因组目录(每样本一个 fasta,可.gz)|"
                   "Genome dir (one fasta per sample, .gz ok)")
@click.option("-o", "--output-dir", default="./genome2tree_output",
              help="输出目录(默认./genome2tree_output)|Output directory")
@click.option("-t", "--threads", default=12, type=int,
              help="线程数(默认12)|Threads (default 12)")
@click.option("--root", default="",
              help="外群物种名(出有根树)|Outgroup species name (rooted tree)")
@click.option("--branch-length", is_flag=True, default=False,
              help="追加 waster_branchlength 枝长计算|Also compute branch lengths")
@click.option("--samples-map", default="",
              callback=lambda c, p, v: _validate_path_exists(v) if v else v,
              help="个体→物种映射文件(个体stem<TAB>物种名)|individual-to-species map")
@click.option("--waster-path", default=None,
              help="waster 路径(默认~/software/ASTER/bin/waster)|waster binary path")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def genome2tree(**kwargs):
    """基因组目录免比对直接构建物种进化树
    |Alignment-free species tree from a directory of genomes (waster).

    示例|Examples: biopytools genome2tree -i genome_dir/ -o results/
    """
    genome2tree_main = _lazy_import_main()
    args = ["genome2tree.py",
            "-i", kwargs["input_dir"],
            "-o", kwargs["output_dir"],
            "-t", str(kwargs["threads"]),
            "--root", kwargs["root"],
            "--log-level", kwargs["log_level"]]
    if kwargs["branch_length"]:
        args.append("--branch-length")
    if kwargs["samples_map"]:
        args.extend(["--samples-map", kwargs["samples_map"]])
    if kwargs["waster_path"]:
        args.extend(["--waster-path", kwargs["waster_path"]])
    original = sys.argv
    sys.argv = args
    try:
        rc = genome2tree_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
