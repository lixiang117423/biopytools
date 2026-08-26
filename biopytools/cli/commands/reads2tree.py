"""reads2tree CLI 包装器|reads2tree Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...reads2tree.main import main as reads2tree_main
        return reads2tree_main
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


@click.command(short_help="fastq 目录免组装直接构建物种树(WASTER)|"
                          "De novo species tree from reads dir (WASTER)")
@click.option("-i", "--input", "input_dir", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="fastq 目录(自动识别双端,可.gz)|FASTQ dir "
                   "(auto-detect paired-end, .gz ok)")
@click.option("-o", "--output-dir", default="./reads2tree_output",
              help="输出目录(默认./reads2tree_output)|Output directory")
@click.option("-t", "--threads", default=12, type=int,
              help="线程数(默认12)|Threads (default 12)")
@click.option("--root", default="",
              help="外群物种名(出有根树)|Outgroup species name (rooted tree)")
@click.option("--branch-length", is_flag=True, default=False,
              help="追加 waster_branchlength 枝长计算|Also compute branch lengths")
@click.option("--samples-map", default="",
              callback=lambda c, p, v: _validate_path_exists(v) if v else v,
              help="个体→物种映射文件(个体名<TAB>物种名)|individual-to-species map")
@click.option("--merge", is_flag=True, default=False,
              help="重叠双端 reads 用 BBMerge 合并(默认 cat 拼接)|"
                   "BBMerge overlapping paired reads (default: cat)")
@click.option("--waster-path", default=None,
              help="waster 路径(默认~/software/ASTER/bin/waster)|waster binary path")
@click.option("--bbmerge-path", default=None,
              help="bbmerge.sh 路径(--merge 时用)|bbmerge.sh path (for --merge)")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
def reads2tree(**kwargs):
    """fastq 目录免组装免比对直接构建物种进化树(WASTER from raw reads)
    |De novo species tree directly from a directory of fastq files (WASTER).

    自动识别双端(_R1/_R2、_1/_2、read1/read2、_R1_001 等),R1+R2 拼接后喂 WASTER;
    重叠双端加 --merge 用 BBMerge 合并。
    |Paired-end fastq auto-detected and concatenated for WASTER; overlapping
    reads can be BBMerge-merged with --merge.

    示例|Examples: biopytools reads2tree -i fastq_dir/ -o results/
    """
    reads2tree_main = _lazy_import_main()
    args = ["reads2tree.py",
            "-i", kwargs["input_dir"],
            "-o", kwargs["output_dir"],
            "-t", str(kwargs["threads"]),
            "--root", kwargs["root"],
            "--log-level", kwargs["log_level"]]
    if kwargs["branch_length"]:
        args.append("--branch-length")
    if kwargs["samples_map"]:
        args.extend(["--samples-map", kwargs["samples_map"]])
    if kwargs["merge"]:
        args.append("--merge")
    if kwargs["waster_path"]:
        args.extend(["--waster-path", kwargs["waster_path"]])
    if kwargs["bbmerge_path"]:
        args.extend(["--bbmerge-path", kwargs["bbmerge_path"]])
    original = sys.argv
    sys.argv = args
    try:
        rc = reads2tree_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
