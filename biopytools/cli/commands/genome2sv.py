"""genome2sv CLI 包装器|genome2sv Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...genome2sv.main import main as genome2sv_main
        return genome2sv_main
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


@click.command(short_help="assembly-to-assembly SV calling (minimap2+svim-asm+SURVIVOR)")
@click.option("-i", "--input", "input_fof", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="样本清单 fof(name<TAB>path)|Sample manifest fof")
@click.option("-r", "--ref", "ref_sample", required=True,
              help="参考样本名(fof 第一列)|Reference sample name")
@click.option("-o", "--output-dir", default="./genome2sv_output",
              help="输出目录|Output directory")
@click.option("-t", "--threads", default=12, type=int,
              help="线程数(每样本,默认12)|Threads per sample (default 12)")
@click.option("--preset", default="asm5",
              type=click.Choice(["asm5", "asm10", "asm20"]),
              help="minimap2 预设(默认asm5)|minimap2 preset (default asm5)")
@click.option("--svim-mode", default="haploid",
              type=click.Choice(["haploid", "diploid"]),
              help="svim-asm 模式(默认haploid)|svim-asm mode (default haploid)")
@click.option("--max-dist", default=1000, type=int,
              help="SURVIVOR 断点最大距离bp(默认1000)|SURVIVOR max breakpoint dist (default 1000)")
@click.option("--min-sv-length", default=50, type=int,
              help="SURVIVOR 最小SV长度bp(默认50)|SURVIVOR min SV length (default 50)")
@click.option("--survivor-type", default=1, type=click.Choice(["0", "1"]),
              help="SV类型一致1/任意0(默认1)|Require same type (default 1)")
@click.option("--survivor-strand", default=1, type=click.Choice(["0", "1"]),
              help="链方向一致1/任意0(默认1)|Require same strand (default 1)")
@click.option("--min-support", default=1, type=int,
              help="SURVIVOR 最小支持调用数(默认1)|SURVIVOR min supporting callers (default 1)")
@click.option("--flank", default=300, type=int,
              help="SV 上下游侧翼长度bp(默认300)|SV flank length bp (default 300)")
@click.option("--log-level", default="INFO", help="日志级别(默认INFO)|Log level (default INFO)")
def genome2sv(**kwargs):
    """fof 清单中参考 vs 其余组装的 SV 调用与合并
    |Assembly-to-assembly SV calling: reference vs other assemblies in a fof.

    示例|Examples: biopytools genome2sv -i samples.fof -r ref -o results/
    """
    genome2sv_main = _lazy_import_main()
    args = ["genome2sv.py",
            "-i", kwargs["input_fof"],
            "-r", kwargs["ref_sample"],
            "-o", kwargs["output_dir"],
            "-t", str(kwargs["threads"]),
            "--preset", kwargs["preset"],
            "--svim-mode", kwargs["svim_mode"],
            "--max-dist", str(kwargs["max_dist"]),
            "--min-sv-length", str(kwargs["min_sv_length"]),
            "--survivor-type", str(kwargs["survivor_type"]),
            "--survivor-strand", str(kwargs["survivor_strand"]),
            "--min-support", str(kwargs["min_support"]),
            "--flank", str(kwargs["flank"]),
            "--log-level", kwargs["log_level"]]
    original = sys.argv
    sys.argv = args
    try:
        rc = genome2sv_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
