"""
func_anno: 蛋白功能注释(IPS+eggnog→GO/KEGG标准表)|Protein functional annotation
"""

import os
import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function."""
    try:
        from ...func_anno import main as func_anno_module
        return func_anno_module.main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


def _is_help_request():
    """检查是否为帮助请求|Check if help request."""
    return any(a in {"-h", "--help"} for a in sys.argv)


def _validate_file(path):
    """验证文件存在|Validate file exists."""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"文件不存在|File not found: {path}")
    return path


@click.command(
    short_help="蛋白功能注释(IPS+eggnog→GO/KEGG表)|Protein annotation (IPS+eggnog→GO/KEGG tables)",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("-i", "--input", required=True,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help="蛋白序列 FASTA|Protein FASTA")
@click.option("-o", "--output-dir", required=True, type=click.Path(),
              help="输出目录|Output directory")
@click.option("-t", "--threads", type=int, default=12, show_default=True,
              help="线程数|Threads")
@click.option("-s", "--sample-name", default=None,
              help="样本名/输出前缀(默认输入文件名)|Sample name (default: input stem)")
@click.option("--ips-result", default=None,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help="复用已有 IPS 结果目录(跳过 IPS)|Reuse existing IPS dir")
@click.option("--eggnog-result", default=None,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help="复用已有 .emapper.annotations(跳过 eggnog)|Reuse existing annotations")
@click.option("--skip-ips", is_flag=True, help="跳过 IPS(只要 GO/KEGG)|Skip IPS")
@click.option("--skip-eggnog", is_flag=True, help="跳过 eggnog|Skip eggnog")
@click.option("--kegg-map", default=None,
              callback=lambda c, p, v: _validate_file(v) if v else None,
              help="外部 KEGG 映射 TSV(补 category)|External KEGG map (fill category)")
@click.option("-m", "--mode", default="mmseqs",
              type=click.Choice(["mmseqs", "diamond", "hmmer"]), show_default=True,
              help="eggnog 搜索模式|eggnog search mode")
@click.option("--data-dir", default=None, help="eggnog DB 目录|eggnog DB dir")
@click.option("--emapper-path", default=None,
              help="emapper.py 路径覆盖|emapper.py path override")
def func_anno(input, output_dir, threads, sample_name, ips_result, eggnog_result,
              skip_ips, skip_eggnog, kegg_map, mode, data_dir, emapper_path):
    """
    蛋白功能注释: interproscan(结构域) + eggnog-mapper(GO/KEGG) → 标准 GO/KEGG 表(衔接下游 R).

    |Protein functional annotation: interproscan (domains) + eggnog-mapper (GO/KEGG)
    → standard GO/KEGG tables for downstream R enrichment.

    示例|Examples: biopytools func-anno -i proteins.fa -o out/ -t 24
    """

    func_anno_main = _lazy_import_main()

    args = ["func_anno.py", "-i", input, "-o", output_dir]
    if threads != 12:
        args.extend(["-t", str(threads)])
    if sample_name:
        args.extend(["-s", sample_name])
    if ips_result:
        args.extend(["--ips-result", ips_result])
    if eggnog_result:
        args.extend(["--eggnog-result", eggnog_result])
    if skip_ips:
        args.append("--skip-ips")
    if skip_eggnog:
        args.append("--skip-eggnog")
    if kegg_map:
        args.extend(["--kegg-map", kegg_map])
    if mode != "mmseqs":
        args.extend(["-m", mode])
    if data_dir:
        args.extend(["--data-dir", data_dir])
    if emapper_path:
        args.extend(["--emapper-path", emapper_path])

    original_argv = sys.argv
    sys.argv = args
    try:
        func_anno_main()
    except SystemExit as e:
        sys.exit(e.code)
    finally:
        sys.argv = original_argv
