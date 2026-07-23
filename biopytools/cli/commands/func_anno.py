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


def _validate_path(path):
    """验证文件或目录存在|Validate file or dir exists."""
    if not _is_help_request() and path and not os.path.exists(path):
        raise click.BadParameter(f"路径不存在|Path not found: {path}")
    return path


@click.command(
    short_help="蛋白功能注释(IPS+eggnog→GO/KEGG表)|Protein annotation (IPS+eggnog→GO/KEGG tables)",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("-i", "--input", required=True,
              callback=lambda c, p, v: _validate_path(v) if v else None,
              help="蛋白序列 FASTA(单文件→by-step) 或目录(多样本→by-sample)"
              "|Protein FASTA (single→by-step) or dir (multi→by-sample)")
@click.option("-o", "--output-dir", required=True, type=click.Path(),
              help="输出目录|Output directory")
@click.option("-t", "--threads", type=int, default=12, show_default=True,
              help="线程数|Threads")
@click.option("-s", "--sample-name", default=None,
              help="样本名/前缀(仅单文件模式生效)|Sample name (single-file mode only)")
@click.option("--by-sample", is_flag=True,
              help="强制 by-sample(单文件也建 sample 子目录, 往同一 -o 多次跑不覆盖)"
              "|Force by-sample layout")
@click.option("--ips-result", default=None,
              callback=lambda c, p, v: _validate_path(v) if v else None,
              help="复用已有 IPS 结果目录(跳过 IPS)|Reuse existing IPS dir")
@click.option("--eggnog-result", default=None,
              callback=lambda c, p, v: _validate_path(v) if v else None,
              help="复用已有 .emapper.annotations(跳过 eggnog)|Reuse existing annotations")
@click.option("--skip-ips", is_flag=True, help="跳过 IPS(只要 GO/KEGG)|Skip IPS")
@click.option("--skip-eggnog", is_flag=True, help="跳过 eggnog|Skip eggnog")
@click.option("--kegg-map", default=None,
              callback=lambda c, p, v: _validate_path(v) if v else None,
              help="外部 KEGG 映射 TSV(补 category)|External KEGG map (fill category)")
@click.option("--kegg-exclude-keywords", default=None,
              help="KEGG 通路 name 黑名单(逗号分隔子串, None=内置植物无关词 cancer/estrogen 等)"
              "|KEGG name blacklist (None=built-in plant-irrelevant)")
@click.option("--kegg-exclude-categories", default="", show_default=True,
              help="排除的 KEGG 分类(逗号分隔, 匹配 category A/B 子串, 如 Human Diseases)"
              "|Exclude KEGG categories (substring match)")
@click.option("-m", "--mode", default="mmseqs",
              type=click.Choice(["mmseqs", "diamond", "hmmer"]), show_default=True,
              help="eggnog 搜索模式|eggnog search mode")
@click.option("--data-dir", default=None, help="eggnog DB 目录|eggnog DB dir")
@click.option("--emapper-path", default=None,
              help="emapper.py 路径覆盖|emapper.py path override")
def func_anno(input, output_dir, threads, sample_name, by_sample, ips_result,
              eggnog_result, skip_ips, skip_eggnog, kegg_map, kegg_exclude_keywords,
              kegg_exclude_categories, mode, data_dir, emapper_path):
    """
    蛋白功能注释: interproscan(结构域) + eggnog-mapper(GO/KEGG) → 标准 GO/KEGG 表(衔接下游 R).

    |Protein functional annotation: interproscan (domains) + eggnog-mapper (GO/KEGG)
    → standard GO/KEGG tables for downstream R enrichment.

    -i 单文件 → by-step(不嵌套); -i 目录 → 多样本 by-sample(每文件一子目录).
    |-i single file → by-step; -i dir → multi-sample by-sample.

    示例|Examples:
      单样本|single: biopytools func-anno -i proteins.fa -o out/ -t 24
      多样本|multi:  biopytools func-anno -i proteins_dir/ -o out/ -t 24
    """

    func_anno_main = _lazy_import_main()

    args = ["func_anno.py", "-i", input, "-o", output_dir]
    if threads != 12:
        args.extend(["-t", str(threads)])
    if sample_name:
        args.extend(["-s", sample_name])
    if by_sample:
        args.append("--by-sample")
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
    if kegg_exclude_keywords is not None:
        args.extend(["--kegg-exclude-keywords", kegg_exclude_keywords])
    if kegg_exclude_categories:
        args.extend(["--kegg-exclude-categories", kegg_exclude_categories])
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
