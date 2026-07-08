"""
eggnog-mapper 功能注释命令|eggnog-mapper functional annotation command
"""

import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function."""
    try:
        from ...eggnog_mapper.main import main as emapper_main
        return emapper_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help="eggNOG功能注释(GO/KEGG/COG/CAZy/Pfam)|eggNOG functional annotation (GO/KEGG/COG/CAZy/Pfam)",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("-i", "--input", required=True,
              help="输入FASTA(蛋白/CDS/基因组)|Input FASTA")
@click.option("-o", "--output-dir", required=True, type=click.Path(),
              help="输出目录|Output directory")
@click.option("--itype", default="proteins",
              type=click.Choice(["proteins", "CDS", "genome", "metagenome"]),
              show_default=True, help="输入类型|Input type")
@click.option("--translate", is_flag=True,
              help="CDS翻译为蛋白|Translate CDS")
@click.option("-m", "--mode", default="mmseqs",
              type=click.Choice(["mmseqs", "diamond", "hmmer", "no_search", "cache"]),
              show_default=True, help="搜索模式|Search mode")
@click.option("--cpu", default=12, show_default=True, help="线程数|Threads")
@click.option("--sensmode", default="sensitive", show_default=True,
              help="灵敏度|Sensitivity")
@click.option("--seed-ortholog-evalue", default=0.001, show_default=True,
              help="seed ortholog E值|evalue")
@click.option("--data-dir", default=None,
              help="DB目录|DB directory (default: ~/database/eggnog)")
@click.option("--prefix", default=None, help="输出前缀|Output prefix")
@click.option("--emapper-path", default=None,
              help="emapper.py路径|emapper.py path override")
@click.option("--resume", is_flag=True, help="续传|Resume")
@click.option("--override", is_flag=True, help="覆盖|Override existing output")
@click.option("--no-format", is_flag=True,
              help="跳过重排版|Skip reformat")
def eggnog_mapper(input, output_dir, itype, translate, mode, cpu, sensmode,
                  seed_ortholog_evalue, data_dir, prefix, emapper_path,
                  resume, override, no_format):
    """
    eggNOG功能注释|eggNOG functional annotation.

    预测蛋白序列的 GO/KEGG/COG/CAZy/Pfam 等功能注释。
    Predict GO/KEGG/COG/CAZy/Pfam functional annotations for protein sequences.

    示例|Examples: biopytools eggnog-mapper -i proteins.faa -o out/
    """

    emapper_main = _lazy_import_main()

    args = ["eggnog_mapper.py"]
    args.extend(["-i", input])
    args.extend(["-o", output_dir])
    if itype != "proteins":
        args.extend(["--itype", itype])
    if translate:
        args.append("--translate")
    if mode != "mmseqs":
        args.extend(["-m", mode])
    if cpu != 12:
        args.extend(["--cpu", str(cpu)])
    if sensmode != "sensitive":
        args.extend(["--sensmode", sensmode])
    if seed_ortholog_evalue != 0.001:
        args.extend(["--seed-ortholog-evalue", str(seed_ortholog_evalue)])
    if data_dir:
        args.extend(["--data-dir", data_dir])
    if prefix:
        args.extend(["--prefix", prefix])
    if emapper_path:
        args.extend(["--emapper-path", emapper_path])
    if resume:
        args.append("--resume")
    if override:
        args.append("--override")
    if no_format:
        args.append("--no-format")

    original_argv = sys.argv
    sys.argv = args
    try:
        emapper_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
