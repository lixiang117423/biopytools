"""疫霉菌基因组注释命令|Oomycete genome annotation command"""

import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...oomycete_anno.main import main as oomycete_main
        return oomycete_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help="疫霉菌基因组注释(T2T Augustus流程)|Oomycete genome annotation (T2T Augustus)",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("-g", "--genome", required=True,
              help="基因组 FASTA|Genome FASTA")
@click.option("-s", "--species", required=True,
              help="Augustus 物种名(简单字母, 如 phytophthora)|Augustus species name")
@click.option("--rnaseq-dirs", multiple=True,
              help="二代 RNA-seq 目录(可多个)|Short RNA-seq dir(s)")
@click.option("--prot-seq", default=None,
              help="同源蛋白文件(Phase2)|Homologous proteins (P2)")
@click.option("--isoseq", default=None,
              help="三代转录本文件(Phase2)|Long-read transcripts (P2)")
@click.option("--read1-pattern", default="_1.clean.fq.gz", show_default=True,
              help="R1 文件后缀模式|R1 suffix pattern")
@click.option("--read2-pattern", default="_2.clean.fq.gz", show_default=True,
              help="R2 文件后缀模式|R2 suffix pattern")
@click.option("--rna-strandness", default="", show_default=True,
              help="链特异性: ''(非链特异性)/FR/RF|Strandness: ''/FR/RF")
@click.option("-o", "--output-dir", default="./oomycete_anno_output", show_default=True,
              help="输出目录|Output directory")
@click.option("-t", "--threads", default=12, show_default=True,
              help="线程数|Threads")
@click.option("--no-soft-masking", "soft_masking", is_flag=True, default=False,
              help="禁用软屏蔽(改用硬屏蔽)|Disable soft masking")
@click.option("--gmes-petap-path", default=None,
              help="GeneMark gmes_petap.pl 路径|GeneMark gmes_petap.pl path")
@click.option("--genemark-perl-env", default=None,
              help="GeneMark perl 提供环境|GeneMark perl provider env")
@click.option("--skip-repeat", is_flag=True, default=False,
              help="跳过重复屏蔽|Skip repeat masking")
@click.option("--skip-rna", is_flag=True, default=False,
              help="跳过 RNA-seq 比对|Skip RNA-seq alignment")
@click.option("--skip-iso", is_flag=True, default=False,
              help="跳过三代转录本(Phase2)|Skip long-read (P2)")
@click.option("--skip-protein", is_flag=True, default=False,
              help="跳过蛋白证据(Phase2)|Skip protein hints (P2)")
@click.option("--skip-ltr", is_flag=True, default=False,
              help="跳过 LTR 注解(Phase2)|Skip LTR annotation (P2)")
def oomycete_anno(genome, species, rnaseq_dirs, prot_seq, isoseq,
                  read1_pattern, read2_pattern, rna_strandness,
                  output_dir, threads, soft_masking,
                  gmes_petap_path, genemark_perl_env,
                  skip_repeat, skip_rna, skip_iso, skip_protein, skip_ltr):
    """
    疫霉菌基因组注释(T2T Augustus流程)|Oomycete genome annotation (T2T Augustus pipeline).

    示例|Examples: biopytools oomycete-anno -g genome.fa -s phytophthora --rnaseq-dirs rna1/ rna2/ -o out/ -t 24
    """
    oomycete_main = _lazy_import_main()

    args = ["oomycete_anno.py",
            "-g", genome,
            "-s", species,
            "--read1-pattern", read1_pattern,
            "--read2-pattern", read2_pattern,
            "--rna-strandness", rna_strandness,
            "-o", output_dir,
            "-t", str(threads)]
    if rnaseq_dirs:
        args.extend(["--rnaseq-dirs"] + list(rnaseq_dirs))
    if prot_seq:
        args.extend(["--prot-seq", prot_seq])
    if isoseq:
        args.extend(["--isoseq", isoseq])
    if soft_masking:
        args.append("--no-soft-masking")
    if gmes_petap_path:
        args.extend(["--gmes-petap-path", gmes_petap_path])
    if genemark_perl_env:
        args.extend(["--genemark-perl-env", genemark_perl_env])
    if skip_repeat:
        args.append("--skip-repeat")
    if skip_rna:
        args.append("--skip-rna")
    if skip_iso:
        args.append("--skip-iso")
    if skip_protein:
        args.append("--skip-protein")
    if skip_ltr:
        args.append("--skip-ltr")

    original_argv = sys.argv
    sys.argv = args
    try:
        oomycete_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
