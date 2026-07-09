"""INDEL分子标记命令|INDEL marker command"""

import sys

import click


def _lazy_import_main():
    """延迟加载主函数|Lazy load main function"""
    try:
        from ...indel_marker.main import main as indel_main
        return indel_main
    except ImportError as e:
        click.echo(f"导入错误|Import Error: {e}", err=True)
        sys.exit(1)


@click.command(
    short_help="抗病/感病INDEL共显性标记开发|R/S INDEL codominant marker development",
    context_settings=dict(help_option_names=["-h", "--help"], max_content_width=120),
)
@click.option("-v", "--vcf", required=True,
              help="多样本合并VCF|Multi-sample VCF")
@click.option("-s", "--samplesheet", required=True,
              help="样品分组TSV(sample_name/group/bam_path)|Samplesheet TSV")
@click.option("-g", "--genome-fasta", required=True,
              help="参考基因组|Reference FASTA")
@click.option("-o", "--output-dir", default="./indel_marker_output", show_default=True,
              help="输出目录|Output directory")
@click.option("-t", "--threads", default=12, show_default=True,
              help="线程数|Threads")
@click.option("--min-indel-size", default=10, show_default=True,
              help="最小INDEL长度|Min INDEL size")
@click.option("--max-indel-size", default=100, show_default=True,
              help="最大INDEL长度|Max INDEL size")
@click.option("--min-quality", default=20.0, show_default=True,
              help="最低QUAL过滤(缺失QUAL保留)|Min QUAL filter (missing QUAL kept)")
@click.option("--max-candidates", default=0, show_default=True,
              help="候选数上限(0=不限)|Candidate cap (0=no limit)")
@click.option("--min-group-consistency", default=0.9, show_default=True,
              help="组内纯合一致比例(1.0=严格)|Within-group consistency")
@click.option("--min-samples-per-group", default=1, show_default=True,
              help="每组最少样品数(默认1)|Min samples per group (default 1)")
@click.option("--min-depth", default=10, show_default=True,
              help="最低覆盖度|Min depth")
@click.option("--deletion-depth-ratio", default=0.3, show_default=True,
              help="deletion骤降阈值|deletion drop threshold")
@click.option("--flank-length", default=300, show_default=True,
              help="侧翼长度|Flank length")
def indel_marker(vcf, samplesheet, genome_fasta, output_dir, threads,
                 min_indel_size, max_indel_size, min_quality, max_candidates,
                 min_group_consistency,
                 min_samples_per_group, min_depth, deletion_depth_ratio,
                 flank_length):
    """
    抗病/感病INDEL共显性标记开发|R/S INDEL codominant marker development.

    示例|Examples: biopytools indel-marker -v v.vcf.gz -s ss.tsv -g ref.fa -o out/
    """
    indel_main = _lazy_import_main()

    args = ["indel_marker.py",
            "-v", vcf, "-s", samplesheet, "-g", genome_fasta,
            "-o", output_dir, "-t", str(threads),
            "--min-indel-size", str(min_indel_size),
            "--max-indel-size", str(max_indel_size),
            "--min-quality", str(min_quality),
            "--max-candidates", str(max_candidates),
            "--min-group-consistency", str(min_group_consistency),
            "--min-samples-per-group", str(min_samples_per_group),
            "--min-depth", str(min_depth),
            "--deletion-depth-ratio", str(deletion_depth_ratio),
            "--flank-length", str(flank_length)]

    original_argv = sys.argv
    sys.argv = args
    try:
        indel_main()
    except SystemExit as e:
        sys.exit(e.code)
    except Exception as e:
        click.echo(f"错误|Error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv
