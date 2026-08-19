"""easyhap CLI 包装器|easyhap Click wrapper"""
import os
import sys

import click


def _lazy_import_main():
    """延迟导入主函数|Lazy import main"""
    try:
        from ...easyhap.main import main as easyhap_main
        return easyhap_main
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


@click.command(short_help="区域单倍型分析(EasyHap)|Regional haplotype analysis (EasyHap)")
@click.option("-i", "--input", "input_vcf", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="phased VCF/VCF.gz/BCF 文件|Phased VCF/VCF.gz/BCF file")
@click.option("--group", required=True,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="样本分组表(样本<TAB>分组, 无表头)|Sample group table")
@click.option("--region", default=None,
              help="单区域 chr:start-end|Single region chr:start-end")
@click.option("--region-file", default=None,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="批量区域文件(TAB 三列)|Batch region file")
@click.option("-o", "--output-dir", default="./easyhap_output",
              help="输出目录(默认./easyhap_output)|Output directory")
@click.option("--mode", default="inbred", type=click.Choice(["inbred", "hybrid"]),
              help="单倍型重建策略(默认inbred)|Mode (default inbred)")
@click.option("--hetero-policy", default="slash",
              type=click.Choice(["slash", "iupac", "missing"]),
              help="inbred 模式杂合编码(默认slash)|Hetero policy (default slash)")
@click.option("--cluster-threshold", default=0.15, type=float,
              help="聚类阈值(默认0.15)|Cluster threshold (default 0.15)")
@click.option("--vcf-backend", default="auto",
              type=click.Choice(["auto", "cyvcf2", "pysam", "plain"]),
              help="VCF 读取后端(默认auto)|VCF backend (default auto)")
@click.option("--no-processed", is_flag=True, default=False,
              help="不写处理表|Skip processed tables")
@click.option("--fisher-groups", default=None,
              help="Fisher 过滤两分组名|Fisher filter groups")
@click.option("--fisher-alpha", default=None, type=float,
              help="Fisher 过滤 P 值阈值|Fisher filter p cutoff")
@click.option("--fisher-adjust", default="none", type=click.Choice(["none", "bh"]),
              help="多重检验校正(默认none)|Fisher adjust (default none)")
@click.option("--plot", is_flag=True, default=False,
              help="生成图形输出|Generate plots")
@click.option("--gff", default=None,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="GFF3/GTF 注释|GFF3/GTF annotation")
@click.option("--traits", default=None,
              callback=lambda c, p, v: _validate_path_exists(v),
              help="性状表|Trait table")
@click.option("--trait-cols", default=None,
              help="性状列名(逗号分隔)|Trait columns")
@click.option("--plot-format", default="pdf",
              help="图格式(默认pdf)|Plot formats (default pdf)")
@click.option("--plot-hap-level", default="hap", type=click.Choice(["hap", "cluster"]),
              help="按单倍型/聚类出图(默认hap)|Plot level (default hap)")
@click.option("--plot-min-count", default=1, type=int,
              help="图中最小类别计数(默认1)|Min count in plots (default 1)")
@click.option("--force", is_flag=True, default=False,
              help="强制重跑已完成区域|Force re-run")
@click.option("--log-level", default="INFO",
              help="日志级别(默认INFO)|Log level (default INFO)")
@click.option("--log-file", default=None,
              help="日志文件(默认99_logs/easyhap.log)|Log file")
def easyhap(**kwargs):
    """phased VCF 区域单倍型分析(上游 EasyHap v1.0)
    |Regional haplotype analysis on phased VCF (upstream EasyHap v1.0).

    示例|Examples: biopytools easyhap -i sample.phased.vcf.gz --group groups.tsv --region Chr1:1-10000 -o out/
    """
    easyhap_main = _lazy_import_main()
    argv_bits = []

    def add(flag, value):
        if value is not None:
            # 增强赋值会把 argv_bits 重绑定为局部变量致 UnboundLocalError
            # | augmented assignment rebinds argv_bits as local
            argv_bits.extend([flag, str(value)])

    add("-i", kwargs["input_vcf"])
    add("--group", kwargs["group"])
    add("--region", kwargs["region"])
    add("--region-file", kwargs["region_file"])
    add("-o", kwargs["output_dir"])
    add("--mode", kwargs["mode"])
    add("--hetero-policy", kwargs["hetero_policy"])
    add("--cluster-threshold", kwargs["cluster_threshold"])
    add("--vcf-backend", kwargs["vcf_backend"])
    if kwargs["no_processed"]:
        argv_bits.append("--no-processed")
    add("--fisher-groups", kwargs["fisher_groups"])
    add("--fisher-alpha", kwargs["fisher_alpha"])
    add("--fisher-adjust", kwargs["fisher_adjust"])
    if kwargs["plot"]:
        argv_bits.append("--plot")
    add("--gff", kwargs["gff"])
    add("--traits", kwargs["traits"])
    add("--trait-cols", kwargs["trait_cols"])
    add("--plot-format", kwargs["plot_format"])
    add("--plot-hap-level", kwargs["plot_hap_level"])
    add("--plot-min-count", kwargs["plot_min_count"])
    if kwargs["force"]:
        argv_bits.append("--force")
    add("--log-level", kwargs["log_level"])
    add("--log-file", kwargs["log_file"])
    args = ["easyhap.py"] + argv_bits
    original = sys.argv
    sys.argv = args
    try:
        rc = easyhap_main()
        sys.exit(rc if rc is not None else 0)
    except SystemExit as e:
        sys.exit(e.code if e.code is not None else 0)
    finally:
        sys.argv = original
