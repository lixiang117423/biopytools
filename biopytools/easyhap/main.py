"""easyhap 命令行入口|easyhap CLI entry"""
import argparse
import sys

from .calculator import EasyHapCalculator
from .config import EasyHapConfig


def parse_arguments(argv=None):
    """解析命令行参数|Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="easyhap",
        description="区域单倍型分析(EasyHap 封装)|Regional haplotype analysis (EasyHap wrapper)")
    parser.add_argument("-i", "--input", dest="input_vcf", required=True,
                        help="phased VCF/VCF.gz/BCF 文件|Phased VCF/VCF.gz/BCF file")
    parser.add_argument("--group", required=True,
                        help="样本分组表(样本<TAB>分组, 无表头)|"
                             "Sample group table (sample<TAB>group, no header)")
    region_group = parser.add_mutually_exclusive_group(required=True)
    region_group.add_argument("--region",
                              help="单区域 chr:start-end|Single region chr:start-end")
    region_group.add_argument("--region-file",
                              help="批量区域文件(TAB 三列 chr start end)|Batch region file")
    parser.add_argument("-o", "--output-dir", default="./easyhap_output",
                        help="输出目录(默认./easyhap_output)|"
                             "Output directory (default ./easyhap_output)")
    parser.add_argument("--mode", default="inbred", choices=["inbred", "hybrid"],
                        help="单倍型重建策略(默认inbred)|"
                             "Haplotype reconstruction mode (default inbred)")
    parser.add_argument("--hetero-policy", default="slash",
                        choices=["slash", "iupac", "missing"],
                        help="inbred 模式杂合编码(默认slash)|"
                             "Hetero encoding in inbred mode (default slash)")
    parser.add_argument("--cluster-threshold", type=float, default=0.15,
                        help="单倍型聚类阈值(默认0.15)|Clustering threshold (default 0.15)")
    parser.add_argument("--vcf-backend", default="auto",
                        choices=["auto", "cyvcf2", "pysam", "plain"],
                        help="VCF 读取后端(默认auto)|VCF reader backend (default auto)")
    parser.add_argument("--no-processed", action="store_true",
                        help="不写 ProcessedVariants/SampleGenotypeTokens 表|"
                             "Skip processed tables")
    parser.add_argument("--fisher-groups",
                        help="Fisher 过滤两分组名(逗号分隔)|"
                             "Two group names for Fisher filter")
    parser.add_argument("--fisher-alpha", type=float,
                        help="Fisher 过滤 P 值阈值(0,1]|Fisher filter p cutoff (0,1]")
    parser.add_argument("--fisher-adjust", default="none", choices=["none", "bh"],
                        help="多重检验校正(默认none)|"
                             "Multiple-testing correction (default none)")
    parser.add_argument("--plot", action="store_true",
                        help="生成图形输出|Generate plots")
    parser.add_argument("--gff",
                        help="GFF3/GTF 注释(基因结构图)|GFF3/GTF annotation for gene plot")
    parser.add_argument("--traits", help="性状表(有表头)|Trait table (with header)")
    parser.add_argument("--trait-cols",
                        help="性状列名(逗号分隔)|Comma-separated trait columns")
    parser.add_argument("--plot-format", default="pdf",
                        help="图格式(逗号分隔 pdf/svg/png, 默认pdf)|"
                             "Plot formats (default pdf)")
    parser.add_argument("--plot-hap-level", default="hap", choices=["hap", "cluster"],
                        help="按单倍型/聚类出图(默认hap)|Plot by hap or cluster (default hap)")
    parser.add_argument("--plot-min-count", type=int, default=1,
                        help="图中最小类别计数(默认1)|Min class count in plots (default 1)")
    parser.add_argument("--force", action="store_true",
                        help="强制重跑已完成区域|Force re-run completed regions")
    parser.add_argument("--log-level", default="INFO",
                        help="日志级别(默认INFO)|Log level (default INFO)")
    parser.add_argument("--log-file",
                        help="日志文件(默认99_logs/easyhap.log)|"
                             "Log file (default 99_logs/easyhap.log)")
    return parser.parse_args(argv)


def main(argv=None):
    """主函数|Main function. Returns exit code."""
    args = parse_arguments(argv)
    config_kwargs = dict(
        input_vcf=args.input_vcf,
        group_file=args.group,
        region=args.region,
        region_file=args.region_file,
        output_dir=args.output_dir,
        mode=args.mode,
        hetero_policy=args.hetero_policy,
        cluster_threshold=args.cluster_threshold,
        vcf_backend=args.vcf_backend,
        no_processed=args.no_processed,
        fisher_groups=args.fisher_groups,
        fisher_alpha=args.fisher_alpha,
        fisher_adjust=args.fisher_adjust,
        plot=args.plot,
        gff=args.gff,
        traits=args.traits,
        trait_cols=args.trait_cols,
        plot_format=args.plot_format,
        plot_hap_level=args.plot_hap_level,
        plot_min_count=args.plot_min_count,
        force=args.force,
        log_level=args.log_level,
        log_file=args.log_file,
    )
    try:
        config = EasyHapConfig(**config_kwargs)
        config.validate()
    except ValueError as e:
        print(f"参数错误|Invalid arguments:\n{e}", file=sys.stderr)
        return 2
    return EasyHapCalculator(config).run()


if __name__ == "__main__":
    raise SystemExit(main())
