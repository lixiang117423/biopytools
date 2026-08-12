"""
选择性扫荡检测命令行入口|Selective Sweep Detection CLI Entry
"""

import argparse
import sys

from .config import SweepMergeConfig
from .merge_stats import SweepMerger
from .pipeline import SweepPipeline
from .utils import SweepModuleLogger


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='选择性扫荡检测:过滤VCF→群体拆分→pi/Tajima\'s D/Fst/RAiSD μ/SweeD CLR→'
                    'composite_score→候选扫荡区域'
                    '|Selective sweep detection: filter VCF -> split pops -> '
                    'pi/TajimaD/Fst/RAiSD mu/SweeD CLR -> composite_score -> candidate regions')

    parser.add_argument('-i', '--input', required=True,
                        help='输入VCF文件|Input VCF file (支持.gz|.gz supported)')
    parser.add_argument('-p', '--pop-info', required=True,
                        help='群体分组文件(样品ID<TAB>分组,无表头)'
                             '|Population info file (sample<TAB>group, no header)')
    parser.add_argument('-o', '--output-dir', default='./selective_sweep_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数(默认12)|Thread count (default 12)')
    parser.add_argument('--win', type=int, default=50000,
                        help='窗口大小(默认50000)|Window size (default 50000)')
    parser.add_argument('--step', type=int, default=50000,
                        help='窗口步长(默认50000)|Window step (default 50000)')
    parser.add_argument('--top-quantile', type=float, default=0.01,
                        help='候选阈值分位数(默认0.01)|Candidate threshold quantile (default 0.01)')
    parser.add_argument('--merge-gap', type=int, default=100000,
                        help='候选窗口合并最大间隔(默认100000)'
                             '|Max gap for merging candidate windows (default 100000)')
    parser.add_argument('--min-maf', type=float, default=0.05,
                        help='过滤MAF阈值(默认0.05)|Filter MAF threshold (default 0.05)')
    parser.add_argument('--max-missing', type=float, default=0.10,
                        help='过滤缺失率阈值(默认0.10)|Filter missing rate threshold (default 0.10)')
    parser.add_argument('--raisd-window', type=int, default=50,
                        help='RAiSD SNP窗口(默认50)|RAiSD SNP window (default 50)')
    parser.add_argument('--raisd-min-samples', type=int, default=15,
                        help='低样本量阈值,低于此值MU分量默认排除(默认15)'
                             '|Low sample threshold; MU excluded below (default 15)')
    parser.add_argument('--include-mu-low-n', action='store_true',
                        help='低样本群体也加入MU分量|Include MU component for low-n pops')
    parser.add_argument('--sweed-grid', type=int, default=10000,
                        help='SweeD CLR网格点数(默认10000)|SweeD CLR grid points (default 10000)')
    parser.add_argument('--sweed-min-samples', type=int, default=15,
                        help='低样本量阈值,低于此值CLR分量默认排除(默认15)'
                             '|Low sample threshold; CLR excluded below (default 15)')
    parser.add_argument('--include-sweed-low-n', action='store_true',
                        help='低样本群体也加入CLR分量|Include CLR component for low-n pops')
    parser.add_argument('--sweed-unfolded', action='store_true',
                        help='使用未折叠SFS(需祖先状态,默认折叠)|Unfolded SFS '
                             '(requires ancestral state; folded by default)')
    parser.add_argument('--xpclr-maxsnps', type=int, default=200,
                        help='XP-CLR窗口最大SNP数(默认200)|XP-CLR max SNPs per window (default 200)')
    parser.add_argument('--xpclr-minsnps', type=int, default=10,
                        help='XP-CLR窗口最小SNP数(默认10)|XP-CLR min SNPs per window (default 10)')
    parser.add_argument('--xpclr-ld', type=float, default=0.95,
                        help='XP-CLR LD加权截断(默认0.95)|XP-CLR LD cutoff (default 0.95)')
    parser.add_argument('--xpclr-min-samples', type=int, default=15,
                        help='低样本量阈值,任一群体低于此值XP-CLR分量默认排除(默认15)'
                             '|Low sample threshold; XP-CLR excluded if either pop below (default 15)')
    parser.add_argument('--include-xpclr-low-n', action='store_true',
                        help='低样本群体对也加入XP-CLR分量|Include XP-CLR component for low-n pairs')
    parser.add_argument('--log-level', default='INFO',
                        help='日志级别(默认INFO)|Log level (default INFO)')
    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        config = SweepMergeConfig(
            input_vcf=args.input,
            pop_info=args.pop_info,
            output_dir=args.output_dir,
            threads=args.threads,
            win=args.win,
            step=args.step,
            top_quantile=args.top_quantile,
            merge_gap=args.merge_gap,
            min_maf=args.min_maf,
            max_missing=args.max_missing,
            raisd_window=args.raisd_window,
            raisd_min_samples=args.raisd_min_samples,
            include_mu_low_n=args.include_mu_low_n,
            sweed_grid=args.sweed_grid,
            sweed_min_samples=args.sweed_min_samples,
            include_sweed_low_n=args.include_sweed_low_n,
            sweed_folded=not args.sweed_unfolded,
            xpclr_maxsnps=args.xpclr_maxsnps,
            xpclr_minsnps=args.xpclr_minsnps,
            xpclr_ld=args.xpclr_ld,
            xpclr_min_samples=args.xpclr_min_samples,
            include_xpclr_low_n=args.include_xpclr_low_n,
            log_level=args.log_level,
        )
        config.validate()

        logger_manager = SweepModuleLogger(config.log_dir, log_level=config.log_level)
        logger = logger_manager.get_logger()
        logger.info(f"输入VCF|Input VCF: {config.input_vcf}")
        logger.info(f"群体分组|Pop info: {config.pop_info}")

        pipeline = SweepPipeline(config, logger)
        pipeline.run()

        merger = SweepMerger(config, logger)
        merger.run()

        logger.info("全部完成|All done")
        sys.exit(0)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
