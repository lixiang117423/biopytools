"""
Pixy主程序模块|Pixy Main Module
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from .config import PixyConfig
from .utils import PixyLogger, PixyChecker, parse_population_file, get_population_summary, check_vcf_has_invariant_sites
from .calculator import PixyCalculator


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='Pixy群体遗传学统计工具 - 计算pi、dxy、fst等|Pixy Population Genetics Statistics Tool - Calculate pi, dxy, fst, etc.',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法|Example usage:
  # 使用窗口计算所有统计量（pixy要求必须指定窗口）|Calculate all statistics with windows (pixy requires window_size)
  %(prog)s -i variants.vcf.gz -p populations.txt -o pixy_output -w 100000

  # 使用BED文件定义窗口|Use BED file to define windows
  %(prog)s -i variants.vcf.gz -p populations.txt -o pixy_output -b windows.bed

  # 只计算pi和dxy（使用窗口）|Calculate only pi and dxy (with windows)
  %(prog)s -i variants.vcf.gz -p populations.txt -o pixy_output -w 100000 --stats pi,dxy

  # 只计算特定染色体（必须指定窗口）|Calculate only specific chromosomes (must specify window)
  %(prog)s -i variants.vcf.gz -p populations.txt -o pixy_output -w 100000 -c chr1,chr2,chr3

注意|Note: pixy要求必须指定窗口大小(-w)、BED文件(-b)或位点文件(-s)|pixy requires window_size (-w), bed_file (-b), or sites_file (-s)
        """
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--vcf-file', required=True,
                         help='输入VCF文件（需用bgzip压缩并建立tabix索引）|Input VCF file (must be bgzip-compressed and tabix-indexed)')
    required.add_argument('-p', '--pop-file', required=True,
                         help='群体文件（两列：样本ID 群体名）|Population file (two columns: sample_id population_name)')
    required.add_argument('-o', '--output-dir', required=True,
                         help='输出目录|Output directory')

    # 统计量选择|Statistics selection
    stats = parser.add_argument_group('统计量选择|Statistics selection')
    stats.add_argument('--stats',
                       default='pi,dxy,fst,watterson,tajima',
                       help='要计算的统计量，逗号分隔（默认: pi,dxy,fst,watterson,tajima）|Statistics to calculate, comma-separated (default: pi,dxy,fst,watterson,tajima)')
    stats.add_argument('--calc-pi', action='store_true', help='计算pi（核苷酸多样性）|Calculate pi (nucleotide diversity)')
    stats.add_argument('--calc-dxy', action='store_true', help='计算dxy（群体间核苷酸差异）|Calculate dxy (nucleotide divergence)')
    stats.add_argument('--calc-fst', action='store_true', help='计算fst（遗传分化系数）|Calculate fst (genetic differentiation)')
    stats.add_argument('--calc-watterson-theta', action='store_true', help='计算Watterson\'s theta|Calculate Watterson\'s theta')
    stats.add_argument('--calc-tajima-d', action='store_true', help='计算Tajima\'s D|Calculate Tajima\'s D')

    # 窗口参数|Window parameters
    windows = parser.add_argument_group('窗口参数|Window parameters')
    windows.add_argument('-w', '--window-size', type=int,
                        help='窗口大小bp（不设置则全基因组计算）|Window size in bp (null for genome-wide)')
    windows.add_argument('-b', '--bed-file',
                        help='BED文件定义窗口（自定义大小窗口）|BED file defining windows (custom-sized windows)')
    windows.add_argument('-s', '--sites-file',
                        help='位点文件（只计算特定位点）|Sites file (calculate only specific sites)')

    # 质控参数|Quality control parameters
    qc = parser.add_argument_group('质控参数|Quality control parameters')
    qc.add_argument('--min-samples', type=int, default=0,
                    help='每个群体最小样本数（默认: 0=不限制）|Minimum samples per population (default: 0=no limit)')
    qc.add_argument('--max-missing', type=float, default=1.0,
                    help='最大缺失率（默认: 1.0=不限制）|Maximum missing rate (default: 1.0=no limit)')
    qc.add_argument('--min-maf', type=float, default=0.0,
                    help='最小等位基因频率（默认: 0.0=不限制）|Minor allele frequency (default: 0.0=no limit)')
    qc.add_argument('--zscore-window', type=int,
                    help='Z-score过滤窗口大小（不设置则不过滤）|Z-score filtering window size (null=no filter)')

    # 染色体参数|Chromosome parameters
    chrom = parser.add_argument_group('染色体参数|Chromosome parameters')
    chrom.add_argument('-c', '--chromosomes',
                      help='指定染色体列表，逗号分隔（不设置则全部）|List of chromosomes, comma-separated (null for all)')

    # 工具路径|Tool paths
    tools = parser.add_argument_group('工具路径|Tool paths')
    tools.add_argument('--pixy-path', default='pixy',
                      help='pixy可执行文件路径（默认: pixy）|pixy executable path (default: pixy)')
    tools.add_argument('--conda-env',
                      default='~/miniforge3/envs/pixy_v.2.0.0',
                      help='conda环境路径|conda environment path')

    # 其他参数|Other parameters
    other = parser.add_argument_group('其他参数|Other parameters')
    other.add_argument('-t', '--threads', type=int, default=12,
                      help='线程数（默认: 12）|Number of threads (default: 12)')
    other.add_argument('--bypass-invariant-check', action='store_true',
                      help='强制绕过不变位点检查（默认自动检测VCF并自动绕过）|Force bypass invariant sites check (default: auto-detect VCF and bypass if needed)')
    other.add_argument('--keep-intermediate', action='store_true',
                      help='保留中间文件|Keep intermediate files')
    other.add_argument('-v', '--verbose', action='store_true',
                      help='详细输出|Verbose output')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 解析--stats参数|Parse --stats parameter
        stats_list = args.stats.split(',')
        calc_pi = 'pi' in stats_list or args.calc_pi
        calc_dxy = 'dxy' in stats_list or args.calc_dxy
        calc_fst = 'fst' in stats_list or args.calc_fst
        calc_watterson_theta = 'watterson' in stats_list or 'watterson_theta' in stats_list or args.calc_watterson_theta
        calc_tajima_d = 'tajima' in stats_list or 'tajima_d' in stats_list or args.calc_tajima_d

        # 如果没有指定任何统计量，使用默认值（所有统计量）|If no statistics specified, use defaults (all statistics)
        if not any([calc_pi, calc_dxy, calc_fst, calc_watterson_theta, calc_tajima_d]):
            calc_pi = calc_dxy = calc_fst = calc_watterson_theta = calc_tajima_d = True

        # 解析染色体列表|Parse chromosome list
        chromosome_list = None
        if args.chromosomes:
            chromosome_list = [c.strip() for c in args.chromosomes.split(',')]

        # 创建配置|Create configuration
        config = PixyConfig(
            vcf_file=args.vcf_file,
            pop_file=args.pop_file,
            output_dir=args.output_dir,
            calc_pi=calc_pi,
            calc_dxy=calc_dxy,
            calc_fst=calc_fst,
            calc_watterson_theta=calc_watterson_theta,
            calc_tajima_d=calc_tajima_d,
            window_size=args.window_size,
            bed_file=args.bed_file if args.bed_file else None,
            sites_file=args.sites_file if args.sites_file else None,
            min_samples=args.min_samples,
            max_missing=args.max_missing,
            min_maf=args.min_maf,
            zscore_window=args.zscore_window,
            chromosomes=chromosome_list,
            pixy_path=args.pixy_path,
            conda_env=args.conda_env,
            threads=args.threads,
            keep_intermediate=args.keep_intermediate,
            verbose=args.verbose,
            bypass_invariant_check=args.bypass_invariant_check
        )

        # 验证配置|Validate configuration
        config.validate()

        # 初始化日志|Initialize logging
        log_file = config.output_path / "pixy.log"
        logger_manager = PixyLogger(log_file)
        logger = logger_manager.get_logger()

        logger.info("=" * 80)
        logger.info("Pixy群体遗传学统计工具|Pixy Population Genetics Statistics Tool")
        logger.info("=" * 80)

        # 检查环境|Check environment
        checker = PixyChecker(logger, config.pixy_path, config.conda_env)
        if not checker.check_pixy():
            logger.error("Pixy环境检查失败|Pixy environment check failed")
            sys.exit(1)

        if not checker.check_dependencies():
            logger.error("Pixy依赖检查失败|Pixy dependencies check failed")
            sys.exit(1)

        # 解析群体文件并显示统计信息|Parse population file and show summary
        pop_map = parse_population_file(config.pop_path)
        pop_summary = get_population_summary(pop_map)

        logger.info(f"样本总数|Total samples: {pop_summary['total_samples']}")
        logger.info(f"群体数|Populations: {pop_summary['total_populations']}")
        logger.info("各群体样本数|Sample count per population:")
        for pop, count in sorted(pop_summary['population_counts'].items()):
            logger.info(f"  {pop}: {count}")

        # 自动检测VCF是否包含不变位点|Auto-detect if VCF contains invariant sites
        has_invariant = check_vcf_has_invariant_sites(config.vcf_path, config.conda_env, logger)

        # 如果VCF不包含不变位点且用户未手动设置bypass，则自动启用bypass
        # If VCF doesn't contain invariant sites and user didn't manually set bypass, auto-enable it
        if not has_invariant and not config.bypass_invariant_check:
            logger.info("自动启用--bypass_invariant_check参数|Auto-enabling --bypass_invariant_check parameter")
            config.bypass_invariant_check = True

        # 创建计算器并运行|Create calculator and run
        calculator = PixyCalculator(config, logger)
        results = calculator.run_all_calculations()

        # 输出结果摘要|Output results summary
        logger.info("=" * 80)
        logger.info("计算结果摘要|Calculation Results Summary")
        logger.info("=" * 80)

        for stat, result in results.items():
            if stat == 'overall_success':
                continue

            status = "[成功|Success]" if result['success'] else "[失败|Failed]"
            logger.info(f"{result['statistic']}: {status}")
            if result['output_file']:
                logger.info(f"  输出文件|Output file: {result['output_file']}")

        # 检查整体成功状态|Check overall success status
        if results['overall_success']:
            logger.info("=" * 80)
            logger.info("所有计算成功完成|All calculations completed successfully!")
            logger.info("=" * 80)
            sys.exit(0)
        else:
            logger.error("=" * 80)
            logger.error("部分计算失败|Some calculations failed")
            logger.error("=" * 80)
            sys.exit(1)

    except KeyboardInterrupt:
        print("\n分析被用户中断|Analysis interrupted by user")
        sys.exit(1)
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行出错|Program execution error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
