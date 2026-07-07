"""
rMVP命令行接口|rMVP Command Line Interface
"""

import sys
import argparse
from pathlib import Path

from .config import RMVPConfig
from .analyzer import RMVPAnalyzer
from .result_parser import RMVPResultParser


def parse_arguments():
    """
    解析命令行参数|Parse command line arguments

    Returns:
        参数对象|Arguments object
    """
    parser = argparse.ArgumentParser(
        description='rMVP GWAS批量分析工具|rMVP Batch GWAS Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  # 基本用法（使用所有3个模型）|Basic usage (all 3 models)
  %(prog)s -i input.vcf.gz -p phenotype.txt -o output

  # 指定模型|Specify models
  %(prog)s -i input.vcf.gz -p phenotype.txt -o output --models GLM MLM

  # 使用conda环境|Use conda environment
  %(prog)s -i input.vcf.gz -p phenotype.txt -o output --r-env ~/miniforge3/envs/rMVP

  # 并行计算|Parallel computing
  %(prog)s -i input.vcf.gz -p phenotype.txt -o output --ncpus 16 --maxLine 50000
        '''
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--vcf',
                        required=True,
                        help='输入VCF文件|Input VCF file')
    parser.add_argument('-p', '--pheno',
                        required=True,
                        help='输入表型文件|Input phenotype file')
    parser.add_argument('-o', '--output',
                        required=True,
                        help='输出目录|Output directory')

    # 可选参数|Optional arguments
    parser.add_argument('--output-prefix',
                        default='RMVP_Result',
                        help='输出前缀|Output prefix (default: RMVP_Result)')

    # 模型选择|Model selection
    parser.add_argument('--models',
                        nargs='+',
                        choices=['GLM', 'MLM', 'FarmCPU'],
                        default=['GLM', 'MLM', 'FarmCPU'],
                        help='分析模型|Analysis models (default: GLM MLM FarmCPU)')

    # R环境|R environment
    parser.add_argument('--r-env',
                        default='~/miniforge3/envs/rMVP',
                        help='R conda环境路径|R conda environment path')
    parser.add_argument('--r-path',
                        help='R可执行文件路径|R executable path')

    # 并行计算|Parallel computing
    parser.add_argument('--ncpus',
                        type=int,
                        default=12,
                        help='CPU核心数|Number of CPU cores (default: 12)')
    parser.add_argument('--maxLine',
                        type=int,
                        default=10000,
                        help='每次读取的SNP数量|Number of SNPs to read at once (default: 10000, smaller uses less memory)')

    # PCA参数|PCA parameters
    parser.add_argument('--n-pc-glm',
                        type=int,
                        default=3,
                        help='GLM模型使用的PC数量|Number of PCs for GLM (default: 3)')
    parser.add_argument('--n-pc-mlm',
                        type=int,
                        default=3,
                        help='MLM模型使用的PC数量|Number of PCs for MLM (default: 3)')
    parser.add_argument('--n-pc-farmcpu',
                        type=int,
                        default=3,
                        help='FarmCPU模型使用的PC数量|Number of PCs for FarmCPU (default: 3)')

    # MLM参数|MLM parameters
    parser.add_argument('--vc-method',
                        choices=['BRENT', 'EMMA', 'HE'],
                        default='BRENT',
                        help='MLM方差组分分析方法|MLM variance component method (default: BRENT)')

    # FarmCPU参数|FarmCPU parameters
    parser.add_argument('--max-loop',
                        type=int,
                        default=10,
                        help='FarmCPU最大迭代次数|FarmCPU max iterations (default: 10)')
    parser.add_argument('--method-bin',
                        choices=['static', 'fast-lmm'],
                        default='static',
                        help='FarmCPU bin方法|FarmCPU bin method (default: static)')

    # 筛选参数|Filtering parameters
    parser.add_argument('--maf',
                        type=float,
                        help='最小等位基因频率阈值|Minor allele frequency threshold')
    parser.add_argument('--miss',
                        type=float,
                        help='缺失率阈值|Missing rate threshold')

    # 输出控制|Output control
    parser.add_argument('--file-type',
                        choices=['jpg', 'pdf', 'tiff'],
                        default='jpg',
                        help='图片格式|Figure format (default: jpg)')
    parser.add_argument('--dpi',
                        type=int,
                        default=300,
                        help='图片分辨率|Figure DPI (default: 300)')
    parser.add_argument('--threshold',
                        type=float,
                        default=0.05,
                        help='显著性阈值|Significance threshold (default: 0.05)')

    # LD去连锁参数|LD pruning parameters（kinship/PCA在去连锁SNP上计算，GWAS用全部SNP）
    parser.add_argument('--ld-pruning',
                        dest='ld_pruning', action='store_true', default=True,
                        help='开启LD去连锁（默认）|Enable LD pruning (default): K/PCA on pruned SNPs')
    parser.add_argument('--no-ld-pruning',
                        dest='ld_pruning', action='store_false',
                        help='关闭LD去连锁，K/PCA用全部SNP|Disable LD pruning, K/PCA use all SNPs')
    parser.add_argument('--ld-window',
                        default='3000kb',
                        help='LD修剪窗口|LD pruning window (e.g. 3000kb 或|or 500, default: 3000kb)')
    parser.add_argument('--ld-step',
                        type=int, default=1,
                        help='LD修剪步长|LD pruning step size (default: 1)')
    parser.add_argument('--ld-r2',
                        type=float, default=0.2,
                        help='LD r2阈值|LD r2 threshold (default: 0.2)')
    parser.add_argument('--plink-path',
                        default=None,
                        help='PLINK可执行文件路径|PLINK executable path (default: conda env Population_genetics)')

    # 日志选项|Logging options
    parser.add_argument('--log-level',
                        choices=['DEBUG', 'INFO', 'WARN', 'ERROR'],
                        default='INFO',
                        help='日志级别|Log level (default: INFO)')
    parser.add_argument('--quiet',
                        action='store_true',
                        help='静默模式，减少输出|Quiet mode, reduce output')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    try:
        # 解析参数|Parse arguments
        args = parse_arguments()

        # 创建输出目录|Create output directory
        output_dir = Path(args.output).expanduser()
        output_dir.mkdir(parents=True, exist_ok=True)

        # 设置日志|Setup logging
        from .utils import RMVPLogger
        log_file = output_dir / f"{args.output_prefix}.log"
        logger_manager = RMVPLogger(log_file, args.log_level)
        logger = logger_manager.get_logger()

        # 验证文件存在性|Validate file existence
        vcf_file = Path(args.vcf).expanduser()
        pheno_file = Path(args.pheno).expanduser()

        if not vcf_file.exists():
            logger.error(f"VCF文件不存在|VCF file does not exist: {vcf_file}")
            sys.exit(1)

        if not pheno_file.exists():
            logger.error(f"表型文件不存在|Phenotype file does not exist: {pheno_file}")
            sys.exit(1)

        # 创建配置|Create config
        config = RMVPConfig(
            vcf_file=str(vcf_file),
            pheno_file=str(pheno_file),
            output_prefix=args.output_prefix,
            output_dir=str(output_dir),
            models=args.models,
            r_env=args.r_env,
            r_path=args.r_path,
            ncpus=args.ncpus,
            maxLine=args.maxLine,
            n_pc_glm=args.n_pc_glm,
            n_pc_mlm=args.n_pc_mlm,
            n_pc_farmcpu=args.n_pc_farmcpu,
            vc_method=args.vc_method,
            max_loop=args.max_loop,
            method_bin=args.method_bin,
            maf=args.maf,
            miss=args.miss,
            file_type=args.file_type,
            dpi=args.dpi,
            threshold=args.threshold,
            ld_pruning=args.ld_pruning,
            ld_window=args.ld_window,
            ld_step=args.ld_step,
            ld_r2=args.ld_r2,
            log_level=args.log_level,
            verbose=not args.quiet,
            **({'plink_path': args.plink_path} if args.plink_path else {})
        )

        # 验证配置|Validate config
        try:
            config.validate()
        except ValueError as e:
            logger.error(f"配置错误|Configuration error: {e}")
            sys.exit(1)

        # 输出分析信息|Output analysis information
        logger.info(f"开始rMVP GWAS分析|Starting rMVP GWAS analysis")
        logger.info(f"VCF文件|VCF file: {vcf_file}")
        logger.info(f"表型文件|Phenotype file: {pheno_file}")
        logger.info(f"输出目录|Output directory: {output_dir}")
        logger.info(f"分析模型|Models: {', '.join(args.models)}")
        logger.info(f"CPU核心数|CPU cores: {args.ncpus}")

        # 创建分析器|Create analyzer
        analyzer = RMVPAnalyzer(config)

        # 运行分析|Run analysis
        success = analyzer.run_analysis()

        if success:
            # 读取表型名称|Read trait names
            trait_names = analyzer._get_trait_names()

            # 解析和整合结果|Parse and integrate results
            logger.info("\n解析和整合结果|Parsing and integrating results")
            result_parser = RMVPResultParser(output_dir, args.output_prefix, logger)

            # 整合结果|Integrate results
            integrated_results = result_parser.parse_and_integrate_results(
                trait_names,
                args.models
            )

            # 保存整合结果|Save integrated results
            if integrated_results:
                result_parser.save_integrated_results(integrated_results)

            # 多表型时按显著性合并结果（>1个表型自动启动）|Merge by significance when multiple traits (>1) detected
            if len(trait_names) > 1:
                logger.info(f"\n检测到{len(trait_names)}个表型，启动结果合并|Detected {len(trait_names)} traits, starting result merge")
                result_parser.merge_results_by_significance(trait_names, args.models)
            else:
                logger.info(f"\n仅{len(trait_names)}个表型，跳过结果合并|Only {len(trait_names)} trait, skipping result merge")

            # 导出 ldblockshow -InGWAS 专用 3 列 TSV(每个 trait×model 一份)
            # |Export ldblockshow-ready 3-col TSV (one per trait×model)
            result_parser.export_ldblockshow_tsv(trait_names, args.models)

            # 收集输出文件|Collect output files
            files = result_parser.collect_output_files(trait_names, args.models)

            # 生成汇总报告|Generate summary report
            result_parser.generate_summary_report(trait_names, args.models)

            # 输出摘要|Output summary
            logger.info("\n分析完成|Analysis completed")
            logger.info(f"  表型数量|Number of traits: {len(trait_names)}")
            logger.info(f"  分析模型|Models: {', '.join(args.models)}")
            logger.info(f"  总分析数|Total analyses: {len(trait_names)} × {len(args.models)} = {len(trait_names) * len(args.models)}")
            logger.info(f"  表格文件|Table files: {len(files['tables'])}")
            logger.info(f"  图片文件|Figure files: {len(files['figures'])}")

            sys.exit(0)
        else:
            logger.error("分析失败|Analysis failed")
            sys.exit(1)

    except KeyboardInterrupt:
        sys.stderr.write("用户中断操作|User interrupted operation\n")
        sys.exit(1)
    except Exception as e:
        sys.stderr.write(f"程序错误|Program error: {e}\n")
        sys.exit(1)


if __name__ == '__main__':
    main()
