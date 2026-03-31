"""
DeepBSA主程序模块|DeepBSA Main Module
"""

import argparse
import sys
from pathlib import Path
from typing import Optional

from .config import DeepBSAConfig
from .utils import DeepBSALogger
from .runner import DeepBSARunner


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='DeepBSA批量分析工具 - 运行多种BSA分析方法|DeepBSA Batch Analysis Tool - Run multiple BSA analysis methods',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例用法|Example usage:
  # 运行所有方法（自动清理VCF注释）|Run all methods (auto clean VCF comments)
  %(prog)s -i filtered_merged_gtx.SNP.vcf

  # 只运行指定方法|Run only specified methods
  %(prog)s -i filtered_merged_gtx.SNP.vcf -m DL,K,ED4

  # 并行运行所有方法（默认）|Run all methods in parallel (default)
  %(prog)s -i filtered_merged_gtx.SNP.vcf

  # 串行运行所有方法|Run all methods serially
  %(prog)s -i filtered_merged_gtx.SNP.vcf --no-parallel

  # 不自动清理注释行|Don't auto clean comment lines
  %(prog)s -i filtered_merged_gtx.SNP.vcf -n

  # 保留清理后的文件|Keep cleaned file
  %(prog)s -i filtered_merged_gtx.SNP.vcf -k

可用方法|Available methods: DL, K, ED4, SNP, SmoothG, SmoothLOD, Ridit
        """
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument('-i', '--input-file', required=True,
                         help='输入文件路径（VCF/CSV格式）|Input file path (VCF/CSV format)')

    # 方法选择|Method selection
    methods = parser.add_argument_group('方法选择|Method selection')
    methods.add_argument('-m', '--methods',
                        help='要运行的方法，逗号分隔（默认：全部）|Methods to run, comma-separated (default: all)')

    # 输出参数|Output parameters
    output = parser.add_argument_group('输出参数|Output parameters')
    output.add_argument('-o', '--output-dir',
                       default='deepbsa_results',
                       help='输出目录（默认：deepbsa_results）|Output directory (default: deepbsa_results)')

    # 运行模式|Running mode
    mode = parser.add_argument_group('运行模式|Running mode')
    mode.add_argument('--no-parallel', dest='parallel', action='store_false',
                     help='串行运行所有方法|Run all methods serially')
    mode.add_argument('-p', '--parallel', dest='parallel', action='store_true',
                     help='并行运行所有方法（默认）|Run all methods in parallel (default)')
    mode.set_defaults(parallel=True)

    # 预处理参数|Preprocessing parameters
    preprocess = parser.add_argument_group('预处理参数|Preprocessing parameters')
    preprocess.add_argument('-n', '--no-auto-clean', action='store_true',
                           help='不自动清理VCF注释行|Don\'t auto clean VCF comment lines')
    preprocess.add_argument('-k', '--keep-clean', action='store_true',
                           help='保留清理后的文件|Keep cleaned file')

    # 工具路径|Tool paths
    tools = parser.add_argument_group('工具路径|Tool paths')
    tools.add_argument('--deepbsa-path',
                      default='~/software/DeepBSA/DeepBSA_linux_v1.4/bin/main.py',
                      help='DeepBSA主程序路径|DeepBSA main script path')
    tools.add_argument('--conda-env',
                      default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/DeepBSA',
                      help='Conda环境路径|Conda environment path')

    # 其他参数|Other parameters
    other = parser.add_argument_group('其他参数|Other parameters')
    other.add_argument('-v', '--verbose', action='store_true',
                      help='详细输出|Verbose output')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 解析方法列表|Parse methods list
        methods_list = None
        if args.methods:
            methods_list = [m.strip() for m in args.methods.split(',')]

        # 创建配置|Create configuration
        config = DeepBSAConfig(
            input_file=args.input_file,
            methods=methods_list,
            output_dir=args.output_dir,
            parallel=args.parallel,
            auto_clean=not args.no_auto_clean,
            keep_clean=args.keep_clean,
            deepbsa_path=args.deepbsa_path,
            conda_env=args.conda_env,
            verbose=args.verbose
        )

        # 验证配置|Validate configuration
        config.validate()

        # 初始化日志|Initialize logging
        log_file = config.output_path / "deepbsa.log"
        logger_manager = DeepBSALogger(log_file)
        logger = logger_manager.get_logger()

        # 创建执行器并运行|Create runner and run
        runner = DeepBSARunner(config, logger)
        results = runner.run()

        # 打印摘要|Print summary
        runner.print_summary(results)

        # 检查整体成功状态|Check overall success status
        success_count = sum(1 for v in results.values() if v)
        total_count = len(results)

        if success_count == total_count:
            logger.info("")
            logger.info("=" * 60)
            logger.info(f"所有方法成功完成|All methods completed successfully! ({success_count}/{total_count})")
            logger.info("=" * 60)
            sys.exit(0)
        else:
            logger.warning("")
            logger.warning("=" * 60)
            logger.warning(f"部分方法失败|Some methods failed ({success_count}/{total_count} 成功|successful)")
            logger.warning("=" * 60)
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
