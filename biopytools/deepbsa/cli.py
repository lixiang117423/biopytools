"""
DeepBSA CLI命令行路由器|DeepBSA CLI Command Router
提供子命令模式：batch, run, merge|Provides subcommands: batch, run, merge
"""

import argparse
import sys
from typing import Optional


def main():
    """DeepBSA CLI主入口|DeepBSA CLI main entry"""
    parser = argparse.ArgumentParser(
        description='DeepBSA BSA分析工具集|DeepBSA BSA Analysis Toolkit',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
推荐工作流|Recommended workflow:
  biopytools deepbsa batch -i input.vcf -o batch_jobs/
  [投递并运行 batch_jobs/run.sh]
  bash batch_jobs/all/merge_results.sh

示例|Examples:
  批量模式|Batch mode:  %(prog)s batch -i input.vcf -o batch_jobs/
  运行模式|Run mode:    %(prog)s run -i input.vcf -o results/
  合并模式|Merge mode:  %(prog)s merge -i results/ -o final/ (高级|Advanced)
        """
    )

    subparsers = parser.add_subparsers(
        dest='command',
        help='子命令|Sub-commands',
        metavar='{batch,run,merge}'
    )

    # batch子命令|batch subcommand
    batch_parser = subparsers.add_parser(
        'batch',
        help='生成批量处理命令|Generate batch processing commands',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    _setup_batch_parser(batch_parser)

    # run子命令|run subcommand
    run_parser = subparsers.add_parser(
        'run',
        help='运行DeepBSA分析方法|Run DeepBSA analysis methods',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    _setup_run_parser(run_parser)

    # merge子命令|merge subcommand
    merge_parser = subparsers.add_parser(
        'merge',
        help='合并分析结果（高级功能）|Merge analysis results (Advanced)',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    _setup_merge_parser(merge_parser)

    # vcf2csv子命令|vcf2csv subcommand
    vcf2csv_parser = subparsers.add_parser(
        'vcf2csv',
        help='VCF转CSV（为DeepBSA准备输入数据）|Convert VCF to CSV for DeepBSA',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    _setup_vcf2csv_parser(vcf2csv_parser)

    # 解析参数|Parse arguments
    args = parser.parse_args()

    # 路由到对应的处理函数|Route to handler
    if args.command == 'batch':
        _handle_batch(args)
    elif args.command == 'run':
        _handle_run(args)
    elif args.command == 'merge':
        _handle_merge(args)
    elif args.command == 'vcf2csv':
        _handle_vcf2csv(args)
    else:
        parser.print_help()
        sys.exit(1)


def _setup_batch_parser(parser: argparse.ArgumentParser):
    """设置batch子命令参数|Setup batch subcommand parser"""
    # 添加描述信息|Add description
    parser.description = """
生成批量处理命令脚本|Generate batch processing command scripts

示例|Examples: biopytools deepbsa batch -i input.vcf -o batch_jobs/
"""

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '-i', '--input-file',
        required=True,
        help='输入文件路径（VCF/CSV格式）|Input file path (VCF/CSV format)'
    )
    required.add_argument(
        '-o', '--output-dir',
        required=True,
        help='批量任务输出目录|Batch jobs output directory'
    )
    required.add_argument(
        '-s', '--script-name',
        default='run_deepbsa_methods.sh',
        help='生成的脚本文件名（默认：run_deepbsa_methods.sh）|Generated script filename (default: run_deepbsa_methods.sh)'
    )

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|Optional arguments')
    optional.add_argument(
        '--threads',
        type=int,
        default=88,
        help='每个方法的线程数（默认：88）|Threads per method (default: 88)'
    )
    optional.add_argument(
        '-m', '--methods',
        help='要运行的方法，逗号分隔（默认：全部）|Methods to run, comma-separated (default: all)'
    )
    optional.add_argument(
        '--smooth-func',
        default='Tri-kernel-smooth',
        choices=['Tri-kernel-smooth', 'LOWESS', 'Moving Average'],
        help='平滑函数（默认：Tri-kernel-smooth）|Smooth function (default: Tri-kernel-smooth)'
    )


def _setup_run_parser(parser: argparse.ArgumentParser):
    # 添加描述信息|Add description
    parser.description = """
运行DeepBSA分析方法|Run DeepBSA analysis methods

示例|Examples: biopytools deepbsa run -i input.vcf -o results/
"""

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '-i', '--input-file',
        required=True,
        help='输入文件路径（VCF/CSV格式）|Input file path (VCF/CSV format)'
    )

    # 输出参数|Output parameters
    output = parser.add_argument_group('输出参数|Output parameters')
    output.add_argument(
        '-o', '--output-dir',
        default='deepbsa_results',
        help='输出目录（默认：deepbsa_results）|Output directory (default: deepbsa_results)'
    )

    # 方法选择|Method selection
    methods = parser.add_argument_group('方法选择|Method selection')
    methods.add_argument(
        '-m', '--methods',
        help='要运行的方法，逗号分隔（默认：全部）|Methods to run, comma-separated (default: all)'
    )

    # 运行模式|Running mode
    mode = parser.add_argument_group('运行模式|Running mode')
    mode.add_argument(
        '--no-parallel',
        dest='parallel',
        action='store_false',
        help='串行运行所有方法|Run all methods serially'
    )
    mode.add_argument(
        '-p', '--parallel',
        dest='parallel',
        action='store_true',
        help='并行运行所有方法（默认）|Run all methods in parallel (default)'
    )
    mode.set_defaults(parallel=True)

    # 预处理参数|Preprocessing parameters
    preprocess = parser.add_argument_group('预处理参数|Preprocessing parameters')
    preprocess.add_argument(
        '-n', '--no-auto-clean',
        action='store_true',
        help='不自动清理VCF注释行|Don\'t auto clean VCF comment lines'
    )
    preprocess.add_argument(
        '-k', '--keep-clean',
        action='store_true',
        help='保留清理后的文件|Keep cleaned file'
    )

    # 工具路径|Tool paths
    tools = parser.add_argument_group('工具路径|Tool paths')
    tools.add_argument(
        '--deepbsa-path',
        default='',
        help='DeepBSA主程序路径（默认：使用内置版本）|DeepBSA main script path (default: builtin version)'
    )
    tools.add_argument(
        '--conda-env',
        default='~/miniforge3/envs/DeepBSA',
        help='Conda环境路径|Conda environment path'
    )

    # 多线程参数|Multithread parameters
    mt = parser.add_argument_group('多线程参数|Multithread parameters')
    mt.add_argument(
        '--threads',
        type=int,
        default=6,
        help='DeepBSA线程数（默认：6）|Number of threads (default: 6)'
    )

    # 平滑参数|Smooth parameters
    smooth = parser.add_argument_group('平滑参数|Smooth parameters')
    smooth.add_argument(
        '--smooth-func',
        default='Tri-kernel-smooth',
        choices=['Tri-kernel-smooth', 'LOWESS', 'Moving Average'],
        help='平滑函数（默认：Tri-kernel-smooth）|Smooth function (default: Tri-kernel-smooth)'
    )

    # 其他参数|Other parameters
    other = parser.add_argument_group('其他参数|Other parameters')
    other.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='详细输出|Verbose output'
    )
    other.add_argument(
        '-f', '--force',
        action='store_true',
        help='强制重新运行所有步骤|Force rerun all steps'
    )
    other.add_argument(
        '--skip-merge',
        action='store_true',
        help='跳过合并步骤|Skip merge step'
    )


def _setup_merge_parser(parser: argparse.ArgumentParser):
    """设置merge子命令参数|Setup merge subcommand parser"""
    parser.description = """
合并DeepBSA运行结果|Merges DeepBSA results into unified format

优先从 npy 文件读取平滑值和原始值（与原始软件完全一致），
如果 npy 不存在则 fallback 到 values.txt 重新计算。

示例|Examples: biopytools deepbsa merge -i ./deepbsa_results -o ./merged_results
"""

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '-i', '--input-dir',
        required=True,
        help='DeepBSA输出目录|DeepBSA output directory'
    )
    required.add_argument(
        '-o', '--output-dir',
        required=True,
        help='合并结果输出目录|Merged results output directory'
    )

    # 其他参数|Other parameters
    other = parser.add_argument_group('可选参数|Optional parameters')
    other.add_argument(
        '-m', '--methods',
        help='要合并的方法，逗号分隔（默认：全部）|Methods to merge, comma-separated (default: all)'
    )
    other.add_argument(
        '--smooth-func',
        default='LOWESS',
        choices=['LOWESS', 'Tri-kernel'],
        help='平滑函数（仅fallback时使用，默认：LOWESS）|Smooth function (fallback only, default: LOWESS)'
    )
    other.add_argument(
        '--smooth-frac',
        type=float,
        default=0.1,
        help='平滑窗口比例（仅fallback时使用，默认：0.1）|Smooth window fraction (fallback only, default: 0.1)'
    )


def _setup_vcf2csv_parser(parser: argparse.ArgumentParser):
    """设置vcf2csv子命令参数|Setup vcf2csv subcommand parser"""
    parser.description = """
将VCF文件转换为DeepBSA可用的CSV格式|Convert VCF file to DeepBSA-compatible CSV format
提取FORMAT字段中的AD（Allele Depth）信息|Extract AD (Allele Depth) from FORMAT field

示例|Examples: biopytools deepbsa vcf2csv -i input.vcf -o ./
"""

    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '-i', '--input-file',
        required=True,
        help='输入VCF文件路径|Input VCF file path'
    )

    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '-o', '--output-file',
        required=True,
        help='输出CSV文件路径|Output CSV file path'
    )


def _handle_batch(args):
    """处理batch命令|Handle batch command"""
    from .batch_generator import BatchGenerator
    from .config import BatchConfig
    from .utils import DeepBSALogger

    try:
        # 解析方法列表|Parse methods list
        methods_list = None
        if args.methods:
            methods_list = [m.strip() for m in args.methods.split(',')]

        # 创建配置|Create configuration
        config = BatchConfig(
            input_file=args.input_file,
            output_dir=args.output_dir,
            methods=methods_list,
            script_name=args.script_name,
            threads=args.threads,
            smooth_func=args.smooth_func
        )
        config.validate()

        # 初始化日志|Initialize logging
        log_file = config.output_path / "batch.log"
        logger_manager = DeepBSALogger(log_file)
        logger = logger_manager.get_logger()

        # 创建生成器并运行|Create generator and run
        generator = BatchGenerator(config, logger)
        success = generator.generate()

        sys.exit(0 if success else 1)

    except KeyboardInterrupt:
        print("\n批量命令生成被用户中断|Batch command generation interrupted by user")
        sys.exit(1)
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行出错|Program execution error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def _handle_run(args):
    """处理run命令|Handle run command"""
    from .config import DeepBSAConfig
    from .utils import DeepBSALogger
    from .runner import DeepBSARunner

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
            verbose=args.verbose,
            force=args.force,
            skip_merge=args.skip_merge,
            threads=args.threads,
            smooth_func=args.smooth_func
        )
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


def _handle_merge(args):
    """处理merge命令|Handle merge command"""
    from pathlib import Path
    from ..merge_deepbsa.main import WindowsDeepbsaMerger
    from ..merge_deepbsa.config import MergeDeepbsaConfig
    from ..merge_deepbsa.utils import MergeDeepbsaLogger

    try:
        # 解析方法列表|Parse methods list
        methods_list = None
        if args.methods:
            methods_list = [m.strip() for m in args.methods.split(',')]

        # 创建配置|Create config
        config = MergeDeepbsaConfig(
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            methods=methods_list,
            smooth_func=args.smooth_func,
            smooth_frac=args.smooth_frac
        )

        # 初始化日志|Initialize logging
        log_file = config.output_path / "merge_deepbsa.log"
        logger_manager = MergeDeepbsaLogger(log_file=str(log_file))
        logger = logger_manager.get_logger()

        # 创建合并器并运行|Create merger and run
        merger = WindowsDeepbsaMerger(config, logger)
        merger.run()

        sys.exit(0)

    except KeyboardInterrupt:
        print("\n合并被用户中断|Merge interrupted by user")
        sys.exit(1)
    except ValueError as e:
        print(f"配置错误|Configuration error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行出错|Program execution error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


def _handle_vcf2csv(args):
    """处理vcf2csv命令|Handle vcf2csv command"""
    from pathlib import Path
    from .vcf2csv import vcf2csv
    from .utils import DeepBSALogger

    try:
        input_path = Path(args.input_file).expanduser().resolve()
        output_path = Path(args.output_file).expanduser().resolve()

        if not input_path.exists():
            print(f"错误|Error: 输入文件不存在|Input file not found: {input_path}")
            sys.exit(1)

        output_path.parent.mkdir(parents=True, exist_ok=True)
        log_file = output_path.parent / "vcf2csv.log"
        logger_manager = DeepBSALogger(log_file)
        logger = logger_manager.get_logger()

        csv_path = vcf2csv(input_path, output_path, logger)

        print(f"\n转换完成|Conversion completed: {csv_path}")
        sys.exit(0)

    except KeyboardInterrupt:
        print("\n转换被用户中断|Conversion interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"程序执行出错|Program execution error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
