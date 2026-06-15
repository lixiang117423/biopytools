"""
重复序列屏蔽主程序模块|Repeat Masking Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import RepeatMaskConfig
from .utils import RepeatMaskLogger, CommandRunner, get_genome_stats
from .core import RepeatMaskPipeline


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="重复序列屏蔽工具|Repeat Masking Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--genome',
                       required=True,
                       help='基因组FASTA文件|Genome FASTA file')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-dir',
                       default='./repeatmask_output',
                       help='输出目录|Output directory (default: ./repeatmask_output)')

    # 分析参数|Analysis parameters
    parser.add_argument('-t', '--threads',
                       type=int,
                       default=12,
                       help='线程数|Number of threads (default: 12)')

    parser.add_argument('-s', '--species',
                       default=None,
                       help='物种名称(用于Dfam/Repbase数据库)|Species name for Dfam/Repbase database')

    # 软件路径|Software paths
    parser.add_argument('--repeatmodeler-path',
                       default='~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatModeler',
                       help='RepeatModeler路径|RepeatModeler path')

    parser.add_argument('--repeatmasker-path',
                       default='~/miniforge3/envs/repeat_identiy/bin/RepeatMasker',
                       help='RepeatMasker路径|RepeatMasker path')

    parser.add_argument('--builddatabase-path',
                       default='BuildDatabase',
                       help='BuildDatabase路径|BuildDatabase path')

    # 流程控制|Pipeline control
    parser.add_argument('--skip-modeler',
                       action='store_true',
                       help='跳过RepeatModeler步骤|Skip RepeatModeler step')

    parser.add_argument('--no-dfam',
                       action='store_true',
                       help='不使用Dfam数据库|Do not use Dfam database')

    parser.add_argument('--masking-mode',
                       choices=['soft', 'hard', 'x'],
                       default='soft',
                       help='屏蔽模式|Masking mode: soft(小写|lowercase), hard(N), x(X) (default: soft)')

    # RepeatModeler参数|RepeatModeler parameters
    parser.add_argument('--use-ltr',
                       action='store_true',
                       default=True,
                       help='运行LTR结构发现|Run LTR structural discovery (default: enabled)')

    parser.add_argument('--no-ltr',
                       action='store_true',
                       help='不运行LTR结构发现|Do not run LTR structural discovery')

    parser.add_argument('--modeler-quick',
                       action='store_true',
                       help='RepeatModeler快速模式|RepeatModeler quick mode')

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建配置|Create configuration
        config = RepeatMaskConfig(
            genome=args.genome,
            output_dir=args.output_dir,
            threads=args.threads,
            species=args.species,
            repeatmodeler_path=args.repeatmodeler_path,
            repeatmasker_path=args.repeatmasker_path,
            builddatabase_path=args.builddatabase_path,
            skip_modeler=args.skip_modeler,
            use_dfam=not args.no_dfam,
            masking_mode=args.masking_mode,
            use_ltr=not args.no_ltr,
            modeler_quick=args.modeler_quick
        )

        # 验证配置|Validate configuration
        config.validate()

        # 初始化日志|Initialize logging
        logger_manager = RepeatMaskLogger(config.output_path)
        logger = logger_manager.get_logger()

        logger.info("开始重复序列屏蔽流程|Starting repeat masking pipeline")
        logger.info(f"输入基因组|Input genome: {config.genome}")
        logger.info(f"输出目录|Output directory: {config.output_dir}")
        logger.info(f"线程数|Threads: {config.threads}")

        if config.species:
            logger.info(f"物种|Species: {config.species}")

        # 获取基因组统计信息|Get genome statistics
        genome_stats = get_genome_stats(config.genome, logger)

        # 初始化命令执行器|Initialize command runner
        cmd_runner = CommandRunner(logger, config.output_path)

        # 创建流程对象|Create pipeline object
        pipeline = RepeatMaskPipeline(config, logger, cmd_runner)

        # 构建从头重复库|Build de novo repeat library
        library_path = pipeline.build_denovo_library()

        # 运行屏蔽|Run masking
        success = pipeline.run_masking(library_path)

        if success:
            logger.info("重复序列屏蔽流程完成|Repeat masking pipeline completed successfully")
            logger.info(f"结果保存在|Results saved in: {config.output_dir}")
            sys.exit(0)
        else:
            logger.error("重复序列屏蔽流程失败|Repeat masking pipeline failed")
            sys.exit(1)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)

    except KeyboardInterrupt:
        print("\n用户中断|User interrupted")
        sys.exit(130)

    except Exception as e:
        print(f"错误|Error: {e}")
        import traceback
        traceback.print_exc()
        sys.exit(1)


if __name__ == '__main__':
    main()
