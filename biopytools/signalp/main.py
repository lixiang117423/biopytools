"""
SignalP 6.0信号肽预测主程序模块|SignalP 6.0 Signal Peptide Prediction Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import SignalPConfig
from .utils import SignalPLogger, SignalPRunner


class SignalPPredictor:
    """SignalP信号肽预测主类|SignalP Signal Peptide Prediction Main Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = SignalPConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = SignalPLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化运行器|Initialize runner
        self.runner = SignalPRunner(self.config, self.logger)

    def run(self):
        """运行SignalP预测|Run SignalP prediction"""
        success = self.runner.run()

        if success:
            self.logger.info("预测流程成功完成|Prediction workflow completed successfully")
            return True
        else:
            self.logger.error("预测流程失败|Prediction workflow failed")
            return False


def parse_arguments():
    """解析命令行参数|Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description="SignalP 6.0信号肽预测工具|SignalP 6.0 Signal Peptide Prediction Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 基本用法|Basic usage
  python -m biopytools.signalp -i proteins.faa -o output_dir

  # 指定生物类型和模式|Specify organism and mode
  python -m biopytools.signalp -i proteins.faa -o output --organism eukarya --mode slow

  # 生成图表|Generate plots
  python -m biopytools.signalp -i proteins.faa -o output --format all
        """
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入FASTA文件(氨基酸序列)|Input FASTA file (amino acid sequences)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='输出目录|Output directory'
    )

    # 可选参数|Optional parameters
    parser.add_argument(
        '--organism',
        default='eukarya',
        choices=['eukarya', 'other', 'euk'],
        help='生物类型|Organism type (default: eukarya)'
    )
    parser.add_argument(
        '--format',
        default='txt',
        choices=['txt', 'png', 'eps', 'all', 'none'],
        help='输出格式|Output format (default: txt)'
    )
    parser.add_argument(
        '--mode',
        default='fast',
        choices=['fast', 'slow', 'slow-sequential'],
        help='预测模式|Prediction mode (default: fast)'
    )
    parser.add_argument(
        '--bsize',
        type=int,
        default=12,
        help='批处理大小|Batch size (default: 12)'
    )
    parser.add_argument(
        '--write-procs',
        type=int,
        default=12,
        help='写入进程数|Number of write processes (default: 12)'
    )
    parser.add_argument(
        '--torch-num-threads',
        type=int,
        default=12,
        help='PyTorch线程数|PyTorch threads (default: 12)'
    )
    parser.add_argument(
        '--signalp-path',
        default='~/miniforge3/envs/signalp6/bin/signalp6',
        help='SignalP程序路径|SignalP program path'
    )
    parser.add_argument(
        '--model-dir',
        default=None,
        help='模型权重目录|Model weights directory'
    )
    parser.add_argument(
        '--skip-resolve',
        action='store_true',
        help='跳过resolve步骤|Skip resolve step'
    )
    parser.add_argument(
        '--cleanup-plots',
        action='store_true',
        default=True,
        help='自动删除plot文件|Auto-delete plot files (default: True)'
    )
    parser.add_argument(
        '--keep-plots',
        action='store_true',
        help='保留plot文件|Keep plot files (overrides --cleanup-plots)'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    # 处理plot文件清理参数|Handle plot cleanup parameter
    cleanup_plots = args.cleanup_plots
    if args.keep_plots:
        cleanup_plots = False

    try:
        # 创建预测器|Create predictor
        predictor = SignalPPredictor(
            fasta_file=args.input,
            output_dir=args.output_dir,
            organism=args.organism,
            format=args.format,
            mode=args.mode,
            bsize=args.bsize,
            write_procs=args.write_procs,
            torch_num_threads=args.torch_num_threads,
            signalp_path=args.signalp_path,
            model_dir=args.model_dir,
            skip_resolve=args.skip_resolve,
            cleanup_plots=cleanup_plots
        )

        # 运行预测|Run prediction
        success = predictor.run()

        if success:
            sys.exit(0)
        else:
            sys.exit(1)

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}")
        sys.exit(1)
    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == '__main__':
    main()
