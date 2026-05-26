"""
DeepLoc 2.1蛋白质亚细胞定位预测主程序模块|DeepLoc 2.1 Protein Subcellular Localization Prediction Main Module
"""

import argparse
import sys
from pathlib import Path

from .config import DeepLocConfig
from .utils import DeepLocLogger, DeepLocRunner


class DeepLocPredictor:
    """DeepLoc亚细胞定位预测主类|DeepLoc Subcellular Localization Prediction Main Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = DeepLocConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = DeepLocLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化运行器|Initialize runner
        self.runner = DeepLocRunner(self.config, self.logger)

    def run(self):
        """运行DeepLoc预测|Run DeepLoc prediction"""
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
        description="DeepLoc 2.1蛋白质亚细胞定位预测工具|DeepLoc 2.1 Protein Subcellular Localization Prediction Tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  # 基本用法|Basic usage
  python -m biopytools.deeploc -i proteins.faa -o output_dir

  # 使用Accurate模型|Use Accurate model
  python -m biopytools.deeploc -i proteins.faa -o output --model Accurate

  # 绘制attention图|Plot attention
  python -m biopytools.deeploc -i proteins.faa -o output --plot

  # 使用GPU加速|Use GPU
  python -m biopytools.deeploc -i proteins.faa -o output --device cuda
        """
    )

    # 必需参数|Required parameters
    parser.add_argument(
        '-i', '--input',
        required=True,
        help='输入FASTA文件(蛋白质序列)|Input FASTA file (protein sequences)'
    )
    parser.add_argument(
        '-o', '--output-dir',
        required=True,
        help='输出目录|Output directory'
    )

    # 可选参数|Optional parameters
    parser.add_argument(
        '--model',
        default='Fast',
        choices=['Fast', 'Accurate'],
        help='预测模型|Prediction model (default: Fast)\n'
             'Fast: ESM1b模型(快速)|Fast: ESM1b model (fast)\n'
             'Accurate: ProtT5模型(高精度)|Accurate: ProtT5 model (accurate)'
    )
    parser.add_argument(
        '--device',
        default='cpu',
        choices=['cpu', 'cuda', 'mps'],
        help='计算设备|Compute device (default: cpu)'
    )
    parser.add_argument(
        '--singularity-image',
        default='~/software/singularity/deeploc2.1_latest.sif',
        help='Singularity镜像路径|Singularity image path'
    )
    parser.add_argument(
        '--database-dir',
        default='~/software/deeploc/database',
        help='数据库目录|Database directory'
    )
    parser.add_argument(
        '--singularity-exec',
        default='~/miniforge3/envs/singularity_v.3.8.7/bin/singularity',
        help='Singularity可执行文件路径|Singularity executable path'
    )
    parser.add_argument(
        '--plot',
        action='store_true',
        help='绘制attention图|Plot attention values'
    )

    return parser.parse_args()


def main():
    """主函数|Main function"""
    args = parse_arguments()

    try:
        # 创建预测器|Create predictor
        predictor = DeepLocPredictor(
            fasta_file=args.input,
            output_dir=args.output_dir,
            model=args.model,
            device=args.device,
            singularity_image=args.singularity_image,
            database_dir=args.database_dir,
            singularity_exec=args.singularity_exec,
            plot=args.plot
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
