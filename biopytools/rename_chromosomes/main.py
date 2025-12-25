"""
Rename Chromosomes Main Module
染色体重命名主模块
"""

import os
import sys
import argparse
from .config import ChromosomeRenameConfig
from .utils import ChromosomeRenameLogger, FastaRenamer


def setup_logger(output_dir: str):
    """设置日志 | Setup logging"""
    logger_manager = ChromosomeRenameLogger(output_dir)
    return logger_manager.get_logger()


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='Rename Chromosomes - 染色体重命名工具',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 基本用法
  %(prog)s -i genome.fa -o renamed.fa -n 20

  # 指定输出目录
  %(prog)s -i input.fa -o output/renamed.fa -n 24
        '''
    )

    # 必需参数 | Required parameters
    required = parser.add_argument_group('required arguments')
    required.add_argument('-i', '--input', required=True,
                         help='[FILE] 输入FASTA文件路径 | Input FASTA file path')
    required.add_argument('-o', '--output', required=True,
                         help='[FILE] 输出FASTA文件路径 | Output FASTA file path')
    required.add_argument('-n', '--number', required=True, type=int,
                         help='[INT] 染色体数量 | Number of chromosomes')

    args = parser.parse_args()

    # 获取输出目录 | Get output directory
    output_dir = os.path.dirname(args.output)
    if not output_dir:
        output_dir = '.'

    # 创建输出目录 | Create output directory
    os.makedirs(output_dir, exist_ok=True)

    # 设置日志 | Setup logging
    logger = setup_logger(output_dir)

    try:
        logger.info("=" * 60)
        logger.info("Chromosome Rename Pipeline")
        logger.info("=" * 60)

        # 初始化配置 | Initialize config
        config = ChromosomeRenameConfig(
            input_file=args.input,
            output_file=args.output,
            chromosome_number=args.number
        )

        config.validate()

        # 初始化重命名器 | Initialize renamer
        renamer = FastaRenamer(logger)

        # 执行重命名 | Perform renaming
        if renamer.rename_sequences(
            config.input_file,
            config.output_file,
            config.chromosome_number
        ):
            logger.info("=" * 60)
            logger.info("Pipeline completed successfully")
            logger.info("=" * 60)
            return 0
        else:
            logger.error("Pipeline failed")
            return 1

    except KeyboardInterrupt:
        logger.warning("用户中断操作")
        return 130
    except ValueError as e:
        logger.error(f"配置错误: {str(e)}")
        return 1
    except Exception as e:
        logger.error(f"Pipeline失败: {str(e)}", exc_info=True)
        return 1


class ChromosomeRenamer:
    """染色体重命名类 | Chromosome Renamer Class"""

    def __init__(self, **kwargs):
        """
        初始化重命名器 | Initialize renamer

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        self.config = ChromosomeRenameConfig(**kwargs)
        self.config.validate()

        # 创建输出目录 | Create output directory
        output_dir = os.path.dirname(self.config.output_file)
        if not output_dir:
            output_dir = '.'
        os.makedirs(output_dir, exist_ok=True)

        # 初始化日志 | Initialize logging
        self.logger_manager = ChromosomeRenameLogger(output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化重命名器 | Initialize renamer
        self.renamer = FastaRenamer(self.logger)

    def rename(self):
        """
        执行重命名 | Perform renaming

        Returns:
            是否成功 | Whether successful
        """
        return self.renamer.rename_sequences(
            self.config.input_file,
            self.config.output_file,
            self.config.chromosome_number
        )


if __name__ == "__main__":
    sys.exit(main())
