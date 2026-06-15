"""
染色体重命名主程序模块|Chromosome Rename Main Module
"""

import argparse
import sys
from .config import ChrRenameConfig
from .utils import ChrRenameLogger, CommandRunner
from .calculator import ChrRenameCalculator


class ChrRenameAnnotator:
    """染色体重命名主类|Main Chromosome Rename Annotator Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = ChrRenameConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = ChrRenameLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化计算器|Initialize calculator
        self.calculator = ChrRenameCalculator(
            self.config,
            self.logger,
            self.cmd_runner
        )

    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            success = self.calculator.run_analysis()
            if not success:
                self.logger.error("分析失败|Analysis failed")
                sys.exit(1)
        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            import traceback
            self.logger.debug(traceback.format_exc())
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='基于minimap2的染色体重命名工具|Chromosome rename tool based on minimap2 alignment',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-r', '--ref', required=True,
                       help='参考基因组FASTA文件|Reference genome FASTA file')
    parser.add_argument('-q', '--query', required=True,
                       help='待重命名的基因组FASTA文件|Query genome FASTA file to rename')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output-dir',
                       default='./chr_rename_output',
                       help='输出目录|Output directory')

    parser.add_argument('-a', '--minimap2-path',
                       default='minimap2',
                       help='minimap2软件路径|minimap2 software path')

    parser.add_argument('-x', '--preset',
                       default='asm5',
                       choices=['asm5', 'asm10', 'asm20'],
                       help='minimap2预设模式|minimap2 preset mode (asm5/asm10/asm20)')

    parser.add_argument('-t', '--threads',
                       type=int, default=12,
                       help='线程数|Number of threads')

    parser.add_argument('-i', '--min-identity',
                       type=float, default=0.9,
                       help='最小序列一致性阈值(0-1)|Minimum identity threshold (0-1)')
    parser.add_argument('-l', '--min-alignment-length',
                       type=int, default=100000,
                       help='最小比对长度(bp)|Minimum alignment length (bp)')

    args = parser.parse_args()

    # 创建重命名器并运行|Create renamer and run
    renamer = ChrRenameAnnotator(
        ref_fasta=args.ref,
        query_fasta=args.query,
        output_dir=args.output_dir,
        minimap2_path=args.minimap2_path,
        preset=args.preset,
        threads=args.threads,
        min_identity=args.min_identity,
        min_alignment_length=args.min_alignment_length
    )

    renamer.run_analysis()


if __name__ == "__main__":
    main()
