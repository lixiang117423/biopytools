"""
VCF2Gene主程序模块|VCF2Gene Main Module
"""

import argparse
import sys
from .config import VCF2GeneConfig
from .utils import VCF2GeneLogger
from .calculator import VariantAnnotator


class VCF2GeneRunner:
    """VCF2Gene运行器|VCF2Gene Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = VCF2GeneConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = VCF2GeneLogger(log_level="INFO")
        self.logger = self.logger_manager.get_logger()

    def run(self):
        """运行注释|Run annotation"""
        try:
            annotator = VariantAnnotator(self.config, self.logger)
            annotator.annotate()
            self.logger.info("程序执行成功|Program executed successfully")
            return True
        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            return False


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='VCF变异基因注释工具|VCF Variant Gene Annotation Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--vcf', required=True,
                       help='输入VCF文件路径|Input VCF file path')
    parser.add_argument('-g', '--gff', required=True,
                       help='输入GFF注释文件路径|Input GFF annotation file path')
    parser.add_argument('-o', '--output', required=True,
                       help='输出结果文件路径|Output result file path')

    # 可选参数|Optional arguments
    parser.add_argument('-t', '--threads',
                       type=int, default=12,
                       help='线程数|Number of threads (reserved for future use)')

    args = parser.parse_args()

    # 创建运行器并运行|Create runner and run
    runner = VCF2GeneRunner(
        vcf_file=args.vcf,
        gff_file=args.gff,
        output_file=args.output,
        threads=args.threads
    )

    success = runner.run()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
