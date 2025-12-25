"""
VCF文件样品名称重命名器主模块 | VCF Sample Name Renamer Main Module
"""

import os
import sys
import time
from .config import VCFRenamerConfig
from .utils import VCFRenamerLogger, VCFRenamerChecker, VCFRenamerProcessor


class VCFRenamer:
    """VCF样品名称重命名器 | VCF Sample Name Renamer"""

    def __init__(self, **kwargs):
        """
        初始化重命名器 | Initialize renamer

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        # 初始化配置 | Initialize configuration
        self.config = VCFRenamerConfig(**kwargs)
        self.config.validate()

        # 获取输出目录 | Get output directory
        output_dir = os.path.dirname(self.config.output_vcf) or "."

        # 初始化日志 | Initialize logging
        self.logger_manager = VCFRenamerLogger(output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化工具类 | Initialize utility classes
        self.checker = VCFRenamerChecker(self.logger)
        self.processor = VCFRenamerProcessor(self.logger)

    def rename(self):
        """
        执行重命名操作 | Perform renaming operation

        Returns:
            是否成功 | Whether successful
        """
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("🧬 开始VCF样品名称重命名 | Starting VCF Sample Name Renaming")
            self.logger.info("=" * 60)
            self.logger.info(f"📋 配置信息 | Configuration:\n{self.config}")

            # 步骤1: 检查bcftools | Step 1: Check bcftools
            if not self.checker.check_bcftools():
                return False

            # 步骤2: 提取样品名称 | Step 2: Extract sample names
            self.logger.info("=" * 60)
            self.logger.info("📋 步骤 1/4: 提取样品名称 | Step 1/4: Extracting sample names")
            self.logger.info("=" * 60)
            samples = self.checker.extract_samples(self.config.input_vcf)
            if not samples:
                self.logger.error("❌ 无法提取样品名称 | Failed to extract sample names")
                return False
            self.config.sample_count = len(samples)
            self.config.old_samples = samples

            # 显示原始样品名称 | Show original sample names
            self.logger.info("📋 原始样品名称 | Original sample names:")
            for i, sample in enumerate(samples, 1):
                self.logger.info(f"  {i}. {sample}")

            # 步骤3: 生成映射表 | Step 3: Generate mapping
            self.logger.info("=" * 60)
            self.logger.info("📝 步骤 2/4: 生成新旧ID映射表 | Step 2/4: Generating sample name mapping")
            self.logger.info("=" * 60)
            if not self.processor.generate_mapping(
                samples,
                self.config.prefix,
                self.config.mapping_file
            ):
                return False

            # 显示新样品名称 | Show new sample names
            self.logger.info("📋 新样品名称 | New sample names:")
            for i in range(1, len(samples) + 1):
                new_name = f"{self.config.prefix}{i}"
                old_name = samples[i - 1]
                self.logger.info(f"  {new_name} <- {old_name}")

            # 步骤4: 重命名样品 | Step 4: Rename samples
            self.logger.info("=" * 60)
            self.logger.info("🔄 步骤 3/4: 使用bcftools重命名样品 | Step 3/4: Renaming samples with bcftools")
            self.logger.info("=" * 60)
            if not self.processor.rename_samples(
                self.config.input_vcf,
                self.config.output_vcf,
                self.config.mapping_file
            ):
                return False

            # 步骤5: 创建索引 | Step 5: Create index
            self.logger.info("=" * 60)
            self.logger.info("📇 步骤 4/4: 创建VCF索引 | Step 4/4: Creating VCF index")
            self.logger.info("=" * 60)
            if not self.processor.index_vcf(self.config.output_vcf):
                return False

            # 完成统计 | Completion statistics
            elapsed_time = time.time() - start_time
            self.logger.info("=" * 60)
            self.logger.info("🎉 重命名完成！| Renaming completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"📊 统计信息 | Statistics:")
            self.logger.info(f"  • 样品数量 | Sample count: {self.config.sample_count}")
            self.logger.info(f"  • 输出文件 | Output file: {self.config.output_vcf}")
            self.logger.info(f"  • 映射文件 | Mapping file: {self.config.mapping_file}")
            self.logger.info(f"  • 总耗时 | Total time: {elapsed_time:.2f} seconds")

            # 如果不保留映射文件，删除它 | If not keeping mapping file, delete it
            if not self.config.keep_mapping:
                self.logger.info(f"🗑️  删除映射文件 | Deleting mapping file: {self.config.mapping_file}")
                os.remove(self.config.mapping_file)

            return True

        except Exception as e:
            self.logger.error(f"❌ 程序执行出错 | Program execution error: {str(e)}", exc_info=True)
            return False


def main():
    """主函数 | Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='🧬 VCF文件样品名称重命名工具 | VCF Sample Name Renamer',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 🚀 基本用法（使用默认前缀S）| Basic usage (default prefix S)
  %(prog)s -i input.vcf.gz -o output.vcf.gz

  # 🏷️ 自定义样品名前缀 | Custom sample name prefix
  %(prog)s -i input.vcf.gz -o output.vcf.gz -p Sample

  # 📋 指定映射文件路径 | Specify mapping file path
  %(prog)s -i input.vcf.gz -o output.vcf.gz -m sample_mapping.txt

  # 🗑️ 不保留映射文件 | Don't keep mapping file
  %(prog)s -i input.vcf.gz -o output.vcf.gz --no-mapping
        '''
    )

    # 必需参数 | Required parameters
    required = parser.add_argument_group('required arguments | 必需参数')
    required.add_argument('-i', '--input', required=True,
                         help='📂 输入VCF文件路径 | Input VCF file path')
    required.add_argument('-o', '--output', required=True,
                         help='📤 输出VCF文件路径 | Output VCF file path')

    # 可选参数 | Optional parameters
    optional = parser.add_argument_group('optional arguments | 可选参数')
    optional.add_argument('-p', '--prefix', default='S',
                         help='🏷️ 新样品名前缀 (默认: S) | New sample name prefix (default: S)')
    optional.add_argument('-m', '--mapping',
                         help='📋 映射文件路径 | Mapping file path')
    optional.add_argument('--no-mapping', action='store_true',
                         help='🗑️ 不保留映射文件 | Don\'t keep mapping file')

    args = parser.parse_args()

    # 创建重命名器并运行 | Create renamer and run
    renamer = VCFRenamer(
        input_vcf=args.input,
        output_vcf=args.output,
        prefix=args.prefix,
        mapping_file=args.mapping,
        keep_mapping=not args.no_mapping
    )

    success = renamer.rename()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()
