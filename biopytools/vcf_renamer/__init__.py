"""
VCF文件样品名称重命名工具包 | VCF Sample Name Renamer Toolkit
功能: 重命名VCF文件中的样品名称，防止软件截断 |
Features: Rename sample names in VCF files to prevent software truncation
作者 | Author: Xiang LI
版本 | Version: 1.0.0
日期 | Date: 2025-12-25

使用示例 | Usage Examples:
    from biopytools.vcf_renamer import VCFRenamer

    # 创建重命名器 | Create renamer
    renamer = VCFRenamer(
        input_vcf="input.vcf.gz",
        output_vcf="output.vcf.gz",
        prefix="S"
    )

    # 执行重命名 | Perform renaming
    renamer.rename()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import VCFRenamer, VCFRenamerConfig

__all__ = ['VCFRenamer', 'VCFRenamerConfig']
