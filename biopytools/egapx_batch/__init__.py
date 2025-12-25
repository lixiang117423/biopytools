"""
EGAPx批量运行工具包 | EGAPx Batch Processing Toolkit
功能: 按染色体批量生成EGAPx运行配置和脚本 |
Features: Generate EGAPx run configs and scripts by chromosome
作者 | Author: Xiang LI
版本 | Version: 1.0.0
日期 | Date: 2025-12-25

使用示例 | Usage Examples:
    from biopytools.egapx_batch import EGAPxBatchGenerator

    # 创建批量生成器 | Create batch generator
    generator = EGAPxBatchGenerator(
        genome="genome.fa",
        yaml_template="template.yaml",
        script_template="template.sh",
        output_dir="output"
    )

    # 执行批量生成 | Perform batch generation
    generator.generate()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import EGAPxBatchGenerator, EGAPxBatchConfig

__all__ = ['EGAPxBatchGenerator', 'EGAPxBatchConfig']
