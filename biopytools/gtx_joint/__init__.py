"""
GTX Joint Calling命令生成工具包 | GTX Joint Calling Command Generator Toolkit
功能: 按染色体或区间生成GTX joint calling命令脚本 |
Features: Generate GTX joint calling command scripts by chromosome or windows
作者 | Author: Xiang LI
版本 | Version: 1.0.0
日期 | Date: 2025-12-25

使用示例 | Usage Examples:
    from biopytools.gtx_joint import GTXJointGenerator

    # 创建命令生成器 | Create command generator
    generator = GTXJointGenerator(
        gtx_exec="/path/to/gtx",
        reference="genome.fa",
        gvcf_dir="./gvcf",
        output_dir="./output"
    )

    # 生成命令脚本 | Generate command script
    generator.generate_commands()

    # 按区间生成 | Generate by windows
    generator.generate_commands(window_size=10000000)
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import GTXJointGenerator, GTXJointConfig

__all__ = ['GTXJointGenerator', 'GTXJointConfig']
