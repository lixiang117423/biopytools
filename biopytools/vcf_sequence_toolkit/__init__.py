"""
序列提取工具包 | Sequence Extraction Toolkit
功能: 从VCF文件和基因组文件中提取特定区间的序列变异信息 | 
Features: Extract sequence variation information from VCF and genome files for specific regions
作者 | Author: Xiang LI  
版本 | Version: v1.0 - 模块化版本 | Modular version
日期 | Date: 2025-07-21

使用示例 | Usage Examples:
    from biopytools.sequence_toolkit import SequenceExtractor, SequenceConfig
    
    # 创建提取器 | Create extractor
    extractor = SequenceExtractor(
        vcf_file="variants.vcf",
        genome_file="genome.fa",
        chrom="chr1",
        start=1000,
        end=1050,
        output_dir="results"
    )
    
    # 运行提取 | Run extraction
    extractor.run_extraction()
"""

__version__ = "1.0.0"
__author__ = "Xiang LI"

from .main import SequenceExtractor, main
from .config import SequenceConfig

__all__ = ['SequenceExtractor', 'SequenceConfig', 'main']

# ===== FILE: setup.py =====
"""
简化的setup.py - 避免与脚本执行冲突
"""
from setuptools import setup, find_packages

# 只在明确调用setup.py时才执行
if __name__ == "__main__" and "setup.py" in __file__:
    setup(
        name="sequence-toolkit",
        version="1.0.0",
        description="从VCF文件和基因组文件中提取特定区间的序列变异信息",
        author="Xiang LI",
        packages=find_packages(),
        install_requires=[
            "pysam>=0.19.0",
            "pandas>=1.3.0",
        ],
        entry_points={
            'console_scripts': [
                'sequence-extractor=sequence_toolkit.main:main',
            ],
        },
        python_requires=">=3.7",
    )

# ===== FILE: sequence_toolkit/__main__.py =====
"""
序列提取工具的模块入口点 | Module entry point for sequence extraction toolkit
支持 python -m sequence_toolkit 调用 | Supports python -m sequence_toolkit calls
"""

import sys
import os

# 避免setup.py干扰
if __name__ == "__main__":
    # 防止被误认为是setup脚本
    if len(sys.argv) > 1 and sys.argv[1] in ['build', 'install', 'sdist', 'bdist_wheel']:
        print("Error: This is not a setup script.")
        print("Usage: python -m sequence_toolkit -h")
        sys.exit(1)
    
    from .main import main
    sys.exit(main())

# ===== FILE: create_clean_entry.py =====
#!/usr/bin/env python3
"""
创建干净的入口脚本，避免setup.py干扰
"""

import os
import sys
import shutil

def create_clean_entry_script(script_name="sequence_extractor", target_dir="/usr/local/bin"):
    """创建干净的入口脚本"""
    
    script_content = '''#!/usr/bin/env python3
"""
序列提取工具 - 干净的入口脚本
避免setup.py和其他配置文件的干扰
"""

import sys
import os

def main():
    # 清理环境，避免setup.py干扰
    os.environ['PYTHONDONTWRITEBYTECODE'] = '1'
    
    # 检查是否是setup相关命令
    if len(sys.argv) > 1 and sys.argv[1] in ['build', 'install', 'sdist', 'bdist_wheel', '--help-commands']:
        print("Error: This is a sequence extraction tool, not a setup script.")
        print("Usage: {} -h".format(os.path.basename(sys.argv[0])))
        return 1
    
    try:
        # 尝试导入并运行
        from biopytools.sequence_toolkit.main import main as sequence_main
        return sequence_main()
    except ImportError as e:
        print(f"Error: Cannot import sequence_toolkit: {e}")
        print("Please ensure the module is properly installed.")
        print("Alternative: python -m sequence_toolkit.main")
        return 1
    except Exception as e:
        print(f"Error: {e}")
        return 1

if __name__ == "__main__":
    sys.exit(main())
'''
    
    # 创建脚本文件
    if os.access(target_dir, os.W_OK):
        script_path = os.path.join(target_dir, script_name)
    else:
        # 如果没有权限写入系统目录，写入当前目录
        script_path = script_name
    
    try:
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        # 设置执行权限
        os.chmod(script_path, 0o755)
        
        print(f"✓ 已创建干净的入口脚本: {script_path}")
        print(f"使用方法: {script_path} -h")
        
        return True
        
    except Exception as e:
        print(f"✗ 创建入口脚本失败: {e}")
        return False

def create_local_entry():
    """在当前目录创建入口脚本"""
    return create_clean_entry_script("sequence_extractor", ".")

if __name__ == "__main__":
    print("创建干净的序列提取工具入口脚本...")
    
    if len(sys.argv) > 1:
        script_name = sys.argv[1]
    else:
        script_name = "sequence_extractor"
    
    success = create_local_entry()
    
    if success:
        print("\\n现在可以使用以下方式运行:")
        print(f"  ./sequence_extractor -h")
        print(f"  python sequence_extractor -h")
        print(f"  python -m sequence_toolkit -h")
    else:
        print("\\n备用方案:")
        print("  python -m sequence_toolkit.main -h")
        print("  python run_sequence_extractor.py -h")
