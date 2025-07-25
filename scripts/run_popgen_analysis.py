#!/usr/bin/env python3
"""
群体遗传分析运行脚本 | Population Genetics Analysis Runner Script
这是一个简化的入口脚本 | Simple entry script

用法 | Usage:
    python run_popgen_analysis.py -v data.vcf.gz -o results -g groups.txt
"""

import sys
import os

# 添加包路径到系统路径
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from biopytools.popgen_toolkit.main import main

if __name__ == "__main__":
    main()
