#!/usr/bin/env python3
"""
VCF筛选模块提取脚本 | VCF Filter Module Extraction Script
从打包文件中提取各个模块到对应目录 | Extract modules from packed file to corresponding directories
"""

import os
import re
from pathlib import Path

def extract_modules():
    """提取模块文件 | Extract module files"""
    
    # 读取当前文件内容 | Read current file content
    current_file = __file__
    with open(current_file, 'r', encoding='utf-8') as f:
        content = f.read()
    
    # 查找所有模块 | Find all modules
    pattern = r'# ===== FILE: ([^=]+) =====\\n(.*?)(?=
