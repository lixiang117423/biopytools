#!/usr/bin/env python3
"""
VCF单体型提取模块提取脚本 v2.0 | VCF Haplotype Extractor Module Extraction Script v2.0
从打包文件中提取各个模块到对应的目录结构中 | Extract modules from package file to corresponding directory structure
"""

import os
import re
from pathlib import Path

def extract_modules_from_file(input_file: str, base_dir: str = "."):
    """
    从包含多个模块的文件中提取各个模块
    Extract individual modules from a file containing multiple modules
    """
    
    # 读取输入文件 | Read input file
    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            content = f.read()
    except FileNotFoundError:
        print(f"错误：找不到文件 {input_file}")
        return False
    except Exception as e:
        print(f"错误：读取文件失败 {e}")
        return False
    
    # 使用正则表达式匹配文件块 | Use regex to match file blocks
    file_pattern = r'# ===== FILE: (.+?) =====\n(.*?)(?=# ===== (?:FILE:|END FILE =====))'
    matches = re.findall(file_pattern, content, re.DOTALL)
    
    if not matches:
        print("错误：未找到有效的文件块")
        return False
    
    base_path = Path(base_dir)
    extracted_files = []
    
    for file_path, file_content in matches:
        # 创建完整的文件路径 | Create full file path
        full_path = base_path / file_path
        
        # 创建目录（如果不存在）| Create directory if it doesn't exist
        full_path.parent.mkdir(parents=True, exist_ok=True)
        
        # 写入文件 | Write file
        try:
            with open(full_path, 'w', encoding='utf-8') as f:
                f.write(file_content.lstrip('\n'))
            
            extracted_files.append(str(full_path))
            print(f"✓ 提取: {full_path}")
            
        except Exception as e:
            print(f"✗ 提取失败 {full_path}: {e}")
            return False
    
    print(f"\n成功提取 {len(extracted_files)} 个文件:")
    for file_path in extracted_files:
        print(f"  - {file_path}")
    
    # 创建空的 __pycache__ 忽略文件 | Create empty __pycache__ ignore file
    pycache_gitignore = base_path / "haplotype_extractor" / "__pycache__" / ".gitignore"
    pycache_gitignore.parent.mkdir(exist_ok=True)
    pycache_gitignore.write_text("*\n!.gitignore\n")
    
    print(f"\n模块提取完成！ | Module extraction completed!")
    print(f"使用方法 | Usage:")
    print(f"  python run_haplotype_extractor.py -v input.vcf -p positions.txt -o output.txt")
    print(f"  或者直接导入模块 | Or import module directly:")
    print(f"  from biopytools.haplotype_extractor import HaplotypeExtractor")
    
    print(f"\n依赖要求 | Dependencies:")
    print(f"  - bcftools (必须安装并在PATH中 | Must be installed and in PATH)")
    print(f"  - Python 3.7+ 标准库 | Python 3.7+ standard library")
    
    print(f"\n位置文件格式示例 | Position file format example:")
    print(f"  # 可以有表头（推荐）| Can have header (recommended)")
    print(f"  CHR\\tPOS")
    print(f"  OV12\\t93635286")
    print(f"  OV12\\t93440535")
    print(f"  ")
    print(f"  # 也可以没有表头 | Can also without header")
    print(f"  OV12\\t93635286")
    print(f"  OV12\\t93440535")
    
    return True

def main():
    """主函数 | Main function"""
    import argparse
    
    parser = argparse.ArgumentParser(description="提取VCF单体型提取模块 v2.0")
    parser.add_argument(
        "-f", "--file", 
        default="haplotype_extractor_package.py",
        help="包含模块的输入文件（默认: haplotype_extractor_package.py）"
    )
    parser.add_argument(
        "-d", "--dir",
        default=".",
        help="输出目录（默认: 当前目录）"
    )
    
    args = parser.parse_args()
    
    success = extract_modules_from_file(args.file, args.dir)
    
    if not success:
        exit(1)

if __name__ == "__main__":
    main()
