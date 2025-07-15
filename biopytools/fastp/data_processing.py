"""
FASTP数据处理模块 | FASTP Data Processing Module
"""

import os
from pathlib import Path
from typing import List, Tuple

class SampleFinder:
    """样本查找器 | Sample Finder"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def find_sample_pairs(self) -> List[Tuple[str, Path, Path]]:
        """查找样本配对文件 | Find sample paired files"""
        sample_pairs = []
        
        # 查找所有read1文件 | Find all read1 files
        read1_pattern = f"*{self.config.read1_suffix}"
        read1_files = list(self.config.input_path.glob(read1_pattern))
        
        if not read1_files:
            self.logger.warning(
                f"在输入目录中未找到匹配模式的文件 | "
                f"No files matching pattern found in input directory: {read1_pattern}"
            )
            return sample_pairs
        
        self.logger.info(f"找到 {len(read1_files)} 个read1文件 | Found {len(read1_files)} read1 files")
        
        for read1_file in read1_files:
            # 提取样本名 | Extract sample name
            sample_name = read1_file.name.replace(self.config.read1_suffix, "")
            
            # 构建read2文件路径 | Construct read2 file path
            read2_file = self.config.input_path / f"{sample_name}{self.config.read2_suffix}"
            
            if read2_file.exists():
                sample_pairs.append((sample_name, read1_file, read2_file))
                self.logger.debug(f"找到配对样本 | Found paired sample: {sample_name}")
            else:
                self.logger.warning(
                    f"找不到样本 {sample_name} 的配对文件 | "
                    f"Cannot find paired file for sample {sample_name}: {read2_file}"
                )
        
        return sample_pairs
    
    def validate_sample_pairs(self, sample_pairs: List[Tuple[str, Path, Path]]) -> bool:
        """验证样本配对 | Validate sample pairs"""
        if not sample_pairs:
            self.logger.error("未找到有效的样本配对文件 | No valid sample pairs found")
            return False
        
        # 检查文件大小 | Check file sizes
        for sample_name, read1_file, read2_file in sample_pairs:
            if read1_file.stat().st_size == 0:
                self.logger.warning(f"样本 {sample_name} 的read1文件为空 | Read1 file is empty for sample {sample_name}")
            
            if read2_file.stat().st_size == 0:
                self.logger.warning(f"样本 {sample_name} 的read2文件为空 | Read2 file is empty for sample {sample_name}")
        
        return True
