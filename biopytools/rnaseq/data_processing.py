"""
RNA-seq数据处理模块 | RNA-seq Data Processing Module
"""

import os
import glob
import re
from typing import List, Dict, Tuple

class FastqPatternParser:
    """FASTQ文件模式解析器 | FASTQ File Pattern Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse_fastq_pattern(self, pattern: str) -> Tuple[str, str, str, str]:
        """解析fastq文件模式以提取前缀、后缀和读取标识符 | Parse fastq file pattern to extract prefix, suffix and read identifiers"""
        if "*" not in pattern:
            raise ValueError("文件模式必须包含*作为样本名称占位符 | File pattern must contain * as sample name placeholder")

        # 分割模式字符串 | Split pattern string
        parts = pattern.split("*")
        if len(parts) != 2:
            raise ValueError("文件模式只能包含一个*占位符 | File pattern can only contain one * placeholder")

        prefix = parts[0]  # *之前的部分 | Part before *
        suffix = parts[1]  # *之后的部分 | Part after *

        # 检测读取指示符 (R1, R2 or 1, 2) | Detect read indicators
        read_indicators = ["R1", "R2", "_1", "_2", ".1", ".2"]
        read1_indicator = None
        read2_indicator = None

        for indicator in read_indicators:
            if indicator in suffix:
                if indicator.endswith("1") or indicator == "R1":
                    read1_indicator = indicator
                    if indicator == "R1":
                        read2_indicator = "R2"
                    elif indicator == "_1":
                        read2_indicator = "_2"
                    elif indicator == ".1":
                        read2_indicator = ".2"
                break

        if not read1_indicator:
            raise ValueError(
                f"无法识别读取指示符 (R1/R2, _1/_2, .1/.2) | Cannot identify read indicator (R1/R2, _1/_2, .1/.2) in pattern '{pattern}'"
            )

        return prefix, suffix, read1_indicator, read2_indicator

class SampleParser:
    """样本解析器 | Sample Parser"""
    
    def __init__(self, logger):
        self.logger = logger
        self.pattern_parser = FastqPatternParser(logger)
    
    def parse_input_samples(self, input_path: str, fastq_pattern: str = None) -> List[Dict]:
        """解析输入样本信息 | Parse input sample information"""
        samples = []

        if os.path.isdir(input_path):
            if fastq_pattern:
                # 使用用户指定的文件模式 | Use user-specified file pattern
                samples = self._parse_with_pattern(input_path, fastq_pattern)
            else:
                # 默认模式：自动搜索常见的fastq文件对 | Default mode: automatically search for common fastq file pairs
                samples = self._parse_with_default_patterns(input_path)
        else:
            # 如果是文件，假设它是样本信息文件 | If it's a file, assume it's a sample information file
            samples = self._parse_sample_file(input_path)

        return samples
    
    def _parse_with_pattern(self, input_path: str, fastq_pattern: str) -> List[Dict]:
        """使用指定模式解析样本 | Parse samples with specified pattern"""
        samples = []
        
        try:
            prefix, suffix, read1_indicator, read2_indicator = self.pattern_parser.parse_fastq_pattern(fastq_pattern)

            # 构建搜索模式 | Build search pattern
            search_pattern = os.path.join(input_path, f"{prefix}*{suffix}")
            fastq_files = glob.glob(search_pattern)

            # 过滤R1文件 | Filter R1 files
            r1_files = [f for f in fastq_files if read1_indicator in os.path.basename(f)]

            for fq1 in r1_files:
                # 构建相应的R2文件路径 | Build corresponding R2 file path
                fq2 = fq1.replace(read1_indicator, read2_indicator)

                if os.path.exists(fq2):
                    # 提取样本名称 | Extract sample name
                    basename = os.path.basename(fq1)
                    sample_name = basename.replace(prefix, "").replace(suffix, "")

                    samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})
                else:
                    self.logger.warning(f"警告：找不到相应的R2文件 | Warning: Cannot find corresponding R2 file: {fq2}")
        
        except Exception as e:
            self.logger.error(f"解析文件模式时出错 | Error parsing file pattern: {e}")
        
        return samples
    
    def _parse_with_default_patterns(self, input_path: str) -> List[Dict]:
        """使用默认模式解析样本 | Parse samples with default patterns"""
        samples = []
        
        # 默认模式列表 | Default pattern list
        patterns = [
            ("*_1.fq.gz", "*_2.fq.gz"),
            ("*_R1.fq.gz", "*_R2.fq.gz"),
            ("*.R1.fastq.gz", "*.R2.fastq.gz"),
            ("*_1.fastq.gz", "*_2.fastq.gz"),
            ("*.1.fq.gz", "*.2.fq.gz"),
            ("*_R1.fq", "*_R2.fq"),
            ("*_1.fastq", "*_2.fastq"),
        ]

        for pattern1, pattern2 in patterns:
            search_path1 = os.path.join(input_path, pattern1)
            fastq_files = glob.glob(search_path1)

            for fq1 in fastq_files:
                # 构建相应的R2文件路径 | Build corresponding R2 file path
                fq2 = fq1.replace(pattern1.replace("*", ""), pattern2.replace("*", ""))

                if os.path.exists(fq2):
                    # 提取样本名称 | Extract sample name
                    basename = os.path.basename(fq1)
                    # 从模式中移除固定部分以获得样本名称 | Remove fixed parts from pattern to get sample name
                    sample_name = basename
                    for to_remove in [pattern1.replace("*", ""), pattern2.replace("*", "")]:
                        if to_remove in sample_name:
                            sample_name = sample_name.replace(to_remove, "")
                            break

                    # 确保样本名称不为空 | Ensure sample name is not empty
                    if not sample_name:
                        sample_name = os.path.splitext(os.path.splitext(basename)[0])[0]

                    samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})

            # 如果找到文件，不要尝试其他模式 | If files found, don't try other patterns
            if samples:
                break
        
        return samples
    
    def _parse_sample_file(self, input_path: str) -> List[Dict]:
        """解析样本信息文件 | Parse sample information file"""
        samples = []
        
        # 格式：sample_name\tfastq1_path\tfastq2_path | Format: sample_name\tfastq1_path\tfastq2_path
        try:
            with open(input_path, "r") as f:
                for line_num, line in enumerate(f, 1):
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        samples.append({"name": parts[0], "fastq1": parts[1], "fastq2": parts[2]})
                    elif line.strip():  # 非空行但格式不正确 | Non-empty line but incorrect format
                        self.logger.warning(f"行 {line_num} 格式不正确，跳过 | Line {line_num} has incorrect format, skipping: {line.strip()}")
        except Exception as e:
            self.logger.error(f"读取样本文件时出错 | Error reading sample file: {e}")
        
        return samples
