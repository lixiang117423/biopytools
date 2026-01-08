"""
FASTP数据处理模块|FASTP Data Processing Module
"""

import os
import re
import fnmatch
from pathlib import Path
from typing import List, Tuple


def extract_sample_name(filename: str, pattern: str) -> str:
    """
    从文件名中提取样本名，支持通配符模式|Extract sample name from filename, supports wildcard patterns

    Args:
        filename: 文件名|Filename
        pattern: 文件名模式（可包含通配符）|Filename pattern (may contain wildcards)

    Returns:
        提取的样本名|Extracted sample name
    """
    # 如果模式不包含通配符，直接替换|If pattern doesn't contain wildcard, replace directly
    if '*' not in pattern and '?' not in pattern:
        return filename.replace(pattern, "")

    # 将通配符模式转换为正则表达式|Convert wildcard pattern to regex
    # 将 * 替换为 (.+) 捕获组以提取样本名|Replace * with (.+) capture group to extract sample name
    regex_pattern = fnmatch.translate(pattern)
    # 修改正则表达式，将通配符部分设为捕获组|Modify regex to make wildcard part a capture group
    regex_pattern = regex_pattern.replace(r'\.\*.*', r'(.+)\..*')
    regex_pattern = regex_pattern.replace(r'\Z', r'\Z')

    match = re.match(regex_pattern, filename)
    if match:
        return match.group(1)

    # 如果正则匹配失败，回退到简单替换|If regex match fails, fallback to simple replacement
    # 去除已知的后缀|Remove known suffix
    name = filename
    # 尝试去除模式中通配符后面的部分|Try to remove part after wildcard in pattern
    if '*' in pattern:
        suffix = pattern.split('*', 1)[1]
        if name.endswith(suffix):
            name = name[:-len(suffix)]
    return name


class SampleFinder:
    """样本查找器|Sample Finder"""

    def __init__(self, config, logger):
        """
        初始化样本查找器|Initialize sample finder

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger

    def find_sample_pairs(self) -> List[Tuple[str, Path, Path]]:
        """
        查找样本配对文件|Find sample paired files

        Returns:
            样本配对列表，每个元素为 (样本名, read1文件, read2文件)|List of sample pairs (sample_name, read1_file, read2_file)
            单末端模式下 read2_file 为 None|In single-end mode, read2_file is None
        """
        sample_pairs = []

        # 单末端模式|Single-end mode
        if self.config.single_end:
            # 查找所有fastq文件|Find all fastq files
            # 如果后缀已经包含通配符，直接使用；否则在前面加*|If suffix contains wildcard, use directly; otherwise add *
            if '*' in self.config.read1_suffix:
                file_pattern = self.config.read1_suffix
            else:
                file_pattern = f"*{self.config.read1_suffix}"
            input_files = list(self.config.input_path.glob(file_pattern))

            if not input_files:
                self.logger.warning(
                    f"在输入目录中未找到匹配模式的文件|"
                    f"No files matching pattern found in input directory: {file_pattern}"
                )
                return sample_pairs

            self.logger.info(f"找到 {len(input_files)} 个单末端文件|Found {len(input_files)} single-end files")

            for input_file in input_files:
                # 提取样本名|Extract sample name
                sample_name = extract_sample_name(input_file.name, self.config.read1_suffix)
                sample_pairs.append((sample_name, input_file, None))
                self.logger.debug(f"找到单末端样本|Found single-end sample: {sample_name}")

        else:
            # 双末端模式|Paired-end mode
            # 查找所有read1文件|Find all read1 files
            # 如果后缀已经包含通配符，直接使用；否则在前面加*|If suffix contains wildcard, use directly; otherwise add *
            if '*' in self.config.read1_suffix:
                read1_pattern = self.config.read1_suffix
            else:
                read1_pattern = f"*{self.config.read1_suffix}"
            read1_files = list(self.config.input_path.glob(read1_pattern))

            if not read1_files:
                self.logger.warning(
                    f"在输入目录中未找到匹配模式的文件|"
                    f"No files matching pattern found in input directory: {read1_pattern}"
                )
                return sample_pairs

            self.logger.info(f"找到 {len(read1_files)} 个read1文件|Found {len(read1_files)} read1 files")

            for read1_file in read1_files:
                # 提取样本名|Extract sample name
                sample_name = extract_sample_name(read1_file.name, self.config.read1_suffix)

                # 构建read2文件路径|Construct read2 file path
                read2_file = self.config.input_path / f"{sample_name}{self.config.read2_suffix}"

                if read2_file.exists():
                    sample_pairs.append((sample_name, read1_file, read2_file))
                    self.logger.debug(f"找到配对样本|Found paired sample: {sample_name}")
                else:
                    self.logger.warning(
                        f"找不到样本 {sample_name} 的配对文件|"
                        f"Cannot find paired file for sample {sample_name}: {read2_file}"
                    )

        return sample_pairs

    def validate_sample_pairs(self, sample_pairs: List[Tuple[str, Path, Path]]) -> bool:
        """
        验证样本配对|Validate sample pairs

        Args:
            sample_pairs: 样本配对列表|List of sample pairs

        Returns:
            验证是否通过|Whether validation passed
        """
        if not sample_pairs:
            self.logger.error("未找到有效的样本文件|No valid sample files found")
            return False

        # 检查文件大小|Check file sizes
        for sample_name, read1_file, read2_file in sample_pairs:
            if read1_file.stat().st_size == 0:
                self.logger.warning(f"样本 {sample_name} 的read1文件为空|Read1 file is empty for sample {sample_name}")

            if read2_file is not None and read2_file.stat().st_size == 0:
                self.logger.warning(f"样本 {sample_name} 的read2文件为空|Read2 file is empty for sample {sample_name}")

        return True
