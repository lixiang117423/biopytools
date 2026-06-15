"""
KMC样本查找模块|KMC Sample Finder Module
基于FASTP模块改写|Adapted from FASTP module
"""

import os
import re
import fnmatch
from pathlib import Path
from typing import List, Tuple, Optional


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


class KMCSampleFinder:
    """KMC样本查找器|KMC Sample Finder"""

    def __init__(self, config, logger):
        """
        初始化样本查找器|Initialize sample finder

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger

    def find_samples(self) -> List[Tuple[str, List[Path]]]:
        """
        查找样本文件|Find sample files

        Returns:
            样本列表，每个元素为 (样本名, 文件路径列表)|List of samples (sample_name, file_paths)
            单末端模式下文件路径列表只有1个文件，双末端有2个文件|Single-end: 1 file, Paired-end: 2 files
        """
        samples = []

        # 文件列表模式|File list mode
        if self.config.input_files:
            self.logger.info("检测到文件列表输入模式|Detected file list input mode")
            for input_file in self.config.input_files:
                file_path = Path(input_file)
                if not file_path.exists():
                    self.logger.warning(f"文件不存在|File does not exist: {input_file}")
                    continue

                # 提取样本名|Extract sample name
                if self.config.sample_names:
                    # 使用用户指定的样本名|Use user-specified sample names
                    idx = self.config.input_files.index(input_file)
                    if idx < len(self.config.sample_names):
                        sample_name = self.config.sample_names[idx]
                    else:
                        sample_name = file_path.stem
                else:
                    sample_name = file_path.stem

                # 检查是否有配对文件|Check for paired file
                file_paths = [file_path]
                if not self.config.single_end:
                    paired_file = self._find_paired_file(file_path)
                    if paired_file:
                        file_paths.append(paired_file)
                        self.logger.info(f"样本|Sample {sample_name}: 找到配对文件|Found paired file")

                samples.append((sample_name, file_paths))

        # 目录模式|Directory mode
        elif self.config.input_path:
            self.logger.info(f"检测到目录输入模式|Detected directory input mode: {self.config.input_path}")
            self.logger.info(f"数据模式|Data mode: {'单末端|Single-end' if self.config.single_end else '双末端|Paired-end'}")
            self.logger.info(f"文件模式|File pattern: *{self.config.read1_suffix}" +
                            (f", *{self.config.read2_suffix}" if not self.config.single_end else ""))

            if self.config.single_end:
                # 单末端模式|Single-end mode
                samples = self._find_single_end_samples()
            else:
                # 双末端模式|Paired-end mode
                samples = self._find_paired_end_samples()

        return samples

    def _find_single_end_samples(self) -> List[Tuple[str, List[Path]]]:
        """查找单末端样本|Find single-end samples"""
        samples = []

        # 查找所有fastq文件|Find all fastq files
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
            return samples

        self.logger.info(f"找到 {len(input_files)} 个单末端文件|Found {len(input_files)} single-end files")

        for input_file in input_files:
            # 提取样本名|Extract sample name
            sample_name = extract_sample_name(input_file.name, self.config.read1_suffix)
            samples.append((sample_name, [input_file]))
            self.logger.debug(f"找到单末端样本|Found single-end sample: {sample_name}")

        return samples

    def _find_paired_end_samples(self) -> List[Tuple[str, List[Path]]]:
        """查找双末端样本配对|Find paired-end sample pairs"""
        samples = []

        # 查找所有read1文件|Find all read1 files
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
            return samples

        self.logger.info(f"找到 {len(read1_files)} 个read1文件|Found {len(read1_files)} read1 files")

        for read1_file in read1_files:
            # 提取样本名|Extract sample name
            sample_name = extract_sample_name(read1_file.name, self.config.read1_suffix)

            # 构建read2文件路径|Construct read2 file path
            read2_file = self.config.input_path / f"{sample_name}{self.config.read2_suffix}"

            if read2_file.exists():
                samples.append((sample_name, [read1_file, read2_file]))
                self.logger.debug(f"找到配对样本|Found paired sample: {sample_name}")
            else:
                self.logger.warning(
                    f"找不到样本 {sample_name} 的配对文件|"
                    f"Cannot find paired file for sample {sample_name}: {read2_file}"
                )

        return samples

    def _find_paired_file(self, read1_file: Path) -> Optional[Path]:
        """
        自动查找配对的Read2文件|Auto-find paired Read2 file

        Args:
            read1_file: Read1文件路径|Read1 file path

        Returns:
            Read2文件路径（如果存在）|Read2 file path (if exists), None otherwise
        """
        # 获取文件所在的目录|Get the directory of the file
        file_dir = read1_file.parent

        # 尝试多种常见的配对文件模式|Try multiple common paired file patterns
        filename = read1_file.name

        # 常见的配对模式|Common pairing patterns
        patterns = [
            filename.replace('_1.', '_2.'),
            filename.replace('.R1.', '.R2.'),
            filename.replace('_R1.', '_R2.'),
            filename.replace('_r1.', '_r2.'),
            filename.replace('.1.', '.2.'),
        ]

        for pattern in patterns:
            paired_file = file_dir / pattern
            if paired_file.exists():
                return paired_file

        return None

    def validate_samples(self, samples: List[Tuple[str, List[Path]]]) -> bool:
        """
        验证样本|Validate samples

        Args:
            samples: 样本列表|List of samples

        Returns:
            验证是否通过|Whether validation passed
        """
        if not samples:
            self.logger.error("未找到有效的样本文件|No valid sample files found")
            return False

        # 检查文件大小|Check file sizes
        for sample_name, file_paths in samples:
            for i, file_path in enumerate(file_paths):
                if file_path.stat().st_size == 0:
                    read_num = f"read{i+1}" if len(file_paths) > 1 else "read"
                    self.logger.warning(f"样本 {sample_name} 的{read_num}文件为空|{read_num.capitalize()} file is empty for sample {sample_name}")

        self.logger.info(f"验证通过|Validation passed: {len(samples)} 个样本|samples")
        return True

    def get_combined_input(self, samples: List[Tuple[str, List[Path]]]) -> List[Tuple[str, str]]:
        """
        获取合并后的输入列表（用于KMC）|Get combined input list (for KMC)

        KMC可以接受多个输入文件，自动处理|KMC can accept multiple input files

        Args:
            samples: 样本列表|List of samples

        Returns:
            (样本名, 输入文件参数)列表|List of (sample_name, input_file_args)
        """
        result = []

        for sample_name, file_paths in samples:
            # KMC支持多个输入文件，用空格分隔|KMC supports multiple input files, space-separated
            input_str = ' '.join(str(f) for f in file_paths)
            result.append((sample_name, input_str))

        return result
