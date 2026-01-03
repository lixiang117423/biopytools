# """
# RNA-seq | RNA-seq Data Processing Module
# """

# import os
# import glob
# import re
# from typing import List, Dict, Tuple

# class FastqPatternParser:
#     """FASTQ | FASTQ File Pattern Parser"""
    
#     def __init__(self, logger):
#         self.logger = logger
    
#     def parse_fastq_pattern(self, pattern: str) -> Tuple[str, str, str, str]:
#         """fastq | Parse fastq file pattern to extract prefix, suffix and read identifiers"""
#         if "*" not in pattern:
#             raise ValueError("* | File pattern must contain * as sample name placeholder")

#         #  | Split pattern string
#         parts = pattern.split("*")
#         if len(parts) != 2:
#             raise ValueError("* | File pattern can only contain one * placeholder")

#         prefix = parts[0]  # * | Part before *
#         suffix = parts[1]  # * | Part after *

#         #  (R1, R2 or 1, 2) | Detect read indicators
#         read_indicators = ["R1", "R2", "_1", "_2", ".1", ".2"]
#         read1_indicator = None
#         read2_indicator = None

#         for indicator in read_indicators:
#             if indicator in suffix:
#                 if indicator.endswith("1") or indicator == "R1":
#                     read1_indicator = indicator
#                     if indicator == "R1":
#                         read2_indicator = "R2"
#                     elif indicator == "_1":
#                         read2_indicator = "_2"
#                     elif indicator == ".1":
#                         read2_indicator = ".2"
#                 break

#         if not read1_indicator:
#             raise ValueError(
#                 f" (R1/R2, _1/_2, .1/.2) | Cannot identify read indicator (R1/R2, _1/_2, .1/.2) in pattern '{pattern}'"
#             )

#         return prefix, suffix, read1_indicator, read2_indicator

# class SampleParser:
#     """ | Sample Parser"""
    
#     def __init__(self, logger):
#         self.logger = logger
#         self.pattern_parser = FastqPatternParser(logger)
    
#     def parse_input_samples(self, input_path: str, fastq_pattern: str = None) -> List[Dict]:
#         """ | Parse input sample information"""
#         samples = []

#         if os.path.isdir(input_path):
#             if fastq_pattern:
#                 #  | Use user-specified file pattern
#                 samples = self._parse_with_pattern(input_path, fastq_pattern)
#             else:
#                 # fastq | Default mode: automatically search for common fastq file pairs
#                 samples = self._parse_with_default_patterns(input_path)
#         else:
#             #  | If it's a file, assume it's a sample information file
#             samples = self._parse_sample_file(input_path)

#         return samples
    
#     def _parse_with_pattern(self, input_path: str, fastq_pattern: str) -> List[Dict]:
#         """ | Parse samples with specified pattern"""
#         samples = []
        
#         try:
#             prefix, suffix, read1_indicator, read2_indicator = self.pattern_parser.parse_fastq_pattern(fastq_pattern)

#             #  | Build search pattern
#             search_pattern = os.path.join(input_path, f"{prefix}*{suffix}")
#             fastq_files = glob.glob(search_pattern)

#             # R1 | Filter R1 files
#             r1_files = [f for f in fastq_files if read1_indicator in os.path.basename(f)]

#             for fq1 in r1_files:
#                 # R2 | Build corresponding R2 file path
#                 fq2 = fq1.replace(read1_indicator, read2_indicator)

#                 if os.path.exists(fq2):
#                     #  | Extract sample name
#                     basename = os.path.basename(fq1)
#                     sample_name = basename.replace(prefix, "").replace(suffix, "")

#                     samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})
#                 else:
#                     self.logger.warning(f"R2 | Warning: Cannot find corresponding R2 file: {fq2}")
        
#         except Exception as e:
#             self.logger.error(f" | Error parsing file pattern: {e}")
        
#         return samples
    
#     def _parse_with_default_patterns(self, input_path: str) -> List[Dict]:
#         """ | Parse samples with default patterns"""
#         samples = []
        
#         #  | Default pattern list
#         patterns = [
#             ("*_1.fq.gz", "*_2.fq.gz"),
#             ("*_R1.fq.gz", "*_R2.fq.gz"),
#             ("*.R1.fastq.gz", "*.R2.fastq.gz"),
#             ("*_1.fastq.gz", "*_2.fastq.gz"),
#             ("*.1.fq.gz", "*.2.fq.gz"),
#             ("*_R1.fq", "*_R2.fq"),
#             ("*_1.fastq", "*_2.fastq"),
#         ]

#         for pattern1, pattern2 in patterns:
#             search_path1 = os.path.join(input_path, pattern1)
#             fastq_files = glob.glob(search_path1)

#             for fq1 in fastq_files:
#                 # R2 | Build corresponding R2 file path
#                 fq2 = fq1.replace(pattern1.replace("*", ""), pattern2.replace("*", ""))

#                 if os.path.exists(fq2):
#                     #  | Extract sample name
#                     basename = os.path.basename(fq1)
#                     #  | Remove fixed parts from pattern to get sample name
#                     sample_name = basename
#                     for to_remove in [pattern1.replace("*", ""), pattern2.replace("*", "")]:
#                         if to_remove in sample_name:
#                             sample_name = sample_name.replace(to_remove, "")
#                             break

#                     #  | Ensure sample name is not empty
#                     if not sample_name:
#                         sample_name = os.path.splitext(os.path.splitext(basename)[0])[0]

#                     samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})

#             #  | If files found, don't try other patterns
#             if samples:
#                 break
        
#         return samples
    
#     def _parse_sample_file(self, input_path: str) -> List[Dict]:
#         """ | Parse sample information file"""
#         samples = []
        
#         # sample_name\tfastq1_path\tfastq2_path | Format: sample_name\tfastq1_path\tfastq2_path
#         try:
#             with open(input_path, "r") as f:
#                 for line_num, line in enumerate(f, 1):
#                     parts = line.strip().split("\t")
#                     if len(parts) >= 3:
#                         samples.append({"name": parts[0], "fastq1": parts[1], "fastq2": parts[2]})
#                     elif line.strip():  #  | Non-empty line but incorrect format
#                         self.logger.warning(f" {line_num}  | Line {line_num} has incorrect format, skipping: {line.strip()}")
#         except Exception as e:
#             self.logger.error(f" | Error reading sample file: {e}")
        
#         return samples

"""
RNA-seq | RNA-seq Data Processing Module
"""

import os
import glob
import re
from typing import List, Dict, Tuple

class FastqPatternParser:
    """FASTQ文件模式解析器|FASTQ File Pattern Parser"""

    def __init__(self, logger):
        self.logger = logger

    def parse_fastq_pattern(self, pattern: str) -> Tuple[str, str, str, str]:
        """解析fastq文件模式以提取前缀、后缀和read标识符|Parse fastq file pattern to extract prefix, suffix and read identifiers"""
        if "*" not in pattern:
            raise ValueError("文件模式必须包含*作为样本名占位符|File pattern must contain * as sample name placeholder")

        # 分割模式字符串|Split pattern string
        parts = pattern.split("*")
        if len(parts) != 2:
            raise ValueError("文件模式只能包含一个*占位符|File pattern can only contain one * placeholder")

        prefix = parts[0]  # *之前的部分|Part before *
        suffix = parts[1]  # *之后的部分|Part after *

        # 检测read标识符 - 扩展支持更多格式|Detect read indicators - extended support for more formats
        read_indicators = [
            # 标准格式|Standard formats
            ("R1", "R2"),
            ("_1", "_2"),
            (".1", ".2"),
            # forward/reverse格式(小写)|forward/reverse formats (lowercase)
            ("_f1", "_r2"),
            ("_f1", "_f2"),  # f1/f2
            (".f1", ".r2"),
            (".f1", ".f2"),
            # forward/reverse格式(大写)|forward/reverse formats (uppercase)
            ("_F1", "_R2"),
            ("_F1", "_F2"),
            (".F1", ".R2"),
            (".F1", ".F2"),
            # 其他常见格式|Other common formats
            ("_read1", "_read2"),
            (".read1", ".read2"),
            ("_pair1", "_pair2"),
            (".pair1", ".pair2"),
        ]

        read1_indicator = None
        read2_indicator = None

        for r1, r2 in read_indicators:
            if r1 in suffix:
                read1_indicator = r1
                read2_indicator = r2
                self.logger.info(f"检测到read标识符|Detected read indicator: {r1} <-> {r2}")
                break

        if not read1_indicator:
            raise ValueError(
                f"无法识别模式中的read标识符|Cannot identify read indicator in pattern '{pattern}'\n"
                f"支持的格式|Supported formats: R1/R2, _1/_2, .1/.2, _f1/_r2, _F1/_R2, _read1/_read2, etc."
            )

        return prefix, suffix, read1_indicator, read2_indicator

class SampleParser:
    """样本解析器|Sample Parser"""

    def __init__(self, logger):
        self.logger = logger
        self.pattern_parser = FastqPatternParser(logger)

    def parse_input_samples(self, input_path: str, fastq_pattern: str = None) -> List[Dict]:
        """解析输入样本信息|Parse input sample information"""
        samples = []

        if os.path.isdir(input_path):
            if fastq_pattern:
                # 使用用户指定的文件模式|Use user-specified file pattern
                self.logger.info(f"使用指定的文件模式|Using specified file pattern: {fastq_pattern}")
                samples = self._parse_with_pattern(input_path, fastq_pattern)
            else:
                # 默认模式：自动搜索常见的fastq文件对|Default mode: automatically search for common fastq file pairs
                self.logger.info("自动检测常见FASTQ文件格式|Auto-detecting common FASTQ file formats")
                samples = self._parse_with_default_patterns(input_path)
        else:
            # 如果是文件，假设它是样本信息文件|If it's a file, assume it's a sample information file
            self.logger.info(f"读取样本信息文件|Reading sample information file: {input_path}")
            samples = self._parse_sample_file(input_path)

        return samples
    
    def _parse_with_pattern(self, input_path: str, fastq_pattern: str) -> List[Dict]:
        """使用指定模式解析样本|Parse samples with specified pattern"""
        samples = []

        try:
            prefix, suffix, read1_indicator, read2_indicator = self.pattern_parser.parse_fastq_pattern(fastq_pattern)

            # 构建搜索模式|Build search pattern
            search_pattern = os.path.join(input_path, fastq_pattern)
            fastq_files = glob.glob(search_pattern)

            self.logger.info(f"搜索模式|Search pattern: {search_pattern}")
            self.logger.info(f"找到 {len(fastq_files)} 个read1文件|Found {len(fastq_files)} read1 files")

            for fq1 in sorted(fastq_files):
                basename = os.path.basename(fq1)

                #
                # 从文件名中提取样本名，移除prefix和suffix
                # 例如："CRR564415_f1.fq.gz" with pattern "*_f1.fq.gz"
                # prefix = "", suffix = "_f1.fq.gz"
                # sample_name = "CRR564415"
                sample_name = basename
                if prefix:
                    sample_name = sample_name.replace(prefix, "", 1)
                if suffix:
                    sample_name = sample_name.replace(suffix, "", 1)

                # 构建read2文件名
                # 将suffix中的read1替换为read2
                # 例如: suffix "_f1.fq.gz" -> read2_suffix "_r2.fq.gz"
                read2_suffix = suffix.replace(read1_indicator, read2_indicator)
                read2_filename = prefix + sample_name + read2_suffix
                fq2 = os.path.join(input_path, read2_filename)

                self.logger.info(f"样本|Sample: {sample_name}")
                self.logger.info(f"  Read1: {basename}")
                self.logger.info(f"  Read2: {read2_filename} (期望|expected)")

                if os.path.exists(fq2):
                    samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})
                    self.logger.info(f"配对成功|Paired successfully")
                else:
                    self.logger.warning(f"找不到配对的read2文件|Cannot find paired read2 file: {fq2}")
                    self.logger.warning(f"跳过该样本|Skipping this sample")

        except Exception as e:
            self.logger.error(f"解析文件模式时出错|Error parsing file pattern: {e}")

        return samples
    
    def _parse_with_default_patterns(self, input_path: str) -> List[Dict]:
        """使用默认模式解析样本|Parse samples with default patterns"""
        samples = []

        # 默认模式列表 - 扩展支持|Default pattern list - extended support
        patterns = [
            # 标准格式|Standard formats
            ("*_1.fq.gz", "*_2.fq.gz"),
            ("*_R1.fq.gz", "*_R2.fq.gz"),
            ("*.R1.fastq.gz", "*.R2.fastq.gz"),
            ("*_1.fastq.gz", "*_2.fastq.gz"),
            ("*.1.fq.gz", "*.2.fq.gz"),
            ("*_R1.fq", "*_R2.fq"),
            ("*_1.fastq", "*_2.fastq"),
            # forward/reverse格式|forward/reverse formats
            ("*_f1.fq.gz", "*_r2.fq.gz"),
            ("*_f1.fq.gz", "*_f2.fq.gz"),
            ("*_F1.fq.gz", "*_R2.fq.gz"),
            ("*.f1.fastq.gz", "*.r2.fastq.gz"),
            ("*_f1.fastq", "*_r2.fastq"),
            # read/pair格式|read/pair formats
            ("*_read1.fq.gz", "*_read2.fq.gz"),
            ("*_pair1.fq.gz", "*_pair2.fq.gz"),
        ]

        for pattern1, pattern2 in patterns:
            # 尝试read1模式|Try read1 pattern
            try:
                # 解析pattern1|Parse pattern1
                prefix, suffix, read1_indicator, read2_indicator = self.pattern_parser.parse_fastq_pattern(pattern1)

                search_path = os.path.join(input_path, pattern1)
                fastq_files = glob.glob(search_path)

                if fastq_files:
                    self.logger.info(f"尝试模式|Trying pattern: {pattern1} / {pattern2}")

                for fq1 in sorted(fastq_files):
                    basename = os.path.basename(fq1)

                    # 提取样本名|Extract sample name
                    sample_name = basename
                    if prefix:
                        sample_name = sample_name.replace(prefix, "", 1)
                    if suffix:
                        sample_name = sample_name.replace(suffix, "", 1)

                    # 构建read2文件名|Build read2 filename
                    read2_suffix = suffix.replace(read1_indicator, read2_indicator)
                    read2_filename = prefix + sample_name + read2_suffix
                    fq2 = os.path.join(input_path, read2_filename)

                    if os.path.exists(fq2):
                        samples.append({"name": sample_name, "fastq1": fq1, "fastq2": fq2})
                        self.logger.info(f"样本配对成功|Sample paired: {sample_name}")

                # 如果找到文件，不再尝试其他模式|If files found, don't try other patterns
                if samples:
                    self.logger.info(f"使用模式|Using pattern: {pattern1} / {pattern2}")
                    break

            except ValueError:
                # 继续尝试下一个模式|Continue to next pattern
                continue

        return samples
    
    def _parse_sample_file(self, input_path: str) -> List[Dict]:
        """解析样本信息文件|Parse sample information file"""
        samples = []

        # 格式: sample_name\tfastq1_path\tfastq2_path|Format: sample_name\tfastq1_path\tfastq2_path
        try:
            with open(input_path, "r") as f:
                for line_num, line in enumerate(f, 1):
                    parts = line.strip().split("\t")
                    if len(parts) >= 3:
                        samples.append({"name": parts[0], "fastq1": parts[1], "fastq2": parts[2]})
                        self.logger.info(f"行 {line_num}|Line {line_num}: {parts[0]}")
                    elif line.strip():  # 非空行但格式不正确|Non-empty line but incorrect format
                        self.logger.warning(f"行 {line_num} 格式不正确，跳过|Line {line_num} has incorrect format, skipping: {line.strip()}")
        except Exception as e:
            self.logger.error(f"读取样本文件时出错|Error reading sample file: {e}")

        return samples