"""
VCF文件修复模块|VCF File Repair Module

用于修复损坏的VCF文件，主要是处理列数不匹配的问题
Repair corrupted VCF files, mainly handling column count mismatches
"""

import gzip
import logging
from pathlib import Path
from typing import Optional, Tuple


class VCFRepairer:
    """VCF文件修复器|VCF File Repairer"""

    def __init__(self, logger: logging.Logger):
        """
        初始化VCF修复器|Initialize VCF repairer

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def check_and_repair(self, input_file: str, output_file: str,
                        auto_repair: bool = False) -> Tuple[bool, Optional[str]]:
        """
        检查并修复VCF文件|Check and repair VCF file

        Args:
            input_file: 输入VCF文件路径|Input VCF file path
            output_file: 输出VCF文件路径|Output VCF file path
            auto_repair: 是否自动修复|Whether to auto-repair

        Returns:
            (是否需要修复, 修复后的文件路径)|(Whether repair needed, repaired file path)
        """
        self.logger.info("检查VCF文件完整性|Checking VCF file integrity")
        self.logger.info(f"输入文件|Input file: {input_file}")

        # 检查文件是否损坏|Check if file is corrupted
        is_corrupted, error_info = self._check_vcf_integrity(input_file)

        if not is_corrupted:
            self.logger.info("VCF文件完整，无需修复|VCF file is intact, no repair needed")
            return False, None

        self.logger.warning(f"检测到VCF文件损坏|Detected VCF file corruption: {error_info}")

        if not auto_repair:
            self.logger.error("文件损坏但未启用自动修复，请使用--repair-vcf参数|File corrupted but auto-repair not enabled, please use --repair-vcf parameter")
            return True, None

        # 执行修复|Perform repair
        self.logger.info("开始修复VCF文件|Starting VCF file repair")
        success = self.repair_vcf(input_file, output_file)

        if success:
            self.logger.info(f"VCF文件修复成功|VCF file repaired successfully: {output_file}")
            return True, output_file
        else:
            self.logger.error("VCF文件修复失败|VCF file repair failed")
            return True, None

    def _check_vcf_integrity(self, vcf_file: str) -> Tuple[bool, Optional[str]]:
        """
        检查VCF文件完整性|Check VCF file integrity

        Args:
            vcf_file: VCF文件路径|VCF file path

        Returns:
            (是否损坏, 错误信息)|(Whether corrupted, error message)
        """
        expected_columns = -1
        line_count = 0

        try:
            with gzip.open(vcf_file, 'rt', encoding='utf-8') as f:
                for line in f:
                    line_count += 1
                    stripped_line = line.strip()

                    # 跳过空行|Skip empty lines
                    if not stripped_line:
                        continue

                    # 处理表头|Process header
                    if stripped_line.startswith("#"):
                        if stripped_line.startswith("#CHROM"):
                            parts = stripped_line.split('\t')
                            expected_columns = len(parts)
                        continue

                    # 检查数据行列数|Check data row column count
                    if expected_columns == -1:
                        return True, "未找到#CHROM表头|#CHROM header not found"

                    parts = stripped_line.split('\t')
                    current_cols = len(parts)

                    if current_cols != expected_columns:
                        error_msg = (f"第{line_count}行列数不匹配|Row {line_count} column count mismatch: "
                                   f"期望{expected_columns}列，实际{current_cols}列|Expected {expected_columns}, got {current_cols}")
                        return True, error_msg

                    # 每100000行报告一次进度|Report progress every 100k lines
                    if line_count % 100000 == 0:
                        self.logger.debug(f"已检查|Checked {line_count} 行|lines")

        except EOFError:
            return True, f"GZIP文件在第{line_count}行意外结束|GZIP file ended unexpectedly at line {line_count}"
        except Exception as e:
            return True, f"读取文件时发生错误|Error reading file: {str(e)}"

        return False, None

    def repair_vcf(self, input_file: str, output_file: str) -> bool:
        """
        修复VCF文件|Repair VCF file

        遇到列数不匹配的行时停止写入，保证输出文件结构完整
        Stop writing when encountering column mismatch, ensuring output file structural integrity

        Args:
            input_file: 输入VCF文件路径|Input VCF file path
            output_file: 输出VCF文件路径|Output VCF file path

        Returns:
            修复是否成功|Whether repair succeeded
        """
        expected_columns = -1
        line_count = 0
        valid_lines = 0
        error_line = -1

        try:
            with gzip.open(input_file, 'rt', encoding='utf-8') as f_in, \
                 gzip.open(output_file, 'wt', encoding='utf-8') as f_out:

                for line in f_in:
                    line_count += 1
                    stripped_line = line.strip()

                    # 跳过空行|Skip empty lines
                    if not stripped_line:
                        continue

                    # 处理表头|Process header
                    if stripped_line.startswith("#"):
                        f_out.write(line)
                        if stripped_line.startswith("#CHROM"):
                            parts = stripped_line.split('\t')
                            expected_columns = len(parts)
                            self.logger.info(f"检测到表头|Header detected: {expected_columns} 列|columns")
                        continue

                    # 检查数据行|Check data row
                    if expected_columns == -1:
                        self.logger.error("在找到#CHROM表头之前遇到数据行|Data row encountered before #CHROM header")
                        return False

                    parts = stripped_line.split('\t')
                    current_cols = len(parts)

                    # 核心检查逻辑|Core check logic
                    if current_cols != expected_columns:
                        self.logger.warning(f"发现损坏行|Found corrupted row at line {line_count}")
                        self.logger.warning(f"期望列数|Expected columns: {expected_columns}, "
                                         f"实际列数|Actual columns: {current_cols}")
                        error_line = line_count
                        break

                    # 写入正常行|Write valid row
                    f_out.write(line)
                    valid_lines += 1

                    # 显示进度|Show progress
                    if valid_lines % 100000 == 0:
                        self.logger.info(f"已处理|Processed {valid_lines} 有效行|valid lines")

        except EOFError:
            self.logger.warning(f"GZIP文件在第{line_count}行意外结束|GZIP file ended unexpectedly at line {line_count}")
        except Exception as e:
            self.logger.error(f"修复过程发生错误|Error during repair: {str(e)}")
            return False

        # 输出统计信息|Output statistics
        self.logger.info("=" * 60)
        self.logger.info("VCF修复完成|VCF repair completed")
        self.logger.info("=" * 60)
        self.logger.info(f"共读取行数|Total lines read: {line_count}")
        self.logger.info(f"有效行数|Valid lines: {valid_lines}")

        if error_line > 0:
            self.logger.warning(f"在第{error_line}行停止写入|Stopped writing at line {error_line}")

        self.logger.info(f"修复后文件|Repaired file: {output_file}")
        self.logger.info("=" * 60)

        return valid_lines > 0
