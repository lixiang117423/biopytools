"""SNP区域基因提取工具函数|SNP Region Gene Extractor Utility Functions"""

import logging
import sys
import subprocess


class SnpRegionLogger:
    """SNP区域基因提取日志管理器|SNP Region Gene Extractor Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'

        level = getattr(logging, log_level.upper(), logging.INFO)

        handlers = [logging.StreamHandler(sys.stdout)]
        if self.log_file:
            handlers.append(logging.FileHandler(self.log_file))

        logging.basicConfig(
            level=level,
            format=log_format,
            datefmt=date_format,
            handlers=handlers
        )

        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def run_command(command, logger):
    """
    运行命令并返回结果|Run command and return result

    参数|Parameters:
        command: 命令列表|Command list
        logger: 日志器|Logger

    返回|Returns:
        tuple: (success: bool, stdout: str, stderr: str)
    """
    try:
        logger.debug(f"执行命令|Executing command: {' '.join(command)}")

        result = subprocess.run(
            command,
            check=True,
            capture_output=True,
            text=True
        )

        return True, result.stdout, result.stderr

    except subprocess.CalledProcessError as e:
        return False, e.stdout, e.stderr

    except FileNotFoundError as e:
        return False, "", str(e)


def parse_snp_file(file_path, logger):
    """
    智能识别并解析SNP文件|Intelligently parse SNP file

    参数|Parameters:
        file_path: SNP文件路径|SNP file path
        logger: 日志器|Logger

    返回|Returns:
        list: SNP列表 [(chrom, pos), ...]
    """
    snp_list = []
    seen = set()

    try:
        with open(file_path, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # 空行|Empty line
                if not line:
                    continue

                # 注释行|Comment line
                if line.startswith('#'):
                    continue

                # 表头行（包含SNP等关键词）|Header line containing keywords
                if any(keyword in line.upper() for keyword in ['SNP', 'POSITION', 'COORDINATE']):
                    logger.debug(f"跳过表头行|Skipping header line: {line}")
                    continue

                # 数据行：必须包含冒号|Data line: must contain colon
                if ':' not in line:
                    logger.debug(f"跳过非数据行（无冒号）|Skipping non-data line (no colon): {line}")
                    continue

                # 解析 Chr:Pos 格式|Parse Chr:Pos format
                try:
                    parts = line.split(':')
                    if len(parts) >= 2:
                        chrom = parts[0].strip()
                        pos = int(parts[1].strip())
                        snp_id = f"{chrom}:{pos}"

                        if snp_id not in seen:
                            seen.add(snp_id)
                            snp_list.append((chrom, pos))
                        else:
                            logger.debug(f"跳过重复SNP|Skipping duplicate SNP: {snp_id}")

                except (ValueError, IndexError) as e:
                    logger.debug(f"跳过格式错误的行|Skipping malformed line: {line} | Error: {e}")
                    continue

        logger.info(f"成功解析SNP文件|Successfully parsed SNP file: {len(snp_list)} 个唯一SNP|unique SNPs")

        return snp_list

    except FileNotFoundError:
        logger.error(f"SNP文件不存在|SNP file not found: {file_path}")
        raise

    except Exception as e:
        logger.error(f"解析SNP文件出错|Error parsing SNP file: {e}")
        raise
