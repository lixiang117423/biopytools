"""
CNCB链接提取器工具函数模块 | CNCB Link Extractor Utility Functions Module
"""

import time
import ftplib
import logging
import os
import stat
from collections import defaultdict
from typing import Dict, List, Tuple, Optional, Set
from pathlib import Path


class CNCBLogger:
    """CNCB日志管理器 | CNCB Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, verbose: bool = False):
        self.log_file = log_file
        self.verbose = verbose
        self.logger = None

        self._setup_logger()

    def _setup_logger(self):
        """设置日志记录器 | Setup logger"""
        self.logger = logging.getLogger("cncb_link_extractor")
        self.logger.setLevel(logging.DEBUG)

        # 清除现有处理器
        self.logger.handlers.clear()

        # 控制台处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO if not self.verbose else logging.DEBUG)

        # 日志格式
        formatter = logging.Formatter('[%(asctime)s] %(message)s',
                                    datefmt='%Y-%m-%d %H:%M:%S')
        console_handler.setFormatter(formatter)
        self.logger.addHandler(console_handler)

        # 文件处理器
        if self.log_file:
            try:
                log_dir = os.path.dirname(os.path.abspath(self.log_file))
                Path(log_dir).mkdir(parents=True, exist_ok=True)

                file_handler = logging.FileHandler(self.log_file, encoding='utf-8')
                file_handler.setLevel(logging.DEBUG)
                file_handler.setFormatter(formatter)
                self.logger.addHandler(file_handler)
            except Exception as e:
                self.logger.warning(f"无法创建日志文件 | Cannot create log file: {e}")

    def get_logger(self):
        """获取日志记录器 | Get logger"""
        return self.logger


class InputFileParser:
    """输入文件解析器 | Input File Parser"""

    @staticmethod
    def read_and_group_by_project(input_file: str) -> Optional[Dict[str, List[str]]]:
        """
        读取输入文件并按项目分组 | Read input file and group by project

        Args:
            input_file: 输入文件路径 | Input file path

        Returns:
            按项目分组的Run ID列表 | Run IDs grouped by project, or None if error
        """
        projects = defaultdict(list)

        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                line_count = 0

                for line_num, line in enumerate(f, 1):
                    line_count += 1
                    line = line.strip()

                    # 跳过空行和注释行
                    if not line or line.startswith('#'):
                        continue

                    # 解析ProjectID和RunID
                    parts = line.split('\t')
                    if len(parts) == 2:
                        project_id, run_id = parts
                        project_id = project_id.strip()
                        run_id = run_id.strip()

                        if project_id and run_id:
                            projects[project_id].append(run_id)
                        else:
                            logging.warning(f"⚠️ 第{line_num}行: 项目ID或Run ID为空 -> '{line}'")
                    else:
                        logging.warning(f"⚠️ 第{line_num}行: 格式不正确，应为两列Tab分隔 -> '{line}'")

            # 对每个项目的Run ID进行排序
            for project_id in projects:
                projects[project_id].sort()
                # 去重
                projects[project_id] = list(set(projects[project_id]))

            logging.info(f"📄 成功解析文件 {input_file}")
            logging.info(f"📊 发现 {len(projects)} 个项目，总计 {sum(len(ids) for ids in projects.values())} 个Run ID")

            return projects

        except FileNotFoundError:
            logging.error(f"❌ 错误：输入文件 {input_file} 未找到！")
            return None
        except Exception as e:
            logging.error(f"❌ 读取输入文件时发生错误 | Error reading input file: {e}")
            return None

    @staticmethod
    def validate_input_file(input_file: str) -> Tuple[bool, str]:
        """
        验证输入文件格式 | Validate input file format

        Args:
            input_file: 输入文件路径 | Input file path

        Returns:
            (是否有效, 错误信息) | (is_valid, error_message)
        """
        if not os.path.exists(input_file):
            return False, f"文件不存在 | File does not exist: {input_file}"

        if not os.path.isfile(input_file):
            return False, f"路径不是文件 | Path is not a file: {input_file}"

        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                valid_lines = 0
                total_lines = 0

                for line_num, line in enumerate(f, 1):
                    total_lines += 1
                    line = line.strip()

                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) == 2 and parts[0].strip() and parts[1].strip():
                        valid_lines += 1
                    else:
                        return False, f"第{line_num}行格式错误 | Line {line_num} format error: '{line}'"

                if valid_lines == 0:
                    return False, "文件中没有有效的数据行 | No valid data lines found in file"

                return True, f"文件验证通过，包含 {valid_lines} 个有效数据行 | File validation passed, contains {valid_lines} valid data lines"

        except UnicodeDecodeError:
            return False, "文件编码错误，请使用UTF-8编码 | File encoding error, please use UTF-8 encoding"
        except Exception as e:
            return False, f"读取文件时发生错误 | Error reading file: {e}"


class FTPConnectionManager:
    """FTP连接管理器 | FTP Connection Manager"""

    def __init__(self, host: str, timeout: int = 60, retry_attempts: int = 3):
        self.host = host
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        self.ftp = None

    def connect(self) -> bool:
        """连接FTP服务器 | Connect to FTP server"""
        for attempt in range(self.retry_attempts):
            try:
                logging.info(f"🌐 正在连接到FTP服务器: {self.host} (尝试 {attempt + 1}/{self.retry_attempts})")
                self.ftp = ftplib.FTP(self.host, timeout=self.timeout)
                self.ftp.login()
                logging.info("✅ FTP连接成功！")
                return True

            except Exception as e:
                logging.warning(f"⚠️ FTP连接失败 (尝试 {attempt + 1}): {e}")
                if attempt < self.retry_attempts - 1:
                    time.sleep(2 ** attempt)  # 指数退避
                else:
                    logging.error(f"❌ FTP连接失败，已达到最大重试次数 | FTP connection failed, max retries reached")
                    return False

    def disconnect(self):
        """断开FTP连接 | Disconnect from FTP server"""
        if self.ftp:
            try:
                self.ftp.quit()
                logging.info("📤 FTP连接已断开")
            except Exception as e:
                logging.warning(f"⚠️ 断开FTP连接时发生错误 | Error disconnecting from FTP: {e}")
            finally:
                self.ftp = None

    def test_connection(self) -> bool:
        """测试FTP连接 | Test FTP connection"""
        if not self.ftp:
            return False

        try:
            self.ftp.voidcmd("NOOP")
            return True
        except Exception:
            return False

    def get_ftp(self):
        """获取FTP对象 | Get FTP object"""
        return self.ftp


class PathCache:
    """路径缓存管理器 | Path Cache Manager"""

    def __init__(self):
        self.cache: Dict[str, str] = {}
        self.hits = 0
        self.misses = 0

    def get(self, key: str) -> Optional[str]:
        """获取缓存路径 | Get cached path"""
        if key in self.cache:
            self.hits += 1
            return self.cache[key]
        self.misses += 1
        return None

    def set(self, key: str, value: str):
        """设置缓存路径 | Set cached path"""
        self.cache[key] = value

    def get_stats(self) -> Dict[str, int]:
        """获取缓存统计信息 | Get cache statistics"""
        total = self.hits + self.misses
        hit_rate = (self.hits / total * 100) if total > 0 else 0
        return {
            "hits": self.hits,
            "misses": self.misses,
            "hit_rate_percent": round(hit_rate, 2)
        }


class FileDownloader:
    """文件下载脚本生成器 | File Download Script Generator"""

    @staticmethod
    def generate_download_script(urls: List[str], script_file: str,
                               make_executable: bool = True) -> bool:
        """
        生成下载脚本 | Generate download script

        Args:
            urls: URL列表 | List of URLs
            script_file: 脚本文件路径 | Script file path
            make_executable: 是否设置执行权限 | Whether to set executable permission

        Returns:
            是否成功 | Success status
        """
        try:
            with open(script_file, 'w', encoding='utf-8') as sh_file:
                # 写入shebang和头部信息
                sh_file.write("#!/bin/bash\n")
                sh_file.write("# Auto-generated download script for CNCB data\n")
                sh_file.write(f"# Generated on: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                sh_file.write(f"# Total URLs: {len(urls)}\n\n")

                # 写入下载命令
                for url in sorted(urls):
                    sh_file.write(f"wget -c '{url}'\n")

                # 添加总结信息
                sh_file.write("\necho 'Download script completed'\n")

            # 设置执行权限
            if make_executable:
                try:
                    st = os.stat(script_file)
                    os.chmod(script_file, st.st_mode | stat.S_IEXEC)
                    logging.info(f"✅ 下载脚本已生成并设置执行权限: {script_file}")
                except Exception as e:
                    logging.warning(f"⚠️ 无法设置脚本执行权限: {e}")
                    logging.info(f"💡 请手动运行: chmod +x {script_file}")
            else:
                logging.info(f"✅ 下载脚本已生成: {script_file}")

            return True

        except Exception as e:
            logging.error(f"❌ 生成下载脚本失败 | Failed to generate download script: {e}")
            return False

    @staticmethod
    def generate_summary_report(projects: Dict[str, List[str]],
                              successful_urls: List[str],
                              failed_ids: List[Tuple[str, str]],
                              output_dir: str) -> str:
        """
        生成总结报告 | Generate summary report

        Args:
            projects: 项目数据 | Project data
            successful_urls: 成功的URL列表 | List of successful URLs
            failed_ids: 失败的ID列表 | List of failed IDs
            output_dir: 输出目录 | Output directory

        Returns:
            报告文件路径 | Report file path
        """
        report_file = os.path.join(output_dir, "CNCB_download_report.txt")

        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("CNCB数据链接提取报告 | CNCB Data Link Extraction Report\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"生成时间 | Generated Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                # 项目统计
                f.write("📊 项目统计 | Project Statistics:\n")
                f.write("-" * 30 + "\n")
                total_projects = len(projects)
                total_ids = sum(len(ids) for ids in projects.values())
                f.write(f"总项目数 | Total Projects: {total_projects}\n")
                f.write(f"总Run ID数 | Total Run IDs: {total_ids}\n\n")

                # 成功统计
                f.write("✅ 成功统计 | Success Statistics:\n")
                f.write("-" * 30 + "\n")
                f.write(f"成功链接数 | Successful URLs: {len(successful_urls)}\n")
                success_rate = (len(successful_urls) / total_ids * 100) if total_ids > 0 else 0
                f.write(f"成功率 | Success Rate: {success_rate:.2f}%\n\n")

                # 失败统计
                f.write("❌ 失败统计 | Failure Statistics:\n")
                f.write("-" * 30 + "\n")
                f.write(f"失败ID数 | Failed IDs: {len(failed_ids)}\n")

                if failed_ids:
                    f.write("\n失败的Run IDs列表 | List of Failed Run IDs:\n")
                    for project, run_id in failed_ids:
                        f.write(f"{project}\t{run_id}\n")

                f.write("\n" + "=" * 60 + "\n")
                f.write("报告生成完成 | Report generation completed\n")

            return report_file

        except Exception as e:
            logging.error(f"❌ 生成报告失败 | Failed to generate report: {e}")
            return ""