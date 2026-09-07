"""
CNCB链接提取器工具函数模块|CNCB Link Extractor Utility Functions Module
"""

import re
import time
import ftplib
import logging
import os
import stat
from collections import Counter, defaultdict
from typing import Dict, List, Tuple, Optional, Set
from pathlib import Path

# ID文件分隔符嗅探:Tab/分号/逗号/空白任意组合(Excel导出常见混排)
# |Delimiter sniffing for ID files: tab/semicolon/comma/whitespace in any combination
_FIELD_SEPARATOR_RE = re.compile(r'[\t;,]+|\s+')


def get_requests():
    """
    惰性导入requests模块,未安装时返回None|Lazily import requests; return None if not installed

    get_link_from_CNCB的HTTP搜索器(ena_searcher/gsa_searcher)共用此函数,
    避免顶层导入让纯CNCB模式也依赖requests
    |Shared by the module's HTTP searchers; avoids a top-level import that would
    make even pure CNCB mode depend on requests
    """
    try:
        import requests
        return requests
    except ImportError:
        return None


def create_retry_session(requests_module, retry_attempts: int):
    """
    创建带重试适配器的HTTP会话(模块内ena/gsa/ncbi搜索器共用)
    |Create an HTTP session with a retry adapter (shared by the module's ena/gsa/ncbi searchers)

    Args:
        requests_module: 已导入的requests模块|The imported requests module
        retry_attempts: 重试次数|Retry count

    Returns:
        挂载了Retry适配器的Session|A Session with a Retry adapter mounted
    """
    from requests.adapters import HTTPAdapter
    from urllib3.util.retry import Retry

    session = requests_module.Session()
    retries = Retry(
        total=retry_attempts,
        backoff_factor=0.5,
        status_forcelist=[500, 502, 503, 504]
    )
    session.mount('https://', HTTPAdapter(max_retries=retries))
    return session


class CNCBLogger:
    """CNCB日志管理器|CNCB Logger Manager"""

    def __init__(self, log_file: Optional[str] = None, verbose: bool = False):
        self.log_file = log_file
        self.verbose = verbose
        self.logger = None

        self._setup_logger()

    def _setup_logger(self):
        """设置日志记录器|Setup logger"""
        self.logger = logging.getLogger("cncb_link_extractor")
        self.logger.setLevel(logging.DEBUG)

        # 清除现有处理器
        self.logger.handlers.clear()

        # 控制台处理器
        console_handler = logging.StreamHandler()
        console_handler.setLevel(logging.INFO if not self.verbose else logging.DEBUG)

        # 日志格式|Log format
        log_format = '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s'
        date_format = '%Y-%m-%d %H:%M:%S'
        formatter = logging.Formatter(log_format, datefmt=date_format)
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
                self.logger.warning(f"无法创建日志文件|Cannot create log file: {e}")

    def get_logger(self):
        """获取日志记录器|Get logger"""
        return self.logger


class InputFileParser:
    """输入文件解析器|Input File Parser"""

    @staticmethod
    def _split_fields(line: str) -> List[str]:
        """
        按Tab/分号/逗号/空白拆分字段并过滤空段|Split fields on tab/semicolon/comma/whitespace and drop empties

        ID文件的分隔符来源不可控(手工整理/Excel导出/其他工具输出),
        逐行自动嗅探而非要求固定分隔符;"A, B"这类混排拆出的空段直接丢弃
        |ID files come from uncontrolled sources (hand-edited, Excel exports,
        other tools), so delimiters are sniffed per line instead of fixed;
        empty segments from mixes like "A, B" are dropped
        """
        return [p.strip() for p in _FIELD_SEPARATOR_RE.split(line) if p.strip()]

    @staticmethod
    def _looks_like_header(fields: List[str]) -> bool:
        """
        两列均为纯字母视为表头行|All-alphabetic first two columns mean a header row

        真实accession(PRJNA/ERR/SRR/CRR等)必含数字,纯字母列只会出现在
        "ProjectID,RunID"这类表头里
        |Real accessions (PRJNA/ERR/SRR/CRR...) always contain digits;
        all-alphabetic columns only appear in headers like "ProjectID,RunID"
        """
        return len(fields) >= 2 and not any(
            any(c.isdigit() for c in f) for f in fields[:2]
        )

    @staticmethod
    def read_and_group_by_project(input_file: str, logger=None) -> Optional[Dict[str, List[str]]]:
        """
        读取输入文件并按项目分组|Read input file and group by project

        每行两列:ProjectID与RunID,分隔符自动识别(Tab/分号/逗号/空格)
        |Two columns per line: ProjectID and RunID, delimiter auto-detected
        (tab/semicolon/comma/space)

        Args:
            input_file: 输入文件路径|Input file path
            logger: 日志记录器|Logger instance

        Returns:
            按项目分组的Run ID列表|Run IDs grouped by project, or None if error
        """
        # 使用传入的logger或默认的
        log = logger or logging.getLogger(__name__)

        projects = defaultdict(list)
        delimiter_rows = Counter()

        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # 跳过空行和注释行
                    if not line or line.startswith('#'):
                        continue

                    # 记录本行用到的分隔符类别,供解析完成后汇报
                    # |Track which delimiter classes this line used, reported after parsing
                    if '\t' in line:
                        delimiter_rows['tab'] += 1
                    if ';' in line:
                        delimiter_rows['semicolon'] += 1
                    if ',' in line:
                        delimiter_rows['comma'] += 1
                    if re.search(r'\s', line):
                        delimiter_rows['space'] += 1

                    fields = InputFileParser._split_fields(line)

                    if InputFileParser._looks_like_header(fields):
                        log.info(f"第{line_num}行检测到表头行,已跳过|Line {line_num} looks like a header row, skipped -> '{line}'")
                        continue

                    if len(fields) == 1:
                        log.warning(f"第{line_num}行|Line {line_num}: 只有一个字段,应为ProjectID和RunID两列|Only one field, expected two columns (ProjectID and Run ID) -> '{line}'")
                        continue

                    if len(fields) > 2:
                        log.warning(f"第{line_num}行|Line {line_num}: 超过两列,取前两列并忽略其余|More than two columns, keeping the first two and ignoring the rest -> '{line}'")

                    projects[fields[0]].append(fields[1])

            # 对每个项目的Run ID去重并排序
            # 原实现先sort再set,set会破坏已排序的顺序
            # |Dedup then sort; the old sort-then-set destroyed the sorted order
            for project_id in projects:
                projects[project_id] = sorted(set(projects[project_id]))

            log.info(f"成功解析文件|Successfully parsed file: {input_file}")
            log.info(f"发现|Found {len(projects)} 个项目|projects，总计|total {sum(len(ids) for ids in projects.values())} 个Run IDs")
            if delimiter_rows:
                stats = ", ".join(f"{name}={count}" for name, count in sorted(delimiter_rows.items()))
                log.info(f"检测到的分隔符(出现行数)|Detected delimiters (rows containing): {stats}")

            return projects

        except FileNotFoundError:
            log.error(f"错误|Error: 输入文件|Input file {input_file} 未找到|not found!")
            return None
        except Exception as e:
            log.error(f"读取输入文件时发生错误|Error reading input file: {e}")
            return None

    @staticmethod
    def validate_input_file(input_file: str) -> Tuple[bool, str]:
        """
        验证输入文件格式|Validate input file format

        校验口径与read_and_group_by_project一致:分隔符自动识别,
        两列及以上即有效(超出两列解析时取前两列),表头行跳过不计数
        |Same rules as read_and_group_by_project: delimiters are auto-detected,
        two or more columns count as valid (extras are dropped at parse time),
        header rows are skipped
        """
        if not os.path.exists(input_file):
            return False, f"文件不存在|File does not exist: {input_file}"

        if not os.path.isfile(input_file):
            return False, f"路径不是文件|Path is not a file: {input_file}"

        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                valid_lines = 0

                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    if not line or line.startswith('#'):
                        continue

                    fields = InputFileParser._split_fields(line)

                    if InputFileParser._looks_like_header(fields):
                        continue

                    if len(fields) >= 2:
                        valid_lines += 1
                    else:
                        return False, f"第{line_num}行格式错误,只解析出一个字段,应为ProjectID和RunID两列|Line {line_num} format error: only one field parsed, expected two columns (ProjectID and Run ID): '{line}'"

                if valid_lines == 0:
                    return False, "文件中没有有效的数据行|No valid data lines found in file"

                return True, f"文件验证通过，包含 {valid_lines} 个有效数据行|File validation passed, contains {valid_lines} valid data lines"

        except UnicodeDecodeError:
            return False, "文件编码错误，请使用UTF-8编码|File encoding error, please use UTF-8 encoding"
        except Exception as e:
            return False, f"读取文件时发生错误|Error reading file: {e}"


class FTPConnectionManager:
    """FTP连接管理器|FTP Connection Manager"""

    def __init__(self, host: str, timeout: int = 60, retry_attempts: int = 3, logger=None):
        self.host = host
        self.timeout = timeout
        self.retry_attempts = retry_attempts
        self.ftp = None
        self.logger = logger or logging.getLogger(__name__)

    def connect(self) -> bool:
        """连接FTP服务器|Connect to FTP server"""
        for attempt in range(self.retry_attempts):
            try:
                self.logger.info(f"正在连接到FTP服务器|Connecting to FTP server: {self.host} (尝试|attempt {attempt + 1}/{self.retry_attempts})")
                self.ftp = ftplib.FTP(self.host, timeout=self.timeout)
                self.ftp.login()
                self.logger.info("FTP连接成功|FTP connection successful")
                return True

            except Exception as e:
                self.logger.warning(f"FTP连接失败|FTP connection failed (尝试|attempt {attempt + 1}): {e}")
                if attempt < self.retry_attempts - 1:
                    time.sleep(2 ** attempt)  # 指数退避
                else:
                    self.logger.error(f"FTP连接失败，已达到最大重试次数|FTP connection failed, max retries reached")
                    return False

    def disconnect(self):
        """断开FTP连接|Disconnect from FTP server"""
        if self.ftp:
            try:
                self.ftp.quit()
                self.logger.info("FTP连接已断开|FTP connection disconnected")
            except Exception as e:
                self.logger.warning(f"断开FTP连接时发生错误|Error disconnecting from FTP: {e}")
            finally:
                self.ftp = None

    def test_connection(self) -> bool:
        """测试FTP连接|Test FTP connection"""
        if not self.ftp:
            return False

        try:
            self.ftp.voidcmd("NOOP")
            return True
        except Exception:
            return False

    def get_ftp(self):
        """获取FTP对象|Get FTP object"""
        return self.ftp


class PathCache:
    """路径缓存管理器|Path Cache Manager"""

    def __init__(self):
        self.cache: Dict[str, str] = {}
        self.hits = 0
        self.misses = 0

    def get(self, key: str) -> Optional[str]:
        """获取缓存路径|Get cached path"""
        if key in self.cache:
            self.hits += 1
            return self.cache[key]
        self.misses += 1
        return None

    def set(self, key: str, value: str):
        """设置缓存路径|Set cached path"""
        self.cache[key] = value

    def get_stats(self) -> Dict[str, int]:
        """获取缓存统计信息|Get cache statistics"""
        total = self.hits + self.misses
        hit_rate = (self.hits / total * 100) if total > 0 else 0
        return {
            "hits": self.hits,
            "misses": self.misses,
            "hit_rate_percent": round(hit_rate, 2)
        }


class FileDownloader:
    """文件下载脚本生成器|File Download Script Generator"""

    @staticmethod
    def generate_download_script(urls: List[str], script_file: str,
                               make_executable: bool = True, logger=None) -> bool:
        """
        生成下载脚本|Generate download script

        Args:
            urls: URL列表|List of URLs
            script_file: 脚本文件路径|Script file path
            make_executable: 是否设置执行权限|Whether to set executable permission
            logger: 日志记录器|Logger instance

        Returns:
            是否成功|Success status
        """
        # 使用传入的logger或默认的
        log = logger or logging.getLogger(__name__)

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
                    os.chmod(script_file, st.st_mode|stat.S_IEXEC)
                    log.info(f"下载脚本已生成并设置执行权限|Download script generated with executable permission: {script_file}")
                except Exception as e:
                    log.warning(f"无法设置脚本执行权限|Cannot set script executable permission: {e}")
                    log.info(f"请手动运行|Please manually run: chmod +x {script_file}")
            else:
                log.info(f"下载脚本已生成|Download script generated: {script_file}")

            return True

        except Exception as e:
            log.error(f" 生成下载脚本失败|Failed to generate download script: {e}")
            return False

    @staticmethod
    def generate_summary_report(projects: Dict[str, List[str]],
                              successful_urls: List[str],
                              failed_ids: List[Tuple[str, str]],
                              output_dir: str,
                              logger=None) -> str:
        """
        生成总结报告|Generate summary report

        Args:
            projects: 项目数据|Project data
            successful_urls: 成功的URL列表|List of successful URLs
            failed_ids: 失败的ID列表|List of failed IDs
            output_dir: 输出目录|Output directory
            logger: 日志记录器|Logger instance

        Returns:
            报告文件路径|Report file path
        """
        # 使用传入的logger或默认的
        log = logger or logging.getLogger(__name__)

        report_file = os.path.join(output_dir, "CNCB_download_report.txt")

        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("CNCB数据链接提取报告|CNCB Data Link Extraction Report\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"生成时间|Generated Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                # 项目统计|Project statistics
                f.write("项目统计|Project Statistics:\n")
                f.write("-" * 30 + "\n")
                total_projects = len(projects)
                total_ids = sum(len(ids) for ids in projects.values())
                f.write(f"总项目数|Total Projects: {total_projects}\n")
                f.write(f"总Run ID数|Total Run IDs: {total_ids}\n\n")

                # 成功统计|Success statistics
                f.write("成功统计|Success Statistics:\n")
                f.write("-" * 30 + "\n")
                f.write(f"成功链接数|Successful URLs: {len(successful_urls)}\n")
                success_rate = (len(successful_urls) / total_ids * 100) if total_ids > 0 else 0
                f.write(f"成功率|Success Rate: {success_rate:.2f}%\n\n")

                # 失败统计|Failure statistics
                f.write("失败统计|Failure Statistics:\n")
                f.write("-" * 30 + "\n")
                f.write(f"失败ID数|Failed IDs: {len(failed_ids)}\n")

                if failed_ids:
                    f.write("\n失败的Run IDs列表|List of Failed Run IDs:\n")
                    for project, run_id in failed_ids:
                        f.write(f"{project}\t{run_id}\n")

                f.write("\n" + "=" * 60 + "\n")
                f.write("报告生成完成|Report generation completed\n")

            return report_file

        except Exception as e:
            log.error(f" 生成报告失败|Failed to generate report: {e}")
            return ""