"""
Reads提取工具函数模块|Reads Extraction Utility Functions Module
"""

import gzip
import logging
import sys


class ExtractReadsLogger:
    """Reads提取日志管理器|Reads Extraction Logger Manager"""

    def __init__(self, log_file=None, log_level="INFO"):
        self.log_file = log_file
        self.setup_logging(log_level)

    def setup_logging(self, log_level):
        """设置日志|Setup logging"""
        # 标准日志格式|Standard log format
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


class FastqReader:
    """FASTQ文件读取器|FASTQ File Reader"""

    def __init__(self, fastq_file, logger):
        self.fastq_file = fastq_file
        self.logger = logger

    def open_file(self):
        """打开FASTQ文件(支持gzip)|Open FASTQ file (support gzip)"""
        if self.fastq_file.endswith('.gz'):
            return gzip.open(self.fastq_file, 'rt')
        else:
            return open(self.fastq_file, 'r')

    def read_reads_by_name(self, read_names):
        """
        根据read名称列表读取reads|Read reads by name list

        Args:
            read_names: set或list，包含需要提取的read名称|set or list of read names to extract

        Returns:
            dict: {read_name: (header, sequence, plus, quality)}
        """
        # 总是创建set的副本，避免修改原始集合|Always create a copy to avoid modifying original set
        read_set = set(read_names)
        found_reads = {}
        total_reads = 0

        # 保存原始目标数量|Save original target count
        target_count = len(read_names)

        self.logger.info(f"开始扫描FASTQ文件|Starting to scan FASTQ file: {self.fastq_file}")
        self.logger.info(f"目标reads数量|Target reads count: {target_count:,}")

        with self.open_file() as f:
            while True:
                header = f.readline().strip()
                if not header:
                    break

                sequence = f.readline().strip()
                plus = f.readline().strip()
                quality = f.readline().strip()

                total_reads += 1

                # 提取read名称(去除@和空格后的部分)|Extract read name (remove @ and parts after space)
                read_name = header[1:].split()[0] if header.startswith('@') else header

                if read_name in read_set:
                    found_reads[read_name] = (header, sequence, plus, quality)
                    read_set.remove(read_name)

                    if len(found_reads) % 10000 == 0:
                        self.logger.info(f"已找到|Found: {len(found_reads):,} reads")

                # 如果所有目标reads都找到了，可以提前退出|Early exit if all reads found
                if not read_set:
                    self.logger.info(f"所有目标reads均已找到|All target reads found")
                    break

                if total_reads % 100000 == 0:
                    self.logger.info(f"已扫描|Scanned: {total_reads:,} reads, 已找到|Found: {len(found_reads):,}")

        self.logger.info(f"扫描完成|Scan completed: 总reads数|Total reads: {total_reads:,}")
        self.logger.info(f"找到目标reads|Found target reads: {len(found_reads):,}/{target_count:,}")

        return found_reads


class MappingFileParser:
    """对应关系文件解析器|Mapping File Parser"""

    def __init__(self, mapping_file, logger):
        self.mapping_file = mapping_file
        self.logger = logger

    def parse(self, contig_filter=None):
        """
        解析对应关系文件|Parse mapping file

        Args:
            contig_filter: 可选，只提取特定contig的reads|Optional, only extract reads for specific contig

        Returns:
            dict: {contig_id: [read_names]} 或 {read_name: contig_id}
        """
        contig_to_reads = {}
        read_to_contig = {}
        total_entries = 0
        skipped = 0

        self.logger.info(f"开始解析对应关系文件|Starting to parse mapping file: {self.mapping_file}")

        with open(self.mapping_file, 'r') as f:
            # 跳过header行(以#开头)|Skip header line (starts with #)
            first_line = f.readline().strip()
            if not first_line.startswith('#'):
                # 如果没有header，seek回到文件开头|If no header, seek back to start
                f.seek(0)
            else:
                self.logger.info(f"检测到header行|Detected header line: {first_line}")

            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    skipped += 1
                    continue

                parts = line.split('\t')
                if len(parts) < 2:
                    self.logger.warning(f"行{line_num}格式错误，跳过|Line {line_num} format error, skipped")
                    skipped += 1
                    continue

                contig_id = parts[0].strip()
                read_name = parts[1].strip()

                # 如果指定了contig过滤|If contig filter specified
                if contig_filter and contig_id != contig_filter:
                    continue

                if contig_id not in contig_to_reads:
                    contig_to_reads[contig_id] = []
                contig_to_reads[contig_id].append(read_name)
                read_to_contig[read_name] = contig_id

                total_entries += 1

                if total_entries % 10000 == 0:
                    self.logger.debug(f"已解析|Parsed: {total_entries:,} entries")

        self.logger.info(f"解析完成|Parse completed: 总条目数|Total entries: {total_entries:,}")
        self.logger.info(f"涉及contig数量|Contigs involved: {len(contig_to_reads):,}")

        if skipped > 0:
            self.logger.warning(f"跳过空行|Skipped empty lines: {skipped}")

        return contig_to_reads, read_to_contig
