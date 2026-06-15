"""
FASTQ Reads提取核心逻辑|FASTQ Reads Extraction Core Logic
"""

import gzip
import os


class FastqReadsExtractor:
    """FASTQ Reads提取器|FASTQ Reads Extractor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def extract(self):
        """执行提取操作|Execute extraction"""
        from .utils import MappingFileParser, FastqReader

        self.logger.info("开始Reads提取流程|Starting reads extraction pipeline")

        # 解析对应关系文件|Parse mapping file
        parser = MappingFileParser(self.config.mapping_file, self.logger)
        contig_to_reads, read_to_contig = parser.parse()

        if not contig_to_reads:
            self.logger.error("未找到任何有效的对应关系|No valid mapping entries found")
            return False

        # 获取所有需要提取的read名称|Get all read names to extract
        all_read_names = set(read_to_contig.keys())
        self.logger.info(f"需要提取的reads总数|Total reads to extract: {len(all_read_names):,}")

        # 从FASTQ文件中读取目标reads|Read target reads from FASTQ file
        reader = FastqReader(self.config.fastq_file, self.logger)
        found_reads = reader.read_reads_by_name(all_read_names)

        if not found_reads:
            self.logger.error("未找到任何目标reads|No target reads found in FASTQ file")
            return False

        # 写入单个输出文件|Write to single output file
        self._write_fastq(found_reads.values(), self.config.output_file)

        self.logger.info(f"Reads提取完成|Reads extraction completed")
        self.logger.info(f"输出文件|Output file: {self.config.output_file}")
        self.logger.info(f"提取reads数|Extracted reads: {len(found_reads):,}")

        return True

    def _write_fastq(self, reads, output_path):
        """
        写入FASTQ文件|Write FASTQ file

        Args:
            reads: iterator or list of (header, sequence, plus, quality)
            output_path: 输出文件路径|Output file path
        """
        # 确保输出目录存在|Ensure output directory exists
        output_dir = os.path.dirname(output_path)
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)

        if output_path.endswith('.gz'):
            open_func = gzip.open
            mode = 'wt'
        else:
            open_func = open
            mode = 'w'

        with open_func(output_path, mode) as f:
            for header, sequence, plus, quality in reads:
                f.write(f"{header}\n{sequence}\n{plus}\n{quality}\n")
