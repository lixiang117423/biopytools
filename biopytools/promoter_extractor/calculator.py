"""
启动子提取器核心计算模块|Promoter Extractor Core Calculation Module
"""

import logging
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Set
from dataclasses import dataclass
from .utils import reverse_complement
from .config import PromoterExtractorConfig


@dataclass
class GeneInfo:
    """基因信息|Gene Information"""
    gene_id: str
    seqid: str
    source: str
    start: int
    end: int
    strand: str
    attributes: Dict[str, str]


class GFF3Parser:
    """GFF3文件解析器|GFF3 File Parser"""

    def __init__(self, logger: logging.Logger):
        """
        初始化GFF3解析器|Initialize GFF3 parser

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger

    def parse_genes(self, gff_file: str, target_genes: Optional[Set[str]] = None) -> Dict[str, GeneInfo]:
        """
        解析GFF3文件中的基因特征|Parse gene features from GFF3 file

        Args:
            gff_file: GFF3文件路径|GFF3 file path
            target_genes: 目标基因ID集合（可选）|Target gene ID set (optional)

        Returns:
            dict: {gene_id: GeneInfo}|Gene information dictionary
        """
        genes = {}

        try:
            with open(gff_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    # 跳过注释行|Skip comment lines
                    if line.startswith('#'):
                        continue

                    line = line.strip()
                    if not line:
                        continue

                    fields = line.split('\t')
                    if len(fields) != 9:
                        self.logger.debug(f"跳过格式错误的行|Skipping malformed line: {line_num}")
                        continue

                    seqid, source, feature, start, end, score, strand, phase, attributes = fields

                    # 只处理gene特征|Only process gene features
                    if feature != 'gene':
                        continue

                    # 解析属性字段|Parse attributes field
                    attr_dict = self._parse_attributes(attributes)
                    gene_id = attr_dict.get('ID') or attr_dict.get('gene_id')

                    if not gene_id:
                        self.logger.warning(f"第{line_num}行|Line {line_num}: 未找到基因ID|Gene ID not found")
                        continue

                    # 如果指定了目标基因列表，跳过不在列表中的基因
                    # If target gene list specified, skip genes not in list
                    if target_genes and gene_id not in target_genes:
                        continue

                    # 解析坐标|Parse coordinates
                    try:
                        start_pos = int(start)
                        end_pos = int(end)
                    except ValueError:
                        self.logger.warning(f"第{line_num}行|Line {line_num}: 坐标格式错误|Invalid coordinate format")
                        continue

                    gene_info = GeneInfo(
                        gene_id=gene_id,
                        seqid=seqid,
                        source=source,
                        start=start_pos,
                        end=end_pos,
                        strand=strand,
                        attributes=attr_dict
                    )

                    genes[gene_id] = gene_info

            self.logger.info(f"从GFF文件中解析出|Parsed from GFF file: {len(genes)} 个基因|genes")

        except Exception as e:
            self.logger.error(f"解析GFF文件失败|Failed to parse GFF file: {e}")
            raise

        return genes

    def _parse_attributes(self, attributes_str: str) -> Dict[str, str]:
        """
        解析GFF3属性字段|Parse GFF3 attributes field

        Args:
            attributes_str: 属性字符串|Attributes string (格式: key1=value1;key2=value2)

        Returns:
            dict: 属性字典|Attributes dictionary
        """
        attr_dict = {}
        for attr in attributes_str.split(';'):
            attr = attr.strip()
            if '=' in attr:
                key, value = attr.split('=', 1)
                attr_dict[key] = value
        return attr_dict


class GenomeSequenceLoader:
    """基因组序列加载器|Genome Sequence Loader"""

    def __init__(self, logger: logging.Logger):
        """
        初始化序列加载器|Initialize sequence loader

        Args:
            logger: 日志对象|Logger object
        """
        self.logger = logger
        self.sequences = {}
        self.seq_lengths = {}

    def load_fasta(self, fasta_file: str) -> bool:
        """
        加载FASTA格式的基因组文件|Load genome file in FASTA format

        Args:
            fasta_file: FASTA文件路径|FASTA file path

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            current_seq_id = None
            current_seq = []

            with open(fasta_file, 'r', encoding='utf-8') as f:
                for line in f:
                    line = line.strip()

                    if line.startswith('>'):
                        # 保存前一个序列|Save previous sequence
                        if current_seq_id is not None:
                            sequence = ''.join(current_seq)
                            self.sequences[current_seq_id] = sequence
                            self.seq_lengths[current_seq_id] = len(sequence)

                        # 开始新序列|Start new sequence
                        current_seq_id = line[1:].split()[0]  # 去掉>和可能的描述
                        current_seq = []
                    else:
                        current_seq.append(line)

                # 保存最后一个序列|Save last sequence
                if current_seq_id is not None:
                    sequence = ''.join(current_seq)
                    self.sequences[current_seq_id] = sequence
                    self.seq_lengths[current_seq_id] = len(sequence)

            self.logger.info(f"加载了|Loaded {len(self.sequences)} 条染色体序列|chromosome sequences")

            return True

        except Exception as e:
            self.logger.error(f"加载FASTA文件失败|Failed to load FASTA file: {e}")
            return False

    def get_sequence(self, seqid: str, start: int, end: int) -> str:
        """
        获取指定区域的序列（坐标从1开始）|Get sequence for specified region (1-based coordinates)

        Args:
            seqid: 序列ID|Sequence ID
            start: 起始位置（包含）|Start position (inclusive, 1-based)
            end: 结束位置（包含）|End position (inclusive, 1-based)

        Returns:
            str: 序列|Sequence
        """
        if seqid not in self.sequences:
            raise KeyError(f"序列ID不存在|Sequence ID not found: {seqid}")

        seq = self.sequences[seqid]
        seq_length = self.seq_lengths[seqid]

        # 转换为0-based坐标|Convert to 0-based coordinates
        start_0based = max(0, start - 1)
        end_0based = min(seq_length, end)

        if start_0based >= end_0based:
            return ""

        return seq[start_0based:end_0based]

    def get_sequence_length(self, seqid: str) -> int:
        """
        获取序列长度|Get sequence length

        Args:
            seqid: 序列ID|Sequence ID

        Returns:
            int: 序列长度|Sequence length
        """
        return self.seq_lengths.get(seqid, 0)


class PromoterCalculator:
    """启动子计算器|Promoter Calculator"""

    def __init__(self, config: PromoterExtractorConfig, logger: logging.Logger):
        """
        初始化启动子计算器|Initialize promoter calculator

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger
        self.stats = {
            'total_genes': 0,
            'extracted': 0,
            'skipped_boundary': 0,
            'skipped_length': 0
        }

    def calculate_promoter_region(self, gene: GeneInfo) -> Optional[Tuple[str, int, int, int]]:
        """
        计算启动子区域|Calculate promoter region

        Args:
            gene: 基因信息|Gene information

        Returns:
            tuple: (seqid, start, end, length) 或 None|or None if invalid
        """
        seqid = gene.seqid
        strand = gene.strand
        promoter_length = self.config.promoter_length

        # 正向链基因|Forward strand gene
        if strand == '+':
            # 启动子在基因上游|Promoter upstream of gene
            promoter_start = gene.start - promoter_length
            promoter_end = gene.start - 1

        # 负向链基因|Reverse strand gene
        elif strand == '-':
            # 启动子在基因下游|Promoter downstream of gene
            promoter_start = gene.end + 1
            promoter_end = gene.end + promoter_length

        else:
            self.logger.warning(f"基因|Gene {gene.gene_id}: 无效的链方向|Invalid strand: {strand}")
            return None

        # 交换确保start < end|Swap to ensure start < end
        if promoter_start > promoter_end:
            promoter_start, promoter_end = promoter_end, promoter_start

        return (seqid, promoter_start, promoter_end)

    def adjust_boundary(self, seqid: str, start: int, end: int,
                       genome_loader: GenomeSequenceLoader) -> Optional[Tuple[int, int, int]]:
        """
        调整边界以适应基因组范围|Adjust boundaries to fit genome range

        Args:
            seqid: 序列ID|Sequence ID
            start: 起始位置|Start position
            end: 结束位置|End position
            genome_loader: 基因组加载器|Genome loader

        Returns:
            tuple: (adjusted_start, adjusted_end, length) 或 None|or None if too short
        """
        seq_length = genome_loader.get_sequence_length(seqid)

        # 调整到边界|Adjust to boundaries
        adjusted_start = max(1, start)
        adjusted_end = min(seq_length, end)

        # 计算实际长度|Calculate actual length
        actual_length = adjusted_end - adjusted_start + 1

        # 检查最小长度要求|Check minimum length requirement
        if actual_length < self.config.min_length:
            return None

        return (adjusted_start, adjusted_end, actual_length)


class PromoterExtractor:
    """启动子提取器主类|Promoter Extractor Main Class"""

    def __init__(self, config: PromoterExtractorConfig, logger: logging.Logger):
        """
        初始化启动子提取器|Initialize promoter extractor

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
        """
        self.config = config
        self.logger = logger
        self.gff_parser = GFF3Parser(logger)
        self.genome_loader = GenomeSequenceLoader(logger)
        self.calculator = PromoterCalculator(config, logger)

    def extract_promoters(self) -> Tuple[List[Tuple[str, str]], List[Tuple[str, str, int, int, str, int]]]:
        """
        提取启动子序列|Extract promoter sequences

        Returns:
            tuple: (fasta_records, bed_records)
                   - fasta_records: [(header, sequence), ...]
                   - bed_records: [(seqid, start, end, name, score, strand), ...]
        """
        # 加载基因组序列|Load genome sequences
        self.logger.info("加载基因组序列|Loading genome sequences...")
        if not self.genome_loader.load_fasta(self.config.genome_file):
            raise RuntimeError("加载基因组文件失败|Failed to load genome file")

        # 解析基因|Parse genes
        self.logger.info("解析GFF文件|Parsing GFF file...")
        target_genes = None
        if self.config.gene_list:
            target_genes = parse_gene_list(self.config.gene_list)
            self.logger.info(f"目标基因数量|Target gene count: {len(target_genes)}")

        genes = self.gff_parser.parse_genes(self.config.gff_file, target_genes)
        self.calculator.stats['total_genes'] = len(genes)

        if not genes:
            self.logger.warning("未找到任何基因|No genes found")
            return [], []

        # 提取启动子|Extract promoters
        self.logger.info("提取启动子序列|Extracting promoter sequences...")
        fasta_records = []
        bed_records = []

        for gene_id, gene in genes.items():
            # 计算启动子区域|Calculate promoter region
            region = self.calculator.calculate_promoter_region(gene)
            if not region:
                continue

            seqid, start, end = region

            # 调整边界|Adjust boundaries
            adjusted = self.calculator.adjust_boundary(
                seqid, start, end, self.genome_loader
            )

            if not adjusted:
                self.calculator.stats['skipped_length'] += 1
                self.logger.debug(f"基因|Gene {gene_id}: 启动子太短|Promoter too short")
                continue

            adj_start, adj_end, length = adjusted

            # 检查是否被截断|Check if truncated
            if length < self.config.promoter_length:
                self.calculator.stats['skipped_boundary'] += 1
                self.logger.debug(f"基因|Gene {gene_id}: 边界截断|Boundary truncated ({length}bp)")

            # 获取序列|Get sequence
            sequence = self.genome_loader.get_sequence(seqid, adj_start, adj_end)

            # 负向链需要反向互补|Reverse strand needs reverse complement
            if gene.strand == '-':
                sequence = reverse_complement(sequence)

            # 生成FASTA记录|Generate FASTA record
            header = f"{gene_id}_promoter"
            fasta_records.append((header, sequence))

            # 生成BED记录|Generate BED record (0-based coordinates)
            bed_start = adj_start - 1  # 转换为0-based|Convert to 0-based
            bed_records.append((seqid, bed_start, adj_end, gene_id, length, gene.strand))

            self.calculator.stats['extracted'] += 1

        self.logger.info(f"成功提取|Successfully extracted: {self.calculator.stats['extracted']} 个启动子|promoters")

        return fasta_records, bed_records


# 导入parse_gene_list函数以便使用|Import parse_gene_list for use
from .utils import parse_gene_list
