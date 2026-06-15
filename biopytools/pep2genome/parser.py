"""
PAF文件解析模块|PAF File Parser Module
解析miniprot输出的PAF格式文件|Parse PAF format files output by miniprot
"""

import os
from dataclasses import dataclass
from typing import List, Dict, Optional
from collections import defaultdict


@dataclass
class PAFRecord:
    """PAF记录|PAF Record

    PAF格式（前12列必选）|PAF format (first 12 columns required):
    1. query_name: 查询序列名称|Query sequence name
    2. query_length: 查询序列长度|Query sequence length
    3. query_start: 查询序列起始位置|Query start position
    4. query_end: 查询序列终止位置|Query end position
    5. strand: 链方向|Strand (+ or -)
    6. target_name: 目标序列名称|Target sequence name
    7. target_length: 目标序列长度|Target sequence length
    8. target_start: 目标序列起始位置|Target start position
    9. target_end: 目标序列终止位置|Target end position
    10. residue_matches: 匹配残基数|Number of residue matches
    11. alignment_block_length: 对齐块长度|Alignment block length
    12. mapping_quality: 映射质量|Mapping quality

    可选Tags（以tab分隔）|Optional tags (tab-separated):
    - AS:i: Alignment score
    - ms:i: Number of matching bases
    - nn:i: Number of N bases
    - tp:A: Alignment type (primary/secondary)
    - cm:i: Number of codons
    - s1:i: Number of splices
    - ...
    """
    query_name: str
    query_length: int
    query_start: int
    query_end: int
    strand: str
    target_name: str
    target_length: int
    target_start: int
    target_end: int
    residue_matches: int
    alignment_block_length: int
    mapping_quality: int

    # 可选字段（从tags解析）|Optional fields (parsed from tags)
    tags: Dict[str, any] = None

    # 原始行|Original line
    raw_line: str = None

    def __post_init__(self):
        if self.tags is None:
            self.tags = {}

    @property
    def identity(self) -> float:
        """计算序列一致性|Calculate sequence identity

        Returns:
            float: 一致性百分比|Identity percentage
        """
        if self.alignment_block_length == 0:
            return 0.0
        return (self.residue_matches / self.alignment_block_length) * 100

    @property
    def query_coverage(self) -> float:
        """计算查询序列覆盖度|Calculate query coverage

        Returns:
            float: 覆盖度百分比|Coverage percentage
        """
        if self.query_length == 0:
            return 0.0
        return ((self.query_end - self.query_start) / self.query_length) * 100

    @property
    def target_coverage(self) -> float:
        """计算目标序列覆盖度|Calculate target coverage

        Returns:
            float: 覆盖度百分比|Coverage percentage
        """
        if self.target_length == 0:
            return 0.0
        return ((self.target_end - self.target_start) / self.target_length) * 100


class PAFParser:
    """PAF文件解析器|PAF File Parser"""

    def __init__(self, logger):
        """初始化|Initialize

        Args:
            logger: 日志器|Logger
        """
        self.logger = logger

    def parse_file(self, paf_file: str) -> List[PAFRecord]:
        """解析PAF文件|Parse PAF file

        Args:
            paf_file: PAF文件路径|PAF file path

        Returns:
            List[PAFRecord]: PAF记录列表|List of PAF records
        """
        self.logger.info(f"解析PAF文件|Parsing PAF file: {paf_file}")

        if not os.path.exists(paf_file):
            raise FileNotFoundError(f"PAF文件不存在|PAF file not found: {paf_file}")

        records = []
        total_lines = 0

        with open(paf_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                total_lines += 1
                try:
                    record = self._parse_line(line)
                    records.append(record)
                except Exception as e:
                    self.logger.warning(f"解析行失败|Failed to parse line: {line[:100]}... | {e}")

        self.logger.info(f"共解析 {total_lines} 行，成功 {len(records)} 条记录|"
                        f"Total {total_lines} lines, successfully parsed {len(records)} records")

        return records

    def _parse_line(self, line: str) -> PAFRecord:
        """解析单行PAF记录|Parse single PAF line

        Args:
            line: PAF行|PAF line

        Returns:
            PAFRecord: PAF记录对象|PAF record object
        """
        fields = line.split('\t')

        if len(fields) < 12:
            raise ValueError(f"PAF行必须至少有12列|PAF line must have at least 12 columns: {len(fields)}")

        # 解析前12列必选字段|Parse first 12 required columns
        record = PAFRecord(
            query_name=fields[0],
            query_length=int(fields[1]),
            query_start=int(fields[2]),
            query_end=int(fields[3]),
            strand=fields[4],
            target_name=fields[5],
            target_length=int(fields[6]),
            target_start=int(fields[7]),
            target_end=int(fields[8]),
            residue_matches=int(fields[9]),
            alignment_block_length=int(fields[10]),
            mapping_quality=int(fields[11]),
            tags={},
            raw_line=line
        )

        # 解析可选tags|Parse optional tags
        for i in range(12, len(fields)):
            tag_field = fields[i]
            if ':' in tag_field:
                parts = tag_field.split(':')
                if len(parts) >= 3:
                    tag_name = parts[0]
                    tag_type = parts[1]
                    tag_value = parts[2]

                    # 根据类型转换值|Convert value based on type
                    if tag_type == 'i':
                        record.tags[tag_name] = int(tag_value)
                    elif tag_type == 'f':
                        record.tags[tag_name] = float(tag_value)
                    elif tag_type == 'A':
                        record.tags[tag_name] = tag_value
                    elif tag_type == 'Z':
                        record.tags[tag_name] = tag_value

        return record

    def group_by_query(self, records: List[PAFRecord]) -> Dict[str, List[PAFRecord]]:
        """按查询序列分组|Group records by query sequence

        Args:
            records: PAF记录列表|List of PAF records

        Returns:
            Dict[str, List[PAFRecord]]: 按查询序列分组的记录|Records grouped by query
        """
        grouped = defaultdict(list)
        for record in records:
            grouped[record.query_name].append(record)
        return dict(grouped)

    def group_by_target(self, records: List[PAFRecord]) -> Dict[str, List[PAFRecord]]:
        """按目标序列分组|Group records by target sequence

        Args:
            records: PAF记录列表|List of PAF records

        Returns:
            Dict[str, List[PAFRecord]]: 按目标序列分组的记录|Records grouped by target
        """
        grouped = defaultdict(list)
        for record in records:
            grouped[record.target_name].append(record)
        return dict(grouped)

    def get_best_alignment(self, records: List[PAFRecord]) -> Optional[PAFRecord]:
        """获取最佳比对（按mapping quality和identity）|Get best alignment

        Args:
            records: PAF记录列表|List of PAF records

        Returns:
            Optional[PAFRecord]: 最佳比对记录|Best alignment record
        """
        if not records:
            return None

        # 按mapping quality和identity排序|Sort by mapping quality and identity
        sorted_records = sorted(
            records,
            key=lambda r: (r.mapping_quality, r.identity),
            reverse=True
        )

        return sorted_records[0]
