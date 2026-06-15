"""
GFF3/BED导出模块|GFF3/BED Export Module
将PAF比对结果导出为GFF3或BED格式|Export PAF alignments to GFF3 or BED format
"""

import os
from typing import List
from .parser import PAFRecord


class GFF3Exporter:
    """GFF3导出器|GFF3 Exporter"""

    def __init__(self, config, logger):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def export_gff3(self, records: List[PAFRecord], output_file: str = None):
        """导出为GFF3格式|Export to GFF3 format

        Args:
            records: PAF记录列表|List of PAF records
            output_file: 输出文件路径|Output file path
        """
        if output_file is None:
            output_file = os.path.join(self.config.output_dir, "alignment.gff3")

        self.logger.info(f"导出GFF3格式|Exporting GFF3 format: {output_file}")

        with open(output_file, 'w') as f:
            # GFF3头部|GFF3 header
            f.write("##gff-version 3\n")
            f.write("##feature-ontology https://github.com/The-Sequence-Ontology/SO-Ontologies/blob/v3.1/Features_and_Properties.md\n")
            f.write("# seqid source type start end score strand phase attributes\n")

            # 写入比对记录|Write alignment records
            for record in records:
                line = self._format_gff3_line(record)
                f.write(line + '\n')

        self.logger.info(f"GFF3文件已保存|GFF3 file saved: {output_file}")

    def _format_gff3_line(self, record: PAFRecord) -> str:
        """格式化单条GFF3记录|Format single GFF3 record

        Args:
            record: PAF记录|PAF record

        Returns:
            str: GFF3格式行|GFF3 format line
        """
        seqid = record.target_name
        source = "miniprot"
        type_ = "protein_match"  # SO:0000753
        start = record.target_start + 1  # GFF3使用1-based坐标|GFF3 uses 1-based coordinates
        end = record.target_end
        score = record.mapping_quality
        strand = record.strand
        phase = "."

        # 属性字段|Attributes field
        attributes = [
            f"ID={record.query_name}",
            f"Name={record.query_name}",
            f"Target={record.query_name} {record.query_start + 1} {record.query_end}",
            f"Identity={record.identity:.2f}",
            f"Query_coverage={record.query_coverage:.2f}",
            f"Target_coverage={record.target_coverage:.2f}"
        ]

        # 添加可选tags|Add optional tags
        if 'AS' in record.tags:
            attributes.append(f"Alignment_score={record.tags['AS']}")
        if 'cm' in record.tags:
            attributes.append(f"Codons={record.tags['cm']}")

        attributes_str = ";".join(attributes)

        fields = [str(x) for x in [seqid, source, type_, start, end, score, strand, phase, attributes_str]]
        return "\t".join(fields)


class BEDExporter:
    """BED导出器|BED Exporter"""

    def __init__(self, config, logger):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def export_bed(self, records: List[PAFRecord], output_file: str = None):
        """导出为BED格式|Export to BED format

        Args:
            records: PAF记录列表|List of PAF records
            output_file: 输出文件路径|Output file path

        Returns:
            str: BED文件路径|BED file path
        """
        if output_file is None:
            output_file = os.path.join(self.config.output_dir, "alignment.bed")

        self.logger.info(f"导出BED格式|Exporting BED format: {output_file}")

        with open(output_file, 'w') as f:
            # BED头部（可选）|BED header (optional)
            # f.write("# chrom chromStart chromEnd name score strand\n")

            # 写入比对记录|Write alignment records
            for record in records:
                line = self._format_bed_line(record)
                f.write(line + '\n')

        self.logger.info(f"BED文件已保存|BED file saved: {output_file}")

        return output_file

    def _format_bed_line(self, record: PAFRecord) -> str:
        """格式化单条BED记录|Format single BED record

        Args:
            record: PAF记录|PAF record

        Returns:
            str: BED格式行|BED format line
        """
        chrom = record.target_name
        chrom_start = record.target_start  # BED使用0-based坐标|BED uses 0-based coordinates
        chrom_end = record.target_end
        name = record.query_name
        score = record.mapping_quality
        strand = record.strand

        # BED6格式|BED6 format
        fields = [chrom, str(chrom_start), str(chrom_end), name, str(score), strand]
        return "\t".join(fields)
