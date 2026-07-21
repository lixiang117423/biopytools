"""
比对可视化生成器|Alignment Visualization Generator
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple
from .text_alignment import TextAlignmentGenerator
from .html_alignment import HTMLAlignmentGenerator

class AlignmentVisualizer:
    """比对可视化生成器(文本和HTML格式)|Alignment visualization generator (text and HTML)"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.text_generator = TextAlignmentGenerator(config, logger)
        self.html_generator = HTMLAlignmentGenerator(config, logger)

    def generate_visualizations(self, blast_results: List[Tuple[str, str, str]]):
        """
        从BLAST结果生成比对可视化|Generate alignment visualizations from BLAST results

        Args:
            blast_results: BLAST结果列表 [(file_name, sample_name, result_file), ...]
                           |BLAST results list [(file_name, sample_name, result_file), ...]
        """
        if self.config.alignment_output == 'none':
            self.logger.info("跳过比对可视化生成|Skipping alignment visualization generation")
            return None

        self.logger.info("=" * 80)
        self.logger.info("生成比对可视化|Generating alignment visualizations")
        self.logger.info("=" * 80)

        # 解析BLAST结果|Parse BLAST results
        alignments_data = self._parse_blast_results(blast_results)

        if not alignments_data:
            self.logger.warning("没有可用于可视化的比对数据|No alignments data available for visualization")
            return None

        # 生成可视化|Generate visualizations
        output_files = {}

        if self.config.alignment_output in ['text', 'both']:
            text_files, text_summary = self.text_generator.generate_alignments(alignments_data)
            output_files['text'] = {'sample_files': text_files, 'summary': text_summary}

        if self.config.alignment_output in ['html', 'both']:
            html_files, html_index = self.html_generator.generate_alignments(alignments_data)
            output_files['html'] = {'sample_files': html_files, 'index': html_index}

        self.logger.info("=" * 80)
        self.logger.info("比对可视化生成完成|Alignment visualization generation completed")
        self.logger.info("=" * 80)

        return output_files

    def _parse_blast_results(self, blast_results: List[Tuple[str, str, str]]) -> Dict:
        """
        从文件解析BLAST结果|Parse BLAST results from files

        Returns:
            Dict: {sample_name: {'file_name': str, 'alignments': [...]}}
        """
        self.logger.info("解析BLAST结果|Parsing BLAST results...")

        alignments_data = {}
        total_parsed = 0
        total_filtered = 0

        for file_name, sample_name, result_file in blast_results:
            if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
                continue

            sample_alignments = []

            try:
                with open(result_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue

                        alignment = self._parse_blast_line(line)
                        if alignment and self._passes_filters(alignment):
                            sample_alignments.append(alignment)
                            total_parsed += 1
                        else:
                            total_filtered += 1

            except Exception as e:
                self.logger.warning(f"读取文件失败|Error reading file {result_file}: {e}")
                continue

            # 限制每个样品的比对数|Limit alignments per sample
            if len(sample_alignments) > self.config.alignment_max_per_sample:
                self.logger.info(f"样品 {sample_name}: 限制为 {self.config.alignment_max_per_sample} 条比对|Sample {sample_name}: limiting to {self.config.alignment_max_per_sample} alignments")
                sample_alignments = sample_alignments[:self.config.alignment_max_per_sample]

            if sample_alignments:
                alignments_data[sample_name] = {
                    'file_name': file_name,
                    'alignments': sample_alignments
                }

        self.logger.info("解析统计|Parsing statistics:")
        self.logger.info(f"  已解析比对数|Total alignments parsed: {total_parsed}")
        self.logger.info(f"  已过滤比对数|Total alignments filtered: {total_filtered}")
        self.logger.info(f"  有比对的样品数|Total samples with alignments: {len(alignments_data)}")

        return alignments_data

    def _parse_blast_line(self, line: str) -> Dict:
        """解析BLAST输出行|Parse BLAST output line"""
        parts = line.split('\t')

        # BLAST输出格式|BLAST output formats:
        # 合并文件(17列)|Merged file (17 columns): Sample qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen coverage qseq sseq
        # 单样品文件(15列)|Single sample file (15 columns): qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qseq sseq (无coverage|no coverage)
        min_fields_merged = 17
        min_fields_single = 15

        if len(parts) < min_fields_single:
            self.logger.debug(f"行列数不足|Line has insufficient columns: {len(parts)} < {min_fields_single}")
            return None

        try:
            # 按列数和内容检测格式|Detect format by checking column count and content
            is_merged_format = len(parts) >= min_fields_merged

            if is_merged_format:
                # 解析合并文件格式(17列,含Sample和coverage)|Parse merged file format (17 columns with Sample and coverage)
                alignment = {
                    'query_id': parts[1],  # qseqid
                    'subject_id': parts[2],  # sseqid
                    'identity': float(parts[3]),  # pident
                    'length': int(parts[4]),  # length
                    'mismatch': int(parts[5]),  # mismatch
                    'gapopen': int(parts[6]),  # gapopen
                    'qstart': int(parts[7]),  # qstart
                    'qend': int(parts[8]),  # qend
                    'sstart': int(parts[9]),  # sstart
                    'send': int(parts[10]),  # send
                    'evalue': parts[11],  # evalue
                    'bitscore': float(parts[12]),  # bitscore
                    'slen': int(parts[13]),  # slen
                    'coverage': float(parts[14]),  # coverage
                    'query_seq': parts[15] if len(parts) > 15 else '',  # qseq
                    'subject_seq': parts[16] if len(parts) > 16 else ''  # sseq
                }
            else:
                # 解析单样品文件格式(15列,无Sample和coverage)|Parse single sample file format (15 columns without Sample and coverage)
                slen = int(parts[12])
                # 从比对信息计算覆盖度|Calculate coverage from alignment info
                try:
                    sstart = int(parts[8])
                    send = int(parts[9])
                    coverage = abs(send - sstart + 1) / slen * 100 if slen > 0 else 0.0
                    coverage = min(coverage, 100.0)
                except Exception:
                    coverage = 100.0

                alignment = {
                    'query_id': parts[0],  # qseqid
                    'subject_id': parts[1],  # sseqid
                    'identity': float(parts[2]),  # pident
                    'length': int(parts[3]),  # length
                    'mismatch': int(parts[4]),  # mismatch
                    'gapopen': int(parts[5]),  # gapopen
                    'qstart': int(parts[6]),  # qstart
                    'qend': int(parts[7]),  # qend
                    'sstart': int(parts[8]),  # sstart
                    'send': int(parts[9]),  # send
                    'evalue': parts[10],  # evalue
                    'bitscore': float(parts[11]),  # bitscore
                    'slen': slen,  # slen
                    'coverage': coverage,  # 计算的覆盖度|calculated coverage
                    'query_seq': parts[13] if len(parts) > 13 else '',  # qseq
                    'subject_seq': parts[14] if len(parts) > 14 else ''  # sseq
                }

            return alignment

        except (ValueError, IndexError) as e:
            self.logger.debug(f"解析BLAST行失败|Failed to parse BLAST line: {e}")
            return None

    def _passes_filters(self, alignment: Dict) -> bool:
        """检查比对是否通过过滤条件|Check if alignment passes the configured filters"""
        # 检查相似度过滤|Check identity filter
        if self.config.alignment_min_identity > 0:
            if alignment['identity'] < self.config.alignment_min_identity:
                return False

        # 检查覆盖度过滤|Check coverage filter
        if self.config.alignment_min_coverage > 0:
            if alignment['coverage'] < self.config.alignment_min_coverage:
                return False

        return True
