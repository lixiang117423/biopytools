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

            # 限制每条查询序列的比对数(按相似度取前N,保证所有query都有展示)
            # |Limit alignments per query (top-N by identity, so every query is shown)
            sample_alignments = self._limit_per_query(sample_alignments)

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
        # 合并文件(18列)|Merged file (18 columns): Sample qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen coverage qseq sseq qlen
        # 单样品文件(16列)|Single sample file (16 columns): qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen qseq sseq qlen (无coverage|no coverage)
        min_fields_merged = 18
        min_fields_single = 16

        if len(parts) < min_fields_single:
            self.logger.debug(f"行列数不足|Line has insufficient columns: {len(parts)} < {min_fields_single}")
            return None

        try:
            # 按列数和内容检测格式|Detect format by checking column count and content
            is_merged_format = len(parts) >= min_fields_merged

            if is_merged_format:
                # 解析合并文件格式(18列,含Sample和coverage)|Parse merged file format (18 columns with Sample and coverage)
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
                    'subject_seq': parts[16] if len(parts) > 16 else '',  # sseq
                    'qlen': int(parts[17]) if len(parts) > 17 else 0  # qlen
                }
            else:
                # 解析单样品文件格式(16列,无Sample和coverage)|Parse single sample file format (16 columns without Sample and coverage)
                slen = int(parts[12])
                # 从比对信息计算查询覆盖度|Calculate query coverage from alignment info
                # 分母必须是查询长度qlen:用目标(整条染色体)作分母恒≈0%|Denominator must be
                # query length qlen: subject (whole chromosome) as denominator is always ~0%
                try:
                    qstart = int(parts[6])
                    qend = int(parts[7])
                    qlen_int = int(parts[15])
                    coverage = abs(qend - qstart + 1) / qlen_int * 100 if qlen_int > 0 else 0.0
                    coverage = min(coverage, 100.0)
                except Exception:
                    coverage = 0.0

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
                    'coverage': coverage,  # 计算的查询覆盖度|calculated query coverage
                    'query_seq': parts[13] if len(parts) > 13 else '',  # qseq
                    'subject_seq': parts[14] if len(parts) > 14 else '',  # sseq
                    'qlen': int(parts[15]) if len(parts) > 15 else 0  # qlen
                }

            return alignment

        except (ValueError, IndexError) as e:
            self.logger.debug(f"解析BLAST行失败|Failed to parse BLAST line: {e}")
            return None

    def _limit_per_query(self, alignments: list) -> list:
        """限制每条查询序列展示的比对数(按相似度取前N)|Limit alignments per query (top-N by identity)

        按query首次出现顺序输出;每条query只保留相似度最高的N条,避免高相似query
        占满样品总额度后,排在后面的query完全没有展示|Output keeps first-appearance
        query order; each query keeps only its top-N alignments so late queries
        are not starved by highly-similar early queries
        """
        max_per_query = self.config.alignment_max_per_query
        if not max_per_query or max_per_query < 1:
            return alignments

        groups = {}
        order = []
        for alignment in alignments:
            query_id = alignment['query_id']
            if query_id not in groups:
                groups[query_id] = []
                order.append(query_id)
            groups[query_id].append(alignment)

        limited = []
        dropped = 0
        for query_id in order:
            members = groups[query_id]
            if len(members) > max_per_query:
                members = sorted(members, key=lambda x: x['identity'], reverse=True)[:max_per_query]
                dropped += len(groups[query_id]) - max_per_query
            limited.extend(members)

        if dropped:
            self.logger.info(f"每条查询序列限制{max_per_query}条比对,共裁剪{dropped}条|"
                             f"Per-query limit {max_per_query}, {dropped} alignments dropped")
        return limited

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
