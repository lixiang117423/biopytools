"""
染色体重命名核心计算逻辑|Chromosome Rename Core Calculation Logic
"""

import os
from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple


class ChrRenameCalculator:
    """染色体重命名计算器|Chromosome Rename Calculator"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 输出文件路径|Output file paths
        self.paf_file = os.path.join(self.config.output_dir, "alignment.paf")
        self.mapping_file = os.path.join(self.config.output_dir, "chromosome_mapping.tsv")
        self.renamed_fasta = os.path.join(self.config.output_dir, "renamed_genome.fa")
        self.summary_file = os.path.join(self.config.output_dir, "rename_summary.txt")

    def run_minimap2(self) -> bool:
        """运行minimap2比对|Run minimap2 alignment"""
        self.logger.info("开始minimap2比对|Starting minimap2 alignment")

        cmd = (
            f"{self.config.minimap2_path} "
            f"-x {self.config.preset} "
            f"-t {self.config.threads} "
            f"{self.config.ref_fasta} "
            f"{self.config.query_fasta} "
            f"-o {self.paf_file}"
        )

        success = self.cmd_runner.run(cmd, "minimap2比对|minimap2 alignment")

        if success and os.path.exists(self.paf_file):
            self.logger.info(f"PAF文件已生成|PAF file generated: {self.paf_file}")
            return True
        else:
            self.logger.error("PAF文件生成失败|Failed to generate PAF file")
            return False

    def parse_paf(self) -> Dict[str, List[Dict]]:
        """
        解析PAF文件（不应用过滤）|Parse PAF file (no filtering)

        Returns:
            Dict: {query_name: [alignments]}
        """
        self.logger.info("解析PAF文件|Parsing PAF file")

        query_alignments = defaultdict(list)

        if not os.path.exists(self.paf_file):
            self.logger.error(f"PAF文件不存在|PAF file not found: {self.paf_file}")
            return {}

        with open(self.paf_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 12:
                    self.logger.warning(f"行{line_num}格式不正确，跳过|Line {line_num} malformed, skipping")
                    continue

                # PAF格式关键字段|PAF format key fields
                # 字段0: query name|字段1: query length|字段2: query start|字段3: query end
                # 字段5: target name|字段6: target length|字段7: target start|字段8: target end
                # 字段9: 残基匹配数|字段10: 比对长度

                query_name = fields[0]
                query_length = int(fields[1])
                query_start = int(fields[2])
                query_end = int(fields[3])
                target_name = fields[5]
                target_length = int(fields[6])
                target_start = int(fields[7])
                target_end = int(fields[8])
                matches = int(fields[9])
                alignment_length = int(fields[10])

                # 计算覆盖度和一致性|Calculate coverage and identity
                coverage = (query_end - query_start) / query_length
                identity = matches / alignment_length if alignment_length > 0 else 0

                # 不在这里过滤，保留所有比对|Don't filter here, keep all alignments
                alignment_info = {
                    'query_name': query_name,
                    'query_length': query_length,
                    'query_start': query_start,
                    'query_end': query_end,
                    'target_name': target_name,
                    'target_length': target_length,
                    'target_start': target_start,
                    'target_end': target_end,
                    'matches': matches,
                    'alignment_length': alignment_length,
                    'coverage': coverage,
                    'identity': identity
                }

                query_alignments[query_name].append(alignment_info)

        self.logger.info(f"共解析{len(query_alignments)}个query序列|Parsed {len(query_alignments)} query sequences")
        return query_alignments

    def _calculate_merged_coverage(self, alignments: List[Dict]) -> float:
        """
        计算合并后所有比对段的整体覆盖度|Calculate overall coverage after merging all alignment segments

        Args:
            alignments: 比对列表|List of alignments

        Returns:
            float: 整体覆盖度|Overall coverage
        """
        if not alignments:
            return 0.0

        qlen = alignments[0]['query_length']

        # 收集所有比对段|Collect all alignment segments
        segments = [(a['query_start'], a['query_end']) for a in alignments]

        # 合并重叠或连续的区间|Merge overlapping or contiguous intervals
        segments.sort()
        merged = []
        for start, end in segments:
            if merged and start <= merged[-1][1]:
                merged[-1] = (merged[-1][0], max(merged[-1][1], end))
            else:
                merged.append([start, end])

        # 计算总覆盖长度|Calculate total covered length
        total_covered = sum(end - start for start, end in merged)

        return total_covered / qlen if qlen > 0 else 0.0

    def build_mapping_relations(self, query_alignments: Dict[str, List[Dict]]) -> Tuple[Dict[str, str], Dict[str, float]]:
        """
        建立染色体映射关系（优先长序列，确保每个target只映射一次）
        Build chromosome mapping relationships (prioritize long sequences, ensure each target is mapped only once)

        合并每个query的所有比对段，计算整体覆盖度，然后进行映射判断
        Merge all alignment segments for each query, calculate overall coverage, then perform mapping

        优先策略：
        1. 按照序列长度降序处理（染色体 > contig）
        2. 使用单一覆盖度阈值（20%）
        3. 每个target只映射一次

        Args:
            query_alignments: {query_name: [alignments]}

        Returns:
            Tuple[Dict, Dict]: ({query_name: target_name}, {query_name: coverage_used})
        """
        self.logger.info("建立染色体映射关系（优先长序列策略）|Building chromosome mapping relationships (prioritize long sequences)")

        # 使用单一覆盖度阈值|Use single coverage threshold
        min_coverage = 0.2

        mapping = {}
        coverage_used = {}
        used_targets = set()  # 记录已使用的target
        multiple_alignments = []

        # 按序列长度降序排序所有query（染色体优先）|Sort all queries by length descending (chromosomes first)
        sorted_queries = sorted(
            query_alignments.keys(),
            key=lambda q: query_alignments[q][0]['query_length'] if query_alignments[q] else 0,
            reverse=True
        )

        self.logger.info(f"共{len(sorted_queries)}个序列，按长度降序处理|Total {len(sorted_queries)} sequences, processing by length descending")
        self.logger.info(f"覆盖度阈值|Coverage threshold: >= {min_coverage:.0%}")

        mapped_count = 0

        for query_name in sorted_queries:
            alignments = query_alignments[query_name]

            if not alignments:
                continue

            # 按target分组|Group by target
            target_to_alignments = {}
            for aln in alignments:
                target = aln['target_name']
                if target not in target_to_alignments:
                    target_to_alignments[target] = []
                target_to_alignments[target].append(aln)

            # 计算每个target的整体覆盖度|Calculate overall coverage for each target
            target_stats = []
            for target, target_alns in target_to_alignments.items():
                # 跳过已使用的target|Skip already used targets
                if target in used_targets:
                    continue

                merged_cov = self._calculate_merged_coverage(target_alns)
                total_aln_len = sum(a['alignment_length'] for a in target_alns)
                weighted_identity = sum(a['identity'] * a['alignment_length'] for a in target_alns) / total_aln_len if total_aln_len > 0 else 0

                # 只过滤覆盖度和比对长度|Filter by coverage and alignment length only
                if merged_cov >= min_coverage and total_aln_len >= self.config.min_alignment_length:
                    target_stats.append({
                        'target': target,
                        'coverage': merged_cov,
                        'identity': weighted_identity,
                        'total_aln_len': total_aln_len
                    })

            if not target_stats:
                self.logger.debug(f"  {query_name}: 无满足条件的target|No eligible targets")
                continue

            # 按覆盖度排序，选择最佳target|Sort by coverage, select best target
            target_stats.sort(key=lambda x: (x['coverage'], x['total_aln_len']), reverse=True)
            best = target_stats[0]

            mapping[query_name] = best['target']
            coverage_used[query_name] = best['coverage']
            used_targets.add(best['target'])
            mapped_count += 1

            query_len = alignments[0]['query_length']
            self.logger.info(f"  {query_name} ({query_len/1e6:.1f}M bp) -> {best['target']} (覆盖度|coverage: {best['coverage']:.2%})")

            # 如果有多个target满足条件，记录|If multiple targets meet criteria, record
            if len(target_stats) > 1:
                multiple_alignments.append({
                    'query': query_name,
                    'query_length': query_len,
                    'best_target': best['target'],
                    'best_coverage': best['coverage'],
                    'all_targets': [t['target'] for t in target_stats]
                })

        # 输出汇总统计|Output summary statistics
        unmapped_final = [q for q in query_alignments.keys() if q not in mapping]
        self.logger.info(f"成功映射|Successfully mapped: {len(mapping)}/{len(query_alignments)} 个序列|sequences")

        if unmapped_final:
            self.logger.warning(f"警告|WARNING: {len(unmapped_final)}个序列未找到有效比对|sequences have no valid alignments")
            for query in unmapped_final[:5]:  # 只显示前5个|Show only first 5
                self.logger.warning(f"  未映射|Unmapped: {query}")
            if len(unmapped_final) > 5:
                self.logger.warning(f"  ... 还有{len(unmapped_final) - 5}个|and {len(unmapped_final) - 5} more")

        if multiple_alignments:
            self.logger.warning(f"警告|WARNING: {len(multiple_alignments)}个序列有多个比对，已选择最佳|sequences have multiple alignments, best selected")
            for item in multiple_alignments[:5]:  # 只显示前5个|Show only first 5
                self.logger.warning(
                    f"  {item['query']} ({item['query_length']/1e6:.1f}M bp) 可能是嵌合/错误组装|may be chimeric/misassembled\n"
                    f"    最佳比对|Best alignment: {item['best_target']} (覆盖度|coverage: {item['best_coverage']:.2%})\n"
                    f"    所有比对目标|All targets: {', '.join(item['all_targets'])}"
                )
            if len(multiple_alignments) > 5:
                self.logger.warning(f"  ... 还有{len(multiple_alignments) - 5}个|and {len(multiple_alignments) - 5} more")

        return mapping, coverage_used

    def detect_one_to_many(self, mapping: Dict[str, str]) -> Dict[str, List[str]]:
        """
        检测一对多映射关系|Detect one-to-many mapping relationships

        Args:
            mapping: {query_name: target_name}

        Returns:
            Dict: {target_name: [query_names]}
        """
        self.logger.info("检测一对多映射关系|Detecting one-to-many mapping relationships")

        target_to_queries = defaultdict(list)
        for query, target in mapping.items():
            target_to_queries[target].append(query)

        # 只保留有多于1个query的target|Only keep targets with more than 1 query
        one_to_many = {t: qs for t, qs in target_to_queries.items() if len(qs) > 1}

        if one_to_many:
            self.logger.info(f"发现{len(one_to_many)}个一对多映射关系|Found {len(one_to_many)} one-to-many mappings")
            for target, queries in one_to_many.items():
                self.logger.info(
                    f"  多个group映射到|Multiple groups mapped to {target}:\n"
                    f"    {', '.join(queries)}"
                )

        return one_to_many

    def write_mapping_table(self, mapping: Dict[str, str], coverage_used: Dict[str, float]):
        """写入映射表文件|Write mapping table file"""
        self.logger.info(f"写入映射表|Writing mapping table: {self.mapping_file}")

        with open(self.mapping_file, 'w') as f:
            f.write("Query_Name\tTarget_Name\tQuery_Length\tTarget_Length\tCoverage\tIdentity\n")
            f.write("Original_Chromosome\tRenamed_Chromosome\tDetails\n")
            f.write("-" * 80 + "\n")

            # 重新解析PAF以获取详细信息|Reparse PAF to get details
            query_alignments = self.parse_paf()

            for query_name in sorted(mapping.keys()):
                target_name = mapping[query_name]

                # 获取比对详情|Get alignment details
                alignments = query_alignments.get(query_name, [])
                if alignments:
                    # 过滤出目标比对|Filter alignments for the target
                    target_alignments = [a for a in alignments if a['target_name'] == target_name]

                    if target_alignments:
                        # 计算统计信息|Calculate statistics
                        qlen = target_alignments[0]['query_length']
                        tlen = target_alignments[0]['target_length']
                        avg_identity = sum(a['identity'] for a in target_alignments) / len(target_alignments)

                        # 使用合并后的覆盖度|Use merged coverage
                        merged_cov = self._calculate_merged_coverage(target_alignments)

                        f.write(
                            f"{query_name}\t{target_name}\t"
                            f"{qlen}\t{tlen}\t"
                            f"{merged_cov:.2%}\t{avg_identity:.2%}\n"
                        )
                    else:
                        f.write(f"{query_name}\t{target_name}\t-\t-\t-\t-\n")
                else:
                    f.write(f"{query_name}\t{target_name}\t-\t-\t-\t-\n")

        self.logger.info("映射表写入完成|Mapping table written successfully")

    def write_summary_report(self, mapping: Dict[str, str], one_to_many: Dict[str, List[str]]):
        """写入汇总报告|Write summary report"""
        self.logger.info(f"写入汇总报告|Writing summary report: {self.summary_file}")

        with open(self.summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("染色体重命名汇总报告|Chromosome Rename Summary Report\n")
            f.write("=" * 80 + "\n\n")

            # 基本信息|Basic information
            f.write("1. 基本信息|Basic Information\n")
            f.write("-" * 80 + "\n")
            f.write(f"参考基因组|Reference genome: {self.config.ref_fasta}\n")
            f.write(f"查询基因组|Query genome: {self.config.query_fasta}\n")
            f.write(f"比对模式|Alignment preset: {self.config.preset}\n")
            f.write(f"最小一致性|Min identity: {self.config.min_identity}\n")
            f.write(f"映射策略|Mapping strategy: 长序列优先|Prioritize long sequences\n")
            f.write(f"  - 按query序列长度降序处理|Process queries by length descending\n")
            f.write(f"  - 每个target只映射一次|Each target mapped only once\n")
            f.write(f"  - 覆盖度阈值|Coverage threshold: >= 20%\n\n")

            # 统计信息|Statistics
            f.write("2. 统计信息|Statistics\n")
            f.write("-" * 80 + "\n")
            f.write(f"成功映射|Successfully mapped: {len(mapping)} 个序列|sequences\n")
            f.write(f"一对多映射|One-to-many mappings: {len(one_to_many)} 个目标|targets\n\n")

            # 一对多映射详情|One-to-many mapping details
            if one_to_many:
                f.write("3. 一对多映射详情|One-to-many Mapping Details\n")
                f.write("-" * 80 + "\n")
                f.write("(可能表示基因组碎片化组装|May indicate fragmented genome assembly)\n\n")
                for target, queries in sorted(one_to_many.items()):
                    f.write(f"{target}:\n")
                    for query in queries:
                        f.write(f"  - {query}\n")
                    f.write("\n")

            f.write("=" * 80 + "\n")
            f.write("报告生成完成|Report generation completed\n")
            f.write("=" * 80 + "\n")

        self.logger.info("汇总报告写入完成|Summary report written successfully")

    def rename_fasta(self, mapping: Dict[str, str]) -> bool:
        """
        重命名FASTA文件中的染色体名，并按参考基因组顺序排列|Rename chromosome names in FASTA file and order by reference genome

        Args:
            mapping: {old_name: new_name}

        Returns:
            bool: 是否成功|Success status
        """
        self.logger.info("开始重命名FASTA序列|Starting FASTA sequence renaming")
        self.logger.info(f"输出文件|Output file: {self.renamed_fasta}")

        if not os.path.exists(self.config.query_fasta):
            self.logger.error(f"查询基因组文件不存在|Query genome file not found: {self.config.query_fasta}")
            return False

        # 读取参考基因组，获取染色体顺序|Read reference genome to get chromosome order
        ref_chromosomes = self._get_ref_chromosome_order()
        self.logger.info(f"参考基因组染色体顺序|Reference chromosome order: {', '.join(ref_chromosomes)}")

        # 读取query基因组所有序列|Read all query genome sequences
        query_sequences = self._read_fasta_sequences(self.config.query_fasta)

        # 构建反向映射: target_name -> [query_names]|Build reverse mapping: target_name -> [query_names]
        target_to_queries = defaultdict(list)
        unmapped_queries = []

        for query_name in query_sequences.keys():
            if query_name in mapping:
                target_name = mapping[query_name]
                target_to_queries[target_name].append(query_name)
            else:
                unmapped_queries.append(query_name)

        # 按照参考基因组顺序输出|Output in reference genome order
        renamed_count = 0
        unmapped_count = 0

        with open(self.renamed_fasta, 'w') as fout:
            # 先输出已映射的序列，按参考基因组顺序|First output mapped sequences in reference order
            for ref_chr in ref_chromosomes:
                if ref_chr in target_to_queries:
                    for query_name in target_to_queries[ref_chr]:
                        self._write_sequence(fout, query_name, query_sequences[query_name], new_name=ref_chr)
                        renamed_count += 1
                        self.logger.debug(f"{query_name} -> {ref_chr}")

            # 最后输出未映射的序列，按长度降序命名为Contig_01, Contig_02...
            # Finally output unmapped sequences, named as Contig_01, Contig_02... by length descending
            if unmapped_queries:
                self.logger.info(f"开始处理{len(unmapped_queries)}个未映射序列|Processing {len(unmapped_queries)} unmapped sequences")

                # 按序列长度降序排序|Sort by sequence length descending
                unmapped_queries.sort(
                    key=lambda q: len(query_sequences[q][1]),  # sequence length
                    reverse=True
                )

                for idx, query_name in enumerate(unmapped_queries, 1):
                    contig_name = f"Contig_{idx:02d}"
                    self._write_sequence(fout, query_name, query_sequences[query_name], new_name=contig_name)
                    unmapped_count += 1
                    seq_len = len(query_sequences[query_name][1])
                    self.logger.info(f"  {query_name} ({seq_len/1e6:.2f}M bp) -> {contig_name}")

        self.logger.info(f"重命名完成|Renaming completed")
        self.logger.info(f"  重命名为染色体|Renamed as chromosomes: {renamed_count} 个序列|sequences")
        if unmapped_count > 0:
            self.logger.info(f"  重命名为Contig|Renamed as Contig: {unmapped_count} 个序列|sequences")

        return True

    def _get_ref_chromosome_order(self) -> List[str]:
        """获取参考基因组的染色体顺序|Get chromosome order from reference genome"""
        chromosomes = []
        with open(self.config.ref_fasta, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    chr_name = line[1:].split()[0]  # 去除>和可能的描述
                    chromosomes.append(chr_name)
        return chromosomes

    def _read_fasta_sequences(self, fasta_file: str) -> Dict[str, Tuple[str, str]]:
        """
        读取FASTA文件所有序列|Read all sequences from FASTA file

        Returns:
            Dict: {sequence_name: (header, sequence)}
        """
        sequences = {}
        current_name = None
        current_header = None
        current_seq = []

        with open(fasta_file, 'r') as f:
            for line in f:
                line = line.rstrip('\n\r')
                if line.startswith('>'):
                    # 保存前一个序列|Save previous sequence
                    if current_name is not None:
                        sequences[current_name] = (current_header, ''.join(current_seq))

                    # 开始新序列|Start new sequence
                    current_header = line
                    current_name = line[1:].split()[0]
                    current_seq = []
                else:
                    if current_name is not None:
                        current_seq.append(line)

            # 保存最后一个序列|Save last sequence
            if current_name is not None:
                sequences[current_name] = (current_header, ''.join(current_seq))

        return sequences

    def _write_sequence(self, fout, query_name: str, sequence_data: Tuple[str, str], new_name: str = None):
        """
        写入序列到文件|Write sequence to file

        Args:
            fout: 输出文件句柄
            query_name: 原始序列名
            sequence_data: (header, sequence)
            new_name: 新名称（如果为None则使用原始名称）
        """
        header, sequence = sequence_data

        if new_name is not None:
            # 使用新名称|Use new name
            fout.write(f">{new_name}\n")
        else:
            # 使用原始名称|Use original name
            fout.write(f"{header}\n")

        # 写入序列，每行80个字符|Write sequence, 80 chars per line
        for i in range(0, len(sequence), 80):
            fout.write(sequence[i:i+80] + '\n')

    def run_analysis(self) -> bool:
        """运行完整分析流程|Run complete analysis pipeline"""
        self.logger.info("开始染色体重命名分析|Starting chromosome rename analysis")

        # 1. 运行minimap2|Run minimap2
        if not self.run_minimap2():
            return False

        # 2. 解析PAF文件|Parse PAF file
        query_alignments = self.parse_paf()
        if not query_alignments:
            self.logger.error("PAF解析失败或无有效比对|PAF parsing failed or no valid alignments")
            return False

        # 3. 建立映射关系（迭代策略）|Build mapping relationships (iterative strategy)
        mapping, coverage_used = self.build_mapping_relations(query_alignments)
        if not mapping:
            self.logger.error("未找到任何有效映射|No valid mappings found")
            return False

        # 4. 检测一对多映射|Detect one-to-many mappings
        one_to_many = self.detect_one_to_many(mapping)

        # 5. 写入映射表|Write mapping table
        self.write_mapping_table(mapping, coverage_used)

        # 6. 写入汇总报告|Write summary report
        self.write_summary_report(mapping, one_to_many)

        # 7. 重命名FASTA文件|Rename FASTA file
        if not self.rename_fasta(mapping):
            return False

        self.logger.info("染色体重命名分析完成|Chromosome rename analysis completed successfully")
        return True
