"""
VCF2Gene核心计算逻辑模块|VCF2Gene Core Calculation Logic Module
"""

import gzip
import bisect
from collections import defaultdict


class GFFParser:
    """GFF文件解析器|GFF File Parser"""

    # 特征类型优先级（数字越小优先级越高，越底层）
    FEATURE_PRIORITY = {
        'exon': 1,
        'cds': 2,
        'five_prime_utr': 3,
        'three_prime_utr': 4,
        'utr': 5,
        'intron': 6,
        'mrna': 7,
        'transcript': 8,
        'gene': 9,
    }

    def __init__(self, logger):
        self.logger = logger
        # 按染色体存储特征区间（仅存储start坐标，用于二分查找）
        self.feature_starts = defaultdict(list)  # {chr: [start1, start2, ...]}
        self.feature_data = defaultdict(list)    # {chr: [(end, feature_type, gene_id), ...]}
        # 按染色体存储基因信息（也使用区间索引）
        self.gene_starts = defaultdict(list)     # {chr: [start1, start2, ...]}
        self.gene_data = defaultdict(list)       # {chr: [(end, gene_id), ...]}
        # 特征ID到基因ID的映射
        self.feature_to_gene = {}

    def parse(self, gff_file):
        """解析GFF文件|Parse GFF file"""
        self.logger.info(f"开始解析GFF文件|Starting to parse GFF file: {gff_file}")

        # 判断是否为压缩文件|Check if file is gzipped
        open_func = gzip.open if gff_file.endswith('.gz') else open

        # 临时存储所有特征，用于后续追溯基因ID
        temp_features = []

        with open_func(gff_file, 'rt') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()

                # 跳过注释行和空行|Skip comment lines and empty lines
                if not line or line.startswith('#'):
                    continue

                fields = line.split('\t')
                if len(fields) < 9:
                    continue

                chrom = fields[0]
                feature_type = fields[2].lower()
                start = int(fields[3])
                end = int(fields[4])
                attributes = fields[8]

                # 解析属性字段|Parse attributes
                attrs_dict = self._parse_attributes(attributes)

                # 提取ID和Parent|Extract ID and Parent
                feature_id = attrs_dict.get('ID')
                parent_ids = attrs_dict.get('Parent', '')

                # 提取gene_id（从多种可能的属性中）
                gene_id = (attrs_dict.get('gene_id') or
                          attrs_dict.get('Gene') or
                          attrs_dict.get('gene'))

                # 如果是gene特征，记录基因信息
                if feature_type == 'gene' and feature_id:
                    # 使用区间索引存储基因（而不是字典）
                    self.gene_starts[chrom].append(start)
                    self.gene_data[chrom].append((end, feature_id))
                    # 基因的gene_id就是自己的ID
                    gene_id = feature_id
                    # 记录基因ID映射
                    self.feature_to_gene[feature_id] = feature_id

                # 临时存储特征信息
                temp_features.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'feature_type': feature_type,
                    'feature_id': feature_id,
                    'parent_ids': parent_ids,
                    'gene_id': gene_id
                })

        # 第二遍：建立完整的Parent-gene关系并添加特征
        # 首先建立ID到特征的快速查找索引
        feature_index = {f['feature_id']: f for f in temp_features if f['feature_id']}

        for feature in temp_features:
            gene_id = feature['gene_id']

            # 如果没有直接的gene_id，从Parent递归追溯
            if not gene_id and feature['parent_ids']:
                gene_id = self._trace_gene_from_parent_recursive(
                    feature['parent_ids'], feature_index, max_depth=5
                )

            # 如果有gene_id，添加到特征列表
            if gene_id:
                # 记录特征ID到基因ID的映射
                if feature['feature_id']:
                    self.feature_to_gene[feature['feature_id']] = gene_id

                # 添加到区间索引（用于快速查找）
                chrom = feature['chrom']
                self.feature_starts[chrom].append(feature['start'])
                self.feature_data[chrom].append(
                    (feature['end'], feature['feature_type'], gene_id)
                )

        # 对特征进行排序（排序后才能使用二分查找）
        for chrom in self.feature_starts:
            # 同时排序starts和data，保持对应关系
            combined = list(zip(self.feature_starts[chrom], self.feature_data[chrom]))
            combined.sort(key=lambda x: x[0])  # 按start排序
            self.feature_starts[chrom] = [x[0] for x in combined]
            self.feature_data[chrom] = [x[1] for x in combined]

        # 对基因也进行排序
        for chrom in self.gene_starts:
            combined = list(zip(self.gene_starts[chrom], self.gene_data[chrom]))
            combined.sort(key=lambda x: x[0])
            self.gene_starts[chrom] = [x[0] for x in combined]
            self.gene_data[chrom] = [x[1] for x in combined]

        gene_count = sum(len(starts) for starts in self.gene_starts.values())
        feature_count = sum(len(starts) for starts in self.feature_starts.values())

        self.logger.info(f"GFF解析完成|GFF parsing completed: {gene_count}个基因, {feature_count}个特征|genes, {feature_count} features")

    def _parse_attributes(self, attributes_str):
        """解析GFF属性字段|Parse GFF attributes field"""
        attrs = {}
        for attr in attributes_str.split(';'):
            attr = attr.strip()
            if '=' in attr:
                key, value = attr.split('=', 1)
                attrs[key.strip()] = value.strip()
        return attrs

    def _trace_gene_from_parent(self, parent_ids, max_depth=5):
        """从Parent ID追溯基因ID（保留用于兼容性）|Trace gene ID from parent ID (kept for compatibility)"""
        return self._trace_gene_from_parent_recursive(parent_ids, {}, max_depth)

    def _trace_gene_from_parent_recursive(self, parent_ids, feature_index, max_depth=5, visited=None):
        """从Parent ID递归追溯基因ID|Recursively trace gene ID from parent ID

        Args:
            parent_ids: 父特征ID（可能包含多个，逗号分隔）
            feature_index: 特征ID到特征信息的索引字典
            max_depth: 最大追溯深度，防止无限循环
            visited: 已访问的ID集合，防止循环引用

        Returns:
            基因ID或None
        """
        if max_depth <= 0:
            return None

        if visited is None:
            visited = set()

        # 处理多个Parent|Handle multiple parents
        parents = parent_ids.split(',') if isinstance(parent_ids, str) else [parent_ids]

        for parent in parents:
            parent = parent.strip()

            # 防止循环引用
            if parent in visited:
                continue
            visited.add(parent)

            # 如果父特征直接映射到基因
            if parent in self.feature_to_gene:
                return self.feature_to_gene[parent]

            # 检查父特征是否是基因ID本身（在gene_starts或gene_data中）
            for chrom_genes in self.gene_data.values():
                for end, gene_id in chrom_genes:
                    if gene_id == parent:
                        return parent

            # 在feature_index中查找父特征，继续追溯
            if parent in feature_index:
                parent_feature = feature_index[parent]
                if parent_feature['parent_ids']:
                    # 递归追溯
                    gene_id = self._trace_gene_from_parent_recursive(
                        parent_feature['parent_ids'],
                        feature_index,
                        max_depth - 1,
                        visited.copy()
                    )
                    if gene_id:
                        return gene_id

        return None



class VariantAnnotator:
    """变异注释器|Variant Annotator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.gff_parser = GFFParser(logger)

    def annotate(self):
        """执行注释|Perform annotation"""
        self.logger.info("开始变异注释流程|Starting variant annotation pipeline")

        # 解析GFF文件|Parse GFF file
        self.gff_parser.parse(self.config.gff_file)

        # 处理VCF文件并流式写入结果|Process VCF file and stream write results
        self.logger.info(f"开始处理VCF文件|Starting to process VCF file: {self.config.vcf_file}")
        variant_count = self._process_vcf_stream()

        self.logger.info(f"注释完成|Annotation completed: {self.config.output_file}")
        self.logger.info(f"共处理|Total processed: {variant_count}个变异|variants")

    def _process_vcf_stream(self):
        """流式处理VCF文件并写入结果|Stream process VCF file and write results"""
        variant_count = 0

        # 判断是否为压缩文件|Check if file is gzipped
        open_func = gzip.open if self.config.vcf_file.endswith('.gz') else open

        with open(self.config.output_file, 'w') as out_f:
            # 写入表头|Write header
            out_f.write('Chr\tPos\tRef\tAlt\tGene\tType\n')

            with open_func(self.config.vcf_file, 'rt') as f:
                for line in f:
                    line = line.strip()

                    # 跳过注释行|Skip comment lines
                    if line.startswith('#'):
                        continue

                    fields = line.split('\t')
                    if len(fields) < 5:
                        continue

                    chrom = fields[0]
                    pos = int(fields[1])
                    ref = fields[3]
                    alt = fields[4]

                    # 处理多个ALT|Handle multiple ALTs
                    alts = alt.split(',')
                    for alt_allele in alts:
                        gene_id, variant_type = self._annotate_variant(chrom, pos)

                        # 直接写入文件|Write directly to file
                        out_f.write(f"{chrom}\t{pos}\t{ref}\t{alt_allele}\t{gene_id}\t{variant_type}\n")
                        variant_count += 1

        return variant_count

    def _annotate_variant(self, chrom, pos):
        """
        注释单个变异|Annotate single variant（使用二分查找优化）

        返回|Returns: (gene_id, variant_type)
            - gene_id: 基因ID或'intergenic'
            - variant_type: 特征类型名称或'intergenic'
        """
        # 使用二分查找快速定位可能的特征区间
        if chrom in self.gff_parser.feature_starts:
            starts = self.gff_parser.feature_starts[chrom]
            data = self.gff_parser.feature_data[chrom]

            # 二分查找：找到第一个start > pos的位置
            idx = bisect.bisect_right(starts, pos)

            # 从idx向前遍历，找到所有start <= pos的特征并检查是否包含pos
            # 由于特征按start排序，一旦start < pos且不包含pos，后面的更不可能包含
            matched_features = []

            for i in range(idx - 1, -1, -1):
                start = starts[i]
                # 如果start已经远小于pos，且已经检查了足够多的特征，就停止
                if i < idx - 100:  # 最多向前检查100个特征
                    if start < pos - 100000:  # start距离pos超过100kb就停止
                        break

                end, feature_type, gene_id = data[i]

                # 只有当start <= pos <= end时才匹配
                if start <= pos <= end:
                    priority = self.gff_parser.FEATURE_PRIORITY.get(feature_type, 999)
                    matched_features.append((priority, feature_type, gene_id))

            if matched_features:
                # 按优先级排序，选择最底层的特征
                matched_features.sort(key=lambda x: x[0])
                _, feature_type, gene_id = matched_features[0]
                return gene_id, feature_type

        # 如果没有匹配到任何特征，检查是否在基因范围内（使用二分查找）
        if chrom in self.gff_parser.gene_starts:
            gene_starts = self.gff_parser.gene_starts[chrom]
            gene_data = self.gff_parser.gene_data[chrom]

            # 二分查找
            idx = bisect.bisect_right(gene_starts, pos)

            # 向前检查是否在某个基因范围内
            for i in range(idx - 1, -1, -1):
                start = gene_starts[i]
                if i < idx - 50:  # 最多检查50个基因
                    if start < pos - 1000000:  # 距离超过1Mb就停止
                        break

                end, gene_id = gene_data[i]
                if start <= pos <= end:
                    return gene_id, 'gene'

        # 基因间区
        return 'intergenic', 'intergenic'
