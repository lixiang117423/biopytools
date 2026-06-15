"""SNP区域基因提取核心处理器|SNP Region Gene Extractor Core Processor"""

import re
import os
from .utils import run_command


class SnpRegionProcessor:
    """SNP区域基因提取处理器|SNP Region Gene Extractor Processor"""

    def __init__(self, config, logger):
        """
        初始化处理器|Initialize processor

        参数|Parameters:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 数据结构|Data structures
        self.gene_index = {}  # {chrom: [gene_dict, ...]}
        self.mrna_gene_map = {}  # {mrna_id: gene_id}

    def parse_gff3(self):
        """
        解析GFF3文件，建立基因索引|Parse GFF3 file and build gene index

        索引结构|Index structure:
            self.gene_index = {
                'Chr01': [
                    {
                        'gene_id': 'Gene001',
                        'start': 1000,
                        'end': 5000,
                        'strand': '+',
                        'mrnas': ['Gene001.1', 'Gene001.2'],
                        'exons': [(1000, 1500), (2000, 2500), ...]
                    },
                    ...
                ],
                ...
            }
            self.mrna_gene_map = {
                'Gene001.1': 'Gene001',
                'Gene001.2': 'Gene001',
                ...
            }
        """
        self.logger.info(f"开始解析GFF3文件|Starting to parse GFF3 file: {self.config.gff_file}")

        gene_data = {}  # 临时存储基因信息 {gene_id: gene_dict}
        mrna_to_gene = {}  # mRNA到基因的映射
        mrna_exons = {}  # mRNA的外显子 {mrna_id: [(start, end), ...]}

        try:
            with open(self.config.gff_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()

                    # 跳过注释行和空行|Skip comments and empty lines
                    if not line or line.startswith('#'):
                        continue

                    parts = line.split('\t')
                    if len(parts) < 9:
                        continue

                    chrom = parts[0]
                    feature_type = parts[2]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attributes = parts[8]

                    # 解析属性字段|Parse attributes
                    attrs_dict = self._parse_attributes(attributes)

                    # 处理gene特征|Process gene feature
                    if feature_type == 'gene':
                        gene_id = attrs_dict.get('ID')
                        if gene_id:
                            gene_data[gene_id] = {
                                'gene_id': gene_id,
                                'chrom': chrom,
                                'start': start,
                                'end': end,
                                'strand': strand,
                                'mrnas': [],
                                'exons': []
                            }

                    # 处理mRNA特征|Process mRNA feature
                    elif feature_type in ['mRNA', 'transcript']:
                        mrna_id = attrs_dict.get('ID')
                        parent_id = attrs_dict.get('Parent')
                        if mrna_id and parent_id:
                            if parent_id in gene_data:
                                gene_data[parent_id]['mrnas'].append(mrna_id)
                                mrna_to_gene[mrna_id] = parent_id
                                mrna_exons[mrna_id] = []

                    # 处理exon特征|Process exon feature
                    elif feature_type == 'exon':
                        parent_id = attrs_dict.get('Parent')
                        if parent_id and parent_id in mrna_exons:
                            mrna_exons[parent_id].append((start, end))

            # 合并外显子信息到基因|Merge exon information to genes
            for gene_id, gene_info in gene_data.items():
                # 收集该基因所有mRNA的外显子|Collect all exons from all mRNAs of this gene
                all_exons = set()
                for mrna_id in gene_info['mrnas']:
                    if mrna_id in mrna_exons:
                        all_exons.update(mrna_exons[mrna_id])

                # 排序外显子|Sort exons
                gene_info['exons'] = sorted(list(all_exons), key=lambda x: x[0])

            # 构建染色体索引|Build chromosome index
            for gene_id, gene_info in gene_data.items():
                chrom_gene = gene_info.copy()
                gene_chrom = gene_info['chrom']  # 从基因信息中获取染色体

                if gene_chrom not in self.gene_index:
                    self.gene_index[gene_chrom] = []

                # 按基因起始位置排序|Sort by gene start position
                self.gene_index[gene_chrom].append(chrom_gene)

            # 对每个染色体的基因按起始位置排序|Sort genes by start position for each chromosome
            for chrom in self.gene_index:
                self.gene_index[chrom].sort(key=lambda x: x['start'])

            # 建立mRNA到基因的映射|Build mRNA to gene mapping
            self.mrna_gene_map = mrna_to_gene

            # 统计信息|Statistics
            total_genes = sum(len(genes) for genes in self.gene_index.values())
            total_mrnas = len(self.mrna_gene_map)

            self.logger.info(f"GFF3解析完成|GFF3 parsing completed: {total_genes} 个基因|genes, {total_mrnas} 个转录本|transcripts")
            self.logger.info(f"覆盖染色体|Covered chromosomes: {len(self.gene_index)}")

        except FileNotFoundError:
            self.logger.error(f"GFF3文件不存在|GFF3 file not found: {self.config.gff_file}")
            raise

        except Exception as e:
            self.logger.error(f"解析GFF3文件出错|Error parsing GFF3 file: {e}")
            raise

    def _parse_attributes(self, attributes_str):
        """
        解析GFF3属性字段|Parse GFF3 attributes field

        参数|Parameters:
            attributes_str: 属性字符串|Attributes string (e.g., "ID=Gene001;Name=Test")

        返回|Returns:
            dict: 属性字典|Attributes dictionary
        """
        attrs = {}
        for attr in attributes_str.split(';'):
            attr = attr.strip()
            if '=' in attr:
                key, value = attr.split('=', 1)
                attrs[key.strip()] = value.strip()
        return attrs

    def get_snp_feature(self, snp_pos, gene_start, gene_end, strand, exon_list):
        """
        判断SNP在基因中的位置特征|Determine SNP feature in gene

        参数|Parameters:
            snp_pos: SNP坐标|SNP position
            gene_start: 基因起始位置|Gene start position
            gene_end: 基因终止位置|Gene end position
            strand: 正负链|Strand (+/-)
            exon_list: 外显子列表|Exon list [(start, end), ...]

        返回|Returns:
            str: 'promoter', 'exon', 'intron', or 'intergenic'
        """
        promoter_length = self.config.promoter

        # 1. 判断启动子|Check promoter
        if strand == '+':
            promoter_start = max(1, gene_start - promoter_length)
            promoter_end = gene_start
        else:  # strand == '-'
            promoter_start = gene_end
            promoter_end = gene_end + promoter_length

        if promoter_start <= snp_pos <= promoter_end:
            return 'promoter'

        # 2. 判断外显子|Check exon
        for exon_start, exon_end in exon_list:
            if exon_start <= snp_pos <= exon_end:
                return 'exon'

        # 3. 判断内含子（在基因范围内但不在外显子）|Check intron
        if gene_start <= snp_pos <= gene_end:
            return 'intron'

        return 'intergenic'

    def calculate_distance(self, snp_pos, gene_start, gene_end, strand, feature, promoter_length):
        """
        计算SNP到基因的距离|Calculate distance from SNP to gene

        参数|Parameters:
            snp_pos: SNP坐标|SNP position
            gene_start: 基因起始位置|Gene start position
            gene_end: 基因终止位置|Gene end position
            strand: 正负链|Strand (+/-)
            feature: SNP特征|SNP feature (promoter/exon/intron/intergenic)
            promoter_length: 启动子长度|Promoter length

        返回|Returns:
            int: 距离值（始终为正或零）|Distance value (always >= 0)
        """
        if feature == 'promoter':
            # 启动子区域：返回到基因边界的距离|Promoter region: distance to gene boundary
            if strand == '+':
                # 正链启动子：SNP在gene_start上游|Positive strand: SNP upstream of gene_start
                return gene_start - snp_pos
            else:  # strand == '-'
                # 负链启动子：SNP在gene_end下游|Negative strand: SNP downstream of gene_end
                return snp_pos - gene_end
        else:
            # exon/intron区域：返回从基因起点的偏移|Exon/intron: offset from gene start
            if strand == '+':
                return snp_pos - gene_start
            else:  # strand == '-'
                # 负链：从基因终点往回算的距离|Negative: distance from gene end
                return gene_end - snp_pos

    def find_genes_for_snp(self, chrom, snp_pos):
        """
        为SNP查找相关基因|Find related genes for SNP

        参数|Parameters:
            chrom: 染色体|Chromosome
            snp_pos: SNP坐标|SNP position

        返回|Returns:
            list: 匹配的基因列表|List of matched genes
                  Each item: {gene_id, mrna_id, strand, feature, distance}
        """
        left = self.config.left
        right = self.config.right
        promoter_length = self.config.promoter

        # 计算查询区间|Calculate query region
        query_start = max(1, snp_pos - left)
        query_end = snp_pos + right

        matched_genes = []

        # 获取该染色体上的基因|Get genes on this chromosome
        if chrom not in self.gene_index:
            return matched_genes

        for gene in self.gene_index[chrom]:
            gene_id = gene['gene_id']
            gene_start = gene['start']
            gene_end = gene['end']
            strand = gene['strand']

            # 计算基因的"有效范围"（包含启动子）|Calculate effective range (including promoter)
            if strand == '+':
                effective_start = max(1, gene_start - promoter_length)
                effective_end = gene_end
            else:  # strand == '-'
                effective_start = gene_start
                effective_end = gene_end + promoter_length

            # 判断是否有重叠|Check if there is overlap
            if not (query_end < effective_start or query_start > effective_end):
                # 有重叠，进一步判断SNP特征|Has overlap, determine SNP feature
                feature = self.get_snp_feature(
                    snp_pos, gene_start, gene_end, strand,
                    gene['exons']
                )

                distance = self.calculate_distance(snp_pos, gene_start, gene_end, strand, feature, promoter_length)

                for mrna_id in gene['mrnas']:
                    matched_genes.append({
                        'gene_id': gene_id,
                        'mrna_id': mrna_id,
                        'strand': strand,
                        'feature': feature,
                        'distance': distance
                    })

        return matched_genes

    def extract_all_sequences(self):
        """
        使用gffread提取所有CDS和蛋白序列|Extract all CDS and protein sequences using gffread

        返回|Returns:
            bool: 成功返回True，失败返回False|True on success, False on failure
        """
        self.logger.info("开始提取所有CDS和蛋白序列|Starting to extract all CDS and protein sequences")

        # 构建gffread命令|Build gffread command
        cmd = [
            self.config.gffread_path,
            '-g', self.config.genome_file,
            '-x', self.config.temp_cds,
            '-y', self.config.temp_protein,
            self.config.gff_file
        ]

        success, stdout, stderr = run_command(cmd, self.logger)

        if success:
            self.logger.info(f"成功提取序列|Successfully extracted sequences")
            self.logger.debug(f"CDS序列文件|CDS output: {self.config.temp_cds}")
            self.logger.debug(f"蛋白序列文件|Protein output: {self.config.temp_protein}")
            return True
        else:
            self.logger.error(f"gffread执行失败|gffread execution failed")
            if stderr:
                self.logger.error(f"错误信息|Error message: {stderr}")
            return False

    def extract_region_sequences(self, mrna_ids):
        """
        使用seqkit提取指定mRNA的序列|Extract sequences for specified mRNAs using seqkit

        参数|Parameters:
            mrna_ids: mRNA ID列表|List of mRNA IDs

        返回|Returns:
            bool: 成功返回True，失败返回False|True on success, False on failure
        """
        self.logger.info(f"开始提取区域基因序列|Starting to extract region gene sequences: {len(mrna_ids)} 个mRNA|mRNAs")

        if not mrna_ids:
            self.logger.warning("没有需要提取的mRNA|No mRNAs to extract")
            return False

        # 提取CDS序列|Extract CDS sequences
        cmd_cds = [
            self.config.seqkit_path,
            'grep',
        ]

        # 为每个mRNA ID添加一个-p参数|Add -p parameter for each mRNA ID
        for mrna_id in mrna_ids:
            cmd_cds.extend(['-p', mrna_id])

        cmd_cds.append(self.config.temp_cds)

        success_cds, stdout_cds, stderr_cds = run_command(cmd_cds, self.logger)

        if success_cds:
            with open(self.config.cds_output, 'w') as f:
                f.write(stdout_cds)
            self.logger.info(f"成功提取CDS序列|Successfully extracted CDS sequences: {self.config.cds_output}")
        else:
            self.logger.error(f"seqkit提取CDS序列失败|seqkit CDS extraction failed")
            if stderr_cds:
                self.logger.error(f"错误信息|Error message: {stderr_cds}")
            return False

        # 提取蛋白序列|Extract protein sequences
        cmd_protein = [
            self.config.seqkit_path,
            'grep',
        ]

        # 为每个mRNA ID添加一个-p参数|Add -p parameter for each mRNA ID
        for mrna_id in mrna_ids:
            cmd_protein.extend(['-p', mrna_id])

        cmd_protein.append(self.config.temp_protein)

        success_protein, stdout_protein, stderr_protein = run_command(cmd_protein, self.logger)

        if success_protein:
            with open(self.config.protein_output, 'w') as f:
                f.write(stdout_protein)
            self.logger.info(f"成功提取蛋白序列|Successfully extracted protein sequences: {self.config.protein_output}")
        else:
            self.logger.error(f"seqkit提取蛋白序列失败|seqkit protein extraction failed")
            if stderr_protein:
                self.logger.error(f"错误信息|Error message: {stderr_protein}")
            return False

        return True

    def cleanup_temp_files(self):
        """清理临时文件|Clean up temporary files"""
        try:
            if os.path.exists(self.config.temp_cds):
                os.remove(self.config.temp_cds)
                self.logger.debug(f"删除临时文件|Deleted temp file: {self.config.temp_cds}")

            if os.path.exists(self.config.temp_protein):
                os.remove(self.config.temp_protein)
                self.logger.debug(f"删除临时文件|Deleted temp file: {self.config.temp_protein}")

        except Exception as e:
            self.logger.warning(f"清理临时文件时出错|Error cleaning temp files: {e}")
