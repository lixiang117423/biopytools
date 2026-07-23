"""
单倍型可视化解析模块|Haplotype Visualization Parser Module

解析GFF和VCF文件|Parse GFF and VCF files
"""

import subprocess
import tempfile
import gzip
import shutil
import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Dict, Optional, Any
from collections import defaultdict

from .utils import parse_gff_attributes


@dataclass
class GFFFeature:
    """GFF特征数据类|GFF Feature Data Class"""
    seqid: str  # 染色体/序列ID|Chromosome/sequence ID
    source: str  # 来源|Source
    feature_type: str  # 特征类型|Feature type (gene, mRNA, CDS, exon, etc.)
    start: int  # 起始位置|Start position (1-based)
    end: int  # 终止位置|End position (1-based)
    score: str  # 分数|Score
    strand: str  # 链|Strand (+, -, or .)
    phase: str  # 相位|Phase (0, 1, 2, or .)
    attributes: Dict[str, str]  # 属性|Attributes
    parent: Optional[str] = None  # 父特征ID|Parent feature ID

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        if isinstance(self.attributes, str):
            self.attributes = parse_gff_attributes(self.attributes)


@dataclass
class Gene:
    """基因数据类|Gene Data Class"""
    gene_id: str  # 基因ID|Gene ID
    gene_name: str  # 基因名称|Gene name
    chrom: str  # 染色体|Chromosome
    start: int  # 起始位置|Start position
    end: int  # 终止位置|End position
    strand: str  # 链|Strand
    mrnas: List['MRNA']  # mRNA转录本列表|mRNA transcript list

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        if self.mrnas is None:
            self.mrnas = []


@dataclass
class MRNA:
    """mRNA数据类|mRNA Data Class"""
    mrna_id: str  # mRNA ID|mRNA ID
    transcript_id: str  # 转录本ID|Transcript ID
    start: int  # 起始位置|Start position
    end: int  # 终止位置|End position
    strand: str  # 链|Strand
    cds_list: List['CDS']  # CDS列表|CDS list
    exon_list: List['Exon']  # 外显子列表|Exon list
    five_prime_utr: Optional['UTR'] = None  # 5' UTR|5' UTR
    three_prime_utr: Optional['UTR'] = None  # 3' UTR|3' UTR

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        if self.cds_list is None:
            self.cds_list = []
        if self.exon_list is None:
            self.exon_list = []


@dataclass
class CDS:
    """CDS数据类|CDS Data Class"""
    start: int  # 起始位置|Start position
    end: int  # 终止位置|End position
    phase: int  # 相位|Phase (0, 1, 2)


@dataclass
class Exon:
    """外显子数据类|Exon Data Class"""
    start: int  # 起始位置|Start position
    end: int  # 终止位置|End position


@dataclass
class UTR:
    """UTR数据类|UTR Data Class"""
    start: int  # 起始位置|Start position
    end: int  # 终止位置|End position
    utr_type: str  # UTR类型|UTR type (five_prime_UTR or three_prime_UTR)


@dataclass
class Variant:
    """变异数据类|Variant Data Class"""
    chrom: str  # 染色体|Chromosome
    pos: int  # 位置|Position (1-based)
    ref: str  # 参考等位基因|Reference allele
    alt: List[str]  # 变异等位基因列表|Alternative allele list
    qual: float  # 质量值|Quality value
    filter: str  # 过滤状态|Filter status
    info: Dict[str, Any]  # INFO字段|INFO field
    format: List[str]  # FORMAT字段列表|FORMAT field list
    samples: Dict[str, Dict[str, Any]]  # 样本基因型数据|Sample genotype data

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        if isinstance(self.format, str):
            self.format = self.format.split(':')
        if isinstance(self.info, str):
            self.info = self._parse_info(self.info)

    @staticmethod
    def _parse_info(info_str: str) -> Dict[str, Any]:
        """解析INFO字段|Parse INFO field"""
        info = {}
        for item in info_str.split(';'):
            item = item.strip()
            if '=' in item:
                key, value = item.split('=', 1)
                info[key] = value
            else:
                info[item] = True
        return info


class GFFParser:
    """GFF文件解析器|GFF File Parser"""

    def __init__(self, logger):
        """初始化GFF解析器|Initialize GFF parser

        Args:
            logger: 日志器|Logger
        """
        self.logger = logger
        self.features: List[GFFFeature] = []
        self.genes: List[Gene] = []

    def parse_file(self, gff_file: Path, chrom: Optional[str] = None,
                   start: Optional[int] = None, end: Optional[int] = None) -> List[Gene]:
        """解析GFF文件（支持区域过滤）|Parse GFF file (with region filtering)

        Args:
            gff_file: GFF文件路径|GFF file path
            chrom: 染色体名称（可选）|Chromosome name (optional)
            start: 起始位置（可选）|Start position (optional)
            end: 终止位置（可选）|End position (optional)

        Returns:
            list: 基因列表|Gene list
        """
        self.logger.info(f"开始解析GFF文件|Start parsing GFF file: {gff_file}")
        if chrom:
            self.logger.info(f"仅解析区域|Parsing region only: {chrom}:{start}-{end}")

        self.features = []
        self.genes = []

        line_count = 0
        filtered_count = 0

        with open(gff_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                line_count += 1

                # 跳过注释行和空行|Skip comment lines and empty lines
                if not line or line.startswith('#'):
                    continue

                # 快速过滤：检查是否在目标区域|Quick filter: check if in target region
                if chrom:
                    fields = line.split('\t')
                    if len(fields) < 5:
                        continue

                    line_chrom = fields[0]
                    line_start = int(fields[3])
                    line_end = int(fields[4])

                    # 跳过不匹配的染色体|Skip non-matching chromosomes
                    if line_chrom != chrom:
                        continue

                    # 跳过不在区域内的特征|Skip features outside region
                    if start and end:
                        if line_end < start or line_start > end:
                            filtered_count += 1
                            continue

                # 解析特征|Parse feature
                try:
                    feature = self._parse_line(line)
                    if feature:
                        self.features.append(feature)
                except Exception as e:
                    self.logger.warning(
                        f"解析第{line_num}行时出错|Error parsing line {line_num}: {e}"
                    )

                # 每处理10000行输出进度|Output progress every 10000 lines
                if line_count % 10000 == 0:
                    self.logger.debug(f"已处理|Processed {line_count} 行，加载|loaded {len(self.features)} 个特征")

        self.logger.info(f"GFF文件读取完成|GFF file reading completed: "
                        f"总行数|Total lines: {line_count}, "
                        f"加载特征|Loaded features: {len(self.features)}, "
                        f"过滤特征|Filtered features: {filtered_count}")

        # 构建基因层次结构|Build gene hierarchy
        self.genes = self._build_gene_hierarchy()

        self.logger.info(f"GFF解析完成，共{len(self.genes)}个基因|"
                        f"GFF parsing completed, total {len(self.genes)} genes")

        return self.genes

    def _parse_line(self, line: str) -> Optional[GFFFeature]:
        """解析GFF行|Parse GFF line

        Args:
            line: GFF行|GFF line

        Returns:
            GFFFeature: GFF特征对象|GFF feature object
        """
        fields = line.split('\t')

        if len(fields) < 9:
            return None

        seqid = fields[0]
        source = fields[1]
        feature_type = fields[2]
        start = int(fields[3])
        end = int(fields[4])
        score = fields[5]
        strand = fields[6]
        phase = fields[7]
        attributes = parse_gff_attributes(fields[8])

        # 获取Parent属性|Get Parent attribute
        parent = attributes.get('Parent') or attributes.get('parent')

        return GFFFeature(
            seqid=seqid,
            source=source,
            feature_type=feature_type,
            start=start,
            end=end,
            score=score,
            strand=strand,
            phase=phase,
            attributes=attributes,
            parent=parent
        )

    def _build_gene_hierarchy(self) -> List[Gene]:
        """构建基因层次结构|Build gene hierarchy

        Returns:
            list: 基因列表|Gene list
        """
        # 按类型分组特征|Group features by type
        features_by_type = defaultdict(list)
        for feature in self.features:
            features_by_type[feature.feature_type].append(feature)

        genes = []

        # 处理gene特征|Process gene features
        for gene_feature in features_by_type.get('gene', []):
            gene_id = gene_feature.attributes.get('ID') or gene_feature.attributes.get('Id')
            gene_name = gene_feature.attributes.get('Name', gene_id)

            gene = Gene(
                gene_id=gene_id,
                gene_name=gene_name,
                chrom=gene_feature.seqid,
                start=gene_feature.start,
                end=gene_feature.end,
                strand=gene_feature.strand,
                mrnas=[]
            )

            # 查找mRNA子特征|Find mRNA sub-features
            for mrna_feature in features_by_type.get('mRNA', []):
                if mrna_feature.parent == gene_id:
                    mrna = self._build_mrna(mrna_feature, features_by_type)
                    gene.mrnas.append(mrna)

            # 如果没有mRNA，尝试使用transcript特征|If no mRNA, try transcript feature
            if not gene.mrnas:
                for transcript_feature in features_by_type.get('transcript', []):
                    if transcript_feature.parent == gene_id:
                        mrna = self._build_mrna(transcript_feature, features_by_type)
                        gene.mrnas.append(mrna)

            genes.append(gene)

        return genes

    def _build_mrna(self, mrna_feature: GFFFeature, features_by_type: Dict) -> MRNA:
        """构建mRNA对象|Build mRNA object

        Args:
            mrna_feature: mRNA特征|mRNA feature
            features_by_type: 按类型分组的特征|Features grouped by type

        Returns:
            MRNA: mRNA对象|mRNA object
        """
        mrna_id = mrna_feature.attributes.get('ID') or mrna_feature.attributes.get('Id')
        transcript_id = mrna_feature.attributes.get('Name') or mrna_id

        mrna = MRNA(
            mrna_id=mrna_id,
            transcript_id=transcript_id,
            start=mrna_feature.start,
            end=mrna_feature.end,
            strand=mrna_feature.strand,
            cds_list=[],
            exon_list=[]
        )

        # 查找CDS子特征|Find CDS sub-features
        for cds_feature in features_by_type.get('CDS', []):
            if cds_feature.parent == mrna_id:
                mrna.cds_list.append(CDS(
                    start=cds_feature.start,
                    end=cds_feature.end,
                    phase=int(cds_feature.phase) if cds_feature.phase != '.' else 0
                ))

        # 查找exon子特征|Find exon sub-features
        for exon_feature in features_by_type.get('exon', []):
            if exon_feature.parent == mrna_id:
                mrna.exon_list.append(Exon(
                    start=exon_feature.start,
                    end=exon_feature.end
                ))

        # 查找UTR子特征|Find UTR sub-features
        for utr_feature in features_by_type.get('five_prime_UTR', []):
            if utr_feature.parent == mrna_id:
                mrna.five_prime_utr = UTR(
                    start=utr_feature.start,
                    end=utr_feature.end,
                    utr_type='five_prime_UTR'
                )

        for utr_feature in features_by_type.get('three_prime_UTR', []):
            if utr_feature.parent == mrna_id:
                mrna.three_prime_utr = UTR(
                    start=utr_feature.start,
                    end=utr_feature.end,
                    utr_type='three_prime_UTR'
                )

        return mrna

    def get_genes_in_region(self, chrom: str, start: int, end: int) -> List[Gene]:
        """获取指定区域的基因|Get genes in specified region

        Args:
            chrom: 染色体|Chromosome
            start: 起始位置|Start position
            end: 终止位置|End position

        Returns:
            list: 基因列表|Gene list
        """
        genes_in_region = []

        for gene in self.genes:
            # 检查染色体|Check chromosome
            if gene.chrom != chrom:
                continue

            # 检查区间重叠|Check interval overlap
            if not (gene.end < start or gene.start > end):
                genes_in_region.append(gene)

        return genes_in_region


class VCFParser:
    """VCF文件解析器|VCF File Parser"""

    def __init__(self, logger):
        """初始化VCF解析器|Initialize VCF parser

        Args:
            logger: 日志器|Logger
        """
        self.logger = logger
        self.variants: List[Variant] = []
        self.sample_names: List[str] = []

    def parse_file(self, vcf_file: Path, chrom: str, start: int, end: int,
                   use_bcftools: bool = True) -> List[Variant]:
        """解析VCF文件指定区域（使用bcftools预筛选）|Parse VCF file for specified region (with bcftools pre-filtering)

        Args:
            vcf_file: VCF文件路径|VCF file path
            chrom: 染色体|Chromosome
            start: 起始位置|Start position
            end: 终止位置|End position
            use_bcftools: 是否使用bcftools预筛选|Whether to use bcftools for pre-filtering

        Returns:
            list: 变异列表|Variant list
        """
        self.logger.info(
            f"开始解析VCF文件区域|Start parsing VCF file region: "
            f"{chrom}:{start}-{end}"
        )

        self.variants = []
        self.sample_names = []

        # 检查是否为压缩文件|Check if compressed file
        is_gzipped = str(vcf_file).endswith('.gz')

        # 使用bcftools预筛选|Use bcftools for pre-filtering
        if use_bcftools:
            actual_vcf_file = self._extract_region_with_bcftools(vcf_file, chrom, start, end)
        else:
            actual_vcf_file = vcf_file

        # 解析VCF文件|Parse VCF file
        try:
            self._parse_vcf_content(actual_vcf_file, chrom, start, end)
        finally:
            # 清理临时文件|Clean up temporary file
            if use_bcftools and actual_vcf_file != vcf_file and actual_vcf_file.exists():
                actual_vcf_file.unlink()
                self.logger.debug(f"清理临时文件|Cleaned temporary file: {actual_vcf_file}")

        self.logger.info(f"VCF解析完成，共{len(self.variants)}个变异|"
                        f"VCF parsing completed, total {len(self.variants)} variants")

        return self.variants

    def _extract_region_with_bcftools(self, vcf_file: Path, chrom: str,
                                     start: int, end: int) -> Path:
        """使用bcftools提取目标区域|Extract target region using bcftools

        Args:
            vcf_file: VCF文件路径|VCF file path
            chrom: 染色体|Chromosome
            start: 起始位置|Start position
            end: 终止位置|End position

        Returns:
            Path: 提取的VCF文件路径|Extracted VCF file path
        """
        self.logger.info(f"使用bcftools提取区域|Extracting region with bcftools: {chrom}:{start}-{end}")

        # 创建临时文件|Create temporary file
        temp_vcf = tempfile.NamedTemporaryFile(
            mode='wb',
            suffix='.vcf',
            delete=False,
            dir=os.path.dirname(str(vcf_file))
        )
        temp_vcf_path = Path(temp_vcf.name)
        temp_vcf.close()

        try:
            # 检查并构建索引|Check and build index
            if not self._check_index_exists(vcf_file):
                self.logger.info(f"VCF索引不存在，正在构建索引|VCF index not found, building index")
                self._build_vcf_index(vcf_file)

            # 构建bcftools命令|Build bcftools command
            region = f"{chrom}:{start}-{end}"
            cmd = [
                'bcftools',
                'view',
                '-r', region,
                '-O', 'v',
                '-o', str(temp_vcf_path),
                str(vcf_file)
            ]

            self.logger.debug(f"执行命令|Executing command: {' '.join(cmd)}")

            # 运行bcftools|Run bcftools
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )

            self.logger.info(f"bcftools提取完成|bcftools extraction completed: {temp_vcf_path}")

            return temp_vcf_path

        except subprocess.CalledProcessError as e:
            self.logger.error(f"bcftools执行失败|bcftools execution failed: {e.stderr}")
            # 如果bcftools失败，返回原文件|If bcftools fails, return original file
            if temp_vcf_path.exists():
                temp_vcf_path.unlink()
            return vcf_file
        except FileNotFoundError:
            self.logger.warning("未找到bcftools，直接读取原VCF文件|bcftools not found, reading original VCF file")
            if temp_vcf_path.exists():
                temp_vcf_path.unlink()
            return vcf_file

    def _check_index_exists(self, vcf_file: Path) -> bool:
        """检查VCF索引是否存在|Check if VCF index exists

        Args:
            vcf_file: VCF文件路径|VCF file path

        Returns:
            bool: 索引是否存在|Whether index exists
        """
        # 检查.csi索引（压缩）|Check .csi index (compressed)
        csi_index = Path(str(vcf_file) + '.csi')
        if csi_index.exists():
            return True

        # 检查.tbi索引（tabix）|Check .tbi index (tabix)
        tbi_index = Path(str(vcf_file) + '.tbi')
        if tbi_index.exists():
            return True

        return False

    def _build_vcf_index(self, vcf_file: Path):
        """构建VCF索引|Build VCF index

        Args:
            vcf_file: VCF文件路径|VCF file path
        """
        self.logger.info(f"开始构建VCF索引|Start building VCF index: {vcf_file}")

        cmd = ['bcftools', 'index', str(vcf_file)]

        try:
            self.logger.debug(f"执行命令|Executing command: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=True
            )
            self.logger.info(f"VCF索引构建完成|VCF index building completed")
        except subprocess.CalledProcessError as e:
            self.logger.warning(f"VCF索引构建失败|VCF index building failed: {e.stderr}")
            raise

    def _parse_vcf_content(self, vcf_file: Path, chrom: str, start: int, end: int):
        """解析VCF文件内容|Parse VCF file content

        Args:
            vcf_file: VCF文件路径|VCF file path
            chrom: 染色体|Chromosome
            start: 起始位置|Start position
            end: 终止位置|End position
        """
        # 检查是否为gzip文件|Check if gzip file
        is_gzipped = str(vcf_file).endswith('.gz')

        # 打开文件|Open file
        if is_gzipped:
            open_func = gzip.open
            mode = 'rt'
        else:
            open_func = open
            mode = 'r'

        with open_func(vcf_file, mode) as f:
            for line in f:
                line = line.strip()

                # 解析header行|Parse header lines
                if line.startswith('##'):
                    continue
                elif line.startswith('#CHROM'):
                    # 解析样本名|Parse sample names
                    fields = line.split('\t')
                    if len(fields) > 9:
                        self.sample_names = fields[9:]
                        self.logger.info(f"检测到{len(self.sample_names)}个样本|"
                                       f"Detected {len(self.sample_names)} samples")
                    continue

                # 解析变异行|Parse variant lines
                variant = self._parse_line(line)
                if variant and variant.chrom == chrom and start <= variant.pos <= end:
                    self.variants.append(variant)

    def _parse_line(self, line: str) -> Optional[Variant]:
        """解析VCF数据行|Parse VCF data line

        Args:
            line: VCF数据行|VCF data line

        Returns:
            Variant: 变异对象|Variant object
        """
        fields = line.split('\t')

        if len(fields) < 8:
            return None

        chrom = fields[0]
        pos = int(fields[1])
        ref = fields[3]
        alt = fields[4].split(',')
        qual = float(fields[5]) if fields[5] != '.' else 0.0
        filter_status = fields[6]
        info = fields[7]

        # 解析FORMAT和样本数据|Parse FORMAT and sample data
        format_fields = []
        samples = {}
        if len(fields) > 8:
            format_fields = fields[8].split(':')
            for i, sample_data in enumerate(fields[9:], start=0):
                sample_name = self.sample_names[i] if i < len(self.sample_names) else f"sample_{i}"
                data_values = sample_data.split(':')

                sample_dict = {}
                for j, key in enumerate(format_fields):
                    if j < len(data_values):
                        value = data_values[j]
                        # 尝试转换为数值|Try to convert to numeric
                        if value == '.':
                            sample_dict[key] = None
                        elif key in ['GT']:
                            # 保留基因型为字符串|Keep genotype as string
                            sample_dict[key] = value
                        else:
                            try:
                                sample_dict[key] = int(value)
                            except ValueError:
                                try:
                                    sample_dict[key] = float(value)
                                except ValueError:
                                    sample_dict[key] = value
                    else:
                        sample_dict[key] = None

                samples[sample_name] = sample_dict

        return Variant(
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            qual=qual,
            filter=filter_status,
            info=info,
            format=format_fields,
            samples=samples
        )

    def get_genotype_matrix(self) -> Dict[str, List[int]]:
        """获取基因型矩阵|Get genotype matrix

        Returns:
            dict: 基因型矩阵字典|Genotype matrix dictionary
                {sample_name: [allele1, allele2, allele1, allele2, ...]}
        """
        genotype_matrix = {}

        for sample in self.sample_names:
            genotype_matrix[sample] = []

        for variant in self.variants:
            for sample in self.sample_names:
                if sample in variant.samples:
                    gt_str = variant.samples[sample].get('GT', './.')
                    gt = self._parse_genotype(gt_str)
                    genotype_matrix[sample].extend(gt)
                else:
                    genotype_matrix[sample].extend([-1, -1])

        return genotype_matrix

    @staticmethod
    def _parse_genotype(gt_str: str) -> List[int]:
        """解析基因型字符串|Parse genotype string

        Args:
            gt_str: 基因型字符串|Genotype string (e.g., "0/1", "1|1", "./.")

        Returns:
            list: 基因型列表|Genotype list [allele1, allele2]
        """
        # 替换分隔符|Replace separator
        gt_str = gt_str.replace('|', '/')

        if gt_str == './.' or gt_str == '.':
            return [-1, -1]

        parts = gt_str.split('/')
        if len(parts) == 2:
            try:
                return [int(parts[0]), int(parts[1])]
            except ValueError:
                return [-1, -1]
        elif len(parts) == 1:
            try:
                return [int(parts[0]), int(parts[0])]
            except ValueError:
                return [-1, -1]
        else:
            return [-1, -1]
