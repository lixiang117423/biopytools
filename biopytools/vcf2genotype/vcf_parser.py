"""
VCF文件解析模块|VCF File Parsing Module
"""

from typing import Iterator, List, Dict, Optional, Any
from collections import Counter
from .utils import FileUtils, GenotypeUtils


class VCFRecord:
    """VCF记录类|VCF Record Class"""

    def __init__(self, chrom: str, pos: int, id_field: str, ref: str, alt: str,
                 qual: str, genotypes: Dict[str, Optional[str]]):
        self.chrom = chrom
        self.pos = pos
        self.id = id_field
        self.ref = ref
        self.alt = alt
        self.qual = qual
        self.genotypes = genotypes


class FastVCFParser:
    """快速VCF解析器（使用cyvcf2）|Fast VCF Parser (using cyvcf2)"""

    def __init__(self, vcf_file: str, logger):
        self.vcf_file = vcf_file
        self.logger = logger
        self.vcf_reader = None

        try:
            import cyvcf2
            self.vcf_reader = cyvcf2.VCF(vcf_file)
            self.logger.info("使用cyvcf2加速解析VCF文件|Using cyvcf2 for fast VCF parsing")
        except ImportError:
            raise ImportError("cyvcf2未安装|cyvcf2 not installed")

    def get_samples(self) -> List[str]:
        """获取样本列表|Get sample list"""
        return list(self.vcf_reader.samples)

    def parse_records(self, target_samples: Optional[List[str]] = None) -> Iterator[Dict[str, Any]]:
        """解析VCF记录，直接yield dict|Parse VCF records, yield dict directly"""
        samples = self.get_samples()

        # 确定目标样本索引和名称|Determine target sample indices and names
        if target_samples:
            sample_indices = [i for i, s in enumerate(samples) if s in target_samples]
            final_samples = [s for s in target_samples if s in samples]
        else:
            sample_indices = list(range(len(samples)))
            final_samples = samples

        for variant in self.vcf_reader:
            # 直接构建行dict，跳过VCFRecord中间对象|Build row dict directly, skip VCFRecord
            row: Dict[str, Any] = {
                'CHROM': variant.CHROM,
                'POS': variant.POS,
                'ID': variant.ID if variant.ID else '.',
                'REF': variant.REF,
                'ALT': ','.join(variant.ALT) if variant.ALT else '.',
                'QUAL': str(variant.QUAL) if variant.QUAL is not None else '.',
            }

            # 提取目标样本基因型|Extract target sample genotypes
            gt_array = variant.genotypes
            genotype_values = []
            for i, sample_idx in enumerate(sample_indices):
                gt = gt_array[sample_idx]
                if gt[0] == -1 or gt[1] == -1:
                    row[final_samples[i]] = './.'
                    genotype_values.append('./.')
                else:
                    separator = '|' if len(gt) > 2 and gt[2] else '/'
                    gt_str = f"{gt[0]}{separator}{gt[1]}"
                    row[final_samples[i]] = gt_str
                    genotype_values.append(gt_str)

            # 统计位点内各基因型数量|Count genotype frequencies at this site
            gt_counts = Counter(genotype_values)
            for gt, count in gt_counts.items():
                row[f'GT_{gt}'] = count

            # 计算当前位点的纯合/杂合比例|Calculate per-site homozygous/heterozygous ratios
            homo_ratio, hetero_ratio = GenotypeUtils.calculate_genotype_stats(genotype_values)
            row['Homozygous_Ratio'] = f"{homo_ratio:.4f}"
            row['Heterozygous_Ratio'] = f"{hetero_ratio:.4f}"

            yield row


class StandardVCFParser:
    """标准VCF解析器（原生Python）|Standard VCF Parser (native Python)"""

    def __init__(self, vcf_file: str, logger):
        self.vcf_file = vcf_file
        self.logger = logger
        self.samples = []
        self.logger.info("使用原生Python解析VCF文件|Using native Python for VCF parsing")

    def _parse_header(self) -> List[str]:
        """解析VCF头部获取样本信息|Parse VCF header to get sample information"""
        with FileUtils.open_file(self.vcf_file) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    if len(fields) > 9:
                        self.samples = fields[9:]  # 样本从第10列开始|Samples start from column 10
                    break
        return self.samples

    def get_samples(self) -> List[str]:
        """获取样本列表|Get sample list"""
        if not self.samples:
            self._parse_header()
        return self.samples

    def parse_records(self, target_samples: Optional[List[str]] = None) -> Iterator[Dict[str, Any]]:
        """解析VCF记录，直接yield dict|Parse VCF records, yield dict directly"""
        samples = self.get_samples()

        # 预构建样本名到索引的映射，O(1)查找|Pre-build sample name to index mapping, O(1) lookup
        sample_to_idx = {name: idx for idx, name in enumerate(samples)}

        if target_samples:
            sample_indices = []
            final_samples = []
            for sample in target_samples:
                if sample in sample_to_idx:
                    sample_indices.append(sample_to_idx[sample])
                    final_samples.append(sample)
                else:
                    self.logger.warning(f"样本 '{sample}' 在VCF文件中未找到|Sample '{sample}' not found in VCF file")
        else:
            sample_indices = list(range(len(samples)))
            final_samples = samples

        with FileUtils.open_file(self.vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue

                fields = line.strip().split('\t')
                if len(fields) < 10:
                    continue

                row: Dict[str, Any] = {
                    'CHROM': fields[0],
                    'POS': int(fields[1]),
                    'ID': fields[2],
                    'REF': fields[3],
                    'ALT': fields[4],
                    'QUAL': fields[5],
                }

                # 提取目标样本基因型|Extract target sample genotypes
                gt_fields = fields[9:]
                genotype_values = []
                for i, sample_idx in enumerate(sample_indices):
                    if sample_idx < len(gt_fields):
                        gt = GenotypeUtils.parse_genotype(gt_fields[sample_idx])
                        gt_display = gt if gt is not None else './.'
                        row[final_samples[i]] = gt_display
                        genotype_values.append(gt_display)
                    else:
                        row[final_samples[i]] = './.'
                        genotype_values.append('./.')

                # 统计位点内各基因型数量|Count genotype frequencies at this site
                gt_counts = Counter(genotype_values)
                for gt, count in gt_counts.items():
                    row[f'GT_{gt}'] = count

                # 计算当前位点的纯合/杂合比例|Calculate per-site homozygous/heterozygous ratios
                homo_ratio, hetero_ratio = GenotypeUtils.calculate_genotype_stats(genotype_values)
                row['Homozygous_Ratio'] = f"{homo_ratio:.4f}"
                row['Heterozygous_Ratio'] = f"{hetero_ratio:.4f}"

                yield row


class VCFParserFactory:
    """VCF解析器工厂|VCF Parser Factory"""

    @staticmethod
    def create_parser(vcf_file: str, logger, use_fast: bool = True):
        """创建VCF解析器|Create VCF parser"""
        if use_fast:
            try:
                return FastVCFParser(vcf_file, logger)
            except ImportError:
                logger.info("cyvcf2未安装，回退到标准解析器|cyvcf2 not installed, falling back to standard parser")

        return StandardVCFParser(vcf_file, logger)
