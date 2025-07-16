"""
VCF文件解析模块 | VCF File Parsing Module
"""

from typing import Iterator, List, Dict, Optional, Any
from .utils import FileUtils, GenotypeUtils

class VCFRecord:
    """VCF记录类 | VCF Record Class"""
    
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
    """快速VCF解析器（使用cyvcf2） | Fast VCF Parser (using cyvcf2)"""
    
    def __init__(self, vcf_file: str, logger):
        self.vcf_file = vcf_file
        self.logger = logger
        self.vcf_reader = None
        
        try:
            import cyvcf2
            self.vcf_reader = cyvcf2.VCF(vcf_file)
            self.logger.info("使用cyvcf2加速解析VCF文件 | Using cyvcf2 for fast VCF parsing")
        except ImportError:
            raise ImportError("cyvcf2未安装 | cyvcf2 not installed")
    
    def get_samples(self) -> List[str]:
        """获取样本列表 | Get sample list"""
        return list(self.vcf_reader.samples)
    
    def parse_records(self, target_samples: Optional[List[str]] = None) -> Iterator[VCFRecord]:
        """解析VCF记录 | Parse VCF records"""
        samples = self.get_samples()
        
        # 确定目标样本索引 | Determine target sample indices
        if target_samples:
            sample_indices = [i for i, s in enumerate(samples) if s in target_samples]
            final_samples = target_samples
        else:
            sample_indices = list(range(len(samples)))
            final_samples = samples
        
        for variant in self.vcf_reader:
            # 提取基因型信息 | Extract genotype information
            genotypes = {}
            gt_array = variant.genotypes
            
            for i, sample_idx in enumerate(sample_indices):
                sample_name = final_samples[i] if target_samples else samples[sample_idx]
                gt = gt_array[sample_idx]
                
                # 转换基因型格式 | Convert genotype format
                if gt[0] == -1 or gt[1] == -1:  # 缺失基因型 | Missing genotype
                    genotypes[sample_name] = None
                else:
                    separator = '|' if len(gt) > 2 and gt[2] else '/'
                    genotypes[sample_name] = f"{gt[0]}{separator}{gt[1]}"
            
            yield VCFRecord(
                chrom=variant.CHROM,
                pos=variant.POS,
                id_field=variant.ID if variant.ID else '.',
                ref=variant.REF,
                alt=','.join(variant.ALT) if variant.ALT else '.',
                qual=str(variant.QUAL) if variant.QUAL is not None else '.',
                genotypes=genotypes
            )

class StandardVCFParser:
    """标准VCF解析器（原生Python） | Standard VCF Parser (native Python)"""
    
    def __init__(self, vcf_file: str, logger):
        self.vcf_file = vcf_file
        self.logger = logger
        self.samples = []
        self.logger.info("使用原生Python解析VCF文件 | Using native Python for VCF parsing")
    
    def _parse_header(self) -> List[str]:
        """解析VCF头部获取样本信息 | Parse VCF header to get sample information"""
        with FileUtils.open_file(self.vcf_file) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    fields = line.strip().split('\t')
                    if len(fields) > 9:
                        self.samples = fields[9:]  # 样本从第10列开始 | Samples start from column 10
                    break
        return self.samples
    
    def get_samples(self) -> List[str]:
        """获取样本列表 | Get sample list"""
        if not self.samples:
            self._parse_header()
        return self.samples
    
    def parse_records(self, target_samples: Optional[List[str]] = None) -> Iterator[VCFRecord]:
        """解析VCF记录 | Parse VCF records"""
        samples = self.get_samples()
        
        # 确定目标样本索引 | Determine target sample indices
        if target_samples:
            sample_indices = []
            final_samples = []
            for sample in target_samples:
                if sample in samples:
                    sample_indices.append(samples.index(sample))
                    final_samples.append(sample)
                else:
                    self.logger.warning(f"样本 '{sample}' 在VCF文件中未找到 | Sample '{sample}' not found in VCF file")
        else:
            sample_indices = list(range(len(samples)))
            final_samples = samples
        
        with FileUtils.open_file(self.vcf_file) as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chrom, pos, id_field, ref, alt, qual = fields[:6]
                format_field = fields[8] if len(fields) > 8 else ""
                
                # 解析基因型 | Parse genotypes
                genotypes = {}
                if len(fields) > 9:
                    gt_fields = fields[9:]
                    
                    for i, sample_idx in enumerate(sample_indices):
                        if sample_idx < len(gt_fields):
                            sample_name = final_samples[i] if target_samples else samples[sample_idx]
                            gt_data = gt_fields[sample_idx]
                            genotypes[sample_name] = GenotypeUtils.parse_genotype(gt_data)
                
                yield VCFRecord(
                    chrom=chrom,
                    pos=int(pos),
                    id_field=id_field,
                    ref=ref,
                    alt=alt,
                    qual=qual,
                    genotypes=genotypes
                )

class VCFParserFactory:
    """VCF解析器工厂 | VCF Parser Factory"""
    
    @staticmethod
    def create_parser(vcf_file: str, logger, use_fast: bool = True):
        """创建VCF解析器 | Create VCF parser"""
        if use_fast:
            try:
                return FastVCFParser(vcf_file, logger)
            except ImportError:
                logger.info("cyvcf2未安装，回退到标准解析器 | cyvcf2 not installed, falling back to standard parser")
        
        return StandardVCFParser(vcf_file, logger)
