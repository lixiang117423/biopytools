"""
序列数据处理模块 | Sequence Data Processing Module
"""

import pysam
from typing import Dict, List, Set, Optional, Tuple
from .utils import SequenceValidator

class GenomeProcessor:
    """基因组处理器 | Genome Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.genome_reader = None
    
    def open_genome(self):
        """打开基因组文件 | Open genome file"""
        try:
            self.genome_reader = pysam.FastaFile(self.config.genome_file)
            self.logger.info(f"基因组文件已打开 | Genome file opened: {self.config.genome_file}")
            return True
        except Exception as e:
            self.logger.error(f"打开基因组文件失败 | Failed to open genome file: {e}")
            return False
    
    def close_genome(self):
        """关闭基因组文件 | Close genome file"""
        if self.genome_reader:
            self.genome_reader.close()
            self.genome_reader = None
    
    def get_reference_sequence(self) -> str:
        """获取参考序列 | Get reference sequence"""
        try:
            # pysam使用0-based坐标 | pysam uses 0-based coordinates
            sequence = self.genome_reader.fetch(
                self.config.chrom, 
                self.config.start, 
                self.config.end + 1
            )
            return sequence.upper()
        except Exception as e:
            self.logger.error(f"获取参考序列失败 | Failed to get reference sequence: {e}")
            return ""

class VCFProcessor:
    """VCF处理器 | VCF Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.vcf_reader = None
        self.validator = SequenceValidator(logger)
    
    def open_vcf(self):
        """打开VCF文件 | Open VCF file"""
        try:
            self.vcf_reader = pysam.VariantFile(self.config.vcf_file)
            self.logger.info(f"VCF文件已打开 | VCF file opened: {self.config.vcf_file}")
            return True
        except Exception as e:
            self.logger.error(f"打开VCF文件失败 | Failed to open VCF file: {e}")
            return False
    
    def close_vcf(self):
        """关闭VCF文件 | Close VCF file"""
        if self.vcf_reader:
            self.vcf_reader.close()
            self.vcf_reader = None
    
    def get_sample_names(self) -> List[str]:
        """获取样品名称列表 | Get sample names list"""
        if not self.vcf_reader:
            return []
        
        sample_names = list(self.vcf_reader.header.samples)
        
        # 应用样品过滤 | Apply sample filtering
        if self.config.sample_list:
            sample_names = [s for s in sample_names if s in self.config.sample_list]
            self.logger.info(f"使用指定样品列表，剩余样品数 | Using specified sample list, remaining samples: {len(sample_names)}")
        
        if self.config.exclude_samples:
            sample_names = [s for s in sample_names if s not in self.config.exclude_samples]
            self.logger.info(f"排除指定样品后，剩余样品数 | After excluding samples, remaining samples: {len(sample_names)}")
        
        return sample_names
    
    def get_variants_in_region(self) -> Dict:
        """获取区间内的变异信息 | Get variants in region"""
        variants = {}
        variant_count = 0
        
        try:
            # 遍历区间内的变异 | Iterate through variants in region
            for record in self.vcf_reader.fetch(
                self.config.chrom, 
                self.config.start, 
                self.config.end + 1
            ):
                pos = record.pos  # VCF位置是1-based | VCF position is 1-based
                
                if self.config.start <= pos <= self.config.end:
                    # 质量过滤 | Quality filtering
                    if self.config.min_qual and record.qual and record.qual < self.config.min_qual:
                        continue
                    
                    variants[pos] = {
                        'ref': record.ref,
                        'alt': [str(alt) for alt in record.alts] if record.alts else [],
                        'qual': record.qual,
                        'samples': {}
                    }
                    
                    # 获取每个样品的基因型 | Get genotype for each sample
                    for sample_name in self.get_sample_names():
                        if sample_name in record.samples:
                            sample = record.samples[sample_name]
                            gt = sample.get('GT', (None, None))
                            
                            if gt is not None and None not in gt:
                                # 处理基因型 | Process genotype
                                alleles = []
                                for allele_idx in gt:
                                    if allele_idx == 0:
                                        alleles.append(record.ref)
                                    elif allele_idx > 0 and allele_idx <= len(record.alts):
                                        alleles.append(str(record.alts[allele_idx-1]))
                                    else:
                                        alleles.append('N')
                                variants[pos]['samples'][sample_name] = alleles
                            else:
                                # 缺失基因型 | Missing genotype
                                variants[pos]['samples'][sample_name] = ['N', 'N']
                    
                    variant_count += 1
            
            self.logger.info(f"在区间内找到 {variant_count} 个变异 | Found {variant_count} variants in region")
            return variants
            
        except Exception as e:
            self.logger.error(f"获取变异信息失败 | Failed to get variants: {e}")
            return {}

class SequenceBuilder:
    """序列构建器 | Sequence Builder"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.validator = SequenceValidator(logger)
    
    def build_sample_sequences(self, reference_seq: str, variants: Dict) -> Dict[str, str]:
        """构建每个样品的序列 | Build sequences for each sample"""
        if not reference_seq:
            self.logger.error("参考序列为空 | Reference sequence is empty")
            return {}
        
        # 获取所有样品名称 | Get all sample names
        sample_names = set()
        for var_info in variants.values():
            sample_names.update(var_info['samples'].keys())
        
        if not sample_names:
            self.logger.warning("未找到样品信息 | No sample information found")
            return {}
        
        sample_sequences = {}
        
        for sample_name in sample_names:
            # 从参考序列开始 | Start with reference sequence
            sequence_list = list(reference_seq)
            variant_applied = 0
            
            # 应用变异 | Apply variants
            for pos, var_info in variants.items():
                if sample_name in var_info['samples']:
                    alleles = var_info['samples'][sample_name]
                    
                    # 选择使用的等位基因 | Select allele to use
                    if self.config.use_first_allele:
                        allele = alleles[0] if alleles else 'N'
                    else:
                        allele = alleles[1] if len(alleles) > 1 else (alleles[0] if alleles else 'N')
                    
                    # 计算相对位置 | Calculate relative position
                    relative_pos = pos - self.config.start
                    
                    # 只处理简单的SNP替换 | Only handle simple SNP substitutions
                    if (len(var_info['ref']) == 1 and len(allele) == 1 and 
                        0 <= relative_pos < len(sequence_list)):
                        sequence_list[relative_pos] = allele
                        variant_applied += 1
            
            sequence = ''.join(sequence_list)
            
            # 验证序列 | Validate sequence
            if self.validator.validate_sequence(sequence, sample_name):
                sample_sequences[sample_name] = sequence
                self.logger.debug(f"样品 {sample_name} 应用了 {variant_applied} 个变异 | Sample {sample_name} applied {variant_applied} variants")
        
        self.logger.info(f"成功构建 {len(sample_sequences)} 个样品的序列 | Successfully built sequences for {len(sample_sequences)} samples")
        return sample_sequences
