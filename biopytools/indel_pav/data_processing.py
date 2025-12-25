"""
VCF数据处理模块 | VCF Data Processing Module
"""

from typing import List, Dict, Tuple, Generator
from .utils import FileHandler
import re

class VCFProcessor:
    """VCF文件处理器 | VCF File Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.samples = []
        self.sample_count = 0
    
    def parse_vcf_header(self) -> List[str]:
        """解析VCF头部信息 | Parse VCF header"""
        self.logger.info("📋 解析VCF头部信息 | Parsing VCF header")
        
        with FileHandler.open_file(self.config.vcf_file) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    # 提取样本名称
                    columns = line.strip().split('\t')
                    self.samples = columns[9:]  # 从第10列开始是样本
                    self.sample_count = len(self.samples)
                    break
        
        if not self.samples:
            raise ValueError("❌ 无法从VCF文件中提取样本信息 | Cannot extract sample information from VCF file")
        
        self.logger.info(f"🧬 找到 {self.sample_count} 个样本 | Found {self.sample_count} samples")
        self.logger.info(f"📊 样本列表 | Sample list: {', '.join(self.samples[:5])}{'...' if len(self.samples) > 5 else ''}")
        
        return self.samples
    
    def is_indel(self, ref: str, alt: str) -> Tuple[bool, str, int]:
        """判断是否为INDEL并返回类型和长度 | Check if variant is INDEL and return type and length"""
        # 排除复杂变异 | Exclude complex variants
        if ',' in alt or alt == '*':
            return False, "", 0
        
        ref_len = len(ref)
        alt_len = len(alt)
        
        # 插入 | Insertion
        if ref_len < alt_len:
            indel_len = alt_len - ref_len
            return True, "insertion", indel_len
        
        # 删除 | Deletion
        elif ref_len > alt_len:
            indel_len = ref_len - alt_len
            return True, "deletion", indel_len
        
        # SNP或其他 | SNP or others
        return False, "", 0
    
    def extract_indel_sequence(self, ref: str, alt: str, variant_type: str) -> str:
        """提取INDEL序列 | Extract INDEL sequence"""
        if variant_type == "insertion":
            # 对于插入，返回插入的序列（去掉共同的前缀）
            common_prefix = 0
            min_len = min(len(ref), len(alt))
            for i in range(min_len):
                if ref[i] == alt[i]:
                    common_prefix += 1
                else:
                    break
            return alt[common_prefix:]
        
        elif variant_type == "deletion":
            # 对于删除，返回删除的序列（去掉共同的前缀）
            common_prefix = 0
            min_len = min(len(ref), len(alt))
            for i in range(min_len):
                if ref[i] == alt[i]:
                    common_prefix += 1
                else:
                    break
            return ref[common_prefix:]
        
        return ""
    
    def parse_genotype(self, gt_info: str) -> int:
        """解析基因型信息 | Parse genotype information"""
        # 提取GT字段（通常是第一个字段）
        gt = gt_info.split(':')[0]
        
        # 处理缺失数据
        if gt in ['./.', '.|.', '.']:
            return 0  # 缺失视为不存在
        
        # 解析基因型
        if '|' in gt:
            alleles = gt.split('|')
        elif '/' in gt:
            alleles = gt.split('/')
        else:
            return 0
        
        # 如果任一等位基因为非参考等位基因（非0），则认为存在该INDEL
        try:
            for allele in alleles:
                if allele != '.' and int(allele) > 0:
                    return 1
            return 0
        except ValueError:
            return 0
    def parse_genotype_for_alt(self, gt_info: str, target_alt: int) -> int:
        """解析针对特定ALT等位基因的基因型 | Parse genotype for specific ALT allele"""
        # 提取GT字段（通常是第一个字段）
        gt = gt_info.split(':')[0]
        
        # 处理缺失数据
        if gt in ['./.', '.|.', '.']:
            return 0  # 缺失视为不存在
        
        # 解析基因型
        if '|' in gt:
            alleles = gt.split('|')
        elif '/' in gt:
            alleles = gt.split('/')
        else:
            return 0
        
        # 检查是否包含目标ALT等位基因
        try:
            for allele in alleles:
                if allele != '.' and int(allele) == target_alt:
                    return 1
            return 0
        except ValueError:
            return 0
    
    def parse_genotype_for_alt(self, gt_info: str, target_alt: int) -> int:
        """解析针对特定ALT等位基因的基因型 | Parse genotype for specific ALT allele"""
        # 提取GT字段（通常是第一个字段）
        gt = gt_info.split(':')[0]
        
        # 处理缺失数据
        if gt in ['./.', '.|.', '.']:
            return 0  # 缺失视为不存在
        
        # 解析基因型
        if '|' in gt:
            alleles = gt.split('|')
        elif '/' in gt:
            alleles = gt.split('/')
        else:
            return 0
        
        # 检查是否包含目标ALT等位基因
        try:
            for allele in alleles:
                if allele != '.' and int(allele) == target_alt:
                    return 1
            return 0
        except ValueError:
            return 0
            
    def filter_by_quality(self, qual: str, info: str) -> bool:
        """根据质量过滤 | Filter by quality"""
        try:
            # 检查质量分数
            if qual != '.' and float(qual) < self.config.min_quality:
                return False
            
            # 检查深度信息（如果有DP字段）
            if 'DP=' in info:
                dp_match = re.search(r'DP=(\d+)', info)
                if dp_match:
                    depth = int(dp_match.group(1))
                    if depth < self.config.min_depth:
                        return False
            
            return True
        except:
            return False
    
    def process_vcf_variants(self) -> Generator[Dict, None, None]:
        """处理VCF变异信息 | Process VCF variants"""
        self.logger.info(f"🔬 开始处理VCF变异，线程数: {self.config.threads} | Starting VCF variant processing, threads: {self.config.threads}")
        
        variant_count = 0
        indel_count = 0
        filtered_count = 0
        
        with FileHandler.open_file(self.config.vcf_file) as f:
            for line_num, line in enumerate(f, 1):
                # 跳过注释行
                if line.startswith('#'):
                    continue
                
                variant_count += 1
                if variant_count % 10000 == 0:
                    self.logger.info(f"📈 已处理 {variant_count} 个变异，发现 {indel_count} 个INDEL | Processed {variant_count} variants, found {indel_count} INDELs")
                
                fields = line.strip().split('\t')
                if len(fields) < 9 + self.sample_count:
                    continue
                
                chrom, pos, var_id, ref, alt, qual, filter_field, info = fields[:8]
                sample_gts = fields[9:]
                
                # 处理多重突变 - 分割多个ALT等位基因
                alt_alleles = alt.split(',')
                
                # 对每个ALT等位基因单独处理
                for alt_idx, single_alt in enumerate(alt_alleles):
                    # 检查是否为INDEL
                    is_indel_var, variant_type, indel_length = self.is_indel(ref, single_alt)
                    if not is_indel_var:
                        continue
                    
                    # 长度过滤
                    if indel_length < self.config.min_length:
                        continue
                    
                    if self.config.max_length is not None and indel_length > self.config.max_length:
                        continue
                    
                    # 质量过滤
                    if not self.filter_by_quality(qual, info):
                        filtered_count += 1
                        continue
                    
                    # 提取INDEL序列
                    indel_sequence = self.extract_indel_sequence(ref, single_alt, variant_type)
                    
                    # 解析样本基因型 - 针对特定ALT等位基因
                    genotypes = []
                    missing_count = 0
                    
                    for gt_info in sample_gts:
                        gt = self.parse_genotype_for_alt(gt_info, alt_idx + 1)  # ALT索引从1开始
                        genotypes.append(gt)
                        if gt == 0:
                            missing_count += 1
                    
                    # 缺失率过滤
                    missing_rate = missing_count / self.sample_count
                    if missing_rate > self.config.max_missing_rate:
                        filtered_count += 1
                        continue
                    
                    indel_count += 1
                    
                    # 计算终止位置
                    end_pos = int(pos) + len(ref) - 1
                    
                    yield {
                        'chrom': chrom,
                        'start': pos,
                        'end': str(end_pos),
                        'ref': ref,  # 确保包含REF字段
                        'sequence': indel_sequence,
                        'type': variant_type,
                        'genotypes': genotypes,
                        'length': indel_length,
                        'missing_rate': missing_rate
                    }
        
        self.logger.info(f"✅ VCF处理完成 | VCF processing completed")
        self.logger.info(f"📊 总变异数: {variant_count} | Total variants: {variant_count}")
        self.logger.info(f"🧬 INDEL数量: {indel_count} | INDEL count: {indel_count}")
        self.logger.info(f"🚫 过滤掉: {filtered_count} | Filtered out: {filtered_count}")
