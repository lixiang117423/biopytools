"""
VCF数据处理模块 | VCF Data Processing Module
"""

import gzip
from typing import List, Tuple, Dict
from collections import OrderedDict

class GenotypeProcessor:
    """基因型处理器 | Genotype Processor"""
    
    @staticmethod
    def convert_gt_to_bases(gt: str, ref: str, alt: str) -> str:
        """
        将基因型转换为碱基组合 | Convert genotype to base combination
        
        Args:
            gt: 基因型字符串 (如 0/1, 1|1) | Genotype string
            ref: 参考等位基因 | Reference allele
            alt: 替代等位基因 | Alternative allele
            
        Returns:
            碱基组合字符串 (如 CA, AA) | Base combination string
        """
        if gt in ('./.', '.|.', '.'):
            return '--'
        
        # 处理多个替代等位基因的情况
        alts = alt.split(',') if ',' in alt else [alt]
        
        # 分离基因型
        separator = '|' if '|' in gt else '/'
        allele_codes = gt.split(separator)
        
        bases = []
        for code in allele_codes:
            if code == '.':
                return '--'
            try:
                idx = int(code)
                if idx == 0:
                    bases.append(ref)
                else:
                    bases.append(alts[idx-1])
            except (IndexError, ValueError):
                return '--'
        
        return ''.join(bases)

class VCFProcessor:
    """VCF文件处理器 | VCF File Processor"""
    
    def __init__(self, logger, command_runner, temp_manager):
        self.logger = logger
        self.command_runner = command_runner
        self.temp_manager = temp_manager
        self.genotype_processor = GenotypeProcessor()
    
    def create_regions_file(self, positions: List[Tuple[str, int]]) -> str:
        """
        创建区域文件用于bcftools查询 | Create regions file for bcftools query
        
        Returns:
            区域文件路径 | Path to regions file
        """
        regions_file = self.temp_manager.create_temp_file(".bed")
        
        with open(regions_file, 'w') as f:
            for chrom, pos in positions:
                # BED格式是0-based，所以pos-1
                f.write(f"{chrom}\t{pos-1}\t{pos}\n")
        
        self.logger.info(f"创建区域文件 | Created regions file: {regions_file}")
        return regions_file
    
    def extract_variants(self, config, positions: List[Tuple[str, int]]) -> str:
        """
        从VCF文件提取指定位点的变异 | Extract variants from VCF file at specified positions
        
        Returns:
            提取后的VCF文件路径 | Path to extracted VCF file
        """
        regions_file = self.create_regions_file(positions)
        extracted_vcf = self.temp_manager.create_temp_file(".vcf")
        
        cmd = [
            config.bcftools_path, 'view',
            '-R', regions_file,
            '-o', extracted_vcf,
            config.vcf_file
        ]
        
        success, stdout, stderr = self.command_runner.run(
            cmd, 
            f"提取 {len(positions)} 个位点的变异信息"
        )
        
        if not success:
            raise RuntimeError(f"VCF提取失败 | VCF extraction failed: {stderr}")
        
        return extracted_vcf
    
    def parse_vcf_to_matrix(self, vcf_path: str, target_positions: List[Tuple[str, int]]) -> Tuple[List[str], Dict]:
        """
        解析VCF文件并转换为矩阵格式 | Parse VCF file and convert to matrix format
        
        Returns:
            (样本列表, 变异数据字典) | (sample list, variant data dictionary)
        """
        self.logger.info("解析VCF文件为矩阵格式 | Parsing VCF file to matrix format")
        
        samples = []
        variant_data = OrderedDict()
        target_positions_set = set(target_positions)
        
        # 打开VCF文件
        open_func = gzip.open if vcf_path.endswith('.gz') else open
        mode = 'rt' if vcf_path.endswith('.gz') else 'r'
        
        with open_func(vcf_path, mode) as f:
            for line in f:
                if line.startswith('#CHROM'):
                    # 获取样本名列表
                    fields = line.strip().split('\t')
                    samples = fields[9:]  # 跳过前9列
                    self.logger.info(f"发现 {len(samples)} 个样本 | Found {len(samples)} samples")
                    continue
                elif line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 10:  # 至少需要有一个样本
                    continue
                
                chrom = fields[0]
                pos = int(fields[1])
                ref = fields[3]
                alt = fields[4]
                
                # 检查是否是目标位点
                if (chrom, pos) not in target_positions_set:
                    continue
                
                # 提取基因型信息
                genotypes = []
                for gt_info in fields[9:]:
                    gt = gt_info.split(':')[0]  # 提取基因型部分
                    base_combo = self.genotype_processor.convert_gt_to_bases(gt, ref, alt)
                    genotypes.append(base_combo)
                
                # 存储变异信息
                variant_key = (chrom, pos)
                variant_data[variant_key] = {
                    'ref': ref,
                    'alt': alt,
                    'genotypes': genotypes
                }
                
                self.logger.debug(f"处理变异 | Processed variant: {chrom}:{pos} {ref}>{alt}")
        
        self.logger.info(f"成功解析 {len(variant_data)} 个变异位点 | Successfully parsed {len(variant_data)} variant sites")
        return samples, variant_data
