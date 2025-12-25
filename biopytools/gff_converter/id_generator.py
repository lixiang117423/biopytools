"""
基因ID生成器模块 | Gene ID Generator Module
"""

from typing import Dict, List, Set
from .utils import ChromosomeExtractor, format_gene_id
from .gff_parser import GFFFeature

class IDGenerator:
    """基因ID生成器 | Gene ID Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.gene_id_mapping = {}  # 原始ID -> 新ID的映射
        self.used_ids = set()      # 已使用的ID集合
    
    def generate_new_ids(self, features: List[GFFFeature]) -> Dict[str, str]:
        """生成新的基因ID | Generate new gene IDs"""
        self.logger.info("🆔 开始生成新的基因ID | Starting to generate new gene IDs")
        
        # 按染色体分组基因 | Group genes by chromosome
        genes_by_chr = self._group_genes_by_chromosome(features)
        
        # 为每个染色体的基因重新编号 | Renumber genes for each chromosome
        for chrom, genes in genes_by_chr.items():
            self._renumber_chromosome_genes(chrom, genes)
        
        self.logger.info(f"✅ 基因ID生成完成 | Gene ID generation completed")
        self.logger.info(f"📊 ID统计 | ID statistics:")
        self.logger.info(f"   - 重新编号基因数 | Renumbered genes: {len(self.gene_id_mapping)}")
        self.logger.info(f"   - 涉及染色体数 | Chromosomes involved: {len(genes_by_chr)}")
        
        return self.gene_id_mapping
    
    def _group_genes_by_chromosome(self, features: List[GFFFeature]) -> Dict[str, List[GFFFeature]]:
        """按染色体分组基因 | Group genes by chromosome"""
        genes_by_chr = {}
        
        for feature in features:
            if feature.type == 'gene':
                chrom = feature.seqid
                if chrom not in genes_by_chr:
                    genes_by_chr[chrom] = []
                genes_by_chr[chrom].append(feature)
        
        # 按起始位置排序 | Sort by start position
        for chrom in genes_by_chr:
            genes_by_chr[chrom].sort(key=lambda x: (x.start, x.end))
            self.logger.info(f"   🧬 {chrom}: {len(genes_by_chr[chrom])} genes")
        
        return genes_by_chr
    
    def _renumber_chromosome_genes(self, chrom: str, genes: List[GFFFeature]):
        """为单个染色体的基因重新编号 | Renumber genes for a single chromosome"""
        # 提取染色体编号 | Extract chromosome number
        chr_num = ChromosomeExtractor.extract_chr_number(chrom)
        
        self.logger.info(f"🔢 为染色体 {chrom} (编号: {chr_num}) 重新编号基因")
        
        current_num = self.config.start_number
        
        for gene in genes:
            # 获取原始基因ID | Get original gene ID
            original_id = gene.get_attribute('ID')
            if not original_id:
                continue
            
            # 去除可能的前缀 | Remove possible prefix
            original_id = original_id.replace('gene-', '').replace('rna-', '')
            
            # 生成新ID | Generate new ID
            new_id = format_gene_id(
                self.config.species_name,
                self.config.species_prefix,
                chr_num,
                current_num
            )
            
            # 确保ID唯一 | Ensure ID uniqueness
            while new_id in self.used_ids:
                current_num += self.config.step
                new_id = format_gene_id(
                    self.config.species_name,
                    self.config.species_prefix,
                    chr_num,
                    current_num
                )
            
            # 记录映射 | Record mapping
            self.gene_id_mapping[original_id] = new_id
            self.used_ids.add(new_id)
            
            current_num += self.config.step
        
        self.logger.info(f"   ✅ {chrom}: 重新编号 {len(genes)} 个基因")
    
    def get_mapped_id(self, original_id: str) -> str:
        """获取映射后的ID | Get mapped ID"""
        # 清理原始ID | Clean original ID
        clean_id = original_id.replace('gene-', '').replace('rna-', '')
        return self.gene_id_mapping.get(clean_id, original_id)

class HierarchicalIDGenerator:
    """层级ID生成器 | Hierarchical ID Generator"""
    
    def __init__(self, id_generator: IDGenerator):
        self.id_generator = id_generator
    
    def generate_child_ids(self, feature: GFFFeature, parent_new_id: str) -> str:
        """生成子特征ID | Generate child feature IDs"""
        if feature.type in ['mRNA', 'transcript', 'lncRNA', 'tRNA', 'rRNA']:
            # 获取转录本编号 | Get transcript number
            transcript_num = self._extract_transcript_number(feature)
            return f"{parent_new_id}.{feature.type}{transcript_num}"
        
        elif feature.type == 'exon':
            # 获取外显子编号 | Get exon number
            exon_num = self._extract_exon_number(feature)
            parent_id = feature.get_attribute('Parent', '')
            if parent_id:
                # 使用父ID作为基础 | Use parent ID as base
                return f"{parent_id}.exon{exon_num}"
        
        elif feature.type in ['CDS', 'five_prime_UTR', 'three_prime_UTR']:
            # CDS和UTR特征 | CDS and UTR features
            parent_id = feature.get_attribute('Parent', '')
            if parent_id:
                feature_num = self._extract_feature_number(feature)
                feature_type_short = {
                    'CDS': 'cds',
                    'five_prime_UTR': '5utr',
                    'three_prime_UTR': '3utr'
                }.get(feature.type, feature.type.lower())
                return f"{parent_id}.{feature_type_short}{feature_num}"
        
        # 默认返回原始ID | Default return original ID
        return feature.get_attribute('ID', '')
    
    def _extract_transcript_number(self, feature: GFFFeature) -> int:
        """提取转录本编号 | Extract transcript number"""
        original_id = feature.get_attribute('ID', '')
        
        # 尝试从ID中提取编号 | Try to extract number from ID
        import re
        match = re.search(r'-R(\d+)$', original_id)
        if match:
            return int(match.group(1))
        
        return 1  # 默认为1
    
    def _extract_exon_number(self, feature: GFFFeature) -> int:
        """提取外显子编号 | Extract exon number"""
        original_id = feature.get_attribute('ID', '')
        
        # 尝试从ID中提取编号 | Try to extract number from ID
        import re
        match = re.search(r'-(\d+)$', original_id)
        if match:
            return int(match.group(1))
        
        return 1  # 默认为1
    
    def _extract_feature_number(self, feature: GFFFeature) -> int:
        """提取特征编号 | Extract feature number"""
        # 简化处理，返回1 | Simplified processing, return 1
        return 1
