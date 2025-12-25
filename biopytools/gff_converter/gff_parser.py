"""
GFF文件解析模块 | GFF File Parsing Module
"""

import re
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from .utils import AttributeParser, ChromosomeExtractor, GFFValidator

@dataclass
class GFFFeature:
    """GFF特征类 | GFF Feature Class"""
    seqid: str          # 序列ID | Sequence ID
    source: str         # 来源 | Source
    type: str           # 特征类型 | Feature type
    start: int          # 起始位置 | Start position
    end: int            # 结束位置 | End position
    score: str          # 得分 | Score
    strand: str         # 链 | Strand
    phase: str          # 相位 | Phase
    attributes: Dict[str, str]  # 属性 | Attributes
    
    def get_attribute(self, key: str, default: str = None) -> str:
        """获取属性值 | Get attribute value"""
        return self.attributes.get(key, default)
    
    def set_attribute(self, key: str, value: str):
        """设置属性值 | Set attribute value"""
        self.attributes[key] = value
    
    def to_gff_line(self) -> str:
        """转换为GFF行 | Convert to GFF line"""
        attr_str = AttributeParser.format_attributes(self.attributes)
        return f"{self.seqid}\t{self.source}\t{self.type}\t{self.start}\t{self.end}\t{self.score}\t{self.strand}\t{self.phase}\t{attr_str}"

class GFFParser:
    """GFF文件解析器 | GFF File Parser"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.features = []
        self.header_lines = []
    
    def parse_file(self) -> List[GFFFeature]:
        """解析GFF文件 | Parse GFF file"""
        self.logger.info(f"🧬 开始解析GFF文件 | Starting to parse GFF file: {self.config.input_file}")
        
        with open(self.config.input_file, 'r') as f:
            line_count = 0
            for line in f:
                line_count += 1
                line = line.strip()
                
                if not line:
                    continue
                
                # 处理头部注释行 | Process header comment lines
                if line.startswith('#'):
                    self.header_lines.append(line)
                    continue
                
                # 验证GFF行格式 | Validate GFF line format
                if not GFFValidator.is_valid_gff_line(line):
                    self.logger.warning(f"⚠️  跳过无效行 | Skipping invalid line {line_count}: {line}")
                    continue
                
                # 解析特征行 | Parse feature line
                try:
                    feature = self._parse_feature_line(line)
                    if feature:
                        self.features.append(feature)
                except Exception as e:
                    self.logger.error(f"❌ 解析第{line_count}行时出错 | Error parsing line {line_count}: {e}")
        
        self.logger.info(f"✅ GFF文件解析完成 | GFF file parsing completed")
        self.logger.info(f"📊 解析统计 | Parsing statistics:")
        self.logger.info(f"   - 头部行数 | Header lines: {len(self.header_lines)}")
        self.logger.info(f"   - 特征总数 | Total features: {len(self.features)}")
        self._log_feature_stats()
        
        return self.features
    
    def _parse_feature_line(self, line: str) -> Optional[GFFFeature]:
        """解析单行特征 | Parse single feature line"""
        fields = line.split('\t')
        
        if len(fields) != 9:
            return None
        
        # 验证坐标 | Validate coordinates
        try:
            start = int(fields[3])
            end = int(fields[4])
            if not GFFValidator.validate_coordinates(start, end):
                self.logger.warning(f"⚠️  无效坐标 | Invalid coordinates: {start}-{end}")
                return None
        except ValueError:
            self.logger.warning(f"⚠️  坐标格式错误 | Invalid coordinate format: {fields[3]}-{fields[4]}")
            return None
        
        # 解析属性 | Parse attributes
        attributes = AttributeParser.parse_attributes(fields[8])
        
        return GFFFeature(
            seqid=fields[0],
            source=fields[1],
            type=fields[2],
            start=start,
            end=end,
            score=fields[5],
            strand=fields[6],
            phase=fields[7],
            attributes=attributes
        )
    
    def _log_feature_stats(self):
        """记录特征统计信息 | Log feature statistics"""
        type_counts = {}
        chromosome_counts = {}
        
        for feature in self.features:
            # 统计特征类型 | Count feature types
            ftype = feature.type
            type_counts[ftype] = type_counts.get(ftype, 0) + 1
            
            # 统计染色体 | Count chromosomes
            chrom = feature.seqid
            chromosome_counts[chrom] = chromosome_counts.get(chrom, 0) + 1
        
        # 输出统计信息 | Output statistics
        self.logger.info("   📈 特征类型分布 | Feature type distribution:")
        for ftype, count in sorted(type_counts.items()):
            self.logger.info(f"      {ftype}: {count}")
        
        self.logger.info("   🧭 染色体分布 | Chromosome distribution:")
        for chrom, count in sorted(chromosome_counts.items()):
            self.logger.info(f"      {chrom}: {count}")
    
    def get_genes_by_chromosome(self) -> Dict[str, List[GFFFeature]]:
        """按染色体分组基因 | Group genes by chromosome"""
        genes_by_chr = {}
        
        for feature in self.features:
            if feature.type == 'gene':
                chrom = feature.seqid
                if chrom not in genes_by_chr:
                    genes_by_chr[chrom] = []
                genes_by_chr[chrom].append(feature)
        
        # 按位置排序 | Sort by position
        for chrom in genes_by_chr:
            genes_by_chr[chrom].sort(key=lambda x: (x.start, x.end))
        
        return genes_by_chr
