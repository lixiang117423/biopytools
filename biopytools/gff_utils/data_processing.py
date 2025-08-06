"""
GFF3数据处理模块 | GFF3 Data Processing Module
"""

import csv
from typing import Dict, List, Tuple, Any
from .utils import AttributeParser, GFFValidator

class GFFParser:
    """GFF3文件解析器 | GFF3 File Parser"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.attribute_parser = AttributeParser(logger)
        self.validator = GFFValidator(logger)
    
    def collect_gene_data(self) -> Dict[str, Dict[str, str]]:
        """
        第一遍扫描：收集所有基因信息 | First pass: collect all gene information
        
        Returns:
            Dict[str, Dict[str, str]]: 基因ID到基因信息的映射
        """
        gene_data = {}
        
        self.logger.info("开始收集基因信息 | Starting to collect gene information")
        
        try:
            with open(self.config.gff3_file, 'r', encoding='utf-8') as infile:
                for line_num, line in enumerate(infile, 1):
                    # 跳过注释行 | Skip comment lines
                    if line.startswith('#'):
                        continue
                    
                    # 验证行格式 | Validate line format
                    if not self.validator.validate_gff_line(line):
                        continue
                    
                    parts = line.strip().split('\t')
                    feature_type = parts[2]
                    
                    # 只处理基因行 | Only process gene lines
                    if feature_type == self.config.gene_type:
                        attributes = self.attribute_parser.parse_attributes(parts[8])
                        gene_id = attributes.get('ID')
                        
                        if gene_id:
                            gene_data[gene_id] = {
                                'chr': parts[0],
                                'start': parts[3],
                                'end': parts[4],
                                'strand': parts[6],
                                'line_num': line_num
                            }
                        else:
                            self.logger.warning(f"基因行缺少ID属性 | Gene line missing ID attribute at line {line_num}")
        
        except FileNotFoundError:
            self.logger.error(f"输入文件未找到 | Input file not found: {self.config.gff3_file}")
            raise
        except Exception as e:
            self.logger.error(f"收集基因信息时发生错误 | Error collecting gene information: {e}")
            raise
        
        self.logger.info(f"收集完成，共发现 {len(gene_data)} 个基因 | Collection completed, found {len(gene_data)} genes")
        return gene_data
    
    def process_transcripts(self, gene_data: Dict[str, Dict[str, str]]) -> List[List[str]]:
        """
        第二遍扫描：处理转录本信息 | Second pass: process transcript information
        
        Args:
            gene_data: 基因信息字典
            
        Returns:
            List[List[str]]: 转录本信息列表
        """
        transcript_data = []
        orphan_transcripts = []
        
        self.logger.info("开始处理转录本信息 | Starting to process transcript information")
        
        try:
            with open(self.config.gff3_file, 'r', encoding='utf-8') as infile:
                for line_num, line in enumerate(infile, 1):
                    # 跳过注释行 | Skip comment lines
                    if line.startswith('#'):
                        continue
                    
                    # 验证行格式 | Validate line format
                    if not self.validator.validate_gff_line(line):
                        continue
                    
                    parts = line.strip().split('\t')
                    feature_type = parts[2]
                    
                    # 只处理转录本行 | Only process transcript lines
                    if feature_type in self.config.transcript_types:
                        attributes = self.attribute_parser.parse_attributes(parts[8])
                        transcript_id = attributes.get('ID')
                        gene_id = attributes.get('Parent')
                        
                        if not transcript_id:
                            self.logger.warning(f"转录本行缺少ID属性 | Transcript line missing ID attribute at line {line_num}")
                            continue
                        
                        if not gene_id:
                            self.logger.warning(f"转录本行缺少Parent属性 | Transcript line missing Parent attribute at line {line_num}")
                            continue
                        
                        # 获取转录本坐标 | Get transcript coordinates
                        transcript_start = parts[3]
                        transcript_end = parts[4]
                        
                        # 查找父基因信息 | Find parent gene information
                        parent_gene_info = gene_data.get(gene_id)
                        
                        if parent_gene_info:
                            # 找到父基因 | Found parent gene
                            chromosome = parent_gene_info['chr']
                            strand = parent_gene_info['strand']
                            gene_start = parent_gene_info['start']
                            gene_end = parent_gene_info['end']
                        else:
                            # 孤儿转录本 | Orphan transcript
                            orphan_transcripts.append(transcript_id)
                            chromosome = parts[0]
                            strand = parts[6]
                            gene_start = 'NA'
                            gene_end = 'NA'
                        
                        # 创建输出行 | Create output row
                        output_row = [
                            gene_id,
                            transcript_id,
                            chromosome,
                            strand,
                            gene_start,
                            gene_end,
                            transcript_start,
                            transcript_end
                        ]
                        transcript_data.append(output_row)
        
        except Exception as e:
            self.logger.error(f"处理转录本信息时发生错误 | Error processing transcript information: {e}")
            raise
        
        # 报告孤儿转录本 | Report orphan transcripts
        if orphan_transcripts:
            self.logger.warning(f"发现 {len(orphan_transcripts)} 个孤儿转录本 | Found {len(orphan_transcripts)} orphan transcripts")
            self.logger.warning(f"孤儿转录本列表 | Orphan transcript list: {', '.join(orphan_transcripts[:10])}")
            if len(orphan_transcripts) > 10:
                self.logger.warning(f"还有 {len(orphan_transcripts) - 10} 个孤儿转录本未列出 | {len(orphan_transcripts) - 10} more orphan transcripts not listed")
        
        self.logger.info(f"处理完成，共处理 {len(transcript_data)} 个转录本 | Processing completed, processed {len(transcript_data)} transcripts")
        return transcript_data
