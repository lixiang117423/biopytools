"""
GFF3基因转录本提取核心模块 | GFF3 Gene Transcript Extraction Core Module
"""

import csv
from typing import List, Dict, Any
from .data_processing import GFFParser

class GFFExtractor:
    """GFF3基因转录本提取器 | GFF3 Gene Transcript Extractor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.gff_parser = GFFParser(config, logger)
    
    def extract_gene_transcript_info(self) -> List[List[str]]:
        """
        提取基因和转录本信息 | Extract gene and transcript information
        
        Returns:
            List[List[str]]: 转录本信息列表
        """
        self.logger.info("开始提取基因和转录本信息 | Starting to extract gene and transcript information")
        self.logger.info(f"输入文件 | Input file: {self.config.gff3_file}")
        self.logger.info(f"基因类型 | Gene type: {self.config.gene_type}")
        self.logger.info(f"转录本类型 | Transcript types: {', '.join(self.config.transcript_types)}")
        
        # 第一遍：收集基因信息 | First pass: collect gene information
        gene_data = self.gff_parser.collect_gene_data()
        
        # 第二遍：处理转录本信息 | Second pass: process transcript information
        transcript_data = self.gff_parser.process_transcripts(gene_data)
        
        self.logger.info("提取完成 | Extraction completed")
        return transcript_data
    
    def write_results(self, transcript_data: List[List[str]]):
        """
        写入结果文件 | Write results file
        
        Args:
            transcript_data: 转录本信息列表
        """
        self.logger.info(f"写入结果文件 | Writing results file: {self.config.output_file}")
        
        try:
            with open(self.config.output_file, 'w', newline='', encoding='utf-8') as outfile:
                writer = csv.writer(outfile, delimiter='\t')
                
                # 写入头部 | Write header
                header = [
                    'Gene_ID', 'Transcript_ID', 'Chromosome', 'Strand',
                    'Gene_Start', 'Gene_End', 'Transcript_Start', 'Transcript_End'
                ]
                writer.writerow(header)
                
                # 写入数据 | Write data
                for row in transcript_data:
                    writer.writerow(row)
            
            self.logger.info(f"结果已保存 | Results saved: {self.config.output_file}")
            
        except Exception as e:
            self.logger.error(f"写入结果文件时发生错误 | Error writing results file: {e}")
            raise
    
    def print_summary(self, transcript_data: List[List[str]]):
        """
        打印统计摘要 | Print summary statistics
        
        Args:
            transcript_data: 转录本信息列表
        """
        if not transcript_data:
            print("没有找到转录本数据 | No transcript data found")
            return
        
        # 统计信息 | Statistics
        total_transcripts = len(transcript_data)
        genes_with_transcripts = set()
        chromosomes = set()
        strands = set()
        orphan_count = 0
        
        for row in transcript_data:
            gene_id, transcript_id, chromosome, strand, gene_start, gene_end, transcript_start, transcript_end = row
            
            genes_with_transcripts.add(gene_id)
            chromosomes.add(chromosome)
            strands.add(strand)
            
            if gene_start == 'NA' or gene_end == 'NA':
                orphan_count += 1
        
        # 打印统计结果 | Print statistics
        print('\n基因转录本提取统计摘要 | Gene Transcript Extraction Summary:')
        print('=' * 60)
        print(f'总转录本数 | Total transcripts: {total_transcripts}')
        print(f'涉及基因数 | Genes involved: {len(genes_with_transcripts)}')
        print(f'染色体数 | Chromosomes: {len(chromosomes)}')
        print(f'链方向 | Strands: {", ".join(sorted(strands))}')
        print(f'孤儿转录本数 | Orphan transcripts: {orphan_count}')
        
        if chromosomes:
            print(f'染色体列表 | Chromosome list: {", ".join(sorted(chromosomes))}')
        
        print('=' * 60)
