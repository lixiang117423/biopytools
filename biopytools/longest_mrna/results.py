"""
最长转录本提取结果处理模块 | Longest mRNA Extraction Results Processing Module
"""

import os
from typing import Dict, Any, List
from .data_processing import TranscriptProcessor

class GeneInfoGenerator:
    """基因信息生成器 | Gene Info Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_gene_info(self, longest_transcripts: Dict[str, Dict[str, Any]], gene_metadata: Dict[str, Dict[str, Any]]):
        """生成基因信息文件 | Generate gene info file"""
        self.logger.info(f"生成基因信息文件 | Generating gene info file: {self.config.gene_info_file}")
        
        try:
            with open(self.config.gene_info_file, 'w', encoding='utf-8') as info_file:
                # 写入头部 | Write header
                info_file.write("mRNA_ID\tgene_ID\tmRNA_start\tmRNA_end\tgene_start\tgene_end\tstrand\tchr\n")
                
                # 写入数据 | Write data
                for gene_id, transcript_info in longest_transcripts.items():
                    gene_data = gene_metadata.get(gene_id, {})
                    transcript_id = transcript_info.get('id', '')
                    
                    info_file.write(
                        f"{transcript_id}\t{gene_id}\t"
                        f"{transcript_info.get('start', '')}\t{transcript_info.get('end', '')}\t"
                        f"{gene_data.get('start', '')}\t{gene_data.get('end', '')}\t"
                        f"{transcript_info.get('strand', '')}\t{transcript_info.get('chrom', '')}\n"
                    )
            
            self.logger.info(f"✓ 基因信息文件生成完成 | Gene info file generated: {self.config.gene_info_file}")
            
        except Exception as e:
            self.logger.error(f"✗ 生成基因信息文件失败 | Failed to generate gene info file: {e}")
            raise

class StatisticsCalculator:
    """统计计算器 | Statistics Calculator"""
    
    def __init__(self, logger):
        self.logger = logger
        self.transcript_processor = TranscriptProcessor(logger)
    
    def calculate_statistics(self, gene_transcripts: Dict[str, List[Dict[str, Any]]], 
                           longest_transcripts: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """计算统计信息 | Calculate statistics"""
        self.logger.info("计算统计信息 | Calculating statistics")
        
        total_genes = len(gene_transcripts)
        multi_isoform_genes = sum(1 for transcripts in gene_transcripts.values() if len(transcripts) > 1)
        
        # 计算平均转录本长度 | Calculate average transcript length
        lengths = []
        for transcript_info in longest_transcripts.values():
            length = self.transcript_processor.calculate_transcript_length(transcript_info)
            lengths.append(length)
        
        avg_length = sum(lengths) / len(lengths) if lengths else 0
        
        stats = {
            'total_genes': total_genes,
            'multi_isoform_genes': multi_isoform_genes,
            'avg_length': avg_length,
            'longest_transcripts_count': len(longest_transcripts)
        }
        
        return stats

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def print_summary(self, stats: Dict[str, Any]):
        """打印总结信息 | Print summary"""
        self.logger.info("=" * 50)
        self.logger.info("提取完成总结 | Extraction Summary")
        self.logger.info("=" * 50)
        self.logger.info(f"总基因数 | Total genes processed: {stats['total_genes']}")
        self.logger.info(f"多转录本基因数 | Genes with multiple transcripts: {stats['multi_isoform_genes']}")
        self.logger.info(f"平均转录本长度 | Average transcript length: {stats['avg_length']:.2f}")
        self.logger.info(f"提取的最长转录本数 | Longest transcripts extracted: {stats['longest_transcripts_count']}")
        self.logger.info(f"输出文件 | Output files:")
        self.logger.info(f"  - 蛋白质序列文件 | Protein sequences: {self.config.output_file}")
        self.logger.info(f"  - 基因信息文件 | Gene info: {self.config.gene_info_file}")
        self.logger.info("=" * 50)
