"""
GFF3结果处理模块 | GFF3 Results Processing Module
"""

import os
from typing import List, Dict, Any

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, transcript_data: List[List[str]]):
        """
        生成总结报告 | Generate summary report
        
        Args:
            transcript_data: 转录本信息列表
        """
        report_file = os.path.join(os.path.dirname(self.config.output_file), "gff_extraction_summary.txt")
        
        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("GFF3基因转录本提取总结报告 | GFF3 Gene Transcript Extraction Summary Report\n")
                f.write("=" * 70 + "\n\n")
                
                # 输入文件信息 | Input file information
                f.write("输入文件信息 | Input File Information:\n")
                f.write(f"  - GFF3文件 | GFF3 file: {self.config.gff3_file}\n")
                f.write(f"  - 基因类型 | Gene type: {self.config.gene_type}\n")
                f.write(f"  - 转录本类型 | Transcript types: {', '.join(self.config.transcript_types)}\n\n")
                
                # 统计信息 | Statistics
                if transcript_data:
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
                    
                    f.write("提取统计 | Extraction Statistics:\n")
                    f.write(f"  - 总转录本数 | Total transcripts: {total_transcripts}\n")
                    f.write(f"  - 涉及基因数 | Genes involved: {len(genes_with_transcripts)}\n")
                    f.write(f"  - 染色体数 | Chromosomes: {len(chromosomes)}\n")
                    f.write(f"  - 孤儿转录本数 | Orphan transcripts: {orphan_count}\n")
                    f.write(f"  - 孤儿转录本比例 | Orphan transcript ratio: {orphan_count/total_transcripts:.2%}\n")
                    
                    if chromosomes:
                        f.write(f"  - 染色体列表 | Chromosome list: {', '.join(sorted(chromosomes))}\n")
                    
                    f.write(f"  - 链方向 | Strands: {', '.join(sorted(strands))}\n")
                    
                    # 基因转录本比例 | Gene-transcript ratio
                    if genes_with_transcripts:
                        avg_transcripts_per_gene = total_transcripts / len(genes_with_transcripts)
                        f.write(f"  - 平均每个基因转录本数 | Average transcripts per gene: {avg_transcripts_per_gene:.2f}\n")
                else:
                    f.write("提取统计 | Extraction Statistics:\n")
                    f.write("  - 未找到转录本数据 | No transcript data found\n")
                
                f.write(f"\n输出文件 | Output file: {self.config.output_file}\n")
                f.write(f"报告生成时间 | Report generated: {self._get_current_time()}\n")
                
        except Exception as e:
            self.logger.error(f"生成总结报告时发生错误 | Error generating summary report: {e}")
            raise
        
        self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
    
    def _get_current_time(self) -> str:
        """获取当前时间字符串 | Get current time string"""
        from datetime import datetime
        return datetime.now().strftime("%Y-%m-%d %H:%M:%S")
