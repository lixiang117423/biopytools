"""
Minimap2结果处理模块 | Minimap2 Results Processing Module
"""

import os
from pathlib import Path

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, unmapped_regions):
        """生成总结报告 | Generate summary report"""
        report_file = os.path.join(self.config.output_dir, "minimap2_summary.txt")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("Minimap2分析总结报告 | Minimap2 Analysis Summary Report\n")
            f.write("=" * 50 + "\n\n")
            
            # 输入文件信息 | Input file information
            f.write("输入文件 | Input Files:\n")
            f.write(f"  - 目标基因组 | Target genome: {self.config.target_genome}\n")
            f.write(f"  - 查询基因组 | Query genome: {self.config.query_genome}\n\n")
            
            # 分析参数 | Analysis parameters
            f.write("分析参数 | Analysis Parameters:\n")
            f.write(f"  - 预设参数 | Preset: {self.config.preset}\n")
            f.write(f"  - 线程数 | Threads: {self.config.threads}\n")
            f.write(f"  - 最小匹配长度 | Minimum match length: {self.config.min_match_length}\n")
            f.write(f"  - 最小未比对长度 | Minimum unmapped length: {self.config.min_unmapped_length}\n")
            f.write(f"  - tp类型过滤 | tp type filter: {self.config.tp_type}\n\n")
            
            # 输出文件信息 | Output file information
            f.write("输出文件 | Output Files:\n")
            f.write(f"  - PAF比对文件 | PAF alignment file: {self.config.paf_file}\n")
            f.write(f"  - BED区间文件 | BED interval file: {self.config.bed_file}\n")
            f.write(f"  - 未比对序列文件 | Unmapped sequence file: {self.config.unmapped_fasta}\n\n")
            
            # 结果统计 | Result statistics
            f.write("结果统计 | Result Statistics:\n")
            if unmapped_regions:
                f.write(f"  - 未比对区间数量 | Number of unmapped regions: {len(unmapped_regions)}\n")
                
                # 按query_name统计 | Statistics by query_name
                query_stats = {}
                total_unmapped_length = 0
                
                for region in unmapped_regions:
                    query_name = region['query_name']
                    length = region['length']
                    
                    if query_name not in query_stats:
                        query_stats[query_name] = {'count': 0, 'total_length': 0}
                    
                    query_stats[query_name]['count'] += 1
                    query_stats[query_name]['total_length'] += length
                    total_unmapped_length += length
                
                f.write(f"  - 未比对序列总长度 | Total unmapped sequence length: {total_unmapped_length:,} bp\n")
                f.write(f"  - 涉及的查询序列数 | Number of query sequences involved: {len(query_stats)}\n\n")
                
                f.write("各查询序列统计 | Statistics by query sequence:\n")
                for query_name, stats in query_stats.items():
                    f.write(f"  - {query_name}: {stats['count']} 个区间, {stats['total_length']:,} bp\n")
            else:
                f.write("  - 未找到符合条件的未比对区间 | No unmapped regions found meeting criteria\n")
        
        self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
