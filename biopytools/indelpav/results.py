"""
PAV结果处理模块|PAV Results Processing Module
"""

import gzip
from pathlib import Path
from typing import List, Dict

class ResultsWriter:
    """结果写入器|Results Writer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def write_pav_results(self, samples: List[str], results: List[Dict]):
        """写入PAV结果|Write PAV results"""
        if not results:
            self.logger.warning("  没有结果可写入|No results to write")
            return
        
        self.logger.info(f" 写入PAV结果到文件|Writing PAV results to file: {self.config.output_file}")
        
        # 决定是否压缩
        if self.config.compress_output or self.config.output_file.endswith('.gz'):
            output_file = self.config.output_file
            if not output_file.endswith('.gz'):
                output_file += '.gz'
            file_handle = gzip.open(output_file, 'wt')
        else:
            file_handle = open(self.config.output_file, 'w')
        
        try:
            # 写入头部
            # header = ['Chromosome', 'Start', 'End', 'INDEL_Sequence', 'INDEL_Type'] + samples
            header = ['Chromosome', 'Start', 'End', 'REF', 'ALT', 'INDEL_Type'] + samples
            file_handle.write('\t'.join(header) + '\n')
            
            # 写入数据
            for result in results:
                # row = [
                #     result['chrom'],
                #     result['start'],
                #     result['end'],
                #     result['sequence'],
                #     result['type']
                # ]
                row = [
                    result['chrom'],
                    result['start'],
                    result['end'],
                    result['ref'],  # 新增REF列
                    result['sequence'],
                    result['type']
                ]
                                
                # 添加样本基因型
                row.extend(map(str, result['genotypes']))
                
                file_handle.write('\t'.join(row) + '\n')
            
            self.logger.info(f" 成功写入 {len(results)} 个INDEL PAV结果|Successfully wrote {len(results)} INDEL PAV results")
            
        finally:
            file_handle.close()
    
    def write_summary_report(self, samples: List[str], results: List[Dict]):
        """写入摘要报告|Write summary report"""
        if not results:
            return
        
        summary_file = Path(self.config.output_file).parent / f"{self.config.base_name}_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("INDEL PAV分析摘要报告|INDEL PAV Analysis Summary Report\n")
            f.write("=" * 60 + "\n\n")
            
            # 输入文件信息
            f.write("输入文件|Input Files:\n")
            f.write(f"  - VCF文件|VCF file: {self.config.vcf_file}\n\n")
            
            # 分析参数
            f.write("分析参数|Analysis Parameters:\n")
            f.write(f"  - 线程数|Threads: {self.config.threads}\n")
            f.write(f"  - 最小INDEL长度|Min INDEL length: {self.config.min_length}\n")
            if self.config.max_length:
                f.write(f"  - 最大INDEL长度|Max INDEL length: {self.config.max_length}\n")
            f.write(f"  - 最小质量分数|Min quality score: {self.config.min_quality}\n")
            f.write(f"  - 最小深度|Min depth: {self.config.min_depth}\n")
            f.write(f"  - 最大缺失率|Max missing rate: {self.config.max_missing_rate}\n")
            f.write(f"  - 包含复杂变异|Include complex: {'是|Yes' if self.config.include_complex else '否|No'}\n\n")
            
            # 结果统计
            f.write("结果统计|Result Statistics:\n")
            f.write(f"  - 样本数量|Sample count: {len(samples)}\n")
            f.write(f"  - INDEL总数|Total INDELs: {len(results)}\n")
            
            # INDEL类型统计
            insertions = sum(1 for r in results if r['type'] == 'insertion')
            deletions = sum(1 for r in results if r['type'] == 'deletion')
            f.write(f"  - 插入|Insertions: {insertions} ({insertions/len(results)*100:.1f}%)\n")
            f.write(f"  - 删除|Deletions: {deletions} ({deletions/len(results)*100:.1f}%)\n")
            
            # 长度统计
            lengths = [r['length'] for r in results]
            f.write(f"  - 平均长度|Average length: {sum(lengths)/len(lengths):.1f} bp\n")
            f.write(f"  - 长度范围|Length range: {min(lengths)}-{max(lengths)} bp\n")
            
            # 染色体分布
            chroms = {}
            for r in results:
                chrom = r['chrom']
                chroms[chrom] = chroms.get(chrom, 0) + 1
            f.write(f"  - 染色体数量|Chromosome count: {len(chroms)}\n")
            
            # 输出文件
            f.write(f"\n输出文件|Output Files:\n")
            f.write(f"  - PAV结果|PAV results: {self.config.output_file}\n")
            f.write(f"  - 摘要报告|Summary report: {summary_file}\n")
        
        self.logger.info(f" 摘要报告已生成|Summary report generated: {summary_file}")
