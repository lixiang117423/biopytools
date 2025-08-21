"""
📊 BLAST结果处理模块
"""

import os
import pandas as pd
from pathlib import Path
from typing import List, Tuple

class ResultsProcessor:
    """📊 BLAST结果处理器"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def process_blast_results(self, blast_results: List[Tuple[str, str, str]]) -> str:
        """处理BLAST结果"""
        self.logger.info("📊 开始处理BLAST结果")
        
        # 创建汇总文件
        summary_file = self.config.output_path / "blast_summary_results.tsv"
        
        # 写入表头
        header = [
            "输入文件", "样品名称", "查询序列ID", "目标序列ID", "序列相似度(%)", 
            "比对长度", "错配数", "Gap数", "查询起始位置", "查询结束位置", 
            "目标起始位置", "目标结束位置", "E-value", "Bit_Score", "目标序列长度", "目标序列覆盖度(%)"
        ]
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write('\t'.join(header) + '\n')
        
        total_hits = 0
        files_with_hits = 0
        
        # 处理每个结果文件
        for file_name, sample_name, result_file in blast_results:
            if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
                continue
            
            self.logger.info(f"📋 处理结果文件: {result_file}")
            
            hit_count = 0
            has_hits = False
            
            with open(result_file, 'r') as rf, open(summary_file, 'a', encoding='utf-8') as sf:
                for line in rf:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 13:
                            try:
                                sstart, send, slen = int(parts[8]), int(parts[9]), int(parts[12])
                                coverage = abs(send - sstart + 1) / slen * 100
                                coverage = min(coverage, 100.0)
                                
                                # 写入汇总文件
                                result_line = [
                                    file_name, sample_name
                                ] + parts[:12] + [parts[12], f"{coverage:.2f}"]
                                
                                sf.write('\t'.join(result_line) + '\n')
                                hit_count += 1
                                has_hits = True
                                
                            except (ValueError, IndexError):
                                self.logger.warning(f"⚠️  跳过格式错误的行: {line.strip()}")
            
            if has_hits:
                files_with_hits += 1
                total_hits += hit_count
                self.logger.info(f"  ✅ 处理完成，发现 {hit_count} 个比对结果")
        
        self.logger.info(f"🎯 汇总统计:")
        self.logger.info(f"  📂 有比对结果的文件数: {files_with_hits}")
        self.logger.info(f"  🎯 总比对命中数: {total_hits}")
        
        return str(summary_file)
    
    def sort_results_by_coverage(self, summary_file: str) -> str:
        """按覆盖度排序结果"""
        self.logger.info("📊 按目标序列覆盖度排序结果")
        
        try:
            df = pd.read_csv(summary_file, sep='\t', encoding='utf-8')
            
            if df.empty:
                self.logger.warning("⚠️  汇总文件为空")
                return summary_file
            
            # 按覆盖度和相似度排序
            df_sorted = df.sort_values(
                by=['目标序列覆盖度(%)', '序列相似度(%)', 'E-value'], 
                ascending=[False, False, True]
            )
            
            # 保存排序后的文件
            sorted_file = summary_file.replace('.tsv', '_sorted.tsv')
            df_sorted.to_csv(sorted_file, sep='\t', index=False, encoding='utf-8')
            
            self.logger.info(f"✅ 排序完成: {sorted_file}")
            return sorted_file
            
        except Exception as e:
            self.logger.error(f"❌ 排序过程中出错: {e}")
            return summary_file
    
    def create_high_quality_results(self, summary_file: str) -> str:
        """创建高质量比对结果"""
        self.logger.info("🌟 筛选高质量比对结果")
        
        try:
            df = pd.read_csv(summary_file, sep='\t', encoding='utf-8')
            
            if df.empty:
                self.logger.warning("⚠️  无数据可筛选")
                return ""
            
            # 筛选高质量结果
            high_quality_df = df[
                (df['E-value'] <= self.config.high_quality_evalue) & 
                (df['序列相似度(%)'] >= self.config.min_identity) &
                (df['目标序列覆盖度(%)'] >= self.config.min_coverage)
            ]
            
            if high_quality_df.empty:
                self.logger.warning("⚠️  未找到符合条件的高质量比对结果")
                return ""
            
            # 保存高质量结果
            high_quality_file = summary_file.replace('.tsv', '_high_quality.tsv')
            high_quality_df.to_csv(high_quality_file, sep='\t', index=False, encoding='utf-8')
            
            self.logger.info(f"🌟 高质量结果文件已创建: {high_quality_file}")
            self.logger.info(f"🌟 高质量比对数量: {len(high_quality_df)}")
            
            return high_quality_file
            
        except Exception as e:
            self.logger.error(f"❌ 创建高质量结果时出错: {e}")
            return ""
