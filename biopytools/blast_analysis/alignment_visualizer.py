"""
🧬 序列比对可视化生成器主模块
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple
from .text_alignment import TextAlignmentGenerator
from .html_alignment import HTMLAlignmentGenerator

class AlignmentVisualizer:
    """🧬 比对可视化生成器（支持文本和HTML）"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.text_generator = TextAlignmentGenerator(config, logger)
        self.html_generator = HTMLAlignmentGenerator(config, logger)
    
    def generate_visualizations(self, blast_results: List[Tuple[str, str, str]]):
        """
        根据配置生成可视化文件
        
        Args:
            blast_results: BLAST结果列表 [(file_name, sample_name, result_file), ...]
        """
        if self.config.alignment_output == 'none':
            self.logger.info("⏭️  跳过比对可视化生成（未启用）")
            return None
        
        self.logger.info("=" * 80)
        self.logger.info("🧬 开始生成序列比对可视化")
        self.logger.info("=" * 80)
        
        # 解析BLAST结果，提取比对数据
        alignments_data = self._parse_blast_results(blast_results)
        
        if not alignments_data:
            self.logger.warning("⚠️  没有符合条件的比对数据，跳过可视化生成")
            return None
        
        # 根据配置生成相应格式
        output_files = {}
        
        if self.config.alignment_output in ['text', 'both']:
            text_files, text_summary = self.text_generator.generate_alignments(alignments_data)
            output_files['text'] = {'sample_files': text_files, 'summary': text_summary}
        
        if self.config.alignment_output in ['html', 'both']:
            html_files, html_index = self.html_generator.generate_alignments(alignments_data)
            output_files['html'] = {'sample_files': html_files, 'index': html_index}
        
        self.logger.info("=" * 80)
        self.logger.info("✅ 比对可视化生成完成！")
        self.logger.info("=" * 80)
        
        return output_files
    
    def _parse_blast_results(self, blast_results: List[Tuple[str, str, str]]) -> Dict:
        """
        解析BLAST结果文件，提取比对数据
        
        Returns:
            Dict: {sample_name: {'file_name': str, 'alignments': [...]}}
        """
        self.logger.info("📊 解析BLAST结果文件...")
        
        alignments_data = {}
        total_parsed = 0
        total_filtered = 0
        
        for file_name, sample_name, result_file in blast_results:
            if not os.path.exists(result_file) or os.path.getsize(result_file) == 0:
                continue
            
            sample_alignments = []
            
            try:
                with open(result_file, 'r', encoding='utf-8') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        
                        alignment = self._parse_blast_line(line)
                        if alignment and self._passes_filters(alignment):
                            sample_alignments.append(alignment)
                            total_parsed += 1
                        else:
                            total_filtered += 1
            
            except Exception as e:
                self.logger.warning(f"⚠️  解析文件失败 {result_file}: {e}")
                continue
            
            # 限制每个样品的比对数量
            if len(sample_alignments) > self.config.alignment_max_per_sample:
                self.logger.info(f"  📌 {sample_name}: 限制显示前 {self.config.alignment_max_per_sample} 个比对")
                sample_alignments = sample_alignments[:self.config.alignment_max_per_sample]
            
            if sample_alignments:
                alignments_data[sample_name] = {
                    'file_name': file_name,
                    'alignments': sample_alignments
                }
        
        self.logger.info(f"📊 解析完成:")
        self.logger.info(f"  ✅ 符合条件的比对: {total_parsed}")
        self.logger.info(f"  🔽 已过滤: {total_filtered}")
        self.logger.info(f"  📂 涉及样品: {len(alignments_data)}")
        
        return alignments_data
    
    def _parse_blast_line(self, line: str) -> Dict:
        """解析BLAST输出的单行"""
        parts = line.split('\t')
        
        # 检查字段数量是否足够
        # 基本字段: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore slen
        # 额外字段: qseq sseq (如果启用了比对可视化)
        min_fields = 13  # 至少需要基本字段
        
        if len(parts) < min_fields:
            return None
        
        try:
            # 解析基本字段
            alignment = {
                'query_id': parts[0],
                'subject_id': parts[1],
                'identity': float(parts[2]),
                'length': int(parts[3]),
                'mismatch': int(parts[4]),
                'gapopen': int(parts[5]),
                'qstart': int(parts[6]),
                'qend': int(parts[7]),
                'sstart': int(parts[8]),
                'send': int(parts[9]),
                'evalue': parts[10],
                'bitscore': float(parts[11]),
                'slen': int(parts[12])
            }
            
            # 计算覆盖度
            coverage = abs(alignment['send'] - alignment['sstart'] + 1) / alignment['slen'] * 100
            alignment['coverage'] = min(coverage, 100.0)
            
            # 尝试提取序列字段（如果存在）
            if len(parts) >= 15:
                # BLAST输出包含序列字段
                alignment['query_seq'] = parts[13]
                alignment['subject_seq'] = parts[14]
            elif self.config.needs_alignment_sequences():
                # 需要序列但BLAST输出中没有，记录警告
                if not hasattr(self, '_sequence_warning_logged'):
                    self.logger.warning("⚠️  BLAST输出中缺少序列字段（qseq/sseq）")
                    self.logger.warning("💡 提示：比对可视化将只显示统计信息，不显示序列对齐")
                    self.logger.warning("💡 建议：重新运行BLAST，或使用 --alignment-output none 禁用可视化")
                    self._sequence_warning_logged = True
                alignment['query_seq'] = ''
                alignment['subject_seq'] = ''
            else:
                # 不需要序列，使用空字符串
                alignment['query_seq'] = ''
                alignment['subject_seq'] = ''
            
            return alignment
            
        except (ValueError, IndexError) as e:
            self.logger.debug(f"解析行失败: {e}")
            return None
    
    def _passes_filters(self, alignment: Dict) -> bool:
        """检查比对是否通过过滤条件"""
        # 相似度过滤
        if self.config.alignment_min_identity > 0:
            if alignment['identity'] < self.config.alignment_min_identity:
                return False
        
        # 覆盖度过滤
        if self.config.alignment_min_coverage > 0:
            if alignment['coverage'] < self.config.alignment_min_coverage:
                return False
        
        return True