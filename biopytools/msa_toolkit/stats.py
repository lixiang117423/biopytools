"""
比对结果统计模块 | Alignment Statistics Module
"""

from Bio import AlignIO
from pathlib import Path

class AlignmentStats:
    """比对统计分析器 | Alignment Statistics Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.stats_file = f"{self.config.output_prefix}_stats.txt"
    
    def calculate_stats(self):
        """计算比对统计信息 | Calculate alignment statistics"""
        try:
            self.logger.info("📊 计算比对统计信息 | Calculating alignment statistics")
            
            # 读取比对结果 | Read alignment
            alignment_file = f"{self.config.output_prefix}.{self.config.output_format}"
            
            # 根据格式读取 | Read based on format
            format_map = {
                'fasta': 'fasta',
                'clustal': 'clustal',
                'phylip': 'phylip-relaxed',
                'nexus': 'nexus'
            }
            
            alignment = AlignIO.read(alignment_file, format_map.get(self.config.output_format, 'fasta'))
            
            # 计算统计信息 | Calculate statistics
            num_sequences = len(alignment)
            alignment_length = alignment.get_alignment_length()
            
            # 计算一致性 | Calculate identity
            identities = []
            for i in range(alignment_length):
                column = alignment[:, i]
                # 计算最常见碱基的频率 | Calculate frequency of most common base
                from collections import Counter
                counts = Counter(column)
                most_common = counts.most_common(1)[0]
                identity = most_common[1] / num_sequences
                identities.append(identity)
            
            avg_identity = sum(identities) / len(identities) * 100
            
            # 保存统计信息 | Save statistics
            with open(self.stats_file, 'w') as f:
                f.write("="*60 + "\n")
                f.write("多序列比对统计报告 | Multiple Sequence Alignment Statistics\n")
                f.write("="*60 + "\n\n")
                f.write(f"📁 输入文件 | Input file: {self.config.input_file}\n")
                f.write(f"🔧 比对方法 | Alignment method: {self.config.method.upper()}\n")
                f.write(f"📊 序列数量 | Number of sequences: {num_sequences}\n")
                f.write(f"📏 比对长度 | Alignment length: {alignment_length}\n")
                f.write(f"🎯 平均一致性 | Average identity: {avg_identity:.2f}%\n")
                f.write(f"🧵 使用线程 | Threads used: {self.config.threads}\n")
                f.write(f"💾 输出文件 | Output file: {alignment_file}\n")
                f.write("\n" + "="*60 + "\n")
            
            self.logger.info(f"✅ 统计报告已生成 | Statistics report generated: {self.stats_file}")
            self.logger.info(f"📊 序列数: {num_sequences} | 比对长度: {alignment_length} | 平均一致性: {avg_identity:.2f}%")
            
            return True
            
        except Exception as e:
            self.logger.warning(f"⚠️  统计计算失败(需要BioPython) | Statistics calculation failed (BioPython required): {e}")
            return False
