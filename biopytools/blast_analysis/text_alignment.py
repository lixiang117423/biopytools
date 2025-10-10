"""
📝 文本格式比对可视化生成器
"""

from pathlib import Path
from typing import Dict, List

class TextAlignmentGenerator:
    """📝 纯文本比对生成器"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.output_dir = config.output_path / config.alignment_output_dir / "text"
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_alignments(self, alignments_data: Dict):
        """生成所有样品的文本比对文件"""
        self.logger.info(f"📝 生成文本格式比对可视化...")
        
        sample_files = []
        all_alignments = []
        
        for sample_name, sample_data in alignments_data.items():
            # 生成单个样品的比对文件
            sample_file = self._write_sample_file(sample_name, sample_data)
            sample_files.append(sample_file)
            all_alignments.extend(sample_data['alignments'])
        
        # 生成汇总文件（可选）
        summary_file = self._write_summary_file(alignments_data, all_alignments)
        
        self.logger.info(f"✅ 文本比对文件已生成: {self.output_dir}")
        return sample_files, summary_file
    
    def _write_sample_file(self, sample_name: str, sample_data: Dict) -> str:
        """写入单个样品的比对文件"""
        output_file = self.output_dir / f"{sample_name}_alignments.txt"
        alignments = sample_data['alignments']
        
        with open(output_file, 'w', encoding='utf-8') as f:
            # 文件头
            f.write("=" * 80 + "\n")
            f.write(f"🧬 BLAST序列比对可视化\n")
            f.write("=" * 80 + "\n")
            f.write(f"样品名称: {sample_name}\n")
            f.write(f"输入文件: {sample_data['file_name']}\n")
            f.write(f"比对数量: {len(alignments)}\n")
            
            if alignments:
                avg_identity = sum(a['identity'] for a in alignments) / len(alignments)
                avg_coverage = sum(a['coverage'] for a in alignments) / len(alignments)
                f.write(f"平均相似度: {avg_identity:.2f}%\n")
                f.write(f"平均覆盖度: {avg_coverage:.2f}%\n")
            
            f.write("=" * 80 + "\n\n")
            
            # 写入每个比对
            for idx, alignment in enumerate(alignments, 1):
                alignment_text = self._format_alignment(idx, alignment)
                f.write(alignment_text)
                f.write("\n")
        
        self.logger.info(f"  ✅ {sample_name}: {len(alignments)} 个比对")
        return str(output_file)
    
    def _write_summary_file(self, alignments_data: Dict, all_alignments: List) -> str:
        """写入汇总文件"""
        output_file = self.output_dir / "all_samples_alignments.txt"
        
        with open(output_file, 'w', encoding='utf-8') as f:
            # 汇总头
            f.write("=" * 80 + "\n")
            f.write(f"🧬 BLAST序列比对汇总报告\n")
            f.write("=" * 80 + "\n")
            f.write(f"总样品数: {len(alignments_data)}\n")
            f.write(f"总比对数: {len(all_alignments)}\n")
            
            if all_alignments:
                avg_identity = sum(a['identity'] for a in all_alignments) / len(all_alignments)
                avg_coverage = sum(a['coverage'] for a in all_alignments) / len(all_alignments)
                f.write(f"平均相似度: {avg_identity:.2f}%\n")
                f.write(f"平均覆盖度: {avg_coverage:.2f}%\n")
            
            f.write("=" * 80 + "\n\n")
            
            # 按样品写入比对
            for sample_name, sample_data in alignments_data.items():
                f.write("\n" + "=" * 80 + "\n")
                f.write(f"样品: {sample_name}\n")
                f.write("=" * 80 + "\n\n")
                
                for idx, alignment in enumerate(sample_data['alignments'], 1):
                    alignment_text = self._format_alignment(idx, alignment)
                    f.write(alignment_text)
                    f.write("\n")
        
        self.logger.info(f"  ✅ 汇总文件: {len(alignments_data)} 个样品")
        return str(output_file)
    
    def _format_alignment(self, index: int, alignment: Dict) -> str:
        """格式化单个比对为文本"""
        lines = []
        
        # 比对头部信息
        lines.append("─" * 80)
        lines.append(f"比对 #{index}: {alignment['query_id']} → {alignment['subject_id']}")
        lines.append("─" * 80)
        lines.append(f"相似度: {alignment['identity']:.2f}% | "
                    f"覆盖度: {alignment['coverage']:.2f}% | "
                    f"E-value: {alignment['evalue']} | "
                    f"Bit Score: {alignment['bitscore']}")
        lines.append(f"比对长度: {alignment['length']} | "
                    f"错配: {alignment['mismatch']} | "
                    f"Gap: {alignment['gapopen']}")
        lines.append("─" * 80)
        lines.append("")
        
        # 检查是否有序列数据
        query_seq = alignment.get('query_seq', '')
        subject_seq = alignment.get('subject_seq', '')
        
        if not query_seq or not subject_seq:
            # 没有序列数据
            lines.append("⚠️  序列数据不可用")
            lines.append("提示: 重新运行BLAST时需要在outfmt中包含 qseq 和 sseq 字段")
            lines.append("")
        else:
            # 格式化序列比对
            # 生成匹配行
            match_line = self._generate_match_line(query_seq, subject_seq)
            
            # 分行显示
            width = self.config.alignment_width
            q_start = alignment['qstart']
            s_start = alignment['sstart']
            
            # 计算实际位置（排除gap）
            q_pos = q_start
            s_pos = s_start
            
            for i in range(0, len(query_seq), width):
                q_segment = query_seq[i:i+width]
                s_segment = subject_seq[i:i+width]
                m_segment = match_line[i:i+width]
                
                # 计算这个片段中非gap字符的数量
                q_bases = sum(1 for c in q_segment if c != '-')
                s_bases = sum(1 for c in s_segment if c != '-')
                
                # 计算片段的结束位置
                q_end_pos = q_pos + q_bases - 1 if q_bases > 0 else q_pos
                s_end_pos = s_pos + s_bases - 1 if s_bases > 0 else s_pos
                
                # 格式化输出（使用10位宽度确保对齐）
                lines.append(f"Query  {q_pos:10d}  {q_segment}  {q_end_pos}")
                lines.append(f"                   {m_segment}")
                lines.append(f"Sbjct  {s_pos:10d}  {s_segment}  {s_end_pos}")
                lines.append("")
                
                # 更新位置（只计算非gap字符）
                q_pos += q_bases
                s_pos += s_bases
            
            # 统计信息
            match_count = match_line.count('|')
            mismatch_count = match_line.count('.')
            gap_count = match_line.count(' ')
            total = len(match_line)
            
            lines.append("统计信息:")
            lines.append(f"  匹配碱基: {match_count} / {total} ({match_count/total*100:.1f}%)")
            lines.append(f"  错配碱基: {mismatch_count} ({mismatch_count/total*100:.1f}%)")
            if gap_count > 0:
                lines.append(f"  Gap碱基: {gap_count} ({gap_count/total*100:.1f}%)")
        
        lines.append("")
        return "\n".join(lines)
    
    def _generate_match_line(self, query_seq: str, subject_seq: str) -> str:
        """生成匹配标记行"""
        match_line = []
        for q, s in zip(query_seq, subject_seq):
            if q == s:
                match_line.append('|')
            elif q == '-' or s == '-':
                match_line.append(' ') 
            else:
                match_line.append('.')
        return ''.join(match_line)