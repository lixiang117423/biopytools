"""
HTML Alignment Visualization Generator
"""

from pathlib import Path
from datetime import datetime
from typing import Dict, List
from .html_templates import (
    get_css_style, get_javascript, 
    get_index_template, get_sample_template
)

class HTMLAlignmentGenerator:
    """HTML alignment visualization generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.output_dir = config.output_path / config.alignment_output_dir / "html"
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_alignments(self, alignments_data: Dict):
        """Generate HTML alignment files for all samples"""
        self.logger.info("Generating HTML format alignment visualizations...")

        sample_files = []
        
        # 生成每个样品的HTML页面
        for sample_name, sample_data in alignments_data.items():
            sample_file = self._generate_sample_page(sample_name, sample_data)
            sample_files.append(sample_file)
        
        # 生成主入口页面
        index_file = self._generate_index_page(alignments_data)

        self.logger.info(f"HTML alignment files generated: {self.output_dir}")
        self.logger.info(f"Open in browser: {index_file}")

        return sample_files, index_file
    
    def _generate_index_page(self, alignments_data: Dict) -> str:
        """生成主入口页面"""
        output_file = self.output_dir / "index.html"
        
        # 计算统计信息
        total_alignments = sum(len(data['alignments']) for data in alignments_data.values())
        all_alignments = [a for data in alignments_data.values() for a in data['alignments']]
        avg_identity = sum(a['identity'] for a in all_alignments) / len(all_alignments) if all_alignments else 0
        avg_coverage = sum(a['coverage'] for a in all_alignments) / len(all_alignments) if all_alignments else 0
        
        # 生成样品列表HTML
        sample_list_html = self._generate_sample_list_html(alignments_data)
        
        # 生成统计HTML
        statistics_html = f"""
        <div class="metrics">
            <div class="metric"><strong>平均相似度:</strong> {avg_identity:.2f}%</div>
            <div class="metric"><strong>平均覆盖度:</strong> {avg_coverage:.2f}%</div>
            <div class="metric"><strong>高质量比对 (≥90%):</strong> {sum(1 for a in all_alignments if a['identity'] >= 90)}</div>
        </div>
        """
        
        # 渲染模板
        template = get_index_template()
        html_content = template.format(
            css_content=get_css_style(self.config.html_theme),
            js_content=get_javascript(),
            analysis_date=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            blast_type=self.config.blast_type.upper(),
            sample_count=len(alignments_data),
            total_alignments=total_alignments,
            sample_list_html=sample_list_html,
            statistics_html=statistics_html
        )
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)
        
        return str(output_file)
    
    def _generate_sample_list_html(self, alignments_data: Dict) -> str:
        """生成样品列表HTML"""
        items = []
        
        for sample_name, sample_data in sorted(alignments_data.items()):
            alignments = sample_data['alignments']
            alignment_count = len(alignments)
            avg_identity = sum(a['identity'] for a in alignments) / alignment_count if alignments else 0
            
            # 相似度颜色标记
            if avg_identity >= 90:
                identity_class = 'identity-high'
            elif avg_identity >= 80:
                identity_class = 'identity-medium'
            else:
                identity_class = 'identity-low'
            
            item_html = f"""
            <div class="sample-item" data-sample-name="{sample_name}" data-identity="{avg_identity:.2f}">
                <div class="sample-header">
                    <h3>{sample_name}</h3>
                    <a href="{sample_name}_alignments.html" class="btn btn-primary">查看详情</a>
                </div>
                <div class="metrics">
                    <div class="metric"><strong>比对数:</strong> {alignment_count}</div>
                    <div class="metric"><strong>平均相似度:</strong> <span class="identity-badge {identity_class}">{avg_identity:.2f}%</span></div>
                    <div class="metric"><strong>输入文件:</strong> {sample_data['file_name']}</div>
                </div>
            </div>
            """
            items.append(item_html)
        
        return '\n'.join(items)
    
    def _generate_sample_page(self, sample_name: str, sample_data: Dict) -> str:
        """生成单个样品的HTML页面"""
        output_file = self.output_dir / f"{sample_name}_alignments.html"
        alignments = sample_data['alignments']
        
        # 计算统计信息
        avg_identity = sum(a['identity'] for a in alignments) / len(alignments) if alignments else 0
        
        # 生成比对列表HTML
        alignments_html = self._generate_alignments_html(alignments)
        
        # 渲染模板
        template = get_sample_template()
        html_content = template.format(
            css_content=get_css_style(self.config.html_theme),
            js_content=get_javascript(),
            sample_name=sample_name,
            file_name=sample_data['file_name'],
            alignment_count=len(alignments),
            avg_identity=avg_identity,
            alignments_html=alignments_html
        )
        
        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        self.logger.info(f"  {sample_name}: {len(alignments)} alignments")
        return str(output_file)
    
    def _generate_alignments_html(self, alignments: List[Dict]) -> str:
        """生成比对列表HTML"""
        items = []
        
        for idx, alignment in enumerate(alignments, 1):
            alignment_id = f"alignment_{idx}"
            
            # 相似度颜色
            identity = alignment['identity']
            if identity >= 95:
                identity_color = '#2E7D32'
            elif identity >= 90:
                identity_color = '#4CAF50'
            elif identity >= 80:
                identity_color = '#FFC107'
            elif identity >= 70:
                identity_color = '#FF9800'
            else:
                identity_color = '#F44336'
            
            # 检查是否有序列数据
            query_seq = alignment.get('query_seq', '')
            subject_seq = alignment.get('subject_seq', '')
            has_sequences = bool(query_seq and subject_seq)
            
            # 格式化序列比对
            if has_sequences:
                formatted_alignment = self._format_alignment_html(alignment)
                
                # 计算统计
                match_count = sum(1 for q, s in zip(query_seq, subject_seq) if q == s)
                mismatch_count = sum(1 for q, s in zip(query_seq, subject_seq) if q != s and q != '-' and s != '-')
                gap_count = sum(1 for q, s in zip(query_seq, subject_seq) if q == '-' or s == '-')
                total = len(query_seq)
                
                stats_html = f"""
                    <div class="alignment-stats">
                        <span>匹配: {match_count} ({match_count/total*100:.1f}%)</span>
                        <span>错配: {mismatch_count} ({mismatch_count/total*100:.1f}%)</span>
                        <span>Gaps: {gap_count} ({gap_count/total*100:.1f}%)</span>
                        <span>比对长度: {alignment['length']}</span>
                    </div>
                    
                    <div style="margin-top: 10px;">
                        <button class="btn btn-secondary" onclick="copySequence('{alignment_id}')">📋 复制序列</button>
                    </div>
                """
            else:
                formatted_alignment = "⚠️  序列数据不可用\n提示: 重新运行BLAST时需要在outfmt中包含 qseq 和 sseq 字段"
                stats_html = f"""
                    <div class="alignment-stats">
                        <span>比对长度: {alignment['length']}</span>
                        <span>错配: {alignment['mismatch']}</span>
                        <span>Gap: {alignment['gapopen']}</span>
                    </div>
                """
            
            item_html = f"""
            <div class="alignment-item" data-alignment-name="{alignment['subject_id']}" data-identity="{identity:.2f}">
                <div class="alignment-header">
                    <div>
                        <h3>比对 #{idx}: {alignment['query_id']} → {alignment['subject_id']}</h3>
                        <div class="metrics">
                            <span class="metric" style="color: {identity_color}; font-weight: bold;">相似度: {identity:.2f}%</span>
                            <span class="metric">覆盖度: {alignment['coverage']:.2f}%</span>
                            <span class="metric">E-value: {alignment['evalue']}</span>
                            <span class="metric">Bit Score: {alignment['bitscore']}</span>
                        </div>
                    </div>
                    <button class="btn btn-primary" onclick="toggleAlignment('{alignment_id}')">展开/折叠</button>
                </div>
                
                <div class="alignment-content" id="{alignment_id}" style="display: none;">
                    <div class="alignment-view">{formatted_alignment}</div>
                    {stats_html}
                </div>
            </div>
            """
            items.append(item_html)
        
        return '\n'.join(items)
    
    def _format_alignment_html(self, alignment: Dict) -> str:
        """格式化单个比对为HTML"""
        query_seq = alignment.get('query_seq', '')
        subject_seq = alignment.get('subject_seq', '')
        
        # 如果没有序列数据，返回提示信息
        if not query_seq or not subject_seq:
            return "⚠️  序列数据不可用\n提示: 重新运行BLAST时需要在outfmt中包含 qseq 和 sseq 字段"
        
        lines = []
        width = self.config.alignment_width
        q_start = alignment['qstart']
        s_start = alignment['sstart']
        
        # 计算实际位置（排除gap）
        q_pos = q_start
        s_pos = s_start
        
        for i in range(0, len(query_seq), width):
            q_segment = query_seq[i:i+width]
            s_segment = subject_seq[i:i+width]
            
            # 生成匹配行并添加颜色
            match_line = []
            for q, s in zip(q_segment, s_segment):
                if q == s:
                    match_line.append('<span class="match">|</span>')
                elif q == '-' or s == '-':
                    match_line.append('<span class="gap"> </span>')
                else:
                    match_line.append('<span class="mismatch">.</span>')
            
            # 计算这个片段中非gap字符的数量
            q_bases = sum(1 for c in q_segment if c != '-')
            s_bases = sum(1 for c in s_segment if c != '-')
            
            # 计算片段的结束位置
            q_end_pos = q_pos + q_bases - 1 if q_bases > 0 else q_pos
            s_end_pos = s_pos + s_bases - 1 if s_bases > 0 else s_pos
            
            # 格式化输出（使用10位宽度确保对齐）
            lines.append(f"Query  {q_pos:10d}  {q_segment}  {q_end_pos}")
            lines.append(f"                   {''.join(match_line)}")
            lines.append(f"Sbjct  {s_pos:10d}  {s_segment}  {s_end_pos}")
            lines.append("")
            
            # 更新位置（只计算非gap字符）
            q_pos += q_bases
            s_pos += s_bases
        
        return '\n'.join(lines)