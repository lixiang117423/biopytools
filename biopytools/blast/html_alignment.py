"""
HTML比对可视化生成器|HTML Alignment Visualization Generator
"""

from pathlib import Path
from datetime import datetime
from typing import Dict, List
from .html_templates import (
    get_css_style, get_javascript,
    get_index_template, get_sample_template,
    get_tab_css, get_merged_javascript, get_merged_template
)

# 合并单文件名(旧模式遗留文件清理时排除自身)|Merged single-file name
# (excluded from stale-file cleanup so it never deletes itself)
MERGED_HTML_FILENAME = "blast_alignments.html"

class HTMLAlignmentGenerator:
    """HTML比对可视化生成器|HTML alignment visualization generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.output_dir = config.output_path / config.alignment_output_dir / "html"
        self.output_dir.mkdir(parents=True, exist_ok=True)
    
    def generate_alignments(self, alignments_data: Dict):
        """生成所有样品的HTML比对文件|Generate HTML alignment files for all samples"""
        self.logger.info("生成HTML格式比对可视化|Generating HTML format alignment visualizations...")

        # 合并模式:输出单个自包含文件,便于分享|Merged mode: single self-contained file for easy sharing
        if self.config.merge_html:
            merged_file = self._generate_merged_page(alignments_data)
            self.logger.info(f"HTML合并文件已生成|Merged HTML file generated: {self.output_dir}")
            self.logger.info(f"浏览器打开|Open in browser: {merged_file}")
            return [], merged_file

        # 旧模式:index页+分样品页|Legacy mode: index page + per-sample pages
        self._remove_stale_merged_file()

        sample_files = []

        for sample_name, sample_data in alignments_data.items():
            sample_file = self._generate_sample_page(sample_name, sample_data)
            sample_files.append(sample_file)

        index_file = self._generate_index_page(alignments_data)

        self.logger.info(f"HTML比对文件已生成|HTML alignment files generated: {self.output_dir}")
        self.logger.info(f"浏览器打开|Open in browser: {index_file}")

        return sample_files, index_file

    def _generate_merged_page(self, alignments_data: Dict) -> str:
        """生成合并单页HTML(概览tab+每样品tab)|Generate merged single-page HTML (overview + per-sample tabs)"""
        # 清理旧模式遗留文件,避免同目录两套HTML并存|Remove stale legacy files first
        self._remove_stale_legacy_files()

        output_file = self.output_dir / MERGED_HTML_FILENAME
        total_alignments = sum(len(data['alignments']) for data in alignments_data.values())

        # tab按钮和样品面板:样品名排序与概览列表一致|tab buttons/panels in the same sorted order as the list
        tab_buttons = ['<button class="tab-btn active" onclick="switchTab(0)">概览|Overview</button>']
        tab_indices = {}
        sample_panels = []
        for i, (sample_name, sample_data) in enumerate(sorted(alignments_data.items()), start=1):
            tab_indices[sample_name] = i
            tab_buttons.append(f'<button class="tab-btn" onclick="switchTab({i})">{sample_name}</button>')
            sample_panels.append(self._generate_sample_panel(sample_name, sample_data, i))

        overview_html = self._generate_merged_overview_html(alignments_data, tab_indices)

        template = get_merged_template()
        html_content = template.format(
            css_content=get_css_style(self.config.html_theme),
            tab_css_content=get_tab_css(),
            js_content=get_merged_javascript(),
            analysis_date=datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            blast_type=self.config.blast_type.upper(),
            sample_count=len(alignments_data),
            total_alignments=total_alignments,
            tab_buttons_html='\n            '.join(tab_buttons),
            overview_html=overview_html,
            sample_panels_html='\n'.join(sample_panels)
        )

        with open(output_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        return str(output_file)

    def _generate_merged_overview_html(self, alignments_data: Dict, tab_indices: Dict) -> str:
        """生成合并页概览tab内容|Generate merged page overview tab content"""
        sample_list_html = self._generate_sample_list_html(alignments_data, tab_indices)

        all_alignments = [a for data in alignments_data.values() for a in data['alignments']]
        avg_identity = sum(a['identity'] for a in all_alignments) / len(all_alignments) if all_alignments else 0
        avg_coverage = sum(a['coverage'] for a in all_alignments) / len(all_alignments) if all_alignments else 0

        statistics_html = f"""
        <div class="metrics">
            <div class="metric"><strong>平均相似度:</strong> {avg_identity:.2f}%</div>
            <div class="metric"><strong>平均覆盖度:</strong> {avg_coverage:.2f}%</div>
            <div class="metric"><strong>高质量比对 (≥90%):</strong> {sum(1 for a in all_alignments if a['identity'] >= 90)}</div>
        </div>
        """

        return f"""
            <section class="search-filter">
                <input type="text" id="search-0" class="filter-search" placeholder="搜索样品名称|Search samples...">
                <select id="identity-filter-0" class="filter-identity">
                    <option value="0">所有相似度</option>
                    <option value="70">≥ 70%</option>
                    <option value="80">≥ 80%</option>
                    <option value="90">≥ 90%</option>
                    <option value="95">≥ 95%</option>
                </select>
                <button onclick="applyFilters()">应用筛选</button>
            </section>

            <section class="sample-list">
                <h2>样品列表</h2>
                {sample_list_html}
            </section>

            <section class="statistics">
                <h2>快速统计|Quick Statistics</h2>
                {statistics_html}
            </section>
        """

    def _generate_sample_panel(self, sample_name: str, sample_data: Dict, tab_index: int) -> str:
        """生成单个样品的tab面板内容|Generate one sample's tab panel content"""
        alignments = sample_data['alignments']
        avg_identity = sum(a['identity'] for a in alignments) / len(alignments) if alignments else 0

        # 比对ID加tab索引前缀,避免跨样品ID冲突|tab-index prefix avoids cross-sample ID collisions
        alignments_html = self._generate_alignments_html(alignments, id_prefix=str(tab_index))

        return f"""
        <div class="tab-panel" id="tab-panel-{tab_index}">
            <section class="sample-detail-header">
                <h2>{sample_name} - Sequence Alignments</h2>
                <div class="detail-metrics">
                    <span>输入文件|Input: {sample_data['file_name']}</span>
                    <span>比对数量|Alignments: {len(alignments)}</span>
                    <span>平均相似度|Avg identity: {avg_identity:.2f}%</span>
                </div>
            </section>

            <section class="search-filter">
                <input type="text" id="search-{tab_index}" class="filter-search" placeholder="搜索目标序列ID|Search subject IDs...">
                <select id="identity-filter-{tab_index}" class="filter-identity">
                    <option value="0">所有相似度</option>
                    <option value="70">≥ 70%</option>
                    <option value="80">≥ 80%</option>
                    <option value="90">≥ 90%</option>
                    <option value="95">≥ 95%</option>
                </select>
                <button onclick="applyFilters()">应用筛选</button>
                <button onclick="expandAll()" class="btn-secondary">展开全部</button>
                <button onclick="collapseAll()" class="btn-secondary">折叠全部</button>
            </section>

            <section class="alignment-list">
                {alignments_html}
            </section>
        </div>
        """

    def _remove_stale_legacy_files(self):
        """清理旧模式遗留的index.html和分样品页(排除合并文件自身)|Remove stale legacy files"""
        if not self.output_dir.exists():
            return
        stale_files = [self.output_dir / "index.html"]
        stale_files.extend(
            p for p in self.output_dir.glob("*_alignments.html")
            if p.name != MERGED_HTML_FILENAME
        )
        for path in stale_files:
            if path.exists():
                path.unlink()
                self.logger.info(f"清理遗留HTML文件|Removed stale HTML file: {path.name}")

    def _remove_stale_merged_file(self):
        """清理合并模式遗留的合并文件|Remove stale merged file"""
        merged_file = self.output_dir / MERGED_HTML_FILENAME
        if merged_file.exists():
            merged_file.unlink()
            self.logger.info(f"清理遗留合并HTML文件|Removed stale merged HTML file: {merged_file.name}")

    def _generate_index_page(self, alignments_data: Dict) -> str:
        """生成主入口页面|Generate index page"""
        output_file = self.output_dir / "index.html"
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
    
    def _generate_sample_list_html(self, alignments_data: Dict, tab_indices: Dict = None) -> str:
        """生成样品列表HTML(合并模式下详情按钮切换tab)|Generate sample list HTML"""
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

            # 合并模式切tab,旧模式跳转分页面|switch tab in merged mode, link in legacy mode
            if tab_indices is not None:
                detail_button = f'<button class="btn btn-primary" onclick="switchTab({tab_indices[sample_name]})">查看详情</button>'
            else:
                detail_button = f'<a href="{sample_name}_alignments.html" class="btn btn-primary">查看详情</a>'

            item_html = f"""
            <div class="sample-item" data-sample-name="{sample_name}" data-identity="{avg_identity:.2f}">
                <div class="sample-header">
                    <h3>{sample_name}</h3>
                    {detail_button}
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

        self.logger.info(f"  样品|Sample {sample_name}: {len(alignments)} 比对|alignments")
        return str(output_file)
    
    def _generate_alignments_html(self, alignments: List[Dict], id_prefix: str = "") -> str:
        """生成比对列表HTML(id_prefix用于合并模式跨样品去重)|Generate alignment list HTML"""
        items = []

        for idx, alignment in enumerate(alignments, 1):
            alignment_id = f"alignment_{id_prefix}_{idx}" if id_prefix else f"alignment_{idx}"
            
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
                        <button class="btn btn-secondary" onclick="copySequence('{alignment_id}')"> 复制序列</button>
                    </div>
                """
            else:
                formatted_alignment = "  序列数据不可用\n提示: 重新运行BLAST时需要在outfmt中包含 qseq 和 sseq 字段"
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
            return "  序列数据不可用\n提示: 重新运行BLAST时需要在outfmt中包含 qseq 和 sseq 字段"
        
        lines = []
        width = self.config.alignment_width
        q_start = alignment['qstart']
        s_start = alignment['sstart']
        
        # 计算实际位置（排除gap）
        q_pos = q_start
        s_pos = s_start
        # 反链判定:sstart>send时subject坐标递减(BLAST反链命中sstart>send)
        # |Reverse strand: sstart>send → subject coords decrement
        s_step = -1 if alignment.get('send', s_start) < s_start else 1

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
            s_end_pos = s_pos + s_step * (s_bases - 1) if s_bases > 0 else s_pos
            
            # 格式化输出（使用10位宽度确保对齐）
            lines.append(f"Query  {q_pos:10d}  {q_segment}  {q_end_pos}")
            lines.append(f"                   {''.join(match_line)}")
            lines.append(f"Sbjct  {s_pos:10d}  {s_segment}  {s_end_pos}")
            lines.append("")
            
            # 更新位置（只计算非gap字符;反链时subject递减）|update (subject decrements on reverse)
            q_pos += q_bases
            s_pos += s_step * s_bases
        
        return '\n'.join(lines)