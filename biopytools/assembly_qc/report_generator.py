"""
报告生成模块|Report Generation Module
"""

import os
from pathlib import Path
from typing import Dict, Any
from datetime import datetime


class ReportGenerator:
    """报告生成器|Report Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_reports(self, results: Dict[str, Any]):
        """生成所有报告|Generate all reports"""
        self.logger.info("生成报告|Generating reports")

        # 生成HTML报告|Generate HTML report
        if self.config.generate_html:
            self._generate_html_report(results)

        # 生成表格|Generate table
        if self.config.generate_table:
            self._generate_publication_table(results)

    def _generate_html_report(self, results: Dict[str, Any]):
        """生成HTML报告|Generate HTML report"""
        self.logger.info("生成HTML报告|Generating HTML report")

        html_file = os.path.join(self.config.output_dir, "assembly_qc_report.html")

        html_content = self._build_html_content(results)

        with open(html_file, 'w', encoding='utf-8') as f:
            f.write(html_content)

        self.logger.info(f"HTML报告已生成|HTML report generated: {html_file}")

    def _build_html_content(self, results: Dict[str, Any]) -> str:
        """构建HTML内容|Build HTML content"""
        genome_stats = results.get('genome_stats', {})
        busco_results = results.get('busco', {})
        lai_results = results.get('lai', {})
        qv_results = results.get('qv', {})
        mapping_results = results.get('mapping', {})
        long_read_mapping_results = results.get('long_read_mapping', {})

        html = f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>基因组组装质量评估报告|Genome Assembly Quality Control Report</title>
    <style>
        body {{
            font-family: Arial, sans-serif;
            margin: 40px;
            background-color: #f5f5f5;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background-color: white;
            padding: 30px;
            border-radius: 8px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.1);
        }}
        h1 {{
            color: #2c3e50;
            border-bottom: 3px solid #3498db;
            padding-bottom: 10px;
        }}
        h2 {{
            color: #34495e;
            margin-top: 30px;
            border-left: 4px solid #3498db;
            padding-left: 10px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin: 20px 0;
        }}
        th, td {{
            border: 1px solid #ddd;
            padding: 12px;
            text-align: left;
        }}
        th {{
            background-color: #3498db;
            color: white;
        }}
        tr:nth-child(even) {{
            background-color: #f2f2f2;
        }}
        .metric {{
            display: inline-block;
            margin: 10px;
            padding: 15px;
            background-color: #ecf0f1;
            border-radius: 5px;
            min-width: 200px;
        }}
        .metric-label {{
            font-size: 14px;
            color: #7f8c8d;
        }}
        .metric-value {{
            font-size: 24px;
            font-weight: bold;
            color: #2c3e50;
            margin-top: 5px;
        }}
        .status-good {{
            color: #27ae60;
        }}
        .status-warning {{
            color: #f39c12;
        }}
        .status-error {{
            color: #e74c3c;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>基因组组装质量评估报告<br>Genome Assembly Quality Control Report</h1>

        <h2>基因组基本信息|Genome Basic Information</h2>
        <table>
            <tr><th>指标|Metric</th><th>值|Value</th></tr>
            <tr><td>样品名称|Sample Name</td><td>{self.config.sample_name}</td></tr>
            <tr><td>总大小|Total Size</td><td>{genome_stats.get('total_size_mb', 'N/A'):.2f} Mb</td></tr>
            <tr><td>Contig数量|Contig Count</td><td>{genome_stats.get('contig_count', 'N/A')}</td></tr>
            <tr><td>Contig N50</td><td>{genome_stats.get('contig_n50_mb', 'N/A'):.2f} Mb</td></tr>
            <tr><td>GC含量|GC Content</td><td>{genome_stats.get('gc_content', 'N/A'):.2f}%</td></tr>
        </table>

        <h2>核心质量评估|Core Quality Evaluation</h2>

        <h3>BUSCO完整性评估|BUSCO Completeness Evaluation</h3>
        {self._format_busco_section(busco_results)}

        <h3>LAI指数评估|LAI Index Evaluation</h3>
        {self._format_lai_section(lai_results)}

        <h2>外部数据验证|External Data Validation</h2>

        <h3>QV质量值|QV Quality Value</h3>
        {self._format_qv_section(qv_results)}

        <h3>Mapping评估|Mapping Evaluation</h3>
        {self._format_mapping_section(mapping_results)}

        <h3>三代数据Mapping评估|Long-read Mapping Evaluation</h3>
        {self._format_long_read_mapping_section(long_read_mapping_results)}

        <h2>综合评级|Overall Quality Rating</h2>
        {self._format_overall_rating(results)}

        <p style="margin-top: 50px; color: #7f8c8d; font-size: 12px;">
            报告生成时间|Report Generated: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}<br>
            BioPyTools Assembly QC Module v{self.__class__.__module__.split('.')[0] if hasattr(self, '__version__') else '1.0.0'}
        </p>
    </div>
</body>
</html>"""
        return html

    def _format_busco_section(self, busco_results: Dict[str, Any]) -> str:
        """格式化BUSCO部分|Format BUSCO section"""
        if not busco_results:
            return '<p class="status-warning">[!] 未评估|Not evaluated</p>'

        complete = busco_results.get('complete', 'N/A')
        if complete != 'N/A':
            status_class = 'status-good' if complete >= 90 else 'status-warning' if complete >= 70 else 'status-error'
            return f"""
        <div class="metric">
            <div class="metric-label">完整度|Completeness</div>
            <div class="metric-value {status_class}">{complete}%</div>
        </div>
        <div class="metric">
            <div class="metric-label">单拷贝|Single Copy</div>
            <div class="metric-value">{busco_results.get('single', 'N/A')}%</div>
        </div>
        <div class="metric">
            <div class="metric-label">碎片化|Fragmented</div>
            <div class="metric-value">{busco_results.get('fragmented', 'N/A')}%</div>
        </div>
        <div class="metric">
            <div class="metric-label">缺失|Missing</div>
            <div class="metric-value">{busco_results.get('missing', 'N/A')}%</div>
        </div>
"""
        else:
            return '<p class="status-warning">[!] 数据不可用|Data not available</p>'

    def _format_lai_section(self, lai_results: Dict[str, Any]) -> str:
        """格式化LAI部分|Format LAI section"""
        if not lai_results:
            return '<p class="status-warning">[!] 未评估|Not evaluated</p>'

        # 检查LAI是否不适用|Check if LAI is not applicable
        if lai_results.get('error'):
            error_msg = lai_results.get('error', 'Unknown error')
            note = lai_results.get('note', '')
            return f'<p class="status-warning">[!] LAI不适用|LAI Not Applicable: {error_msg}</p>'

        lai_score = lai_results.get('lai_score', 'N/A')
        if lai_score != 'N/A':
            status_class = 'status-good' if lai_score >= 15 else 'status-warning' if lai_score >= 10 else 'status-error'
            return f'<div class="metric"><div class="metric-label">LAI指数|LAI Index</div><div class="metric-value {status_class}">{lai_score}</div></div>'
        else:
            return '<p class="status-warning">[!] 数据不可用|Data not available</p>'

    def _format_lai_for_table(self, lai_results: Dict[str, Any]) -> str:
        """格式化LAI值用于表格输出|Format LAI value for table output"""
        if not lai_results:
            return 'N/A'

        # 检查LAI是否不适用|Check if LAI is not applicable
        if lai_results.get('error'):
            return 'N/A'

        lai_score = lai_results.get('lai_score', 'N/A')
        return f"{lai_score}" if lai_score != 'N/A' else 'N/A'

    def _format_qv_section(self, qv_results: Dict[str, Any]) -> str:
        """格式化QV部分|Format QV section"""
        if not qv_results:
            return '<p class="status-warning">[!] 未评估（未提供reads数据）|Not evaluated (No reads data provided)</p>'

        html = ""
        # 支持多种数据类型的QV值和错误率|Support QV values and error rates for multiple data types
        if 'ngs_qv_value' in qv_results:
            error_rate = qv_results.get('ngs_error_rate', 'N/A')
            error_str = f"{error_rate:.2e}" if error_rate != 'N/A' else 'N/A'
            html += f'<p>NGS数据QV值|NGS QV Value: <strong>{qv_results["ngs_qv_value"]}</strong> (错误率|Error rate: {error_str})</p>'
        if 'long_read_qv_value' in qv_results:
            error_rate = qv_results.get('long_read_error_rate', 'N/A')
            error_str = f"{error_rate:.2e}" if error_rate != 'N/A' else 'N/A'
            html += f'<p>三代数据QV值|Long-read QV Value: <strong>{qv_results["long_read_qv_value"]}</strong> (错误率|Error rate: {error_str})</p>'
        # 如果只有一个qv_value（向后兼容）|If only qv_value exists (backward compatibility)
        if 'qv_value' in qv_results and 'ngs_qv_value' not in qv_results and 'long_read_qv_value' not in qv_results:
            html += f'<p>QV值|QV Value: <strong>{qv_results["qv_value"]}</strong></p>'

        return html if html else '<p>N/A</p>'

    def _format_mapping_section(self, mapping_results: Dict[str, Any]) -> str:
        """格式化Mapping部分|Format mapping section"""
        if not mapping_results:
            return '<p class="status-warning">[!] 未评估（未提供reads数据）|Not evaluated (No reads data provided)</p>'

        html = ""
        if 'total_samples' in mapping_results:
            html += f"<p>样品数量|Number of Samples: {mapping_results['total_samples']}</p>"

        mapping_rate = mapping_results.get('overall_mapping_rate', mapping_results.get('mapping_rate', 'N/A'))
        mean_coverage = mapping_results.get('mean_coverage', 'N/A')

        if mapping_rate != 'N/A':
            html += f"""
        <div class="metric">
            <div class="metric-label">比对率|Mapping Rate</div>
            <div class="metric-value status-good">{mapping_rate:.2f}%</div>
        </div>
"""
        if mean_coverage != 'N/A':
            html += f"""
        <div class="metric">
            <div class="metric-label">平均覆盖度|Mean Coverage</div>
            <div class="metric-value">{mean_coverage:.2f}X</div>
        </div>
"""
        return html if html else '<p>数据不可用|Data not available</p>'

    def _format_long_read_mapping_section(self, long_read_mapping_results: Dict[str, Any]) -> str:
        """格式化三代数据Mapping部分|Format long-read mapping section"""
        if not long_read_mapping_results:
            return '<p class="status-warning">[!] 未评估（未提供三代数据）|Not evaluated (No long-read data provided)</p>'

        html = ""
        if 'total_samples' in long_read_mapping_results:
            html += f"<p>样品数量|Number of Samples: {long_read_mapping_results['total_samples']}</p>"

        total_reads = long_read_mapping_results.get('total_reads', 'N/A')
        mapped_reads = long_read_mapping_results.get('mapped_reads', 'N/A')
        overall_mapping_rate = long_read_mapping_results.get('overall_mapping_rate', 'N/A')

        if total_reads != 'N/A':
            html += f"""
        <div class="metric">
            <div class="metric-label">总reads数|Total Reads</div>
            <div class="metric-value">{total_reads:,}</div>
        </div>
"""
        if mapped_reads != 'N/A':
            html += f"""
        <div class="metric">
            <div class="metric-label">比对reads数|Mapped Reads</div>
            <div class="metric-value">{mapped_reads:,}</div>
        </div>
"""
        if overall_mapping_rate != 'N/A':
            html += f"""
        <div class="metric">
            <div class="metric-label">比对率|Mapping Rate</div>
            <div class="metric-value status-good">{overall_mapping_rate:.2f}%</div>
        </div>
"""
        return html if html else '<p>数据不可用|Data not available</p>'

    def _format_overall_rating(self, results: Dict[str, Any]) -> str:
        """格式化综合评级|Format overall rating"""
        # 简化的评级逻辑|Simplified rating logic
        busco_results = results.get('busco', {})
        complete = busco_results.get('complete', 0)

        if complete >= 95:
            rating = "[5] 优秀|Excellent"
            color = "status-good"
        elif complete >= 90:
            rating = "[4] 良好|Good"
            color = "status-good"
        elif complete >= 80:
            rating = "[3] 一般|Fair"
            color = "status-warning"
        elif complete >= 70:
            rating = "[2] 较差|Poor"
            color = "status-warning"
        else:
            rating = "[1] 很差|Very Poor"
            color = "status-error"

        return f'<p class="{color}" style="font-size: 24px; font-weight: bold;">{rating}</p>'

    def _generate_publication_table(self, results: Dict[str, Any]):
        """生成发表用表格|Generate publication table"""
        self.logger.info("生成发表用表格|Generating publication table")

        genome_stats = results.get('genome_stats', {})
        busco_results = results.get('busco', {})
        lai_results = results.get('lai', {})
        qv_results = results.get('qv', {})
        mapping_results = results.get('mapping', {})
        long_read_mapping_results = results.get('long_read_mapping', {})

        # 准备表格数据|Prepare table data
        table_data = {
            'sample': self.config.sample_name,
            'size': f"{genome_stats.get('total_size_mb', 'N/A'):.2f}" if genome_stats.get('total_size_mb') != 'N/A' else 'N/A',
            'n50': f"{genome_stats.get('contig_n50_mb', 'N/A'):.2f}" if genome_stats.get('contig_n50_mb') != 'N/A' else 'N/A',
            'gc': f"{genome_stats.get('gc_content', 'N/A'):.1f}" if genome_stats.get('gc_content') != 'N/A' else 'N/A',
            'busco': f"{busco_results.get('complete', 'N/A')}" if busco_results else 'N/A',
            'lai': self._format_lai_for_table(lai_results),
            'qv_ngs': f"{qv_results.get('ngs_qv_value', 'N/A')}" if qv_results and 'ngs_qv_value' in qv_results else 'N/A',
            'qv_long': f"{qv_results.get('long_read_qv_value', 'N/A')}" if qv_results and 'long_read_qv_value' in qv_results else 'N/A',
            'qv_ngs_error': f"{qv_results.get('ngs_error_rate', 'N/A'):.2e}" if qv_results and 'ngs_error_rate' in qv_results else 'N/A',
            'qv_long_error': f"{qv_results.get('long_read_error_rate', 'N/A'):.2e}" if qv_results and 'long_read_error_rate' in qv_results else 'N/A',
            'ngs_mapping_rate': f"{mapping_results.get('overall_mapping_rate', mapping_results.get('mapping_rate', 'N/A')):.2f}" if mapping_results and mapping_results.get('overall_mapping_rate', 'N/A') != 'N/A' else 'N/A',
            'ngs_mapping_cov': self._extract_coverage_fraction(mapping_results),
            'long_mapping_rate': f"{long_read_mapping_results.get('overall_mapping_rate', 'N/A'):.2f}" if long_read_mapping_results and long_read_mapping_results.get('overall_mapping_rate', 'N/A') != 'N/A' else 'N/A',
            'long_mapping_cov': self._extract_coverage_fraction(long_read_mapping_results),
        }

        # 生成TSV|Generate TSV
        if self.config.table_format in ['tsv', 'both']:
            tsv_file = os.path.join(self.config.output_dir, "assembly_qc_table.tsv")
            self._write_tsv(table_data, tsv_file)
            self.logger.info(f"TSV表格已生成|TSV table generated: {tsv_file}")

        # 生成Excel|Generate Excel
        if self.config.table_format in ['xlsx', 'both']:
            try:
                import pandas as pd
                xlsx_file = os.path.join(self.config.output_dir, "assembly_qc_table.xlsx")

                df = pd.DataFrame([table_data])
                df.to_excel(xlsx_file, index=False, engine='openpyxl')

                self.logger.info(f"Excel表格已生成|Excel table generated: {xlsx_file}")
            except ImportError:
                self.logger.warning("pandas或openpyxl未安装，无法生成Excel表格|pandas or openpyxl not installed, cannot generate Excel table")

    def _extract_coverage_fraction(self, results: Dict[str, Any]) -> str:
        """
        从结果中提取覆盖比例|Extract coverage fraction from results
        覆盖比例存储在各个样品的coverage_fraction字段中|Coverage fraction stored in each sample's coverage_fraction field
        """
        if not results:
            return 'N/A'

        # 检查是否有samples字段（多样品情况）|Check if has samples field (multi-sample case)
        if 'samples' in results:
            samples = results['samples']
            if not samples:
                return 'N/A'

            # 计算所有样品的平均覆盖比例|Calculate average coverage fraction across all samples
            coverages = []
            for sample_name, sample_data in samples.items():
                if 'coverage_fraction' in sample_data:
                    coverages.append(sample_data['coverage_fraction'])

            if coverages:
                avg_coverage = sum(coverages) / len(coverages)
                return f"{avg_coverage:.2f}"
            else:
                return 'N/A'

        # 单样品情况|Single sample case
        elif 'coverage_fraction' in results:
            return f"{results['coverage_fraction']:.2f}"

        return 'N/A'

    def _write_tsv(self, data: dict, filename: str):
        """写入TSV文件|Write TSV file"""
        with open(filename, 'w', encoding='utf-8') as f:
            # 标题行|Header row
            f.write('\t'.join(data.keys()) + '\n')
            # 数据行|Data row
            f.write('\t'.join(str(v) for v in data.values()) + '\n')
