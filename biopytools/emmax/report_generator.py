"""
📋 GWAS报告生成模块 | GWAS Report Generation Module
"""

import json
import os
from datetime import datetime
from pathlib import Path
from typing import Dict, List, Optional

import numpy as np
import pandas as pd


class ReportGenerator:
    """报告生成器 | Report Generator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def generate_excel_report(
        self, results_file: Path, significant_snps_file: Path
    ) -> bool:
        """生成Excel报告 | Generate Excel report"""
        self.logger.info("📊 生成Excel报告 | Generating Excel report")

        try:
            # 创建Excel写入器 | Create Excel writer
            with pd.ExcelWriter(self.config.excel_report, engine="openpyxl") as writer:
                # 读取数据 | Read data
                results_df = pd.read_csv(results_file, sep="\t")
                significant_snps_df = pd.read_csv(significant_snps_file, sep="\t")

                # 写入主要结果 | Write main results
                results_df.to_excel(writer, sheet_name="All_Results", index=False)

                # 写入显著结果 | Write significant results
                significant_snps_df.to_excel(
                    writer, sheet_name="Significant_SNPs", index=False
                )

                # 创建汇总表 | Create summary table
                summary_data = {
                    "Metric": [
                        "Total SNPs",
                        "Significant SNPs",
                        "Significance Rate (%)",
                        "Most Significant SNP",
                        "Min p-value",
                        "Analysis Date",
                    ],
                    "Value": [
                        len(results_df),
                        len(significant_snps_df),
                        f"{len(significant_snps_df) / len(results_df) * 100:.4f}",
                        significant_snps_df.iloc[0]["snp_id"]
                        if len(significant_snps_df) > 0
                        else "N/A",
                        f"{results_df['p_value'].min():.2e}",
                        datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
                    ],
                }

                summary_df = pd.DataFrame(summary_data)
                summary_df.to_excel(writer, sheet_name="Summary", index=False)

            self.logger.info(
                f"✅ Excel报告已生成 | Excel report generated: {self.config.excel_report}"
            )
            return True

        except Exception as e:
            self.logger.error(
                f"❌ Excel报告生成失败 | Excel report generation failed: {e}"
            )
            return False

    def generate_html_report(
        self, results_file: Path, significant_snps_file: Path
    ) -> bool:
        """生成HTML报告 | Generate HTML report"""
        self.logger.info("🌐 生成HTML报告 | Generating HTML report")

        try:
            # 读取数据 | Read data
            results_df = pd.read_csv(results_file, sep="\t")
            significant_snps_df = pd.read_csv(significant_snps_file, sep="\t")

            # 生成HTML内容 | Generate HTML content
            html_content = f"""
<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GWAS分析报告 | GWAS Analysis Report</title>
    <style>
        body {{ font-family: Arial, sans-serif; margin: 40px; }}
        .header {{ text-align: center; margin-bottom: 40px; }}
        .section {{ margin-bottom: 30px; }}
        .table {{ border-collapse: collapse; width: 100%; }}
        .table th, .table td {{ border: 1px solid #ddd; padding: 8px; text-align: left; }}
        .table th {{ background-color: #f2f2f2; }}
        .summary {{ background-color: #f9f9f9; padding: 20px; border-radius: 5px; }}
        .image {{ text-align: center; margin: 20px 0; }}
        .image img {{ max-width: 100%; height: auto; }}
    </style>
</head>
<body>
    <div class="header">
        <h1>🧬 GWAS分析报告 | GWAS Analysis Report</h1>
        <p>生成时间 | Generated: {datetime.now().strftime("%Y-%m-%d %H:%M:%S")}</p>
    </div>
    
    <div class="section summary">
        <h2>📊 分析概要 | Analysis Summary</h2>
        <table class="table">
            <tr><th>指标 | Metric</th><th>值 | Value</th></tr>
            <tr><td>总SNP数 | Total SNPs</td><td>{len(results_df):,}</td></tr>
            <tr><td>显著SNP数 | Significant SNPs</td><td>{len(significant_snps_df):,}</td></tr>
            <tr><td>显著率 | Significance Rate</td><td>{len(significant_snps_df) / len(results_df) * 100:.4f}%</td></tr>
            <tr><td>最小p值 | Min p-value</td><td>{results_df["p_value"].min():.2e}</td></tr>
            <tr><td>最显著SNP | Most Significant SNP</td><td>{significant_snps_df.iloc[0]["snp_id"] if len(significant_snps_df) > 0 else "N/A"}</td></tr>
        </table>
    </div>
    
    <div class="section">
        <h2>📈 曼哈顿图 | Manhattan Plot</h2>
        <div class="image">
            <img src="{os.path.basename(self.config.manhattan_plot)}" alt="Manhattan Plot">
        </div>
    </div>
    
    <div class="section">
        <h2>📊 QQ图 | QQ Plot</h2>
        <div class="image">
            <img src="{os.path.basename(self.config.qq_plot)}" alt="QQ Plot">
        </div>
    </div>
    
    <div class="section">
        <h2>🎯 显著SNP列表 | Significant SNPs List</h2>
        <table class="table">
            <tr>
                <th>排名 | Rank</th>
                <th>SNP ID</th>
                <th>染色体 | Chrom</th>
                <th>位置 | Position</th>
                <th>β值 | Beta</th>
                <th>标准误 | SE</th>
                <th>p值 | p-value</th>
                <th>-log10(p)</th>
            </tr>
"""

            # 添加显著SNP表格 | Add significant SNPs table
            for _, row in significant_snps_df.head(20).iterrows():  # 只显示前20个
                html_content += f"""
            <tr>
                <td>{row.get("rank", "")}</td>
                <td>{row["snp_id"]}</td>
                <td>{row["chrom"]}</td>
                <td>{row["bp"]}</td>
                <td>{row["beta"]:.6f}</td>
                <td>{row["se"]:.6f}</td>
                <td>{row["p_value"]:.2e}</td>
                <td>{row["neg_log10_p"]:.2f}</td>
            </tr>
"""

            html_content += (
                """
        </table>
        <p><em>注：仅显示前20个最显著的SNP | Note: Only top 20 most significant SNPs are shown</em></p>
    </div>
    
    <div class="section">
        <h2>📁 输出文件 | Output Files</h2>
        <ul>
            <li>曼哈顿图 | Manhattan Plot: """
                + os.path.basename(self.config.manhattan_plot)
                + """</li>
            <li>QQ图 | QQ Plot: """
                + os.path.basename(self.config.qq_plot)
                + """</li>
            <li>显著SNP | Significant SNPs: """
                + os.path.basename(self.config.significant_snps)
                + """</li>
            <li>Excel报告 | Excel Report: """
                + os.path.basename(self.config.excel_report)
                + """</li>
        </ul>
    </div>
    
    <div class="section">
        <h2>⚙️ 分析参数 | Analysis Parameters</h2>
        <ul>
            <li>MAF阈值 | MAF threshold: """
                + str(self.config.maf_threshold)
                + """</li>
            <li>缺失率阈值 | Missing threshold: """
                + str(self.config.missing_threshold)
                + """</li>
            <li>最小深度 | Min depth: """
                + str(self.config.depth_min)
                + """</li>
            <li>最大深度 | Max depth: """
                + str(self.config.depth_max)
                + """</li>
            <li>最小质量值 | Min quality: """
                + str(self.config.qual_min)
                + """</li>
            <li>显著性阈值 | Significance threshold: """
                + str(self.config.p_value_threshold)
                + """</li>
        </ul>
    </div>
    
    <footer style="text-align: center; margin-top: 50px; padding-top: 20px; border-top: 1px solid #ddd;">
        <p>🧬 AutoGWAS v1.0.0 | Generated by MiniMax Agent</p>
    </footer>
</body>
</html>
"""
            )

            # 保存HTML文件 | Save HTML file
            with open(self.config.html_report, "w", encoding="utf-8") as f:
                f.write(html_content)

            self.logger.info(
                f"✅ HTML报告已生成 | HTML report generated: {self.config.html_report}"
            )
            return True

        except Exception as e:
            self.logger.error(
                f"❌ HTML报告生成失败 | HTML report generation failed: {e}"
            )
            return False
