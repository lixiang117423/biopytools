"""
DeepBSA结果合并模块|DeepBSA Results Merger Module
"""

import logging
import shutil
from pathlib import Path
from typing import List, Dict
import pandas as pd
from PIL import Image
import re

from .plot_data_calculator import extract_and_calculate_plot_data


class DeepBSAMerger:
    """DeepBSA结果合并类|DeepBSA Results Merger Class"""

    def __init__(self, output_dir: Path, logger):
        """初始化合并器|Initialize merger

        Args:
            output_dir: DeepBSA输出目录|DeepBSA output directory
            logger: 日志器|Logger
        """
        self.output_dir = output_dir
        self.logger = logger
        self.merged_dir = output_dir / "merged_results"
        self.merged_dir.mkdir(exist_ok=True)

        # 各方法结果存放在each/目录下|Individual method results are in 'each/' directory
        self.each_dir = output_dir / "each"

    def merge_csv_results(self, methods: List[str]) -> Path:
        """合并所有方法的CSV结果|Merge CSV results from all methods

        Args:
            methods: 方法列表|Method list

        Returns:
            Path: 合并后的CSV文件路径|Merged CSV file path
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("合并CSV结果|Merging CSV Results")
        self.logger.info("=" * 60)

        all_results = []

        for method in methods:
            method_dir = self.each_dir / method
            if not method_dir.exists():
                self.logger.warning(f"方法目录不存在|Method directory not found: {method_dir}")
                continue

            # 查找CSV文件|Find CSV files
            csv_files = list(method_dir.glob("Results/variant/*.csv"))

            if not csv_files:
                self.logger.warning(f"未找到CSV文件|No CSV files found for {method}")
                continue

            for csv_file in csv_files:
                self.logger.info(f"处理|Processing: {csv_file}")

                try:
                    # 读取CSV|Read CSV
                    df = pd.read_csv(csv_file)

                    # 只保留有结果的行（Value列不是"-"）|Only keep rows with results (Value != "-")
                    df_filtered = df[df['Value'] != '-'].copy()

                    if len(df_filtered) > 0:
                        # 添加方法名列，SNP方法显示为δSNP|Add method column, SNP shows as δSNP
                        method_name = "δSNP" if method == "SNP" else method
                        df_filtered['Method'] = method_name
                        df_filtered['Source_File'] = csv_file.name

                        all_results.append(df_filtered)
                        self.logger.info(f"  找到|Found {len(df_filtered)} 个有结果的QTL|QTLs with results")
                    else:
                        self.logger.warning(f"  没有有结果的行|No rows with results")

                except Exception as e:
                    self.logger.error(f"  读取CSV失败|Failed to read CSV: {e}")

        if not all_results:
            self.logger.error("没有找到任何结果|No results found")
            return None

        # 合并所有结果|Merge all results
        merged_df = pd.concat(all_results, ignore_index=True)

        # 重新排列列：Method, QTL, Chr, Left, Peak, Right, Value, Source_File
        merged_df = merged_df[['Method', 'QTL', 'Chr', 'Left', 'Peak', 'Right', 'Value', 'Source_File']]

        # 保存合并的CSV|Save merged CSV
        output_csv = self.merged_dir / "merged_results.csv"
        merged_df.to_csv(output_csv, index=False)

        self.logger.info("")
        self.logger.info(f"合并结果保存到|Merged results saved to: {output_csv}")
        self.logger.info(f"总计|Total: {len(merged_df)} 个QTL|QTLs")
        self.logger.info("")
        self.logger.info("各方法统计|Statistics by method:")
        for method in merged_df['Method'].unique():
            count = len(merged_df[merged_df['Method'] == method])
            self.logger.info(f"  {method}: {count} 个QTL|QTLs")

        return output_csv

    def merge_png_images(self, methods: List[str]) -> Path:
        """合并所有方法的PNG图片为PDF文件（每页一个图）|Merge PNG images from all methods into PDF (one image per page)

        Args:
            methods: 方法列表|Method list

        Returns:
            Path: 合并后的PDF文件路径|Merged PDF file path
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("合并PNG图片为PDF|Merging PNG Images into PDF")
        self.logger.info("=" * 60)

        # 收集所有PNG图片|Collect all PNG images
        all_images = []

        for method in methods:
            method_dir = self.each_dir / method
            if not method_dir.exists():
                continue

            # 查找PNG文件|Find PNG files
            png_files = list(method_dir.glob("Results/variant/*.png"))

            if not png_files:
                self.logger.warning(f"未找到PNG文件|No PNG files found for {method}")
                continue

            for png_file in png_files:
                # 读取图片|Read image
                img = Image.open(png_file)
                # 转换为RGB模式（处理RGBA等格式）|Convert to RGB mode (handle RGBA, etc.)
                if img.mode != 'RGB':
                    img = img.convert('RGB')

                method_name = "δSNP" if method == "SNP" else method
                all_images.append((method_name, img))
                self.logger.info(f"  {method}: {png_file.name}")

        if not all_images:
            self.logger.error("没有找到任何PNG图片|No PNG images found")
            return None

        # 保存为PDF（每页一个图）|Save as PDF (one image per page)
        output_file = self.merged_dir / "all_methods_merged.pdf"

        # 将第一张图片保存为PDF|Save first image as PDF
        if all_images:
            first_img = all_images[0][1]
            first_img.save(output_file, save_all=True, append_images=[img for _, img in all_images[1:]])

        self.logger.info("")
        self.logger.info(f"PDF保存到|PDF saved to: {output_file}")
        self.logger.info(f"总计|Total: {len(all_images)} 页|pages")

        # 显示每页的方法名|Show method names for each page
        self.logger.info("")
        self.logger.info("PDF页面内容|PDF pages content:")
        for idx, (method_name, _) in enumerate(all_images, 1):
            self.logger.info(f"  第{idx}页|Page {idx}: {method_name}")

        return output_file

    def create_summary_report(self, methods: List[str]) -> Path:
        """创建汇总报告|Create summary report

        Args:
            methods: 方法列表|Method list

        Returns:
            Path: 报告文件路径|Report file path
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("创建汇总报告|Creating Summary Report")
        self.logger.info("=" * 60)

        report_file = self.merged_dir / "summary_report.txt"

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("DeepBSA分析汇总报告|DeepBSA Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"输出目录|Output directory: {self.output_dir}\n")
            f.write(f"分析方法|Methods: {', '.join(methods)}\n\n")

            # 统计每个方法的结果
            f.write("-" * 80 + "\n")
            f.write("各方法QTL统计|QTL Statistics by Method\n")
            f.write("-" * 80 + "\n\n")

            for method in methods:
                method_dir = self.each_dir / method
                if not method_dir.exists():
                    f.write(f"{method}: 目录不存在|Directory not found\n\n")
                    continue

                csv_files = list(method_dir.glob("Results/variant/*.csv"))
                if not csv_files:
                    f.write(f"{method}: 未找到结果|No results found\n\n")
                    continue

                for csv_file in csv_files:
                    try:
                        df = pd.read_csv(csv_file)
                        df_valid = df[df['Value'] != '-']

                        # SNP方法显示为δSNP|SNP method shows as δSNP
                        method_name = "δSNP" if method == "SNP" else method

                        f.write(f"{method_name} ({csv_file.name}):\n")
                        f.write(f"  总QTL数|Total QTLs: {len(df)}\n")
                        f.write(f"  有效QTL数|Valid QTLs: {len(df_valid)}\n")

                        if len(df_valid) > 0:
                            f.write(f"  最高值|Max Value: {df_valid['Value'].max()}\n")
                            f.write(f"  最低值|Min Value: {df_valid['Value'].min()}\n")

                        f.write("\n")

                    except Exception as e:
                        f.write(f"{method}: 读取失败|Failed to read: {e}\n\n")

        self.logger.info(f"报告保存到|Report saved to: {report_file}")

        return report_file

    def extract_plot_data(self, methods: List[str]) -> Path:
        """提取并计算绘图数据（包含原始值、平滑值、阈值）|Extract and calculate plot data (raw, smooth, threshold)

        在Python中完成所有计算：
        - 提取原始数据（values.txt）
        - 计算平滑曲线（LOWESS或Tri-cube核）
        - 计算阈值（中位数+3倍标准差）
        - R只需负责可视化

        Args:
            methods: 方法列表|Method list

        Returns:
            Path: 绘图数据表路径|Plot data table path
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("提取并计算绘图数据|Extract and Calculate Plot Data")
        self.logger.info("=" * 60)
        self.logger.info("在Python中完成所有计算，R只负责可视化|All calculations in Python, R only for visualization")

        # 提取原始数据并计算平滑曲线和阈值
        # Extract raw data and calculate smooth curves and thresholds
        df = extract_and_calculate_plot_data(
            each_dir=self.each_dir,
            methods=methods,
            logger=self.logger,
            smooth_func="LOWESS",    # 使用LOWESS平滑|Use LOWESS smoothing
            smooth_frac=0.1          # 窗口大小10%|Window size 10%
        )

        if df is None:
            self.logger.error("没有找到任何绘图数据|No plot data found")
            return None

        # 保存为CSV|Save as CSV
        output_file = self.merged_dir / "plot_data_for_R.csv"
        df.to_csv(output_file, index=False)

        self.logger.info("")
        self.logger.info(f"绘图数据表保存到|Plot data table saved to: {output_file}")
        self.logger.info(f"总计|Total: {len(df)} 个数据点|data points")
        self.logger.info("")
        self.logger.info("数据列|Data columns:")
        self.logger.info(f"  - Method: 方法|Method ({len(df['Method'].unique())} 个方法)")
        self.logger.info(f"  - Chromosome: 染色体|Chromosome ({len(df['Chromosome'].unique())} 个染色体)")
        self.logger.info(f"  - Position: 位置(bp)|Position(bp)")
        self.logger.info(f"  - Position_Mb: 位置(Mb)|Position(Mb)")
        self.logger.info(f"  - Raw_Value: 原始值|Raw value")
        self.logger.info(f"  - Smooth_Value: 平滑值|Smoothed value")
        self.logger.info(f"  - Threshold: 阈值|Threshold")
        self.logger.info(f"  - Above_Threshold: 是否超过阈值|Above threshold or not")
        self.logger.info("")
        self.logger.info("各方法阈值|Thresholds by method:")
        for method in sorted(df['Method'].unique()):
            threshold = df[df['Method'] == method]['Threshold'].iloc[0]
            count_above = df[df['Method'] == method]['Above_Threshold'].sum()
            self.logger.info(f"  {method}: 阈值|threshold={threshold:.4f}, 超阈值点|above threshold={count_above}")

        return output_file

    def run(self, methods: List[str]) -> Dict[str, Path]:
        """运行合并流程|Run merge process

        Args:
            methods: 方法列表|Method list

        Returns:
            dict: 合并结果文件路径|Merged result file paths
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("开始合并DeepBSA结果|Starting DeepBSA Results Merge")
        self.logger.info("=" * 60)

        results = {}

        # 合并CSV|Merge CSV
        csv_file = self.merge_csv_results(methods)
        if csv_file:
            results['csv'] = csv_file

        # 合并PNG为PDF（每页一个图）|Merge PNG into PDF (one image per page)
        merged_pdf = self.merge_png_images(methods)
        if merged_pdf:
            results['pdf'] = merged_pdf

        # 提取绘图数据|Extract plot data for R
        plot_data_file = self.extract_plot_data(methods)
        if plot_data_file:
            results['plot_data'] = plot_data_file

        # 创建报告|Create report
        report_file = self.create_summary_report(methods)
        if report_file:
            results['report'] = report_file

        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("合并完成|Merge Completed")
        self.logger.info("=" * 60)
        self.logger.info(f"合并结果保存在|Merged results saved in: {self.merged_dir}")
        self.logger.info(f"各方法结果保存在|Individual results saved in: {self.each_dir}")

        return results
