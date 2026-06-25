"""
DeepBSA结果合并模块|DeepBSA Results Merger Module
"""

import logging
import shutil
from pathlib import Path
from typing import List, Dict, Optional
import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont
import re

from .plot_data_calculator import extract_and_calculate_plot_data


def add_method_label_to_image(img, method_name: str, label_height: int = 60) -> Image.Image:
    """在图片顶部添加方法名标签|Add method name label at the top of image

    Args:
        img: 原始图片|Original image
        method_name: 方法名|Method name
        label_height: 标签区域高度|Height of label area

    Returns:
        Image: 添加了标签的图片|Image with label added
    """
    # 获取图片尺寸|Get image dimensions
    width, height = img.size

    # 创建新图片（增加标签区域高度）|Create new image (add label area height)
    new_img = Image.new('RGB', (width, height + label_height), color='white')
    new_img.paste(img, (0, label_height))

    # 在标签区域绘制文字|Draw text in label area
    draw = ImageDraw.Draw(new_img)

    # 尝试使用系统字体，如果失败则使用默认字体
    # Try to use system font, fallback to default if failed
    try:
        # 尝试加载中文字体|Try to load Chinese font
        font = ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans-Bold.ttf", 32)
    except Exception:
        try:
            font = ImageFont.truetype("/System/Library/Fonts/Helvetica.ttc", 32)
        except Exception:
            # 使用默认字体|Use default font
            font = ImageFont.load_default()

    # 计算文字位置（居中）|Calculate text position (center)
    text_bbox = draw.textbbox((0, 0), method_name, font=font)
    text_width = text_bbox[2] - text_bbox[0]
    text_x = (width - text_width) // 2
    text_y = (label_height - 40) // 2

    # 绘制方法名|Draw method name
    draw.text((text_x, text_y), method_name, fill='black', font=font)

    return new_img


class DeepBSAMerger:
    """DeepBSA结果合并类|DeepBSA Results Merger Class"""

    # 可用方法列表|Available methods list
    AVAILABLE_METHODS = ["DL", "K", "ED4", "SNP", "SmoothG", "SmoothLOD", "Ridit"]

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

        # 检测路径结构|Detect path structure
        # batch模式: output/{method}/Results/
        # 普通模式: output/methods/{method}/Results/
        methods_dir = output_dir / "methods"
        if methods_dir.exists() and any(methods_dir.iterdir()):
            # 普通run模式: 使用methods/目录
            # Normal run mode: use methods/ directory
            self.methods_dir = methods_dir
            self.is_batch_mode = False
        else:
            # batch模式: 直接使用output目录
            # Batch mode: use output directory directly
            self.methods_dir = output_dir
            self.is_batch_mode = True

        if self.is_batch_mode:
            self.logger.debug("检测到batch模式路径结构|Detected batch mode path structure")
        else:
            self.logger.debug("检测到普通run模式路径结构|Detected normal run mode path structure")

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
            method_dir = self.methods_dir / method
            if not method_dir.exists():
                self.logger.warning(f"方法目录不存在|Method directory not found: {method_dir}")
                continue

            # 查找CSV文件|Find CSV files (递归查找，因为子目录名基于输入文件名|Recursive search, subdir name based on input filename)
            csv_files = list(method_dir.glob("Results/**/*.csv"))

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
                        # 添加方法名列，SNP方法显示为ΔSNP|Add method column, SNP shows as ΔSNP
                        method_name = "ΔSNP" if method == "SNP" else method
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

        # 保存合并的Excel|Save merged Excel
        output_excel = self.merged_dir / "merged_results.xlsx"
        merged_df.to_excel(output_excel, index=False, engine='openpyxl')

        self.logger.info("")
        self.logger.info(f"合并结果保存到|Merged results saved to:")
        self.logger.info(f"  Excel: {output_excel}")
        self.logger.info(f"总计|Total: {len(merged_df)} 个QTL|QTLs")
        self.logger.info("")
        self.logger.info("各方法统计|Statistics by method:")
        for method in merged_df['Method'].unique():
            count = len(merged_df[merged_df['Method'] == method])
            self.logger.info(f"  {method}: {count} 个QTL|QTLs")

        return output_excel

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
            method_dir = self.methods_dir / method
            if not method_dir.exists():
                continue

            # 查找PNG文件|Find PNG files (递归查找，因为子目录名基于输入文件名|Recursive search, subdir name based on input filename)
            png_files = list(method_dir.glob("Results/**/*.png"))

            if not png_files:
                self.logger.warning(f"未找到PNG文件|No PNG files found for {method}")
                continue

            for png_file in png_files:
                # 读取图片|Read image
                img = Image.open(png_file)
                # 转换为RGB模式（处理RGBA等格式）|Convert to RGB mode (handle RGBA, etc.)
                if img.mode != 'RGB':
                    img = img.convert('RGB')

                method_name = "ΔSNP" if method == "SNP" else method

                # 在图片顶部添加方法名标签|Add method name label at the top of image
                img_with_label = add_method_label_to_image(img, method_name)

                all_images.append((method_name, img_with_label))
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
                method_dir = self.methods_dir / method
                if not method_dir.exists():
                    f.write(f"{method}: 目录不存在|Directory not found\n\n")
                    continue

                csv_files = list(method_dir.glob("Results/**/*.csv"))
                if not csv_files:
                    f.write(f"{method}: 未找到结果|No results found\n\n")
                    continue

                for csv_file in csv_files:
                    try:
                        df = pd.read_csv(csv_file)
                        df_valid = df[df['Value'] != '-']

                        # SNP方法显示为ΔSNP|SNP method shows as ΔSNP
                        method_name = "ΔSNP" if method == "SNP" else method

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

    def _check_npy_available(self, methods: List[str]) -> bool:
        """检查是否有可用的 npy 文件|Check if npy files are available

        Args:
            methods: 方法列表|Method list

        Returns:
            bool: 是否有 npy 文件|Whether npy files exist
        """
        for method in methods:
            method_dir = self.methods_dir / method
            npy_dir = method_dir / "Results" / "deepbsa_input"
            if npy_dir.exists():
                npy_files = list(npy_dir.glob("all_data_for_plot_*.npy"))
                if npy_files:
                    return True
        return False

    def _extract_plot_data_from_npy(self, methods: List[str]) -> Optional[pd.DataFrame]:
        """从 npy 文件读取绘图数据（批量模式，各方法在子目录中）|Read plot data from npy files (batch mode)

        批量模式路径|Batch mode path: {methods_dir}/{method}/Results/deepbsa_input/*.npy

        Args:
            methods: 方法列表|Method list

        Returns:
            pd.DataFrame: 绘图数据|Plot data
        """
        from biopytools.deepbsa.plot_data_calculator import (
            PlotDataCalculator, _read_chrom_positions_from_values
        )

        calculator = PlotDataCalculator(self.logger)
        all_rows = []

        for method in methods:
            method_dir = self.methods_dir / method
            npy_dir = method_dir / "Results" / "deepbsa_input"

            raw_path = npy_dir / f"all_data_for_plot_{method}.npy"
            smooth_path = npy_dir / f"smooth_data_for_plot_{method}.npy"

            if not raw_path.exists() or not smooth_path.exists():
                continue

            # 查找 values.txt|Find values.txt
            txt_files = list(npy_dir.glob(f"{method} values.txt"))
            if not txt_files:
                self.logger.warning(f"values.txt不存在|values.txt not found for {method}")
                continue

            values_path = txt_files[0]

            self.logger.info(f"处理|Processing: {method}")
            self.logger.info(f"  npy目录|npy dir: {npy_dir}")

            try:
                raw = np.load(raw_path, allow_pickle=True)
                smooth = np.load(smooth_path, allow_pickle=True)

                chrom_names, chrom_positions = _read_chrom_positions_from_values(values_path, self.logger)

                if len(raw) != len(smooth) or len(raw) != len(chrom_names):
                    self.logger.error(f"  染色体数量不匹配|Chromosome count mismatch for {method}")
                    continue

                # 计算阈值|Calculate threshold
                smooth_concat = np.concatenate([np.array(smooth[i], dtype=np.float64) for i in range(len(raw))])
                threshold = calculator.calculate_threshold(smooth_concat)

                display_name = "\u0394SNP" if method == "SNP" else method
                self.logger.info(f"  阈值|Threshold: {threshold:.4f}")

                above_count = 0
                for chrom_idx in range(len(raw)):
                    for pos, raw_val, smooth_val in zip(
                        chrom_positions[chrom_idx], raw[chrom_idx], smooth[chrom_idx]
                    ):
                        is_above = float(smooth_val) >= threshold
                        if is_above:
                            above_count += 1
                        all_rows.append({
                            'Method': display_name,
                            'Chromosome': chrom_names[chrom_idx],
                            'Position': int(pos),
                            'Position_Mb': pos / 1e6,
                            'Raw_Value': float(raw_val),
                            'Smooth_Value': float(smooth_val),
                            'Threshold': threshold,
                            'Above_Threshold': is_above
                        })

                self.logger.info(f"  数据点|Data points: {len(raw)} 条染色体|chromosomes, "
                                f"{sum(len(chrom_positions[i]) for i in range(len(raw)))} 个数据点|data points")
                self.logger.info(f"  超阈值点|Points above threshold: {above_count}")

            except Exception as e:
                self.logger.error(f"  处理失败|Processing failed for {method}: {e}")

        if not all_rows:
            return None

        return pd.DataFrame(all_rows)

    def extract_plot_data(self, methods: List[str]) -> Path:
        """提取并计算绘图数据（包含原始值、平滑值、阈值）|Extract and calculate plot data (raw, smooth, threshold)

        优先从 npy 文件读取（与原始软件完全一致），如果 npy 不存在则 fallback 到 values.txt 重新计算。

        Args:
            methods: 方法列表|Method list

        Returns:
            Path: 绘图数据表路径|Plot data table path
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("提取并计算绘图数据|Extract and Calculate Plot Data")
        self.logger.info("=" * 60)

        # 优先尝试从 npy 文件读取（与原始软件完全一致）
        # Try reading from npy files first (exact reproduction)
        npy_available = self._check_npy_available(methods)
        if npy_available:
            self.logger.info("从npy文件读取绘图数据|Reading plot data from npy files")
            df = self._extract_plot_data_from_npy(methods)
        else:
            self.logger.info("npy文件不存在，使用values.txt重新计算|npy files not found, recalculating from values.txt")
            df = extract_and_calculate_plot_data(
                each_dir=self.methods_dir,
                methods=methods,
                logger=self.logger,
                smooth_func="LOWESS",
                smooth_frac=0.1
            )

        if df is None:
            self.logger.error("没有找到任何绘图数据|No plot data found")
            return None

        # 保存为CSV|Save as CSV
        output_csv = self.merged_dir / "plot_data_for_R.csv"
        df.to_csv(output_csv, index=False)

        self.logger.info("")
        self.logger.info(f"绘图数据表保存到|Plot data table saved to:")
        self.logger.info(f"  CSV: {output_csv}")
        self.logger.info(f"总计|Total: {len(df)} 个数据点|data points")
        self.logger.info("")
        self.logger.info("数据列|Data columns:")
        self.logger.info(f"  - Method: 方法|Method ({len(df['Method'].unique())} 个方法|methods)")
        self.logger.info(f"  - Chromosome: 染色体|Chromosome ({len(df['Chromosome'].unique())} 个染色体|chromosomes)")
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

        return output_csv

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
        self.logger.info(f"各方法结果保存在|Individual results saved in: {self.methods_dir}")

        return results
