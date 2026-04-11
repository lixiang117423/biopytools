"""
合并Windows版DeepBSA结果 - 主模块|Merge Windows DeepBSA Results - Main Module

合并Windows版DeepBSA的运行结果，输出：
1. merged_results.xlsx - 各方法QTL合并表
2. plot_data_for_R.csv - 绘图数据表（含平滑值和阈值）
"""

import pandas as pd
import numpy as np
from pathlib import Path
from typing import List, Dict, Optional
import re

from .config import MergeDeepbsaConfig
from .utils import MergeDeepbsaLogger


# Windows版DeepBSA的已知方法列表
KNOWN_METHODS = ["DL", "K", "ED4", "SNP", "SmoothG", "SmoothLOD", "Ridit"]


def extract_method_from_csv_filename(filename: str) -> Optional[str]:
    """从CSV文件名提取方法名|Extract method name from CSV filename

    Windows版DeepBSA CSV文件名格式|Windows DeepBSA CSV filename format:
        0-{method}-Tri-kernel-smooth-0.1-{value}.csv
        0-ΔSNP-Tri-kernel-smooth-0.1-{value}.csv  (ΔSNP特殊处理)

    Args:
        filename: 文件名（不含路径）|Filename without path

    Returns:
        方法名或None|Method name or None
    """
    # 去掉扩展名|Remove extension
    name = Path(filename).stem

    # 匹配模式|Match pattern: 0-{method}-Tri-kernel-smooth-...
    match = re.match(r'^\d+-(.+)-Tri-kernel-smooth-', name)
    if match:
        method = match.group(1)
        # ΔSNP特殊处理|Special handling for ΔSNP
        if method == "ΔSNP":
            return "SNP"  # 内部用SNP，输出时显示为ΔSNP
        return method

    return None


def extract_method_from_values_filename(filename: str) -> Optional[str]:
    """从values.txt文件名提取方法名|Extract method name from values.txt filename

    Windows版DeepBSA values.txt文件名格式|Windows DeepBSA values.txt filename format:
        {method} values.txt
        ΔSNP values.txt

    Args:
        filename: 文件名（不含路径）|Filename without path

    Returns:
        方法名或None|Method name or None
    """
    name = filename.replace(" values.txt", "").strip()
    if name in KNOWN_METHODS:
        return name
    if name == "ΔSNP":
        return "SNP"
    return name if name else None


class WindowsDeepbsaMerger:
    """Windows版DeepBSA结果合并器|Windows DeepBSA Results Merger"""

    def __init__(self, config: MergeDeepbsaConfig, logger):
        """初始化合并器|Initialize merger

        Args:
            config: 配置对象|Config object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger
        self.methods_list = config.methods

    def _detect_methods(self) -> List[str]:
        """自动检测输入目录中的方法|Auto-detect methods in input directory

        Returns:
            检测到的方法列表|List of detected methods
        """
        methods = set()

        # 从CSV文件名检测|Detect from CSV filenames
        for csv_file in self.config.input_path.glob("*.csv"):
            method = extract_method_from_csv_filename(csv_file.name)
            if method:
                methods.add(method)

        # 从values.txt文件名检测|Detect from values.txt filenames
        for values_file in self.config.input_path.glob("* values.txt"):
            method = extract_method_from_values_filename(values_file.name)
            if method:
                methods.add(method)

        # 从npy文件名检测|Detect from npy filenames
        for npy_file in self.config.input_path.glob("all_data_for_plot_*.npy"):
            name = npy_file.stem
            # 格式|Format: all_data_for_plot_{method}
            method_name = name.replace("all_data_for_plot_", "")
            # ΔSNP -> SNP (内部用SNP)|ΔSNP -> SNP (use SNP internally)
            if method_name == "\u0394SNP":
                method_name = "SNP"
            if method_name:
                methods.add(method_name)

        return sorted(methods)

    def merge_csv_results(self) -> Optional[pd.DataFrame]:
        """合并所有方法的QTL CSV结果|Merge QTL CSV results from all methods

        Returns:
            合并后的DataFrame|Merged DataFrame
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("合并QTL结果|Merging QTL Results")
        self.logger.info("=" * 60)

        all_results = []

        for csv_file in sorted(self.config.input_path.glob("*.csv")):
            method = extract_method_from_csv_filename(csv_file.name)
            if method is None:
                self.logger.debug(f"跳过非标准CSV文件|Skipping non-standard CSV: {csv_file.name}")
                continue

            if self.methods_list and method not in self.methods_list:
                continue

            self.logger.info(f"处理|Processing: {csv_file.name} (方法|method: {method})")

            try:
                df = pd.read_csv(csv_file)

                if 'Value' not in df.columns:
                    self.logger.warning(f"  CSV缺少Value列|CSV missing Value column: {csv_file.name}")
                    continue

                # 只保留有结果的行|Only keep rows with results
                df_filtered = df[df['Value'] != '-'].copy()

                if len(df_filtered) > 0:
                    # SNP方法显示为ΔSNP|SNP method shows as ΔSNP
                    method_display = "ΔSNP" if method == "SNP" else method
                    df_filtered['Method'] = method_display
                    df_filtered['Source_File'] = csv_file.name
                    all_results.append(df_filtered)
                    self.logger.info(f"  找到|Found {len(df_filtered)} 个QTL")
                else:
                    self.logger.warning(f"  没有有结果的行|No rows with results")

            except Exception as e:
                self.logger.error(f"  读取CSV失败|Failed to read CSV: {e}")

        if not all_results:
            self.logger.error("没有找到任何QTL结果|No QTL results found")
            return None

        merged_df = pd.concat(all_results, ignore_index=True)

        # 确保列顺序|Ensure column order
        expected_cols = ['Method', 'QTL', 'Chr', 'Left', 'Peak', 'Right', 'Value', 'Source_File']
        available_cols = [c for c in expected_cols if c in merged_df.columns]
        merged_df = merged_df[available_cols]

        # 保存为Excel|Save as Excel
        output_file = self.config.output_path / "merged_results.xlsx"
        merged_df.to_excel(output_file, index=False, engine='openpyxl')

        self.logger.info("")
        self.logger.info(f"合并结果保存到|Merged results saved to: {output_file}")
        self.logger.info(f"总计|Total: {len(merged_df)} 个QTL")
        self.logger.info("")
        self.logger.info("各方法统计|Statistics by method:")
        for method in merged_df['Method'].unique():
            count = len(merged_df[merged_df['Method'] == method])
            self.logger.info(f"  {method}: {count} 个QTL")

        return merged_df

    def extract_plot_data(self) -> Optional[pd.DataFrame]:
        """提取绘图数据并计算平滑值和阈值|Extract plot data and calculate smooth values and thresholds

        优先从 npy 文件读取（与原始软件完全一致），如果 npy 不存在则 fallback 到 values.txt 重新计算。

        Returns:
            包含绘图数据的DataFrame|DataFrame with plot data
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("提取并计算绘图数据|Extract and Calculate Plot Data")
        self.logger.info("=" * 60)

        # 优先尝试从 npy 文件读取（与原始软件完全一致）
        # Try reading from npy files first (exact reproduction)
        npy_available = self._check_npy_available()
        if npy_available:
            self.logger.info("从npy文件读取绘图数据|Reading plot data from npy files")
            result_df = self._extract_plot_data_from_npy()
        else:
            self.logger.info("npy文件不存在，使用values.txt重新计算|npy files not found, recalculating from values.txt")
            result_df = self._extract_plot_data_from_values()

        if result_df is None:
            self.logger.error("没有找到任何绘图数据|No plot data found")
            return None

        # 保存为CSV|Save as CSV
        output_file = self.config.output_path / "plot_data_for_R.csv"
        result_df.to_csv(output_file, index=False)

        self.logger.info("")
        self.logger.info(f"绘图数据保存到|Plot data saved to: {output_file}")

        return result_df

    def _check_npy_available(self) -> bool:
        """检查输入目录中是否有 npy 文件|Check if npy files exist in input directory"""
        npy_files = list(self.config.input_path.glob("all_data_for_plot_*.npy"))
        return len(npy_files) > 0

    def _extract_plot_data_from_npy(self) -> Optional[pd.DataFrame]:
        """从 npy 文件读取绘图数据|Read plot data from npy files"""
        from biopytools.deepbsa.plot_data_calculator import extract_plot_data_from_npy

        return extract_plot_data_from_npy(
            input_dir=self.config.input_path,
            methods=self.methods_list,
            logger=self.logger
        )

    def _extract_plot_data_from_values(self) -> Optional[pd.DataFrame]:
        """从 values.txt 读取原始数据并重新计算|Read raw data from values.txt and recalculate"""
        from biopytools.deepbsa.plot_data_calculator import PlotDataCalculator

        calculator = PlotDataCalculator(self.logger)

        all_raw_data = []

        for values_file in sorted(self.config.input_path.glob("* values.txt")):
            method = extract_method_from_values_filename(values_file.name)
            if method is None:
                continue
            if self.methods_list and method not in self.methods_list:
                continue

            self.logger.info(f"读取|Reading: {values_file.name} (方法|method: {method})")

            try:
                data_list = []
                with open(values_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if not line:
                            continue
                        parts = line.split('\t')
                        if len(parts) == 3:
                            chromosome, position, value = parts
                            data_list.append({
                                'Method': method,
                                'Chromosome': chromosome,
                                'Position': int(position),
                                'Value': float(value)
                            })

                if data_list:
                    all_raw_data.extend(data_list)
                    self.logger.info(f"  提取|Extracted {len(data_list)} 个数据点|data points")

            except Exception as e:
                self.logger.error(f"  读取文件失败|Failed to read file: {e}")

        if not all_raw_data:
            return None

        df = pd.DataFrame(all_raw_data)
        df['Method'] = df['Method'].apply(lambda x: '\u0394SNP' if x == 'SNP' else x)

        self.logger.info(f"共提取|Total extracted: {len(df)} 个数据点|data points")
        self.logger.info(f"方法|Methods: {', '.join(sorted(df['Method'].unique()))}")

        result_df = calculator.process_all_methods(
            df,
            smooth_func=self.config.smooth_func,
            smooth_frac=self.config.smooth_frac
        )

        return result_df

    def run(self):
        """运行合并流程|Run merge process"""
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("合并Windows版DeepBSA结果|Merging Windows DeepBSA Results")
        self.logger.info("=" * 60)
        self.logger.info(f"输入目录|Input directory: {self.config.input_path}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_path}")

        # 验证配置|Validate config
        self.config.validate()

        # 自动检测方法（如果未指定）|Auto-detect methods if not specified
        if not self.methods_list:
            self.methods_list = self._detect_methods()
            self.logger.info(f"自动检测方法|Auto-detected methods: {self.methods_list}")

        results = {}

        # 合并QTL结果|Merge QTL results
        csv_df = self.merge_csv_results()
        if csv_df is not None:
            results['csv'] = True

        # 提取绘图数据|Extract plot data
        plot_df = self.extract_plot_data()
        if plot_df is not None:
            results['plot_data'] = True

        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("合并完成|Merge Completed")
        self.logger.info("=" * 60)
        self.logger.info(f"结果保存在|Results saved in: {self.config.output_path}")

        return results


def main():
    """主入口函数|Main entry function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='合并Windows版DeepBSA结果|Merge Windows DeepBSA Results',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i ./jicaiBSA_Visualize_Results -o ./merged_results
  %(prog)s -i ./results -o ./merged -m DL,K,ED4,SNP
"""
    )

    required = parser.add_argument_group('必需参数|Required arguments')
    required.add_argument(
        '-i', '--input-dir',
        required=True,
        help='Windows版DeepBSA输出目录|Windows DeepBSA output directory'
    )
    required.add_argument(
        '-o', '--output-dir',
        required=True,
        help='合并结果输出目录|Merged results output directory'
    )

    optional = parser.add_argument_group('可选参数|Optional arguments')
    optional.add_argument(
        '-m', '--methods',
        help='要合并的方法，逗号分隔（默认：全部）|Methods to merge, comma-separated (default: all)'
    )
    optional.add_argument(
        '--smooth-func',
        default='LOWESS',
        choices=['LOWESS', 'Tri-kernel'],
        help='平滑函数|Smoothing function (default: LOWESS)'
    )
    optional.add_argument(
        '--smooth-frac',
        type=float,
        default=0.1,
        help='平滑窗口比例|Smoothing window fraction (default: 0.1)'
    )

    args = parser.parse_args()

    # 解析方法列表|Parse methods list
    methods_list = None
    if args.methods:
        methods_list = [m.strip() for m in args.methods.split(',')]

    # 创建配置|Create config
    config = MergeDeepbsaConfig(
        input_dir=args.input_dir,
        output_dir=args.output_dir,
        methods=methods_list,
        smooth_func=args.smooth_func,
        smooth_frac=args.smooth_frac
    )

    # 创建日志|Create logger
    log_file = config.output_path / "merge_deepbsa.log"
    logger_mgr = MergeDeepbsaLogger(log_file=str(log_file))
    logger = logger_mgr.get_logger()

    # 执行合并|Run merge
    merger = WindowsDeepbsaMerger(config, logger)
    merger.run()


if __name__ == '__main__':
    main()
