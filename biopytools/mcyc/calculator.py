"""
甲烷循环基因丰度矩阵计算器 | Methane Cycle Gene Abundance Matrix Calculator
"""

import os
import time
import numpy as np
import pandas as pd
from typing import Tuple, Optional


class MatrixCalculator:
    """矩阵计算器类 | Matrix Calculator Class"""

    def __init__(self, config):
        """
        初始化矩阵计算器 | Initialize matrix calculator

        Args:
            config: MCycConfig配置对象 | MCycConfig configuration object
        """
        self.config = config
        self.raw_matrix = None

    def print_log(self, msg: str, level: str = "INFO"):
        """打印日志 | Print log"""
        icon_map = {
            "INFO": "ℹ️ ", "WARN": "⚠️ ", "ERROR": "❌ ",
            "SUCCESS": "✅ ", "PROCESS": "⏳ ", "TOOL": "🛠️ "
        }
        icon = icon_map.get(level, "ℹ️ ")
        print(f"{icon} [{time.strftime('%H:%M:%S')}] {msg}")

    def load_raw_matrix(self) -> bool:
        """
        加载原始矩阵 | Load raw matrix

        Returns:
            bool: 是否加载成功 | Whether loading was successful
        """
        try:
            if not os.path.exists(self.config.raw_output):
                self.print_log(f"原始矩阵文件不存在 | Raw matrix file not found: {self.config.raw_output}", "ERROR")
                return False

            self.print_log("正在加载原始矩阵... | Loading raw matrix...", "PROCESS")

            # 使用comment='#'跳过注释行 | Use comment='#' to skip comment lines
            self.raw_matrix = pd.read_csv(
                self.config.raw_output,
                sep='\t',
                index_col=0,
                comment='#'
            )

            if self.raw_matrix.empty:
                self.print_log("读取到的矩阵为空！| Read matrix is empty!", "ERROR")
                return False

            self.print_log(f"成功读取矩阵：{self.raw_matrix.shape[0]} 个基因, {self.raw_matrix.shape[1]} 个样品 | "
                         f"Successfully loaded matrix: {self.raw_matrix.shape[0]} genes, {self.raw_matrix.shape[1]} samples",
                         "SUCCESS")

            return True

        except Exception as e:
            self.print_log(f"矩阵加载失败 | Matrix loading failed: {e}", "ERROR")
            return False

    def calculate_raw_counts(self) -> bool:
        """
        计算原始计数矩阵 | Calculate raw counts matrix

        Returns:
            bool: 是否计算成功 | Whether calculation was successful
        """
        try:
            self.print_log("输出原始计数矩阵... | Outputting raw counts matrix...", "PROCESS")

            # 保存原始计数矩阵 | Save raw counts matrix
            self.raw_matrix.to_csv(self.config.raw_counts_output, sep='\t')

            self.print_log(f"原始计数矩阵已保存到 | Raw counts matrix saved to: {self.config.raw_counts_output}", "SUCCESS")
            return True

        except Exception as e:
            self.print_log(f"原始计数计算失败 | Raw counts calculation failed: {e}", "ERROR")
            return False

    def calculate_tpm(self) -> bool:
        """
        计算TPM (CPM) 矩阵 | Calculate TPM (CPM) matrix

        Returns:
            bool: 是否计算成功 | Whether calculation was successful
        """
        try:
            self.print_log("计算TPM (CPM)... | Calculating TPM (CPM)...", "PROCESS")

            # 计算库大小 | Calculate library sizes
            lib_sizes = self.raw_matrix.sum(axis=0)
            # 避免除零错误 | Avoid division by zero
            lib_sizes[lib_sizes == 0] = 1

            # 计算TPM (CPM): (Reads / LibrarySize) * 1e6
            tpm_matrix = self.raw_matrix.div(lib_sizes, axis=1) * 1e6

            # 保存TPM矩阵 | Save TPM matrix
            tpm_matrix.to_csv(self.config.tpm_output, sep='\t')

            self.print_log(f"TPM矩阵已保存到 | TPM matrix saved to: {self.config.tpm_output}", "SUCCESS")
            return True

        except Exception as e:
            self.print_log(f"TPM计算失败 | TPM calculation failed: {e}", "ERROR")
            return False

    def calculate_clr(self) -> bool:
        """
        计算CLR矩阵 | Calculate CLR matrix

        Returns:
            bool: 是否计算成功 | Whether calculation was successful
        """
        try:
            self.print_log("计算CLR... | Calculating CLR...", "PROCESS")

            def clr_transform(col):
                """CLR变换函数 | CLR transformation function"""
                x = col + 1  # 添加伪计数避免log(0) | Add pseudocount to avoid log(0)
                gm = np.exp(np.mean(np.log(x)))
                if gm == 0:
                    return np.log(x)  # 极其罕见情况 | Very rare case
                return np.log(x / gm)

            # 对每列进行CLR变换 | Apply CLR transformation to each column
            clr_matrix = self.raw_matrix.apply(clr_transform, axis=0)

            # 保存CLR矩阵 | Save CLR matrix
            clr_matrix.to_csv(self.config.clr_output, sep='\t')

            self.print_log(f"CLR矩阵已保存到 | CLR matrix saved to: {self.config.clr_output}", "SUCCESS")
            return True

        except Exception as e:
            self.print_log(f"CLR计算失败 | CLR calculation failed: {e}", "ERROR")
            return False

    def run_all_calculations(self) -> bool:
        """
        运行所有矩阵计算 | Run all matrix calculations

        Returns:
            bool: 是否全部计算成功 | Whether all calculations were successful
        """
        self.print_log("开始矩阵计算... | Starting matrix calculations...", "PROCESS")

        # 1. 加载原始矩阵 | Load raw matrix
        if not self.load_raw_matrix():
            return False

        # 2. 计算原始计数 | Calculate raw counts
        if not self.calculate_raw_counts():
            return False

        # 3. 计算TPM | Calculate TPM
        if not self.calculate_tpm():
            return False

        # 4. 计算CLR | Calculate CLR
        if not self.calculate_clr():
            return False

        self.print_log("所有矩阵计算完成！| All matrix calculations completed!", "SUCCESS")
        return True

    def get_matrix_info(self) -> Tuple[int, int]:
        """
        获取矩阵信息 | Get matrix information

        Returns:
            Tuple[int, int]: (基因数量, 样品数量) | (Gene count, sample count)
        """
        if self.raw_matrix is not None:
            return self.raw_matrix.shape
        return (0, 0)

    def get_output_summary(self) -> dict:
        """
        获取输出文件摘要 | Get output file summary

        Returns:
            dict: 输出文件信息 | Output file information
        """
        outputs = {
            "raw_counts": {"file": self.config.raw_counts_output, "description": "原始计数 | Raw Counts"},
            "tpm": {"file": self.config.tpm_output, "description": "相对丰度 (TPM) | Relative Abundance (TPM)"},
            "clr": {"file": self.config.clr_output, "description": "CLR变换 | CLR Transformed"}
        }

        summary = {}
        for key, info in outputs.items():
            if os.path.exists(info["file"]):
                summary[key] = {
                    "exists": True,
                    "path": info["file"],
                    "description": info["description"]
                }
            else:
                summary[key] = {
                    "exists": False,
                    "path": info["file"],
                    "description": info["description"]
                }

        return summary