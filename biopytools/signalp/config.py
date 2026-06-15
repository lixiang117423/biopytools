"""
SignalP 6.0信号肽预测配置管理模块|SignalP 6.0 Signal Peptide Prediction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class SignalPConfig:
    """SignalP信号肽预测配置类|SignalP Signal Peptide Prediction Configuration Class"""

    # 必需参数|Required parameters
    fasta_file: str
    output_dir: str

    # 可选参数|Optional parameters
    organism: str = "eukarya"  # eukarya, other, euk
    format: str = "txt"  # txt, png, eps, all, none
    mode: str = "fast"  # fast, slow, slow-sequential
    bsize: int = 12
    write_procs: int = 12
    torch_num_threads: int = 12
    skip_resolve: bool = False
    model_dir: Optional[str] = None
    cleanup_plots: bool = True  # 自动删除plot文件|Auto-delete plot files

    # SignalP程序路径|SignalP program path
    signalp_path: str = "~/miniforge3/envs/signalp6/bin/signalp6"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.fasta_file = os.path.normpath(os.path.abspath(self.fasta_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        # 展开~符号并标准化路径|Expand ~ symbol and normalize path
        self.signalp_path = os.path.normpath(expand_path(self.signalp_path))

        # 处理organism参数别名|Handle organism parameter aliases
        if self.organism.lower() in ["euk", "eukarya"]:
            self.organism = "eukarya"

        # 验证format参数|Validate format parameter
        valid_formats = ["txt", "png", "eps", "all", "none"]
        if self.format not in valid_formats:
            raise ValueError(f"无效的format参数|Invalid format parameter: {self.format} "
                           f"(必须是|must be one of: {', '.join(valid_formats)})")

        # 验证mode参数|Validate mode parameter
        valid_modes = ["fast", "slow", "slow-sequential"]
        if self.mode not in valid_modes:
            raise ValueError(f"无效的mode参数|Invalid mode parameter: {self.mode} "
                           f"(必须是|must be one of: {', '.join(valid_modes)})")

        # 处理model_dir|Handle model_dir
        if self.model_dir:
            # 展开~符号并标准化路径|Expand ~ symbol and normalize path
            self.model_dir = os.path.normpath(expand_path(self.model_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.fasta_file):
            errors.append(f"FASTA文件不存在|FASTA file not found: {self.fasta_file}")

        # 检查FASTA文件扩展名|Check FASTA file extension
        valid_extensions = [".fa", ".fasta", ".faa", ".ffn", ".fna"]
        file_ext = os.path.splitext(self.fasta_file)[1].lower()
        if file_ext and file_ext not in valid_extensions:
            errors.append(f"输入文件应为FASTA格式|Input file should be in FASTA format: {self.fasta_file}")

        # 检查SignalP程序|Check SignalP program
        if not os.path.exists(self.signalp_path):
            errors.append(f"SignalP程序不存在|SignalP program not found: {self.signalp_path}")

        # 检查model_dir|Check model_dir
        if self.model_dir and not os.path.exists(self.model_dir):
            errors.append(f"模型目录不存在|Model directory not found: {self.model_dir}")

        # 检查数值参数|Check numeric parameters
        if self.bsize <= 0:
            errors.append(f"批处理大小必须为正数|Batch size must be positive: {self.bsize}")

        if self.write_procs <= 0:
            errors.append(f"写入进程数必须为正数|Write processes must be positive: {self.write_procs}")

        if self.torch_num_threads <= 0:
            errors.append(f"PyTorch线程数必须为正数|PyTorch threads must be positive: {self.torch_num_threads}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
