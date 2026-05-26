"""
DeepLoc 2.1蛋白质亚细胞定位预测配置管理模块|DeepLoc 2.1 Protein Subcellular Localization Prediction Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class DeepLocConfig:
    """DeepLoc亚细胞定位预测配置类|DeepLoc Subcellular Localization Prediction Configuration Class"""

    # 必需参数|Required parameters
    fasta_file: str
    output_dir: str

    # 可选参数|Optional parameters
    model: str = "Fast"  # Accurate, Fast
    plot: bool = False  # 是否绘制attention图
    device: str = "cpu"  # cpu, cuda, mps

    # 路径配置|Path configuration
    singularity_image: str = "~/software/singularity/deeploc2.1_latest.sif"
    database_dir: str = "~/software/deeploc/database"
    singularity_exec: str = "~/miniforge3/envs/singularity_v.3.8.7/bin/singularity"

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.fasta_file = os.path.normpath(os.path.abspath(self.fasta_file))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
        self.singularity_image = os.path.normpath(os.path.abspath(self.singularity_image))
        self.database_dir = os.path.normpath(os.path.abspath(self.database_dir))
        self.singularity_exec = os.path.normpath(os.path.abspath(self.singularity_exec))

        # 验证model参数|Validate model parameter
        valid_models = ["Accurate", "Fast"]
        if self.model not in valid_models:
            raise ValueError(f"无效的model参数|Invalid model parameter: {self.model} "
                           f"(必须是|must be one of: {', '.join(valid_models)})")

        # 验证device参数|Validate device parameter
        valid_devices = ["cpu", "cuda", "mps"]
        if self.device not in valid_devices:
            raise ValueError(f"无效的device参数|Invalid device parameter: {self.device} "
                           f"(必须是|must be one of: {', '.join(valid_devices)})")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not os.path.exists(self.fasta_file):
            errors.append(f"FASTA文件不存在|FASTA file not found: {self.fasta_file}")

        # 检查FASTA文件扩展名|Check FASTA file extension
        valid_extensions = [".fa", ".faa", ".fasta", ".ffn", ".fna"]
        file_ext = os.path.splitext(self.fasta_file)[1].lower()
        if file_ext and file_ext not in valid_extensions:
            errors.append(f"输入文件应为FASTA格式|Input file should be in FASTA format: {self.fasta_file}")

        # 检查singularity镜像|Check singularity image
        if not os.path.exists(self.singularity_image):
            errors.append(f"Singularity镜像不存在|Singularity image not found: {self.singularity_image}")

        # 检查数据库目录|Check database directory
        if not os.path.exists(self.database_dir):
            errors.append(f"数据库目录不存在|Database directory not found: {self.database_dir}")

        # 检查singularity可执行文件|Check singularity executable
        if not os.path.exists(self.singularity_exec):
            errors.append(f"Singularity可执行文件不存在|Singularity executable not found: {self.singularity_exec}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
