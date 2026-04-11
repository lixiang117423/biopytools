"""
合并Windows版DeepBSA结果 - 配置模块|Merge Windows DeepBSA Results - Config Module
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List


@dataclass
class MergeDeepbsaConfig:
    """合并Windows版DeepBSA结果配置类|Config for merging Windows DeepBSA results"""

    # 必需参数|Required parameters
    input_dir: str
    output_dir: str

    # 可选参数|Optional parameters
    methods: Optional[List[str]] = None
    smooth_func: str = "LOWESS"
    smooth_frac: float = 0.1

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_path = Path(self.input_dir)
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        if not self.input_path.exists():
            errors.append(f"输入目录不存在|Input directory not found: {self.input_dir}")
        if not self.input_path.is_dir():
            errors.append(f"输入路径不是目录|Input path is not a directory: {self.input_dir}")

        csv_files = list(self.input_path.glob("*.csv"))
        values_files = list(self.input_path.glob("* values.txt"))
        npy_files = list(self.input_path.glob("all_data_for_plot_*.npy"))
        if not csv_files and not values_files and not npy_files:
            errors.append(
                f"输入目录中未找到CSV、values.txt或npy文件"
                f"|No CSV, values.txt, or npy files found in: {self.input_dir}"
            )

        if errors:
            raise ValueError("\n".join(errors))
