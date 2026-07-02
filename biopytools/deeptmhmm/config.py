"""
DeepTMHMM 1.0跨膜螺旋/信号肽预测配置类|DeepTMHMM 1.0 Transmembrane Helix & Signal Peptide Prediction Configuration
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

from ..common.paths import get_tool_path, expand_path


@dataclass
class DeeptmhmmConfig:
    """DeepTMHMM配置类|DeepTMHMM Configuration Class"""

    # 必需参数|Required parameters
    input_file: str
    output_dir: str

    # 可选参数|Optional parameters
    output_prefix: Optional[str] = None

    # 环境与工具路径|Environment and tool paths
    conda_env: str = 'deeptmhmm_v.1.0'
    deeptmhmm_dir: str = '~/software/deeptmhmm/DeepTMHMM-Academic-License-v1.0'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 展开用户/环境输入路径|Expand user/env input paths
        self.input_file = expand_path(self.input_file)
        self.output_dir = expand_path(self.output_dir)

        # 工具目录支持环境变量/配置文件覆盖|Tool dir supports env var/config override
        self.deeptmhmm_dir = get_tool_path(
            'deeptmhmm',
            self.deeptmhmm_dir,
            'DEEPTMHMM_DIR'
        )

        # conda环境名支持环境变量覆盖|conda env name supports env var override
        self.conda_env = os.environ.get('DEEPTMHMM_ENV', self.conda_env)

        # 派生路径|Derived paths
        self.predict_py = os.path.join(self.deeptmhmm_dir, 'predict.py')
        self.python_bin = expand_path(f'~/miniforge3/envs/{self.conda_env}/bin/python')

        # 输出前缀默认取输入文件名|Output prefix defaults to input filename
        if self.output_prefix is None:
            self.output_prefix = os.path.splitext(os.path.basename(self.input_file))[0]

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not self.input_file:
            errors.append("输入文件不能为空|Input file cannot be empty")
        elif not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        # 检查DeepTMHMM安装目录|Check DeepTMHMM install directory
        if not os.path.isdir(self.deeptmhmm_dir):
            errors.append(f"DeepTMHMM目录不存在|DeepTMHMM directory not found: {self.deeptmhmm_dir}")
        elif not os.path.exists(self.predict_py):
            errors.append(f"predict.py不存在|predict.py not found: {self.predict_py}")

        # 检查conda环境python|Check conda env python
        if not os.path.exists(self.python_bin):
            errors.append(
                f"conda环境python不存在|conda env python not found: {self.python_bin} "
                f"(环境|env: {self.conda_env})"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
