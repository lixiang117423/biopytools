"""
重复序列屏蔽配置管理模块|Repeat Masking Configuration Management Module
"""

import os
from ..common.paths import expand_path
import shutil
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class RepeatMaskConfig:
    """重复序列屏蔽配置类|Repeat Masking Configuration Class"""

    # 必需参数|Required parameters
    genome: str

    # 输出配置|Output configuration
    output_dir: str = './repeatmask_output'

    # 分析参数|Analysis parameters
    threads: int = 12
    species: Optional[str] = None  # 物种名称，用于Dfam/Repbase数据库|Species name for Dfam/Repbase database

    # 软件路径配置|Software path configuration
    repeatmodeler_path: str = '~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatModeler'
    repeatmasker_path: str = '~/miniforge3/envs/repeat_identiy/bin/RepeatMasker'
    builddatabase_path: str = 'BuildDatabase'

    # 流程控制|Pipeline control
    skip_modeler: bool = False  # 跳过RepeatModeler步骤|Skip RepeatModeler step
    use_dfam: bool = True  # 是否使用Dfam数据库|Whether to use Dfam database
    masking_mode: str = 'soft'  # 屏蔽模式: soft(小写)|hard(N)|x(X)|Masking mode

    # RepeatModeler参数|RepeatModeler parameters
    use_ltr: bool = True  # 是否运行LTR结构发现|Whether to run LTR structural discovery
    modeler_quick: bool = False  # RepeatModeler快速模式|RepeatModeler quick mode

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.genome = os.path.normpath(os.path.abspath(self.genome))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 展开工具路径|Expand tool paths
        self.repeatmodeler_path = expand_path(self.repeatmodeler_path)
        self.repeatmasker_path = expand_path(self.repeatmasker_path)

        # 从基因组文件名生成base_name|Generate base_name from genome filename
        genome_name = Path(self.genome).stem
        if genome_name.endswith('.fasta') or genome_name.endswith('.fa'):
            genome_name = Path(genome_name).stem
        self.base_name = f"{genome_name}_repeat"

        # 验证masking_mode|Validate masking_mode
        valid_modes = ['soft', 'hard', 'x']
        if self.masking_mode not in valid_modes:
            raise ValueError(f"无效的屏蔽模式|Invalid masking mode: {self.masking_mode} " +
                           f"(必须是|must be one of {valid_modes})")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome file
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file does not exist: {self.genome}")

        # 检查线程数|Check thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        # 检查软件路径|Check software paths
        if not shutil.which(self.builddatabase_path.split()[0]):
            errors.append(f"BuildDatabase未找到|BuildDatabase not found: {self.builddatabase_path}")

        if not shutil.which(self.repeatmodeler_path.split()[0]):
            errors.append(f"RepeatModeler未找到|RepeatModeler not found: {self.repeatmodeler_path}")

        if not shutil.which(self.repeatmasker_path.split()[0]):
            errors.append(f"RepeatMasker未找到|RepeatMasker not found: {self.repeatmasker_path}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
