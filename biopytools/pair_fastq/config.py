"""
FASTQ配对修复配置管理模块|FASTQ Pair Fixing Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class FastqPairConfig:
    """FASTQ配对修复配置类|FASTQ Pair Fixing Configuration Class"""

    # 必需参数|Required parameters
    input_dir: str
    output_dir: str

    # 文件名模式|Filename patterns
    suffix1: str = "_1.fq.gz"
    suffix2: str = "_2.fq.gz"

    # 处理参数|Processing parameters
    threads: int = 12
    tool: str = "repair"  # 工具选择|Tool selection: seqkit or repair (默认repair)

    # 工具路径|Tool paths
    seqkit_bin: str = "seqkit"
    repair_sh: str = "repair.sh"  # repair.sh脚本名（通过conda run调用）
    repair_conda_env: str = "bbmap_v.39.81"  # conda环境名称

    # repair.sh特定参数|repair.sh specific parameters
    repair_memory: str = "300g"  # repair.sh的内存参数|Memory for repair.sh (default: 300g)
    save_singletons: bool = True  # 是否保存未配对的reads|Whether to save singleton reads

    # 高级选项|Advanced options
    dry_run: bool = False
    verbose: bool = False

    # 日志选项|Logging options
    log_file: Optional[str] = None
    log_level: str = "INFO"  # DEBUG, INFO, WARNING, ERROR, CRITICAL

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.input_dir = os.path.normpath(os.path.abspath(self.input_dir))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 创建输出路径对象|Create output path object
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 展开工具路径|Expand tool path
        self.seqkit_bin = os.path.normpath(expand_path(self.seqkit_bin))
        # repair_sh只保存脚本名，不需要展开路径
        # repair_sh: 只保存脚本名，通过conda run调用

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入目录|Check input directory
        if not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        if not os.path.isdir(self.input_dir):
            errors.append(f"输入路径不是目录|Input path is not a directory: {self.input_dir}")

        # 检查目录是否为空|Check if directory is empty
        if os.path.exists(self.input_dir):
            if not any(os.scandir(self.input_dir)):
                errors.append(f"输入目录为空|Input directory is empty: {self.input_dir}")

        # 检查线程数|Check thread count
        if self.threads < 1:
            errors.append("线程数必须大于0|Threads must be greater than 0")

        # 检查后缀格式|Check suffix format
        if not self.suffix1 or not self.suffix2:
            errors.append("文件后缀不能为空|File suffix cannot be empty")

        # 检查工具可用性|Check tool availability
        if self.tool == "seqkit":
            if not self._check_seqkit():
                errors.append(f"seqkit工具未找到或不可执行|seqkit tool not found or not executable: {self.seqkit_bin}")
        elif self.tool == "repair":
            if not self._check_conda():
                errors.append(f"conda环境未找到或不可用|conda environment not found or unavailable: {self.repair_conda_env}")
        else:
            errors.append(f"无效的工具选择|Invalid tool selection: {self.tool} (必须为seqkit或repair|must be seqkit or repair)")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def _check_seqkit(self):
        """检查seqkit是否可用|Check if seqkit is available"""
        try:
            import subprocess
            result = subprocess.run(
                ["which", self.seqkit_bin],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except Exception:
            return False

    def _check_repair(self):
        """检查repair.sh是否可用|Check if repair.sh is available"""
        try:
            import subprocess
            result = subprocess.run(
                ["which", self.repair_sh],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except Exception:
            return False

    def _check_conda(self):
        """检查conda环境是否可用|Check if conda environment is available"""
        try:
            import subprocess
            result = subprocess.run(
                ["conda", "run", "-n", self.repair_conda_env, "--no-capture-output", "which", self.repair_sh],
                capture_output=True,
                text=True
            )
            return result.returncode == 0
        except Exception:
            return False
