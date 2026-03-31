"""
DeepBSA配置管理模块|DeepBSA Configuration Management Module
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class DeepBSAConfig:
    """DeepBSA配置类|DeepBSA Configuration Class"""

    # 必需参数|Required parameters
    input_file: str = ""  # 输入VCF/CSV文件路径|Input VCF/CSV file path

    # 方法选择|Method selection
    methods: Optional[List[str]] = None  # 要运行的方法列表（None=全部）|Methods to run (None=all)

    # 输出参数|Output parameters
    output_dir: str = "deepbsa_results"  # 输出目录|Output directory

    # 运行模式|Running mode
    parallel: bool = True  # 并行运行所有方法（默认）|Run all methods in parallel (default)

    # 预处理参数|Preprocessing parameters
    auto_clean: bool = False  # 不清理VCF，让DeepBSA自己处理VCF到CSV的转换|Don't clean VCF, let DeepBSA handle VCF to CSV conversion
    keep_clean: bool = False  # 保留清理后的文件|Keep cleaned file

    # 工具路径|Tool paths
    deepbsa_path: str = "~/software/DeepBSA/DeepBSA_linux_v1.4/bin/main.py"  # DeepBSA主程序路径|DeepBSA main script path
    conda_env: str = "/share/org/YZWL/yzwl_lixg/miniforge3/envs/DeepBSA"  # conda环境路径|conda environment path

    # 其他参数|Other parameters
    verbose: bool = False  # 详细输出|Verbose output

    # 可用方法列表|Available methods list
    AVAILABLE_METHODS = ["DL", "K", "ED4", "SNP", "SmoothG", "SmoothLOD", "Ridit"]

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 路径展开|Path expansion
        if self.input_file:
            self.input_path = Path(self.input_file).expanduser().resolve()
        else:
            self.input_path = None

        self.output_path = Path(self.output_dir).expanduser().resolve()
        # 创建输出目录|Create output directory
        self.output_path.mkdir(parents=True, exist_ok=True)
        self.output_dir = str(self.output_path)

        self.deepbsa_script = Path(self.deepbsa_path).expanduser().resolve()
        self.conda_env_path = Path(self.conda_env).expanduser().resolve()

        # 从conda环境路径提取环境名|Extract environment name from conda environment path
        self.conda_env_name = self.conda_env_path.name

        # 方法列表|Methods list
        if self.methods is None:
            self.methods_list = self.AVAILABLE_METHODS
        else:
            self.methods_list = self.methods

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input file
        if not self.input_file:
            errors.append("输入文件不能为空|Input file cannot be empty")
        elif not self.input_path.exists():
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        # 检查DeepBSA脚本|Check DeepBSA script
        if not self.deepbsa_path:
            errors.append("DeepBSA脚本路径不能为空|DeepBSA script path cannot be empty")
        elif not self.deepbsa_script.exists():
            errors.append(f"DeepBSA脚本不存在|DeepBSA script does not exist: {self.deepbsa_path}")

        # 检查conda环境|Check conda environment
        if not self.conda_env:
            errors.append("Conda环境路径不能为空|Conda environment path cannot be empty")
        elif not self.conda_env_path.exists():
            errors.append(f"Conda环境不存在|Conda environment does not exist: {self.conda_env}")

        # 检查方法|Check methods
        if self.methods_list:
            invalid_methods = [m for m in self.methods_list if m not in self.AVAILABLE_METHODS]
            if invalid_methods:
                errors.append(f"无效的方法|Invalid methods: {invalid_methods}")
                errors.append(f"可用方法|Available methods: {', '.join(self.AVAILABLE_METHODS)}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
