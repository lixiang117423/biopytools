"""
DeepBSA配置管理模块|DeepBSA Configuration Management Module
"""

from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
import os
from ..common.paths import expand_path


@dataclass
class BatchConfig:
    """DeepBSA批量任务配置类|DeepBSA Batch Task Configuration Class"""

    # 必需参数|Required parameters
    input_file: str  # 输入VCF/CSV文件路径|Input VCF/CSV file path
    output_dir: str  # 批量任务输出目录|Batch tasks output directory
    script_name: str = "run_deepbsa_methods.sh"  # 生成的脚本文件名|Generated script filename

    # 可选参数|Optional parameters
    methods: Optional[List[str]] = None  # 要运行的方法列表（None=全部）|Methods to run (None=all)
    threads: int = 88  # 每个方法的线程数|Threads per method
    smooth_func: str = "Tri-kernel-smooth"  # 平滑函数|Smooth function

    # 可用方法列表|Available methods list
    AVAILABLE_METHODS = ["DL", "K", "ED4", "SNP", "SmoothG", "SmoothLOD", "Ridit"]

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 路径展开|Path expansion
        self.input_path = Path(self.input_file).expanduser().resolve()
        self.output_path = Path(self.output_dir).expanduser().resolve()

        # 创建输出目录|Create output directory
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 脚本文件路径|Script file path
        self.script_path = self.output_path / self.script_name

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

        # 检查方法|Check methods
        if self.methods_list:
            invalid_methods = [m for m in self.methods_list if m not in self.AVAILABLE_METHODS]
            if invalid_methods:
                errors.append(f"无效的方法|Invalid methods: {invalid_methods}")
                errors.append(f"可用方法|Available methods: {', '.join(self.AVAILABLE_METHODS)}")

        if errors:
            raise ValueError("\n".join(errors))

        return True


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
    deepbsa_path: str = ""  # DeepBSA主程序路径（空字符串=自动检测内置版本）|DeepBSA main script path (empty string=auto-detect builtin version)
    use_builtin: bool = True  # 是否使用内置DeepBSA版本|Whether to use builtin DeepBSA version
    conda_env: str = "~/miniforge3/envs/DeepBSA"  # conda环境路径|conda environment path

    # 其他参数|Other parameters
    verbose: bool = False  # 详细输出|Verbose output
    force: bool = False  # 强制重新运行所有步骤（不跳过已完成）|Force rerun all steps (don't skip completed)
    skip_merge: bool = False  # 跳过合并步骤|Skip merge step
    threads: int = 6  # DeepBSA线程数（默认6，充分利用64核集群）|Number of threads (default: 6, fully utilize 64-core cluster)
    smooth_func: str = "Tri-kernel-smooth"  # 平滑函数|Smooth function

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

        # DeepBSA路径解析|DeepBSA path resolution
        # 优先级：1. 显式指定的deepbsa_path  2. 内置版本（如果use_builtin=True）  3. 环境变量
        # Priority: 1. Explicitly specified deepbsa_path  2. Builtin version (if use_builtin=True)  3. Environment variable
        if self.deepbsa_path and not self.use_builtin:
            # 用户显式指定了路径且不使用内置版本
            # User explicitly specified path and doesn't want to use builtin version
            self.deepbsa_script = Path(self.deepbsa_path).expanduser().resolve()
            self._using_builtin = False
        else:
            # 尝试使用内置版本|Try to use builtin version
            builtin_script = Path(__file__).parent / "deepbsa_builtin" / "main.py"
            if self.use_builtin and builtin_script.exists():
                self.deepbsa_script = builtin_script
                self._using_builtin = True
            elif self.deepbsa_path:
                # 回退到用户指定的路径
                # Fallback to user-specified path
                self.deepbsa_script = Path(self.deepbsa_path).expanduser().resolve()
                self._using_builtin = False
            else:
                # 尝试环境变量或默认路径
                # Try environment variable or default path
                env_path = os.environ.get('DEEPBSA_PATH', '~/software/DeepBSA/DeepBSA_multithread/bin/main.py')
                self.deepbsa_script = Path(env_path).expanduser().resolve()
                self._using_builtin = False

        # 展开conda环境路径|Expand conda env path
        self.conda_env = expand_path(self.conda_env)
        self.conda_env_path = Path(self.conda_env).resolve()

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
        elif self.input_path and not self.input_path.exists():
            errors.append(f"输入文件不存在|Input file does not exist: {self.input_file}")

        # 检查DeepBSA脚本|Check DeepBSA script
        if not self.deepbsa_script.exists():
            errors.append(f"DeepBSA脚本不存在|DeepBSA script does not exist: {self.deepbsa_script}")

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
