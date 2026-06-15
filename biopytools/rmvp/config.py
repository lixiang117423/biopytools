"""
rMVP配置类|rMVP Configuration Class
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List


@dataclass
class RMVPConfig:
    """rMVP GWAS分析配置类|rMVP GWAS Analysis Configuration Class"""

    # 输入文件|Input files
    vcf_file: str = ""
    pheno_file: str = ""
    output_prefix: str = "RMVP_Result"
    output_dir: str = "."

    # 分析模型|Analysis models
    models: List[str] = field(default_factory=lambda: ["GLM", "MLM", "FarmCPU"])

    # R环境|R environment
    r_env: Optional[str] = 'rMVP'  # conda环境名称或路径|conda env name or path
    r_path: Optional[str] = None  # R可执行文件路径|R executable path（运行时设置）|Set at runtime
    r_type: Optional[str] = None  # R类型标识: "conda" 或 "direct"|"conda" or "direct"

    # 运行时检测到的环境信息|Runtime detected environment info
    r_env_name: Optional[str] = None  # 检测到的conda环境名称|Detected conda env name
    r_env_type: Optional[str] = None  # 检测到的R类型|Detected R type: "conda" or "direct"

    # 并行计算|Parallel computing
    ncpus: int = 12  # CPU核心数
    maxLine: int = 10000  # 每次读取的SNP数量（影响内存）

    # PCA参数|PCA parameters
    n_pc_glm: int = 3  # GLM模型使用的PC数量
    n_pc_mlm: int = 3  # MLM模型使用的PC数量
    n_pc_farmcpu: int = 3  # FarmCPU模型使用的PC数量

    # MLM参数|MLM parameters
    vc_method: str = "BRENT"  # 方差组分分析方法: "BRENT", "EMMA", "HE"

    # FarmCPU参数|FarmCPU parameters
    max_loop: int = 10  # 最大迭代次数
    method_bin: str = "static"  # bin方法: "static" or "FaST-LLM"

    # 筛选参数|Filtering parameters
    maf: Optional[float] = None  # 最小等位基因频率阈值
    miss: Optional[float] = None  # 缺失率阈值

    # 输出控制|Output control
    file_output: List[str] = field(default_factory=lambda: ["pmap", "pmap.signal", "plot", "log"])
    file_type: str = "jpg"  # 图片格式: "jpg", "pdf", "tiff"
    dpi: int = 300  # 图片分辨率
    threshold: float = 0.05  # 显著性阈值

    # 日志|Logging
    log_level: str = "INFO"
    verbose: bool = True  # 是否输出详细信息

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.vcf_file = str(Path(self.vcf_file).expanduser()) if self.vcf_file else ""
        self.pheno_file = str(Path(self.pheno_file).expanduser()) if self.pheno_file else ""
        self.output_dir = Path(self.output_dir).expanduser()

        # 处理r_env：如果是路径，提取环境名称|Handle r_env: if it's a path, extract env name
        if self.r_env:
            if '/' in self.r_env or '\\' in self.r_env:
                # 是路径，提取环境名称|It's a path, extract env name
                env_path = Path(self.r_env)
                self.r_env = env_path.name  # 提取最后一部分作为环境名称|Extract last part as env name

        # 展开r_path（如果提供）|Expand r_path (if provided)
        if self.r_path:
            self.r_path = str(Path(self.r_path).expanduser())

        # 规范化模型名称|Normalize model names
        # rMVP的模型名称：GLM, MLM, FarmCPU（注意大小写）|rMVP model names: GLM, MLM, FarmCPU (case-sensitive)
        normalized_models = []
        for m in self.models:
            m_upper = m.upper()
            if m_upper == "FARMCPU":
                normalized_models.append("FarmCPU")  # 保持正确的大小写|Keep correct case
            else:
                normalized_models.append(m_upper)
        self.models = normalized_models

        self.vc_method = self.vc_method.upper()
        self.method_bin = self.method_bin.lower()
        self.log_level = self.log_level.upper()

        # 验证模型参数|Validate model parameters
        valid_models = {"GLM", "MLM", "FarmCPU"}
        for model in self.models:
            if model not in valid_models:
                raise ValueError(f"无效的模型|Invalid model: {model}. 必须是|Must be one of {valid_models}")

        # 确保输出目录存在|Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        if not self.vcf_file:
            raise ValueError("VCF文件不能为空|VCF file cannot be empty")

        if not Path(self.vcf_file).exists():
            raise ValueError(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        if not self.pheno_file:
            raise ValueError("表型文件不能为空|Phenotype file cannot be empty")

        if not Path(self.pheno_file).exists():
            raise ValueError(f"表型文件不存在|Phenotype file does not exist: {self.pheno_file}")

        if self.vc_method not in ["BRENT", "EMMA", "HE"]:
            raise ValueError("vc_method必须是|vc_method must be one of: BRENT, EMMA, HE")

        if self.method_bin not in ["static", "fast-lmm"]:
            raise ValueError("method_bin必须是|method_bin must be one of: static, fast-lmm")
