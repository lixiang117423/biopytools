"""
VCF2PCA配置管理模块|VCF2PCA Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

# 路径管理工具|Path management utilities
try:
    from common.paths import get_tool_path, expand_path
except ImportError:
    # 如果common模块不可用，使用简化版本|Fallback if common module unavailable
    def get_tool_path(tool_name, default_path, env_var=None):
        """获取工具路径|Get tool path"""
        if env_var and os.environ.get(env_var):
            return os.path.expandvars(os.path.expanduser(os.environ[env_var]))
        return os.path.expandvars(os.path.expanduser(default_path))

    def expand_path(path):
        """展开路径|Expand path"""
        return os.path.expandvars(os.path.expanduser(path))


@dataclass
class VCF2PCAConfig:
    """VCF2PCA配置类|VCF2PCA Configuration Class"""

    # 必需参数|Required parameters
    vcf_file: str
    output_dir: str

    # 后端选择|Backend selection
    backend: str = 'v2p'  # 'v2p' (VCF2PCACluster) 或 'plink'

    # PCA参数|PCA parameters
    components: int = 10

    # 质控参数|Quality control parameters (仅PLINK后端|PLINK backend only)
    maf: float = 0.05
    missing_rate: float = 0.1
    hwe_pvalue: float = 1e-6
    skip_qc: bool = False

    # 聚类参数|Clustering parameters (仅V2P后端|V2P backend only)
    cluster: bool = False
    cluster_method: str = 'kmeans'  # 'kmeans', 'dbscan', 'em'
    cluster_k: int = 3  # K-means的K值|K value for K-means

    # 可视化参数|Visualization parameters
    plot: bool = False
    plot_2d: bool = True  # 2D图|2D plot
    plot_3d: bool = False  # 3D图|3D plot

    # 工具路径|Tool paths
    vcf2pca_path: str = get_tool_path(
        'VCF2PCACluster',
        '~/software/VCF2PCACluster-1.42/bin/VCF2PCACluster',
        'VCF2PCA_PATH'
    )
    plink_path: str = get_tool_path(
        'plink',
        'plink',
        'PLINK_PATH'
    )
    bcftools_path: str = get_tool_path(
        'bcftools',
        'bcftools',
        'BCFTOOLS_PATH'
    )

    # 其他参数|Other parameters
    threads: int = 12
    sample_info_file: Optional[str] = None

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        self.vcf_file = expand_path(self.vcf_file)
        self.output_dir = expand_path(self.output_dir)
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 展开工具路径|Expand tool paths
        self.vcf2pca_path = expand_path(self.vcf2pca_path)
        self.plink_path = expand_path(self.plink_path)
        self.bcftools_path = expand_path(self.bcftools_path)

        # 展开样本信息文件|Expand sample info file
        if self.sample_info_file:
            self.sample_info_file = expand_path(self.sample_info_file)

        # 验证backend参数|Validate backend parameter
        if self.backend not in ['v2p', 'plink']:
            raise ValueError(
                f"不支持的backend|Unsupported backend: {self.backend}. "
                f"必须是'must be 'v2p' 或 'plink'"
            )

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查VCF文件|Check VCF file
        if not os.path.exists(self.vcf_file):
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")

        # 检查样本信息文件|Check sample info file
        if self.sample_info_file and not os.path.exists(self.sample_info_file):
            errors.append(
                f"样本信息文件不存在|Sample info file does not exist: "
                f"{self.sample_info_file}"
            )

        # 检查参数范围|Check parameter ranges
        if self.components <= 0:
            errors.append(
                f"主成分数量必须为正整数|Number of components must be positive: "
                f"{self.components}"
            )

        if not 0 < self.maf < 0.5:
            errors.append(f"MAF阈值必须在0-0.5之间|MAF threshold must be between 0-0.5: {self.maf}")

        if not 0 < self.missing_rate < 1:
            errors.append(
                f"缺失率阈值必须在0-1之间|Missing rate threshold must be between 0-1: "
                f"{self.missing_rate}"
            )

        # 检查工具路径|Check tool paths
        if self.backend == 'v2p':
            if not os.path.exists(self.vcf2pca_path):
                errors.append(
                    f"VCF2PCACluster不存在|VCF2PCACluster not found: {self.vcf2pca_path}"
                )
        elif self.backend == 'plink':
            # PLINK路径检查在运行时进行|PLINK path check at runtime
            pass

        if errors:
            raise ValueError("\n".join(errors))

        return True
