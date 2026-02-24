"""
iSeq下载配置管理模块|iSeq Download Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional


@dataclass
class ISeqConfig:
    """iSeq下载配置类|iSeq Download Configuration Class"""

    # 必需参数|Required parameters
    accession: str  # 项目/样本/实验ID|Project/Sample/Experiment ID

    # 路径配置|Path configuration
    iseq_path: str = '/share/org/YZWL/yzwl_lixg/miniforge3/envs/iseq_v.1.9.8/bin/iseq'  # iSeq软件路径|iSeq software path
    conda_env: str = 'iseq_v.1.9.8'  # conda环境名|conda environment name
    output_dir: str = './iseq_output'  # 输出目录|Output directory

    # 下载选项|Download options
    metadata_only: bool = False  # 仅下载元数据|Only download metadata
    gzip: bool = True  # 下载gzip格式FASTQ|Download FASTQ in gzip format
    fastq: bool = False  # 转换为FASTQ格式|Convert to FASTQ format
    merge: Optional[str] = None  # 合并选项(ex/sa/st)|Merge option (ex/sa/st)

    # 性能参数|Performance parameters
    threads: int = 16  # 线程数|Number of threads (for conversion/compression)
    parallel: int = 10  # 并行连接数|Number of parallel connections
    speed: Optional[int] = None  # 下载速度限制(MB/s)|Download speed limit (MB/s)

    # 数据库选项|Database options
    database: str = 'ena'  # 数据库选择(ena/sra)|Database selection (ena/sra)
    protocol: str = 'ftp'  # 协议选择(ftp/https)|Protocol selection (ftp/https)

    # 高级选项|Advanced options
    use_aspera: bool = False  # 使用Aspera下载|Use Aspera for download
    skip_md5: bool = False  # 跳过MD5校验|Skip MD5 check
    quiet: bool = False  # 静默模式|Quiet mode (suppress progress bars)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.iseq_path = os.path.normpath(os.path.abspath(self.iseq_path))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 验证merge参数|Validate merge parameter
        if self.merge is not None:
            valid_merge_values = ['ex', 'sa', 'st']
            if self.merge not in valid_merge_values:
                raise ValueError(
                    f"无效的merge参数|Invalid merge parameter: {self.merge}. "
                    f"必须是|Must be one of {valid_merge_values}"
                )

        # 验证database参数|Validate database parameter
        valid_databases = ['ena', 'sra']
        if self.database not in valid_databases:
            raise ValueError(
                f"无效的database参数|Invalid database parameter: {self.database}. "
                f"必须是|Must be one of {valid_databases}"
            )

        # 验证protocol参数|Validate protocol parameter
        valid_protocols = ['ftp', 'https']
        if self.protocol not in valid_protocols:
            raise ValueError(
                f"无效的protocol参数|Invalid protocol parameter: {self.protocol}. "
                f"必须是|Must be one of {valid_protocols}"
            )

        # gzip和fastq不能同时为True|gzip and fastq cannot both be True
        if self.gzip and self.fastq:
            raise ValueError("gzip和fastq参数不能同时为True|gzip and fastq cannot both be True")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查accession是否为空|Check if accession is empty
        if not self.accession or self.accession.strip() == "":
            errors.append("accession不能为空|accession cannot be empty")

        # 检查iseq_path是否存在|Check if iSeq path exists
        if not os.path.exists(self.iseq_path):
            errors.append(f"iSeq路径不存在|iSeq path does not exist: {self.iseq_path}")

        # 检查threads是否为正数|Check if threads is positive
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 检查parallel是否为正数|Check if parallel is positive
        if self.parallel <= 0:
            errors.append(f"并行连接数必须为正数|Parallel connections must be positive: {self.parallel}")

        # 检查speed是否为正数|Check if speed is positive
        if self.speed is not None and self.speed <= 0:
            errors.append(f"下载速度必须为正数|Download speed must be positive: {self.speed}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
