"""
RNA-Bloom组装配置管理模块|RNA-Bloom Assembly Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List
from ..common.paths import get_tool_path, expand_path


@dataclass
class RNABloomConfig:
    """RNA-Bloom组装配置类|RNA-Bloom Assembly Configuration Class"""

    # 输入文件|Input files
    left_reads: Optional[str] = None  # 左端reads|Left reads (paired-end)
    right_reads: Optional[str] = None  # 右端reads|Right reads (paired-end)
    single_end_forward: Optional[str] = None  # 单端正向reads|Single-end forward reads
    single_end_reverse: Optional[str] = None  # 单端反向reads|Single-end reverse reads
    long_reads: Optional[str] = None  # 长reads|Long reads
    cell_list: Optional[str] = None  # 单细胞列表文件|Single-cell list file

    # 输出目录|Output directory
    output_dir: str = './rnabloom_output'

    # 处理参数|Processing parameters
    threads: int = 12  # 线程数|Number of threads

    # RNA-Bloom路径配置|RNA-Bloom path configuration
    rnabloom_path: str = "rnabloom"  # 默认使用工具名，让系统从PATH查找|Default uses tool name, let system find from PATH

    # Bloom filter配置|Bloom filter configuration
    memory_gb: Optional[float] = None  # Bloom filter总大小(GB)|Total Bloom filter size in GB
    false_positive_rate: Optional[float] = None  # 假阳性率|False positive rate
    num_kmers: Optional[int] = None  # 唯一k-mer数量|Number of unique kmers

    # 数据类型配置|Data type configuration
    stranded: bool = False  # 链特异性|Strand-specific
    revcomp_left: bool = False  # 反向互补左端reads|Reverse-complement left reads
    revcomp_right: bool = False  # 反向互补右端reads|Reverse-complement right reads
    is_pacbio: bool = False  # 是否为PacBio数据|Whether data is from PacBio

    # 参考引导组装|Reference-guided assembly
    reference_transcripts: Optional[str] = None  # 参考转录本文件|Reference transcript file

    # 输出选项|Output options
    min_length: int = 200  # 最小转录本长度|Minimum transcript length
    write_uracil: bool = False  # 输出尿嘧啶(U)而非胸腺嘧啶(T)|Write uracil (U) instead of thymine (T)
    export_non_redundant: bool = True  # 导出去冗余转录本|Export non-redundant transcripts

    # 步骤控制|Step control
    stage: Optional[int] = None  # 停止阶段|Stop at stage (1-3)

    # 跳过检查|Skip checks
    skip_auto_bloom_filter_size: bool = True  # 跳过自动Bloom filter大小设置|Skip auto Bloom filter size

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.output_dir = expand_path(self.output_dir)

        # 展开输入文件路径|Expand input file paths
        if self.left_reads:
            self.left_reads = expand_path(self.left_reads)
        if self.right_reads:
            self.right_reads = expand_path(self.right_reads)
        if self.single_end_forward:
            self.single_end_forward = expand_path(self.single_end_forward)
        if self.single_end_reverse:
            self.single_end_reverse = expand_path(self.single_end_reverse)
        if self.long_reads:
            self.long_reads = expand_path(self.long_reads)
        if self.cell_list:
            self.cell_list = expand_path(self.cell_list)
        if self.reference_transcripts:
            self.reference_transcripts = expand_path(self.reference_transcripts)

        # 只有在rnabloom_path不是默认值时才展开|Only expand rnabloom_path if not default
        if self.rnabloom_path and self.rnabloom_path != "rnabloom":
            self.rnabloom_path = expand_path(self.rnabloom_path)

        # 如果没有指定任何reads，报错|If no reads specified, raise error
        if not any([
            self.left_reads, self.right_reads,
            self.single_end_forward, self.single_end_reverse,
            self.long_reads, self.cell_list
        ]):
            raise ValueError(
                "必须指定至少一种输入reads类型|Must specify at least one input read type: "
                "--left/--right (paired-end), --sef/--ser (single-end), --long (long reads), or --cell-list (single-cell)"
            )

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件是否存在|Check if input files exist
        def check_file(file_path, desc):
            if file_path and not os.path.exists(file_path):
                errors.append(f"{desc}不存在|{desc} does not exist: {file_path}")

        check_file(self.left_reads, "左端reads文件|Left reads file")
        check_file(self.right_reads, "右端reads文件|Right reads file")
        check_file(self.single_end_forward, "单端正向reads文件|Single-end forward reads file")
        check_file(self.single_end_reverse, "单端反向reads文件|Single-end reverse reads file")
        check_file(self.long_reads, "长reads文件|Long reads file")
        check_file(self.cell_list, "单细胞列表文件|Single-cell list file")
        check_file(self.reference_transcripts, "参考转录本文件|Reference transcript file")

        # 检查配对的完整性|Check pairing completeness
        if (self.left_reads and not self.right_reads) or (self.right_reads and not self.left_reads):
            errors.append(
                "左端和右端reads必须同时指定|Both left and right reads must be specified together"
            )

        # 检查rnabloom工具（仅当用户提供自定义路径时）|Check rnabloom tool (only when user provides custom path)
        if self.rnabloom_path != "rnabloom" and not os.path.exists(self.rnabloom_path):
            errors.append(f"RNA-Bloom工具不存在|RNA-Bloom tool does not exist: {self.rnabloom_path}")

        # 检查参数合法性|Check parameter validity
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive integer: {self.threads}")

        if self.memory_gb is not None and self.memory_gb <= 0:
            errors.append(f"内存大小必须为正数|Memory size must be positive: {self.memory_gb}")

        if self.false_positive_rate is not None and (self.false_positive_rate <= 0 or self.false_positive_rate >= 1):
            errors.append(f"假阳性率必须在0-1之间|False positive rate must be between 0 and 1: {self.false_positive_rate}")

        if self.min_length < 0:
            errors.append(f"最小长度必须>=0|Minimum length must be >= 0: {self.min_length}")

        if self.stage is not None and self.stage not in [1, 2, 3]:
            errors.append(f"停止阶段必须是1-3|Stage must be 1-3: {self.stage}")

        # 检查互斥参数|Check mutually exclusive parameters
        if self.cell_list and any([self.left_reads, self.right_reads, self.long_reads]):
            errors.append(
                "单细胞模式(--cell-list)与其他输入参数互斥|"
                "Single-cell mode (--cell-list) is mutually exclusive with other input parameters"
            )

        if self.long_reads and self.reference_transcripts:
            errors.append(
                "参考引导组装(--ref)不支持长reads数据|"
                "Reference-guided assembly (--ref) is not supported for long-read data"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_assembly_mode(self) -> str:
        """获取组装模式|Get assembly mode

        Returns:
            str: 组装模式描述|Assembly mode description
        """
        if self.cell_list:
            return "单细胞混合组装|Single-cell pooled assembly"
        elif self.long_reads:
            if any([self.single_end_forward, self.single_end_reverse]):
                return "长reads + 短reads组装|Long-read + short-read assembly"
            else:
                return "长reads组装|Long-read assembly"
        elif any([self.left_reads, self.right_reads, self.single_end_forward, self.single_end_reverse]):
            return "短reads组装|Short-read assembly"
        else:
            return "未知模式|Unknown mode"

    def get_input_files(self) -> List[str]:
        """获取所有输入文件列表|Get list of all input files

        Returns:
            List[str]: 输入文件路径列表|List of input file paths
        """
        files = []
        if self.left_reads:
            files.append(self.left_reads)
        if self.right_reads:
            files.append(self.right_reads)
        if self.single_end_forward:
            files.append(self.single_end_forward)
        if self.single_end_reverse:
            files.append(self.single_end_reverse)
        if self.long_reads:
            files.append(self.long_reads)
        if self.cell_list:
            files.append(self.cell_list)
        if self.reference_transcripts:
            files.append(self.reference_transcripts)
        return files
