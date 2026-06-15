"""
PanMAN配置管理模块|PanMAN Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path


@dataclass
class PanMANConfig:
    """PanMAN配置类|PanMAN Configuration Class"""

    # 模式选择|Mode selection
    mode: str  # "build", "extract", or "generate_pangraph"

    # PanGraph生成模式参数|PanGraph generation mode parameters
    fasta_file: Optional[str] = None      # 输入FASTA文件 (用于生成PanGraph)

    # 构建模式参数|Build mode parameters
    pangraph_file: Optional[str] = None   # PanGraph JSON
    gfa_file: Optional[str] = None        # GFA文件
    msa_file: Optional[str] = None        # MSA文件 (FASTA格式)
    newick_file: Optional[str] = None     # Newick树文件

    # 提取模式参数|Extract mode parameters
    panman_file: Optional[str] = None     # 输入PanMAN文件

    # 通用参数|Common parameters
    output_prefix: str = "output"
    output_dir: str = "./panman_output"
    reference: Optional[str] = None       # 参考序列名称

    # 提取选项|Extraction options
    extract_summary: bool = False         # 提取摘要
    extract_fasta: bool = False           # 提取FASTA
    extract_msa: bool = False             # 提取MSA
    extract_vcf: bool = False             # 提取VCF
    extract_gfa: bool = False             # 提取GFA
    extract_newick: bool = False          # 提取Newick
    extract_extended_newick: bool = False # 提取扩展Newick
    extract_maf: bool = False             # 提取MAF
    extract_aa: bool = False              # 提取氨基酸翻译
    extract_subnet: bool = False          # 提取子网络
    extract_annotate: bool = False        # 注释节点
    extract_reroot: bool = False          # 重新扎根
    extract_create_network: bool = False  # 创建网络
    extract_print_mutations: bool = False # 打印突变

    # 高级选项|Advanced options
    conda_env: str = "panman_v.0.1.4"     # Conda环境名
    conda_base: Optional[str] = None      # Conda基础路径 (默认从环境变量CONDA_BASE读取)
    backend: str = "conda"                # 后端选择 (conda/docker/singularity) - 改为conda默认
    threads: int = 12
    pangraph_path: Optional[str] = None   # PanGraph可执行文件路径
    pangraph_sif: Optional[str] = None    # PanGraph Singularity SIF镜像路径
    sif_image: Optional[str] = None       # PanMAN Singularity SIF镜像路径
    singularity_path: Optional[str] = None # Singularity可执行文件路径

    # 新增高级参数|New advanced parameters
    input_file: Optional[str] = None      # 输入文件路径 (用于subnet/annotate/create-network)
    acr_method: str = "fitch"             # ACR方法 (fitch/mppa)
    tree_id: Optional[str] = None         # 树ID (用于VCF提取)
    range_query_index: Optional[str] = None  # 范围查询index参数
    range_start: Optional[int] = None     # 范围查询起始坐标
    range_end: Optional[int] = None       # 范围查询结束坐标

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        if self.fasta_file:
            self.fasta_file = os.path.normpath(os.path.abspath(self.fasta_file))
        if self.pangraph_file:
            self.pangraph_file = os.path.normpath(os.path.abspath(self.pangraph_file))
        if self.gfa_file:
            self.gfa_file = os.path.normpath(os.path.abspath(self.gfa_file))
        if self.msa_file:
            self.msa_file = os.path.normpath(os.path.abspath(self.msa_file))
        if self.newick_file:
            self.newick_file = os.path.normpath(os.path.abspath(self.newick_file))
        if self.panman_file:
            self.panman_file = os.path.normpath(os.path.abspath(self.panman_file))
        if self.reference:
            self.reference = self.reference.strip()
        if self.pangraph_path:
            self.pangraph_path = os.path.normpath(os.path.abspath(self.pangraph_path))
        if self.pangraph_sif:
            self.pangraph_sif = os.path.normpath(os.path.abspath(self.pangraph_sif))
        if self.sif_image:
            self.sif_image = os.path.normpath(os.path.abspath(self.sif_image))
        if self.singularity_path:
            self.singularity_path = os.path.normpath(os.path.abspath(self.singularity_path))

        # 设置默认 PanGraph SIF 镜像路径|Set default PanGraph SIF image path
        if not self.pangraph_sif:
            default_pangraph_sif = expand_path("~/software/singularity/pangraph_0.7.3.sif")
            if os.path.exists(default_pangraph_sif):
                self.pangraph_sif = default_pangraph_sif

        # 设置默认 PanMAN SIF 镜像路径|Set default PanMAN SIF image path
        if self.backend == "singularity" and not self.sif_image:
            default_sif = expand_path("~/software/singularity/panman_latest.sif")
            if os.path.exists(default_sif):
                self.sif_image = default_sif

        # 设置默认 singularity 路径|Set default singularity path
        if self.backend == "singularity" and not self.singularity_path:
            default_singularity = expand_path("~/miniforge3/envs/singularity_v.3.8.7/bin/singularity")
            if os.path.exists(default_singularity):
                self.singularity_path = default_singularity

        # 输出目录标准化|Normalize output directory
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 验证模式|Validate mode
        if self.mode not in ["build", "extract", "generate_pangraph"]:
            errors.append(
                f"无效的模式|Invalid mode: {self.mode} "
                "(必须为 'build', 'extract' 或 'generate_pangraph'|"
                "must be 'build', 'extract' or 'generate_pangraph')"
            )

        # PanGraph生成模式验证|PanGraph generation mode validation
        if self.mode == "generate_pangraph":
            if not self.fasta_file:
                errors.append(
                    "PanGraph生成模式需要FASTA文件|"
                    "PanGraph generation mode requires FASTA file"
                )
            elif not os.path.exists(self.fasta_file):
                errors.append(f"FASTA文件不存在|FASTA file does not exist: {self.fasta_file}")

        # 构建模式验证|Build mode validation
        if self.mode == "build":
            # 检查输入文件|Check input files
            input_files = [
                (self.pangraph_file, "PanGraph文件|PanGraph file"),
                (self.gfa_file, "GFA文件|GFA file"),
                (self.msa_file, "MSA文件|MSA file")
            ]

            has_input = False
            for file_path, file_desc in input_files:
                if file_path:
                    has_input = True
                    if not os.path.exists(file_path):
                        errors.append(f"{file_desc}不存在|does not exist: {file_path}")

            if not has_input:
                errors.append(
                    "构建模式需要至少一个输入文件 (Pangraph/GFA/MSA)|"
                    "Build mode requires at least one input file (Pangraph/GFA/MSA)"
                )

            # 检查Newick文件|Check Newick file
            if not self.newick_file:
                errors.append(
                    "构建模式需要Newick树文件|"
                    "Build mode requires Newick tree file"
                )
            elif not os.path.exists(self.newick_file):
                errors.append(f"Newick文件不存在|Newick file does not exist: {self.newick_file}")

        # 提取模式验证|Extract mode validation
        elif self.mode == "extract":
            if not self.panman_file:
                errors.append(
                    "提取模式需要PanMAN文件|"
                    "Extract mode requires PanMAN file"
                )
            elif not os.path.exists(self.panman_file):
                errors.append(f"PanMAN文件不存在|PanMAN file does not exist: {self.panman_file}")

            # 检查是否至少选择了一个提取选项|Check if at least one extraction option is selected
            extract_options = [
                self.extract_summary, self.extract_fasta, self.extract_msa,
                self.extract_vcf, self.extract_gfa, self.extract_newick,
                self.extract_extended_newick, self.extract_maf, self.extract_aa,
                self.extract_subnet, self.extract_annotate, self.extract_reroot,
                self.extract_create_network, self.extract_print_mutations
            ]
            if not any(extract_options):
                errors.append(
                    "提取模式需要至少选择一个提取选项|"
                    "Extract mode requires at least one extraction option"
                )

            # VCF提取需要参考序列|VCF extraction requires reference
            if self.extract_vcf and not self.reference:
                errors.append(
                    "VCF提取需要指定参考序列名称|"
                    "VCF extraction requires reference sequence name"
                )

            # 子网络提取需要输入文件|Subnet extraction requires input file
            if self.extract_subnet and not self.input_file:
                errors.append(
                    "子网络提取需要指定输入文件|"
                    "Subnet extraction requires input file with node list"
                )

            # 注释需要输入文件|Annotation requires input file
            if self.extract_annotate and not self.input_file:
                errors.append(
                    "节点注释需要指定输入文件|"
                    "Node annotation requires input file with annotations"
                )

            # 重新扎根需要参考序列|Reroot requires reference
            if self.extract_reroot and not self.reference:
                errors.append(
                    "重新扎根需要指定参考序列名称|"
                    "Reroot requires reference sequence name"
                )

            # 创建网络需要输入文件|Create network requires input file
            if self.extract_create_network and not self.input_file:
                errors.append(
                    "创建网络需要指定输入文件|"
                    "Create network requires input file with PanMAN list"
                )

            # 范围查询需要所有参数|Range query requires all parameters
            if self.range_query_index and (self.range_start is None or self.range_end is None):
                errors.append(
                    "范围查询需要指定起始和结束坐标|"
                    "Range query requires both start and end coordinates"
                )

            # ACR方法验证|ACR method validation
            if self.acr_method not in ["fitch", "mppa"]:
                errors.append(
                    f"不支持的ACR方法|Unsupported ACR method: {self.acr_method} "
                    "(仅支持 'fitch' 或 'mppa'|only 'fitch' or 'mppa' supported)"
                )

        # 验证线程数|Validate thread count
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        # 验证后端|Validate backend
        if self.backend not in ["conda", "docker", "singularity"]:
            errors.append(
                f"不支持的后端|Unsupported backend: {self.backend} "
                "(仅支持 'conda', 'docker' 或 'singularity'|"
                "only 'conda', 'docker' or 'singularity' supported)"
            )

        # Singularity 后端验证|Singularity backend validation
        if self.backend == "singularity":
            if not self.sif_image:
                errors.append(
                    "Singularity后端需要指定SIF镜像路径|"
                    "Singularity backend requires SIF image path"
                )
            elif not os.path.exists(self.sif_image):
                errors.append(f"SIF镜像不存在|SIF image does not exist: {self.sif_image}")

            if not self.singularity_path:
                errors.append(
                    "Singularity后端需要指定singularity可执行文件路径|"
                    "Singularity backend requires singularity executable path"
                )
            elif not os.path.exists(self.singularity_path):
                errors.append(f"Singularity可执行文件不存在|Singularity executable does not exist: {self.singularity_path}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_input_type(self):
        """获取输入类型|Get input type"""
        if self.pangraph_file:
            return "pangraph"
        elif self.gfa_file:
            return "gfa"
        elif self.msa_file:
            return "msa"
        else:
            return None

    def get_input_file(self):
        """获取输入文件路径|Get input file path"""
        if self.pangraph_file:
            return self.pangraph_file
        elif self.gfa_file:
            return self.gfa_file
        elif self.msa_file:
            return self.msa_file
        else:
            return None
