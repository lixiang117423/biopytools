"""
Cactus配置类|Cactus Configuration Class
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List


def expand_path(path: str) -> str:
    """
    展开路径中的~和环境变量|Expand ~ and environment variables in path

    Args:
        path: 输入路径|Input path

    Returns:
        展开后的绝对路径|Expanded absolute path
    """
    return str(Path(path).expanduser().resolve())


@dataclass
class CactusConfig:
    """Cactus泛基因组分析配置类|Cactus Pangenome Analysis Configuration Class"""

    # 必需参数|Required parameters
    seqfile: str = ""  # 序列文件|Sequence file (seqFile.txt)
    output_dir: str = ""  # 输出目录|Output directory
    reference: str = ""  # 参考基因组名称|Reference genome name (must be first in seqfile)

    # Singularity相关|Singularity related
    singularity_path: str = "~/miniforge3/envs/singularity_v.3.8.7/bin/singularity"  # Singularity可执行文件路径
    cactus_sif: str = "~/software/singularity/cactus_v3.1.4.sif"  # Cactus SIF镜像路径

    # 流程参数|Pipeline parameters
    jobstore: str = "cactus-jobstore"  # Toil jobstore目录|Toil jobstore directory
    out_name: str = "cactus_output"  # 输出文件前缀|Output file prefix
    cleanup: bool = True  # 是否清理jobstore|Whether to cleanup jobstore

    # 输出格式|Output formats
    # 注意：HAL是默认输出（{out_name}.full.hal），不需要在列表中指定
    # Note: HAL is default output ({out_name}.full.hal), no need to specify in list
    output_formats: List[str] = field(default_factory=lambda: ["gfa", "gbz"])  # 输出格式列表

    # 性能参数|Performance parameters
    threads: int = 12  # CPU核心数|Number of CPU cores
    max_memory: str = "100G"  # 最大内存|Maximum memory
    batch_system: Optional[str] = None  # 批处理系统|Batch system (e.g., "slurm")

    # 绑定目录|Bind directories
    bind_paths: Optional[List[str]] = None  # 需要绑定的目录路径|Directories to bind
    work_dir: Optional[str] = None  # 工作目录（用于临时文件）|Work directory (for temporary files)

    # 日志|Logging
    log_level: str = "INFO"
    verbose: bool = True

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开所有路径|Expand all paths
        if self.seqfile:
            self.seqfile = expand_path(self.seqfile)
        if self.output_dir:
            self.output_dir = expand_path(self.output_dir)

        # Singularity路径展开|Expand Singularity paths
        self.singularity_path = expand_path(self.singularity_path)
        self.cactus_sif = expand_path(self.cactus_sif)

        # 规范化输出格式|Normalize output formats
        self.output_formats = [fmt.lower() for fmt in self.output_formats]
        # 注意：HAL是默认输出（{out_name}.full.hal），不需要指定|Note: HAL is default output, no need to specify
        valid_formats = {"gfa", "gbz", "odgi", "hal", "vg", "vcf", "xg"}
        for fmt in self.output_formats:
            if fmt not in valid_formats:
                raise ValueError(f"无效的输出格式|Invalid output format: {fmt}. 必须是|Must be one of {valid_formats}")

        self.log_level = self.log_level.upper()

        # 确保输出目录存在|Ensure output directory exists
        if self.output_dir:
            Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # 设置默认绑定路径|Set default bind paths
        if self.bind_paths is None:
            # 默认绑定输入文件目录、输出目录和临时目录|Bind input, output and temp directories by default
            self.bind_paths = []

            # 绑定/tmp目录（Toil需要）|Bind /tmp directory (required by Toil)
            self.bind_paths.append("/tmp")

        # 设置workDir|Set workDir
        if self.work_dir is None:
            # 默认使用输出目录下的work子目录|Default to work subdirectory in output dir
            self.work_dir = str(Path(self.output_dir) / "work")
        else:
            self.work_dir = expand_path(self.work_dir)

        # 确保workDir存在|Ensure workDir exists
        Path(self.work_dir).mkdir(parents=True, exist_ok=True)

        # 将workDir添加到绑定路径|Add workDir to bind paths
        if self.work_dir not in self.bind_paths:
            self.bind_paths.append(self.work_dir)

            # 绑定输入文件目录|Bind input file directory
            if self.seqfile:
                seqfile_dir = str(Path(self.seqfile).parent)
                if seqfile_dir not in self.bind_paths:
                    self.bind_paths.append(seqfile_dir)

            # 绑定输出目录|Bind output directory
            if self.output_dir:
                output_dir = self.output_dir
                if output_dir not in self.bind_paths:
                    self.bind_paths.append(output_dir)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        if not self.seqfile:
            raise ValueError("序列文件不能为空|Sequence file cannot be empty")

        if not Path(self.seqfile).exists():
            raise ValueError(f"序列文件不存在|Sequence file does not exist: {self.seqfile}")

        if not self.output_dir:
            raise ValueError("输出目录不能为空|Output directory cannot be empty")

        if not self.reference:
            raise ValueError("参考基因组名称不能为空|Reference genome name cannot be empty")

        # 检查Singularity是否存在|Check if Singularity exists
        if not Path(self.singularity_path).exists():
            raise ValueError(f"Singularity不存在|Singularity does not exist: {self.singularity_path}")

        # 检查Cactus SIF是否存在|Check if Cactus SIF exists
        if not Path(self.cactus_sif).exists():
            raise ValueError(f"Cactus SIF不存在|Cactus SIF does not exist: {self.cactus_sif}")

        # 验证序列文件格式|Validate sequence file format
        self._validate_seqfile()

        return True

    def _validate_seqfile(self):
        """
        验证序列文件格式|Validate sequence file format

        序列文件格式|Sequence file format:
            两列格式：样本名 + FASTA路径（制表符或空格分隔）|Two columns: sample_name + FASTA path (tab or space separated)
            第一个样本必须是参考基因组|First sample must be reference genome
            支持相对路径和绝对路径|Support relative and absolute paths
        """
        try:
            with open(self.seqfile, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]

            if not lines:
                raise ValueError(f"序列文件为空|Sequence file is empty: {self.seqfile}")

            # 解析第一行，提取样本名和路径|Parse first line, extract sample name and path
            first_line = lines[0]
            parts = first_line.split()  # 按空白字符分割|Split by whitespace

            if len(parts) < 2:
                raise ValueError(f"序列文件格式错误：第一行应包含两列（样本名 + 路径）|Invalid seqfile format: first line should have 2 columns (sample_name + path)")

            # 第二个部分是路径（样本名可能有空格，所以我们取最后一个作为路径）|Second part is path
            # 但实际上第一列是样本名，后续是路径|First column is sample name, rest is path
            sample_name = parts[0]
            genome_path = ' '.join(parts[1:]) if len(parts) > 2 else parts[1]

            # 检查第一个文件是否存在|Check if first file exists
            if not Path(genome_path).exists():
                # 尝试相对于seqfile的路径|Try path relative to seqfile
                seqfile_dir = Path(self.seqfile).parent
                relative_path = seqfile_dir / genome_path
                if not relative_path.exists():
                    raise ValueError(f"参考基因组文件不存在|Reference genome file does not exist: {genome_path}")

            # 检查至少有2个基因组|Check at least 2 genomes
            if len(lines) < 2:
                raise ValueError(f"序列文件至少需要包含2个基因组（参考+查询）|Sequence file must contain at least 2 genomes (reference + query)")

        except Exception as e:
            raise ValueError(f"序列文件验证失败|Sequence file validation failed: {e}")

    def get_jobstore_path(self) -> str:
        """获取jobstore完整路径|Get full jobstore path"""
        return str(Path(self.output_dir) / self.jobstore)

    def get_bind_args(self) -> List[str]:
        """
        获取Singularity绑定参数|Get Singularity bind arguments

        Returns:
            绑定参数列表|List of bind arguments
        """
        bind_args = []

        if self.bind_paths:
            for bind_path in self.bind_paths:
                # 展开路径|Expand path
                expanded_path = expand_path(bind_path)
                # Singularity绑定格式|Singularity bind format: --bind host_path:container_path
                # 如果容器内路径相同，可以简化为|If container path is same, can simplify to: --bind host_path
                bind_args.append(f"--bind {expanded_path}")

        return bind_args

    def __repr__(self):
        """配置的字符串表示|String representation of configuration"""
        return (
            f"CactusConfig(\n"
            f"  seqfile={self.seqfile!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  reference={self.reference!r},\n"
            f"  singularity_path={self.singularity_path!r},\n"
            f"  cactus_sif={self.cactus_sif!r},\n"
            f"  jobstore={self.jobstore!r},\n"
            f"  out_name={self.out_name!r},\n"
            f"  cleanup={self.cleanup!r},\n"
            f"  output_formats={self.output_formats!r},\n"
            f"  threads={self.threads!r},\n"
            f"  max_memory={self.max_memory!r},\n"
            f"  batch_system={self.batch_system!r},\n"
            f"  bind_paths={self.bind_paths!r},\n"
            f"  log_level={self.log_level!r}\n"
            f")"
        )
