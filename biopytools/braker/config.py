"""
BRAKER3基因组注释配置管理模块|BRAKER3 Genome Annotation Configuration Management Module
"""

import os
from ..common.paths import expand_path, resolve_legacy_path, get_samtools_path
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


def contains_non_ascii(text: str) -> bool:
    """
    检测字符串是否包含非ASCII字符（包括中文等）|Detect if string contains non-ASCII characters

    Args:
        text: 待检测字符串|String to check

    Returns:
        bool: 如果包含非ASCII字符返回True，否则返回False|True if contains non-ASCII, False otherwise
    """
    # 使用正则表达式检测非ASCII字符（Unicode字符范围>127）
    # Use regex to detect non-ASCII characters (Unicode range > 127)
    return bool(re.search(r'[^\x00-\x7F]', text))


def smart_normalize_path(path: str, base_dir: str = None) -> str:
    """
    智能路径标准化：如果路径包含非ASCII字符则保持相对路径，否则转换为绝对路径
    Smart path normalization: Keep relative path if contains non-ASCII, else convert to absolute

    Args:
        path: 输入路径|Input path
        base_dir: 基础目录（用于相对路径处理）|Base directory for relative path resolution

    Returns:
        str: 标准化后的路径|Normalized path
    """
    if not path:
        return path

    # 如果路径包含非ASCII字符，保持为相对路径（如果本来就是）
    # If path contains non-ASCII characters, keep as relative (if it was relative)
    is_relative = not os.path.isabs(path)

    if is_relative:
        # 检查完整路径是否包含非ASCII字符
        # Check if full path contains non-ASCII characters
        if base_dir:
            full_path = os.path.normpath(os.path.join(base_dir, path))
        else:
            full_path = os.path.normpath(os.path.abspath(path))

        if contains_non_ascii(full_path):
            # 包含非ASCII字符，保持相对路径
            # Contains non-ASCII, keep relative path
            return os.path.normpath(path)
        else:
            # 不包含非ASCII字符，可以安全转换为绝对路径
            # No non-ASCII, safe to convert to absolute path
            return full_path
    else:
        # 本来就是绝对路径，检查是否包含非ASCII
        # Already absolute, check for non-ASCII
        if contains_non_ascii(path):
            # 绝对路径包含非ASCII，尝试转换为相对路径
            # Absolute path contains non-ASCII, try to convert to relative
            if base_dir:
                try:
                    rel_path = os.path.relpath(path, base_dir)
                    # 检查相对路径是否不包含上层目录引用（更安全）
                    # Check if relative path doesn't contain parent references (safer)
                    if not rel_path.startswith('..'):
                        return rel_path
                except ValueError:
                    # 在不同驱动器上无法转换（Windows）
                    # Cannot convert on different drives (Windows)
                    pass
            # 无法转换，返回原路径（可能会失败，但已尽力）
            # Cannot convert, return original (may fail, but tried our best)
            return path
        else:
            # 绝对路径且不含非ASCII，标准化后返回
            # Absolute and no non-ASCII, normalize and return
            return os.path.normpath(path)


def get_safe_absolute_path(path: str, base_dir: str = None) -> str:
    """
    获取安全的绝对路径（如果包含非ASCII，通过符号链接规避）
    Get safe absolute path (use symlink to workaround if contains non-ASCII)

    Args:
        path: 输入路径|Input path
        base_dir: 基础目录|Base directory

    Returns:
        str: 安全的绝对路径（不含非ASCII）|Safe absolute path (no non-ASCII)
    """
    import hashlib

    # 标准化路径|Normalize path
    if base_dir and not os.path.isabs(path):
        full_path = os.path.normpath(os.path.join(base_dir, path))
    else:
        full_path = os.path.normpath(path)

    # 如果已经是绝对路径且不含非ASCII，直接返回
    # If already absolute and no non-ASCII, return directly
    if os.path.isabs(full_path) and not contains_non_ascii(full_path):
        return full_path

    # 转换为绝对路径|Convert to absolute path
    abs_path = os.path.abspath(full_path)

    # 如果绝对路径不含非ASCII，返回
    # If absolute path has no non-ASCII, return
    if not contains_non_ascii(abs_path):
        return abs_path

    # 包含非ASCII，需要创建符号链接到用户临时目录
    # Contains non-ASCII, need to create symlink to user temp directory
    path_hash = hashlib.md5(abs_path.encode('utf-8')).hexdigest()[:12]

    # 使用 ~/tmp 而不是系统 /tmp
    # Use ~/tmp instead of system /tmp
    user_home = os.path.expanduser("~")
    user_temp_base = os.path.join(user_home, "tmp")
    os.makedirs(user_temp_base, exist_ok=True)

    safe_link_name = f"braker_link_{path_hash}"
    safe_link_path = os.path.join(user_temp_base, safe_link_name)

    # 创建符号链接|Create symlink
    if os.path.islink(safe_link_path):
        # 符号链接已存在，检查目标|Symlink exists, check target
        existing_target = os.readlink(safe_link_path)
        if existing_target != abs_path:
            os.remove(safe_link_path)
            os.symlink(abs_path, safe_link_path)
    elif os.path.exists(safe_link_path):
        # 不是符号链接，删除|Not symlink, remove
        import shutil
        shutil.rmtree(safe_link_path)
        os.symlink(abs_path, safe_link_path)
    else:
        # 创建新符号链接|Create new symlink
        os.symlink(abs_path, safe_link_path)

    return safe_link_path


@dataclass
class BrakerConfig:
    """BRAKER3注释配置类|BRAKER3 Annotation Configuration Class"""

    # ===== 必需参数|Required parameters =====
    genome: str  # 基因组FASTA文件|Genome FASTA file
    species: str  # 物种名称|Species name

    # ===== 输入数据|Input data =====
    prot_seq: Optional[str] = None  # 近缘物种蛋白质序列或文件夹|Protein sequences file or directory
    isoseq: Optional[str] = None  # 三代全长转录本文件夹|Long-read transcript directory
    rnaseq_dirs: Optional[List[str]] = None  # 二代RNA-seq数据目录列表|List of short-read RNA-seq directories
    rnaseq_sets: Optional[List[str]] = None  # RNA-seq集合ID列表|List of RNA-seq set IDs (已废弃|Deprecated)

    # ===== 文件识别模式|File identification patterns =====
    read1_pattern: str = "_1.clean.fq.gz"  # R1文件模式|R1 file pattern (fastp处理后默认后缀)
    read2_pattern: str = "_2.clean.fq.gz"  # R2文件模式|R2 file pattern

    # ===== Singularity镜像配置|Singularity image configuration =====
    use_singularity: bool = True  # 使用Singularity镜像|Use Singularity image
    singularity_image: str = "~/software/singularity/braker3_devel.sif"

    # ===== 镜像内工具路径|Tool paths inside image =====
    braker_in_image: str = "/opt/BRAKER/scripts/braker.pl"

    # ===== 宿主机工具路径|Host machine tool paths =====
    singularity_bin: str = "~/miniforge3/envs/singularity_v.3.8.7/bin/singularity"
    repeatmodeler_bin: str = "~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/RepeatModeler"
    repeatmasker_bin: str = "~/miniforge3/envs/repeat_identiy/bin/RepeatMasker"
    build_database_bin: str = "~/miniforge3/envs/repeatmodeler_v.2.0.7/bin/BuildDatabase"
    minimap2_bin: str = "~/miniforge3/envs/Genome_dedup/bin/minimap2"
    hisat2_bin: str = "~/miniforge3/envs/RNA_Seq/bin/hisat2"
    hisat2_build_bin: str = "~/miniforge3/envs/RNA_Seq/bin/hisat2-build"

    # ===== repeat_refine 工具(repeat库过滤+证据还原)|repeat_refine tools =====
    hmmscan_bin: str = "~/miniforge3/envs/braker_v.3.0.8/bin/hmmscan"
    mmseqs_bin: str = "~/miniforge3/envs/eggnog-mapper_v.2.1.15/bin/mmseqs"
    miniprot_bin: str = "~/miniforge3/envs/braker_v.3.0.8/bin/miniprot"
    samtools_bin: str = ""  # 空则用 common.paths.get_samtools_path() 兜底|Empty => fallback

    # ===== 输出配置|Output configuration =====
    output_dir: str = "./braker_output"  # 输出目录|Output directory

    # ===== 流程参数|Pipeline parameters =====
    threads: int = 12  # 线程数|Number of threads
    use_fungus: bool = True  # 使用真菌模式|Use fungus mode (suitable for oomycetes)
    soft_masking: bool = True  # 软屏蔽重复序列|Soft masking of repeats

    # ===== 步骤控制|Step control =====
    skip_repeat: bool = False  # 跳过重复序列屏蔽|Skip repeat masking
    skip_long_reads: bool = False  # 跳过三代转录本处理|Skip long-read processing
    skip_short_reads: bool = False  # 跳过二代RNA-seq处理|Skip short-read processing
    skip_repeat_filter: bool = False  # 跳过repeat库过滤(方案1)|Skip repeat library filtering
    skip_rescue: bool = True  # 跳过证据还原(默认关闭,filter库级已处理假重复)|Skip rescue (default off)

    # ===== BRAKER3特定参数|BRAKER3 specific parameters =====
    busco_lineage: Optional[str] = None  # BUSCO谱系|BUSCO lineage
    utr: bool = False  # 预测UTR|Predict UTR
    training_genes: Optional[str] = None  # 训练基因集|Training gene set
    use_existing: bool = False  # 使用已有参数|Use existing parameters

    # ===== repeat_refine 参数(repeat库过滤+证据还原)|repeat_refine params =====
    pfam_db: str = "~/database/eggnog/pfam/Pfam-A.hmm"  # Pfam-A HMM 库|Pfam-A HMM DB
    te_domain_evalue: float = 1e-5  # TE domain hmmscan E-value 阈值|TE domain E-value cutoff
    filter_min_orf_len: int = 30  # 过滤用最小 ORF 长度(aa)|Min ORF length (aa) for filter
    prot_homology_evalue: float = 1e-5  # prot_seq 同源 E-value 阈值|Protein homology E-value
    prot_homology_pident: float = 50.0  # prot_seq 同源 identity 阈值(%)|Protein homology identity
    rescue_min_cds_len: int = 100  # rescue 蛋白证据最小覆盖长度(bp)|Min CDS overlap (bp)
    rescue_min_identity: float = 70  # rescue 蛋白最小 identity(%)|Min protein identity (%)
    rescue_min_depth: int = 5  # rescue RNA-seq 最小覆盖度|Min RNA-seq depth

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 记录当前工作目录，用于智能路径处理|Record current working directory for smart path handling
        self.working_dir = os.getcwd()

        # 展开所有含~和环境变量的路径|Expand all paths with ~ and env vars
        self.output_dir = expand_path(self.output_dir)
        self.singularity_image = expand_path(self.singularity_image)
        self.repeatmodeler_bin = expand_path(self.repeatmodeler_bin)
        self.repeatmasker_bin = expand_path(self.repeatmasker_bin)
        self.singularity_bin = expand_path(self.singularity_bin)
        self.build_database_bin = expand_path(self.build_database_bin)
        self.minimap2_bin = expand_path(self.minimap2_bin)
        self.hisat2_bin = expand_path(self.hisat2_bin)
        self.hisat2_build_bin = expand_path(self.hisat2_build_bin)
        self.hmmscan_bin = expand_path(self.hmmscan_bin)
        self.mmseqs_bin = expand_path(self.mmseqs_bin)
        self.miniprot_bin = expand_path(self.miniprot_bin)
        self.pfam_db = expand_path(self.pfam_db)
        # samtools: 用户未指定则用 common.paths 兜底|Fallback to common.paths if unset
        if self.samtools_bin:
            self.samtools_bin = expand_path(self.samtools_bin)
        else:
            self.samtools_bin = get_samtools_path()

        # 智能路径标准化：检测非ASCII字符，自动决定使用相对或绝对路径
        # Smart path normalization: Detect non-ASCII characters, auto decide relative or absolute
        self.genome = smart_normalize_path(self.genome, self.working_dir)
        self.output_dir = smart_normalize_path(self.output_dir, self.working_dir)
        self.singularity_image = smart_normalize_path(self.singularity_image, self.working_dir)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 智能标准化可选输入文件路径|Smart normalize optional input file paths
        if self.prot_seq:
            self.prot_seq = smart_normalize_path(self.prot_seq, self.working_dir)
        if self.isoseq:
            self.isoseq = smart_normalize_path(self.isoseq, self.working_dir)
        if self.rnaseq_dirs:
            self.rnaseq_dirs = [smart_normalize_path(d, self.working_dir) for d in self.rnaseq_dirs]
        if self.training_genes:
            self.training_genes = smart_normalize_path(self.training_genes, self.working_dir)

        # 定义子目录|Define subdirectories
        # 优先下划线规范名，回退点号老名用于断点续传|Prefer underscore, fall back to legacy dot name
        self.repeat_dir = resolve_legacy_path(self.output_dir, "01_repeat_masking")
        self.long_reads_dir = resolve_legacy_path(self.output_dir, "02_long_reads")
        self.short_reads_dir = resolve_legacy_path(self.output_dir, "03_short_reads")
        self.braker_dir = resolve_legacy_path(self.output_dir, "04_braker_annotation")
        self.log_dir = os.path.join(self.output_dir, "logs")

        # 创建子目录|Create subdirectories
        for dir_path in [self.repeat_dir, self.long_reads_dir,
                         self.short_reads_dir, self.braker_dir,
                         self.log_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

        # 为BRAKER3创建安全的工作目录路径（不含非ASCII）
        # Create safe working directory path for BRAKER3 (no non-ASCII)
        # 如果braker_dir包含中文，在~/tmp创建真实目录而不是符号链接
        # If braker_dir contains Chinese, create real directory in ~/tmp instead of symlink
        if contains_non_ascii(self.braker_dir if os.path.isabs(self.braker_dir) else os.path.join(self.working_dir, self.braker_dir)):
            # 使用真实目录而不是符号链接
            # Use real directory instead of symlink
            import hashlib
            import tempfile
            import shutil
            # 哈希包含输出目录和基因组路径，确保不同基因组使用不同工作目录
            # Hash includes output dir and genome path to ensure different genomes use different dirs
            path_hash = hashlib.md5(f"{self.braker_dir}_{self.genome}".encode('utf-8')).hexdigest()[:12]
            user_home = os.path.expanduser("~")
            user_tmp = os.path.join(user_home, "tmp")
            os.makedirs(user_tmp, exist_ok=True)

            # 创建真实的braker输出目录
            # Create real braker output directory
            safe_dir_name = f"braker_run_{path_hash}"
            self.braker_safe_dir = os.path.join(user_tmp, safe_dir_name)

            # 如果目录不存在则创建（不清理已有目录，避免丢失上次运行结果）
            # Create if not exists (don't clean existing to avoid losing previous results)
            os.makedirs(self.braker_safe_dir, exist_ok=True)

            # 标记需要将结果复制回原始目录
            # Mark that results need to be copied back to original directory
            self.braker_needs_copyback = True
        else:
            # 不含中文，直接使用原目录
            # No Chinese, use original directory directly
            self.braker_safe_dir = self.braker_dir
            self.braker_needs_copyback = False

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")

        # 检查可选输入文件|Check optional input files
        if self.prot_seq and not os.path.exists(self.prot_seq):
            errors.append(f"蛋白质文件不存在|Protein file not found: {self.prot_seq}")

        if self.isoseq and not os.path.exists(self.isoseq):
            errors.append(f"Iso-Seq文件不存在|Iso-Seq file not found: {self.isoseq}")

        if self.rnaseq_dirs:
            for i, rnaseq_dir in enumerate(self.rnaseq_dirs):
                if not os.path.exists(rnaseq_dir):
                    errors.append(f"RNA-seq目录不存在|RNA-seq directory not found [{i}]: {rnaseq_dir}")

        # 检查Singularity镜像|Check Singularity image
        if self.use_singularity and not os.path.exists(self.singularity_image):
            errors.append(f"Singularity镜像不存在|Singularity image not found: {self.singularity_image}")

        # 检查宿主机工具|Check host machine tools
        if not self.skip_repeat:
            if not os.path.exists(self.repeatmodeler_bin):
                errors.append(f"RepeatModeler不存在|RepeatModeler not found: {self.repeatmodeler_bin}")
            if not os.path.exists(self.repeatmasker_bin):
                errors.append(f"RepeatMasker不存在|RepeatMasker not found: {self.repeatmasker_bin}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        # 检查输入数据完整性|Check input data completeness
        has_protein = self.prot_seq is not None
        has_long_reads = self.isoseq is not None
        has_short_reads = self.rnaseq_dirs is not None and len(self.rnaseq_dirs) > 0

        if not has_protein and not has_long_reads and not has_short_reads:
            errors.append(
                "必须提供至少一种证据数据(蛋白质/三代转录本/二代RNA-seq)|"
                "Must provide at least one evidence type (proteins/long-reads/short-reads)"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
