"""
基因组组装配置管理模块 |Genome Assembly Configuration Management Module
"""

import os
import glob
from dataclasses import dataclass
from pathlib import Path

@dataclass
class AssemblyConfig:
    """组装配置类|Assembly Configuration Class"""

    # 必需参数|Required parameters
    hifi_data: str  # HiFi数据文件|HiFi data file

    # 可选参数|Optional parameters
    hic_r1: str = None     # Hi-C R1文件|Hi-C R1 file
    hic_r2: str = None     # Hi-C R2文件|Hi-C R2 file

    # 基本参数|Basic parameters
    prefix: str = "genome_sample"  # 样本前缀|Sample prefix
    threads: int = 88             # 线程数|Number of threads
    genome_size: str = "1.45g"    # 基因组大小|Genome size estimate
    n_hap: int = 2                # 倍性|Ploidy

    # Hifiasm 参数|Hifiasm parameters
    purge_level: int = None       # purge level (-l): 0=no purging, 1=light, 2/3=aggressive [default: 3 for unzip]
    hom_cov: int = None           # homozygous read coverage (--hom-cov) [default: auto]

    # 路径参数|Path parameters
    base_dir: str = "./assembly_output"  # 基础输出目录|Base output directory

    # 执行控制|Execution control
    resume: bool = True  # 断点续传（默认启用）|Resume from previous run (enabled by default)

    # NGS polish参数|NGS polish parameters
    ngs_data: str = None      # NGS二代数据目录|NGS second-generation data directory
    high_cov: float = 95.0    # 高质量contig覆盖度阈值|High quality contig coverage threshold
    medium_cov_min: float = 30.0  # 中等质量contig最小覆盖度|Medium quality contig minimum coverage
    ngs_pattern: str = "_1.clean.fq.gz"  # NGS文件匹配模式|NGS file matching pattern

    # Purge_Dups去冗余参数|Purge_Dups deduplication parameters
    enable_purge_dups: bool = True  # 是否启用去冗余（默认启用）|Whether to enable deduplication (enabled by default)
    purge_dups_path: str = '~/miniforge3/envs/purge_dups_v.1.2.6'  # Purge_Dups软件路径|Purge_Dups software path
    purge_dups_threads: int = None  # 去冗余线程数（默认使用assembly的threads）|Deduplication threads (default: use assembly threads)
    purge_dups_read_type: str = 'hifi'  # 去冗余reads类型|Deduplication reads type

    # 内部属性|Internal attributes
    work_dir: str = None
    raw_dir: str = None
    fasta_dir: str = None
    log_dir: str = None
    stat_dir: str = None
    has_hic: bool = False  # 是否有Hi-C数据|Whether has Hi-C data
    has_ngs: bool = False  # 是否有NGS数据|Whether has NGS data
    ngs_polish_dir: str = None  # NGS polish输出目录|NGS polish output directory
    purge_dups_dir: str = None  # Purge_Dups去冗余输出目录|Purge_Dups deduplication output directory
    has_purge_dups: bool = False  # 是否启用去冗余|Whether to enable deduplication

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.base_dir = os.path.normpath(os.path.abspath(self.base_dir))
        self.work_dir = os.path.join(self.base_dir, self.prefix)

        # 将输入文件路径转换为绝对路径|Convert input file paths to absolute paths
        self.hifi_data = os.path.normpath(os.path.abspath(self.hifi_data))
        if self.hic_r1:
            self.hic_r1 = os.path.normpath(os.path.abspath(self.hic_r1))
        if self.hic_r2:
            self.hic_r2 = os.path.normpath(os.path.abspath(self.hic_r2))

        # 检查是否有Hi-C数据|Check if has Hi-C data
        self.has_hic = (self.hic_r1 is not None and self.hic_r2 is not None)

        # 检查是否有NGS数据|Check if has NGS data
        self.has_ngs = (self.ngs_data is not None)
        if self.has_ngs:
            self.ngs_data = os.path.normpath(os.path.abspath(self.ngs_data))

        # 检查是否启用Purge_Dups|Check if Purge_Dups is enabled
        self.has_purge_dups = self.enable_purge_dups
        if self.purge_dups_threads is None:
            self.purge_dups_threads = self.threads  # 默认使用assembly的线程数|Use assembly threads by default

        # 定义子目录|Define subdirectories
        # 目录编号始终固定，避免重复|Fixed directory numbering to avoid conflicts
        self.raw_dir = os.path.join(self.work_dir, "01.raw_output")
        self.fasta_dir = os.path.join(self.work_dir, "02.fasta")
        self.ngs_polish_dir = os.path.join(self.work_dir, "03.ngs_polish")  # 仅在有NGS时创建|Only created when has NGS
        self.purge_dups_dir = os.path.join(self.work_dir, "04.purge_dups")  # 去冗余目录|Deduplication directory
        self.stat_dir = os.path.join(self.work_dir, "05.statistics")  # 统计目录（编号后移）|Statistics directory (number shifted)
        self.log_dir = os.path.join(self.work_dir, "06.logs")  # 日志目录（编号后移）|Logs directory (number shifted)

        # 创建目录结构|Create directory structure
        dirs_to_create = [self.work_dir, self.raw_dir, self.fasta_dir, self.stat_dir, self.log_dir]
        if self.has_ngs:
            dirs_to_create.append(self.ngs_polish_dir)
        if self.has_purge_dups:
            dirs_to_create.append(self.purge_dups_dir)

        for dir_path in dirs_to_create:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
    
    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查输入文件|Check input files
        if not os.path.exists(self.hifi_data):
            errors.append(f" HiFi文件不存在|HiFi file not found: {self.hifi_data}")

        # 仅当有Hi-C数据时检查Hi-C文件|Check Hi-C files only when Hi-C data is provided
        if self.has_hic:
            if not os.path.exists(self.hic_r1):
                errors.append(f" Hi-C R1文件不存在|Hi-C R1 file not found: {self.hic_r1}")

            if not os.path.exists(self.hic_r2):
                errors.append(f" Hi-C R2文件不存在|Hi-C R2 file not found: {self.hic_r2}")

        # 检查NGS数据目录|Check NGS data directory
        if self.has_ngs:
            if not os.path.exists(self.ngs_data):
                errors.append(f" NGS数据目录不存在|NGS data directory not found: {self.ngs_data}")
            if not os.path.isdir(self.ngs_data):
                errors.append(f" NGS数据路径不是目录|NGS data path is not a directory: {self.ngs_data}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f" 线程数必须为正整数|Threads must be positive: {self.threads}")

        if self.n_hap <= 0:
            errors.append(f" 倍性必须为正整数|Ploidy must be positive: {self.n_hap}")

        if not self.genome_size.lower().endswith(('g', 'm', 'k')):
            errors.append(f" 基因组大小格式错误|Genome size format error: {self.genome_size}")

        if self.has_ngs and (self.high_cov <= 0 or self.high_cov > 100):
            errors.append(f" 高质量覆盖度阈值必须在0-100之间|High coverage threshold must be between 0-100: {self.high_cov}")

        if self.has_ngs and (self.medium_cov_min < 0 or self.medium_cov_min >= self.high_cov):
            errors.append(f" 中等质量覆盖度最小值必须在0-{self.high_cov}之间|Medium coverage minimum must be between 0-{self.high_cov}: {self.medium_cov_min}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_completed_steps(self) -> dict:
        """
        检查各步骤完成状态|Check completion status of each step

        Returns:
            dict: 各步骤完成状态|Completion status of each step
        """
        steps = {
            'assembly': False,      # Hifiasm组装完成|Hifiasm assembly completed
            'fasta_converted': False,  # GFA转FASTA完成|GFA to FASTA converted
            'reads_mapped': False,  # Contig-reads映射已生成|Contig-reads mapping generated
        }

        # 如果有NGS数据，添加NGS polish步骤检查|Add NGS polish steps if NGS data provided
        if self.has_ngs:
            steps.update({
                'bwa_alignment': False,      # BWA比对完成|BWA alignment completed
                'coverage_filter': False,    # 覆盖度过滤完成|Coverage filter completed
                'filtered_reads': False,     # 筛选reads完成|Filtered reads extracted
                'reassembly': False,         # 重新组装完成|Reassembly completed
            })

        # 根据是否有Hi-C数据定义GFA文件|Define GFA files based on Hi-C data availability
        if self.has_hic:
            primary_gfa = f"{self.prefix}.hic.p_ctg.gfa"
            primary_fasta = f"{self.prefix}.primary.fa"
            mapping_file = f"{self.prefix}.p_ctg.contig_reads.tsv"
            # 尝试.bp.前缀的文件名|Try .bp. prefixed filename
            primary_gfa_bp = f"{self.prefix}.bp.hic.p_ctg.gfa"
        else:
            primary_gfa = f"{self.prefix}.p_ctg.gfa"
            primary_fasta = f"{self.prefix}.primary.fa"
            mapping_file = f"{self.prefix}.p_ctg.contig_reads.tsv"
            # 尝试.bp.前缀的文件名|Try .bp. prefixed filename
            primary_gfa_bp = f"{self.prefix}.bp.p_ctg.gfa"

        # 检查组装是否完成|Check if assembly is completed
        # 先检查标准文件名，如果不存在则尝试.bp.前缀的文件名
        # Check standard filename first, if not exist try .bp. prefixed filename
        gfa_path = os.path.join(self.raw_dir, primary_gfa)
        gfa_path_bp = os.path.join(self.raw_dir, primary_gfa_bp)

        if os.path.exists(gfa_path):
            file_size = os.path.getsize(gfa_path)
            if file_size > 0:
                steps['assembly'] = True
        elif os.path.exists(gfa_path_bp):
            file_size = os.path.getsize(gfa_path_bp)
            if file_size > 0:
                steps['assembly'] = True

        # 检查FASTA转换是否完成|Check if FASTA conversion is completed
        fasta_path = os.path.join(self.fasta_dir, primary_fasta)
        if os.path.exists(fasta_path):
            file_size = os.path.getsize(fasta_path)
            if file_size > 0:
                steps['fasta_converted'] = True

        # 检查contig-reads映射是否已生成|Check if contig-reads mapping is generated
        mapping_path = os.path.join(self.fasta_dir, mapping_file)
        if os.path.exists(mapping_path):
            file_size = os.path.getsize(mapping_path)
            if file_size > 0:
                steps['reads_mapped'] = True

        # 如果有NGS数据，检查NGS polish各步骤|Check NGS polish steps if NGS data provided
        if self.has_ngs:
            # 检查BWA比对是否完成|Check if BWA alignment completed
            bam_file = os.path.join(self.ngs_polish_dir, "01.bwa_alignment", "bam", f"{self.prefix}.bam")
            if os.path.exists(bam_file) and os.path.getsize(bam_file) > 0:
                steps['bwa_alignment'] = True

            # 检查coverage filter是否完成|Check if coverage filter completed
            high_quality_list = os.path.join(self.ngs_polish_dir, "02.coverage_filter",
                                            f"{self.prefix}_high_quality.list")
            if os.path.exists(high_quality_list) and os.path.getsize(high_quality_list) > 0:
                steps['coverage_filter'] = True

            # 检查筛选的reads是否已提取|Check if filtered reads extracted
            filtered_reads = os.path.join(self.ngs_polish_dir, "03.filtered_reads",
                                         f"{self.prefix}_high_quality_reads.fq.gz")
            if os.path.exists(filtered_reads) and os.path.getsize(filtered_reads) > 0:
                steps['filtered_reads'] = True

            # 检查重新组装是否完成|Check if reassembly completed
            polished_genome = os.path.join(self.ngs_polish_dir, f"{self.prefix}.polished.fa")
            if os.path.exists(polished_genome) and os.path.getsize(polished_genome) > 0:
                steps['reassembly'] = True

        # 如果启用了Purge_Dups，添加去冗余步骤检查|Add Purge_Dups steps if enabled
        if self.has_purge_dups:
            steps.update({
                'purge_dups_coverage': False,  # 覆盖度计算完成|Coverage calculation completed
                'purge_dups_cutoffs': False,   # 阈值计算完成|Cutoffs calculation completed
                'purge_dups_split_align': False,  # 分割和自比对完成|Split and alignment completed
                'purge_dups_purge': False,  # 去冗余完成|Deduplication completed
                'purge_dups_seqs': False,  # 序列提取完成|Sequence extraction completed
            })

            # 检查Purge_Dups输出文件|Check Purge_Dups output files
            from pathlib import Path
            purge_output = Path(self.purge_dups_dir)

            # 检查最终输出文件|Check final output files
            purged_fa = purge_output / "sequences" / f"{self.prefix}_purged.purge.fa"
            if purged_fa.exists() and purged_fa.stat().st_size > 0:
                steps['purge_dups_seqs'] = True
                steps['purge_dups_coverage'] = True
                steps['purge_dups_cutoffs'] = True
                steps['purge_dups_split_align'] = True
                steps['purge_dups_purge'] = True

        return steps
