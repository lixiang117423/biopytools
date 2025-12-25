"""
HapHiC配置管理模块 | HapHiC Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, List


@dataclass
class HapHiCConfig:
    """HapHiC配置类 | HapHiC Configuration Class"""

    # 必需输入 | Required inputs
    asm_file: str
    hic_file: str
    nchrs: int
    hic_file_type: str = "bam"  # bam 或 fastq

    # 工具路径 | Tool paths
    haphic_bin: str = "haphic"
    bwa_bin: str = "bwa"
    samtools_bin: str = "samtools"

    # BWA参数 | BWA parameters
    read1_pattern: str = "_R1.fastq.gz"
    read2_pattern: str = "_R2.fastq.gz"

    # 输出配置 | Output configuration
    output_dir: Optional[str] = None
    prefix: Optional[str] = None

    # Hi-C数据处理参数 | Hi-C data processing parameters
    mapq_threshold: int = 1
    edit_distance: int = 3
    min_re_sites: int = 25

    # 内部变量 | Internal variables
    bam_file: Optional[str] = None  # 用于兼容性，由hic_file生成
    hic2_file: Optional[str] = None  # 第二个Hi-C FASTQ文件

    # 聚类参数 | Clustering parameters
    min_inflation: float = 1.0
    max_inflation: float = 3.0
    inflation_step: float = 0.2
    nx: int = 80
    min_group_len: int = 0

    # 排序和定向参数 | Ordering and orientation parameters
    processes: int = 8
    fast_sorting: bool = True
    allhic_optimization: bool = False

    # 组装校正参数 | Assembly correction parameters
    correct_nrounds: int = 2
    correct_min_coverage: float = 10.0
    correct_resolution: Optional[int] = None  # 自定义分辨率，None表示自动优化

    # 单倍型分相参数 | Haplotype phasing parameters
    remove_allelic_links: Optional[int] = None
    phasing_weight: float = 1.0
    gfa_files: Optional[str] = None

    # 可视化参数 | Visualization parameters
    generate_plots: bool = False
    bin_size: int = 500
    min_len: float = 1.0
    separate_plots: bool = False

    # 性能参数 | Performance parameters
    threads: int = 8
    memory_limit: Optional[str] = None

    # 断点续传参数 | Resume parameters
    force_rerun: bool = False  # 强制重新运行所有步骤

    # 高级选项 | Advanced options
    quick_view: bool = False
    re_sites: str = "GATC"
    skip_clustering: bool = False

    # 输出格式选项 | Output format options
    output_agp: bool = True
    output_fasta: bool = True
    output_juicebox: bool = True

    # Juicebox配置 | Juicebox configuration
    generate_juicebox: bool = True
    matlock_bin: str = "matlock"
    three_d_dna_dir: str = "/share/org/YZWL/yzwl_lixg/software/3d-dna"
    agp2assembly_script: str = "/share/org/YZWL/yzwl_lixg/software/3d-dna/utils/agp2assembly.py"
    asm_visualizer_script: str = "/share/org/YZWL/yzwl_lixg/software/3d-dna/visualize/run-assembly-visualizer.sh"

    # BWA比对配置 | BWA alignment configuration
    samblaster_bin: str = "samblaster"
    haphic_filter_bam_bin: str = "/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/filter_bam"
    use_samblaster: bool = True
    use_haphic_filter: bool = True

    # Pipeline模式新增参数 | Pipeline mode additional parameters
    density_lower: str = "0.2X"
    density_upper: str = "1.9X"
    read_depth_upper: str = "1.5X"
    topN: int = 10
    rank_sum_hard_cutoff: int = 0
    rank_sum_upper: str = "1.5X"
    concordance_ratio_cutoff: float = 0.2
    nwindows: int = 50
    max_read_pairs: int = 200
    min_read_pairs: int = 20
    remove_concentrated_links: bool = False
    bin_size_kbp: int = -1
    flank: int = 500
    expansion: int = 2
    max_iter: int = 200
    pruning: float = 0.0001
    max_ctg_len: float = 10000  # Kbp
    min_density_ratio: float = 4.0
    ambiguous_cutoff: float = 0.6
    reassign_nrounds: int = 5
    no_additional_rescue: bool = False
    skip_fast_sort: bool = False
    flanking_region: int = 0
    density_cal_method: str = "multiplication"
    confidence_cutoff: float = 1.0
    skip_allhic: bool = False
    mutprob: float = 0.2
    ngen: int = 5000
    npop: int = 100
    seed: int = 42
    Ns: int = 100
    max_width: int = 60
    sort_by_input: bool = False
    keep_letter_case: bool = False

    # 日志配置 | Logging configuration
    verbose: bool = False
    log_file: Optional[str] = None

    # 测试模式配置 | Test mode configuration
    dry_run: bool = False

    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 标准化路径 | Normalize paths
        self.asm_file = os.path.normpath(os.path.abspath(self.asm_file))

        # 只有在bam_file不为None时才处理 | Process bam_file only if not None
        if self.bam_file is not None:
            self.bam_file = os.path.normpath(os.path.abspath(self.bam_file))

        # 只有在haphic_bin不是默认值且不是None时才处理 | Process haphic_bin only if not default and not None
        if self.haphic_bin and self.haphic_bin != "haphic":
            self.haphic_bin = os.path.normpath(os.path.abspath(self.haphic_bin))

        # 设置输出目录 | Set output directory
        if not self.output_dir:
            self.output_dir = os.path.dirname(os.path.abspath(self.asm_file))
        else:
            self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 设置前缀 | Set prefix
        if not self.prefix:
            self.prefix = os.path.splitext(os.path.basename(self.asm_file))[0]

        # 创建输出目录 | Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # 设置默认日志文件 | Set default log file
        if not self.log_file:
            self.log_file = os.path.join(self.output_dir, f"{self.prefix}_haphic.log")

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []

        # 检查必需文件 | Check required files
        if not os.path.exists(self.asm_file):
            errors.append(f"基因组组装文件不存在 | Assembly file does not exist: {self.asm_file}")

        if not os.path.exists(self.hic_file):
            errors.append(f"Hi-C文件不存在 | Hi-C file does not exist: {self.hic_file}")

        # 检查文件格式 | Check file formats
        if not self.asm_file.endswith(('.fa', '.fasta', '.FA', '.FASTA')):
            errors.append("基因组文件必须是FASTA格式 | Assembly file must be in FASTA format")

        if self.hic_file_type == "bam":
            if not self.hic_file.endswith(('.bam', '.BAM')):
                errors.append("Hi-C文件必须是BAM格式 | Hi-C file must be in BAM format when --hic-file-type is bam")
        elif self.hic_file_type == "fastq":
            if os.path.isfile(self.hic_file) and not self.hic_file.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')):
                errors.append("Hi-C文件必须是FASTQ格式 | Hi-C file must be in FASTQ format when --hic-file-type is fastq")

        # 检查HapHiC工具 | Check HapHiC tool
        if self.haphic_bin != "haphic" and not os.path.exists(self.haphic_bin):
            errors.append(f"HapHiC可执行文件不存在 | HapHiC executable does not exist: {self.haphic_bin}")

        # 检查BWA工具 (如果需要) | Check BWA tool (if needed)
        if self.hic_file_type == "fastq":
            if self.bwa_bin != "bwa" and not os.path.exists(self.bwa_bin):
                errors.append(f"BWA可执行文件不存在 | BWA executable does not exist: {self.bwa_bin}")
            if self.samtools_bin != "samtools" and not os.path.exists(self.samtools_bin):
                errors.append(f"Samtools可执行文件不存在 | Samtools executable does not exist: {self.samtools_bin}")

        # 检查数值参数 | Check numeric parameters
        if self.nchrs <= 0:
            errors.append("染色体数必须大于0 | Number of chromosomes must be greater than 0")

        if self.mapq_threshold < 0:
            errors.append("MAPQ阈值不能为负数 | MAPQ threshold cannot be negative")

        if self.edit_distance < 0:
            errors.append("编辑距离不能为负数 | Edit distance cannot be negative")

        if self.min_inflation <= 0:
            errors.append("最小膨胀参数必须大于0 | Min inflation must be greater than 0")

        if self.max_inflation <= self.min_inflation:
            errors.append("最大膨胀参数必须大于最小膨胀参数 | Max inflation must be greater than min inflation")

        if self.threads < 1:
            errors.append("线程数必须大于等于1 | Number of threads must be >= 1")

        if self.processes < 1:
            errors.append("进程数必须大于等于1 | Number of processes must be >= 1")

        # 检查高级选项 | Check advanced options
        if self.quick_view and self.nchrs < 1:
            # Quick view mode doesn't care about nchrs, but needs at least 1 for parameter validation
            pass

        if self.remove_allelic_links and self.remove_allelic_links <= 0:
            errors.append("等位基因链接移除参数必须大于0 | Remove allelic links parameter must be greater than 0")

        if self.phasing_weight < 0 or self.phasing_weight > 1:
            errors.append("分相权重必须在0-1之间 | Phasing weight must be between 0 and 1")

        # 检查RE位点参数 | Check restriction site parameters
        if not self.re_sites:
            errors.append("限制性内切酶位点不能为空 | Restriction enzyme site cannot be empty")

        if self.bin_size <= 0:
            errors.append("装箱大小必须大于0 | Bin size must be greater than 0")

        if self.min_len <= 0:
            errors.append("最小长度必须大于0 | Min length must be greater than 0")

        # 检查输出目录可写性 | Check output directory writability
        try:
            test_file = os.path.join(self.output_dir, '.haphic_test')
            with open(test_file, 'w') as f:
                f.write('test')
            os.remove(test_file)
        except (PermissionError, OSError) as e:
            errors.append(f"输出目录不可写 | Output directory not writable: {e}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_restriction_sites(self):
        """获取限制性内切酶位点列表 | Get restriction enzyme sites list"""
        return [site.strip() for site in self.re_sites.split(',') if site.strip()]

    def get_output_files(self):
        """获取输出文件路径 | Get output file paths"""
        base_path = os.path.join(self.output_dir, self.prefix)

        return {
            "corrected_asm": f"{base_path}_corrected.asm.fa",
            "scaffolds_agp": f"{base_path}_scaffolds.agp",
            "scaffolds_raw_agp": f"{base_path}_scaffolds.raw.agp",
            "scaffolds_fasta": f"{base_path}_scaffolds.fa",
            "juicebox_script": f"{base_path}_juicebox.sh",
            "contact_map": f"{base_path}_contact_map.pdf",
            "hic_file": os.path.join(self.output_dir, "06.juicebox", f"{self.prefix}.hic"),
            "assembly_file": os.path.join(self.output_dir, "06.juicebox", f"{self.prefix}.assembly"),
            "log_file": self.log_file
        }

    def get_command_prefix(self):
        """获取命令前缀 | Get command prefix"""
        return f"{self.haphic_bin}"

    def get_common_options(self):
        """获取通用选项 | Get common options"""
        options = []

        if self.verbose:
            options.append("--verbose")

        if self.threads != 8:
            options.extend(["--threads", str(self.threads)])

        return options

    def get_pipeline_options(self):
        """获取pipeline选项 | Get pipeline options"""
        options = self.get_common_options()

        # 限制性酶位点
        if self.re_sites != "GATC":
            options.extend(["--RE", self.re_sites])

        # 组装校正
        if self.correct_nrounds > 0:
            options.extend(["--correct_nrounds", str(self.correct_nrounds)])

        # 等位基因处理
        if self.remove_allelic_links:
            options.extend(["--remove_allelic_links", str(self.remove_allelic_links)])

        # 分相权重
        if self.phasing_weight != 1.0:
            options.extend(["--phasing_weight", str(self.phasing_weight)])

        # GFA文件
        if self.gfa_files:
            options.extend(["--gfa", self.gfa_files])

        # 快速查看模式
        if self.quick_view:
            options.append("--quick_view")

        # 聚类参数
        if self.min_inflation != 1.0:
            options.extend(["--min_inflation", str(self.min_inflation)])

        if self.max_inflation != 3.0:
            options.extend(["--max_inflation", str(self.max_inflation)])

        if self.inflation_step != 0.2:
            options.extend(["--inflation_step", str(self.inflation_step)])

        if self.nx != 80:
            options.extend(["--Nx", str(self.nx)])

        if self.min_group_len != 0:
            options.extend(["--min_group_len", str(self.min_group_len)])

        # 其他选项
        if self.processes != 8:
            options.extend(["--processes", str(self.processes)])

        return options

    def get_summary(self):
        """获取配置摘要 | Get configuration summary"""
        return {
            "assembly_file": self.asm_file,
            "hic_file": self.hic_file,
            "hic_file_type": self.hic_file_type,
            "nchrs": self.nchrs,
            "output_dir": self.output_dir,
            "prefix": self.prefix,
            "threads": self.threads,
            "processes": self.processes,
            "quick_view": self.quick_view,
            "restriction_sites": self.get_restriction_sites(),
            "correct_nrounds": self.correct_nrounds,
            "remove_allelic_links": self.remove_allelic_links,
            "phasing_weight": self.phasing_weight,
            "gfa_files": self.gfa_files,
            "haphic_bin": self.haphic_bin,
            "bwa_bin": self.bwa_bin if self.hic_file_type == "fastq" else "N/A",
            "samtools_bin": self.samtools_bin if self.hic_file_type == "fastq" else "N/A"
        }