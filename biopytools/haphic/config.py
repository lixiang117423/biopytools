"""
HapHiC配置管理模块|HapHiC Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List
from ..common.paths import expand_path, resolve_legacy_path, get_tool_path


@dataclass
class HapHiCConfig:
    """HapHiC配置类|HapHiC Configuration Class"""

    # 必需输入|Required inputs
    asm_file: str
    hic_file: str
    nchrs: int
    hic_file_type: str = "bam"  # bam 或 fastq

    # 工具路径|Tool paths
    # 注意：默认值使用conda环境路径，遵循开发规范第13章
    # Note: Default values use conda environment paths, following development guide Chapter 13
    haphic_bin: str = field(
        default_factory=lambda: get_tool_path(
            'haphic', '~/miniforge3/envs/haphic/bin/haphic', 'HAPHIC_PATH'
        )
    )
    bwa_bin: str = field(
        default_factory=lambda: get_tool_path(
            'bwa', '~/miniforge3/envs/Population_genetics/bin/bwa', 'BWA_PATH'
        )
    )
    samtools_bin: str = field(
        default_factory=lambda: get_tool_path(
            'samtools', '~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools', 'SAMTOOLS_PATH'
        )
    )

    # BWA参数|BWA parameters
    read1_pattern: str = "_R1.fastq.gz"
    read2_pattern: str = "_R2.fastq.gz"

    # 输出配置|Output configuration
    output_dir: Optional[str] = None
    prefix: Optional[str] = None

    # Hi-C数据处理参数|Hi-C data processing parameters
    mapq_threshold: int = 1
    edit_distance: int = 3
    re_site_cutoff: int = 5  # step1过滤RE位点阈值|Step1 RE site filtering threshold (upstream default: 5)
    min_re_sites: int = 25  # step2重分配RE位点阈值|Step2 reassignment RE site threshold (upstream default: 25)

    # 内部变量|Internal variables
    bam_file: Optional[str] = None  # 用于兼容性，由hic_file生成
    hic2_file: Optional[str] = None  # 第二个Hi-C FASTQ文件

    # 聚类参数|Clustering parameters
    min_inflation: float = 1.1  # upstream default: 1.1
    max_inflation: float = 3.0
    inflation_step: float = 0.1  # upstream default: 0.1
    nx: int = 80
    min_group_len: float = 5.0  # upstream default: 5 (Mbp)

    # 排序和定向参数|Ordering and orientation parameters
    processes: int = 8
    fast_sorting: bool = True
    allhic_optimization: bool = False

    # 组装校正参数|Assembly correction parameters
    correct_nrounds: int = 2  # 模块默认启用组装校正并重跑HapHiC|Module default: enable correction and re-run
    correct_min_coverage: float = 10.0  # 模块自用，自动优化分辨率|Module-specific, auto-optimize resolution
    correct_resolution: Optional[int] = None  # 自定义分辨率，None表示自动优化
    median_cov_ratio: float = 0.2  # 覆盖率截断乘数|Coverage cutoff multiplier
    region_len_ratio: float = 0.1  # 高覆盖区域长度比阈值|Large high-coverage region length ratio
    min_region_cutoff: int = 5000  # 高覆盖区域最小长度(bp)|Min length for high-coverage regions

    # 单倍型分相参数|Haplotype phasing parameters
    remove_allelic_links: Optional[int] = 0  # upstream default: 0 (disabled)
    phasing_weight: float = 1.0
    gfa_files: Optional[str] = None

    # 可视化参数|Visualization parameters
    generate_plots: bool = False
    bin_size: int = 500
    min_len: float = 1.0
    separate_plots: bool = False

    # 性能参数|Performance parameters
    threads: int = 64
    memory_limit: Optional[str] = "300G"

    # 断点续传参数|Resume parameters
    force_rerun: bool = False  # 强制重新运行所有步骤

    # 高级选项|Advanced options
    quick_view: bool = False
    re_sites: str = "GATC"
    skip_clustering: bool = False

    # 输出格式选项|Output format options
    output_agp: bool = True
    output_fasta: bool = True
    output_juicebox: bool = True

    # Juicebox配置|Juicebox configuration
    generate_juicebox: bool = True
    matlock_bin: str = field(
        default_factory=lambda: get_tool_path(
            'matlock', '~/miniforge3/envs/juicer_v.1.6/bin/matlock', 'MATLOCK_PATH'
        )
    )
    three_d_dna_dir: str = field(
        default_factory=lambda: get_tool_path(
            'three_d_dna', '~/software/3d-dna', 'THREE_D_DNA_DIR'
        )
    )
    agp2assembly_script: str = field(
        default_factory=lambda: get_tool_path(
            'agp2assembly', '~/software/3d-dna/utils/agp2assembly.py', 'AGP2ASSEMBLY_PATH'
        )
    )
    asm_visualizer_script: str = field(
        default_factory=lambda: get_tool_path(
            'asm_visualizer', '~/software/3d-dna/visualize/run-assembly-visualizer.sh', 'ASM_VISUALIZER_PATH'
        )
    )

    # BWA比对配置|BWA alignment configuration
    samblaster_bin: str = field(
        default_factory=lambda: get_tool_path(
            'samblaster', '~/miniforge3/envs/haphic/bin/samblaster', 'SAMBLASTER_PATH'
        )
    )
    haphic_filter_bam_bin: str = field(
        default_factory=lambda: get_tool_path(
            'filter_bam', '~/miniforge3/envs/haphic/bin/filter_bam', 'FILTER_BAM_PATH'
        )
    )
    use_samblaster: bool = True
    use_haphic_filter: bool = True

    # Pipeline模式新增参数|Pipeline mode additional parameters
    aln_format: str = "auto"  # 比对文件格式|Alignment file format (upstream default: auto)
    normalize_by_nlinks: bool = False  # 按链接数归一化|Normalize by number of links
    dense_matrix: bool = False  # 使用稠密矩阵|Use dense matrix
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
    bin_size_kbp: int = -1  # 聚类分箱大小(kbp)|Clustering bin size (kbp), -1=auto
    flank: int = 500  # 邻接矩阵侧翼区域(kbp)|Adjacency matrix flank region (kbp)
    expansion: int = 2
    max_iter: int = 200
    pruning: float = 0.0001
    max_ctg_len: float = 10000  # Kbp
    min_links: int = 25  # 重分配最小链接数|Min links for reassignment
    min_link_density: float = 0.0001  # 重分配最小链接密度|Min link density for reassignment
    min_density_ratio: float = 4.0
    ambiguous_cutoff: float = 0.6
    reassign_nrounds: int = 5
    no_additional_rescue: bool = False
    skip_fast_sort: bool = False
    flanking_region: int = 0
    density_cal_method: str = "multiplication"
    confidence_cutoff: float = 1.0
    skip_allhic: bool = False
    skip_ga: bool = False  # 跳过ALLHiC遗传算法|Skip ALLHiC genetic algorithm
    mutprob: float = 0.2
    ngen: int = 5000
    npop: int = 100
    seed: int = 42
    Ns: int = 100
    max_width: int = 60
    sort_by_input: bool = False
    keep_letter_case: bool = False

    # 超长读长参数|Ultra-long read parameters
    ul_file: Optional[str] = None  # 超长读长BAM文件|Ultra-long read BAM file
    min_ul_mapq: int = 30
    min_ul_alignment_length: int = 10000  # bp
    max_distance_to_end: int = 100  # bp
    max_overlap_ratio: float = 0.5
    max_gap_len: int = 10000  # bp
    min_ul_support: int = 2

    # 日志配置|Logging configuration
    verbose: bool = False
    log_file: Optional[str] = None

    # 测试模式配置|Test mode configuration
    dry_run: bool = False

    def _expand_tool_path(self, path: str) -> str:
        """
        展开工具路径|Expand tool path

        对于命令名（不包含/），保持原样让系统通过PATH查找
        对于包含/的路径，展开~并转为绝对路径

        Args:
            path: 工具路径或命令名|Tool path or command name

        Returns:
            展开后的路径|Expanded path
        """
        # 如果是纯命令名（不包含路径分隔符），保持原样
        # If it's a pure command name (no path separator), keep as is
        if '/' not in path:
            return path

        # 包含路径分隔符，展开~并转为绝对路径
        # Contains path separator, expand ~ and convert to absolute path
        expanded = expand_path(path)
        return os.path.normpath(os.path.abspath(expanded))

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 标准化路径|Normalize paths
        self.asm_file = os.path.normpath(os.path.abspath(self.asm_file))

        # 只有在bam_file不为None时才处理|Process bam_file only if not None
        if self.bam_file is not None:
            self.bam_file = os.path.normpath(os.path.abspath(self.bam_file))

        # 处理工具路径|Process tool paths - 使用expand_path展开~符号（遵循开发规范第11.3.1节）
        # 对于命令名，保持原样让系统通过PATH查找
        # For command names, keep as is to let system find via PATH
        if self.haphic_bin:
            self.haphic_bin = self._expand_tool_path(self.haphic_bin)

        if self.bwa_bin:
            self.bwa_bin = self._expand_tool_path(self.bwa_bin)

        if self.samtools_bin:
            self.samtools_bin = self._expand_tool_path(self.samtools_bin)

        if self.matlock_bin:
            self.matlock_bin = self._expand_tool_path(self.matlock_bin)

        if self.samblaster_bin:
            self.samblaster_bin = self._expand_tool_path(self.samblaster_bin)

        if self.haphic_filter_bam_bin:
            self.haphic_filter_bam_bin = self._expand_tool_path(self.haphic_filter_bam_bin)

        # 展开3d-dna相关脚本路径|Expand 3d-dna script paths
        if self.three_d_dna_dir:
            self.three_d_dna_dir = os.path.normpath(os.path.abspath(expand_path(self.three_d_dna_dir)))

        if self.agp2assembly_script:
            self.agp2assembly_script = os.path.normpath(os.path.abspath(expand_path(self.agp2assembly_script)))

        if self.asm_visualizer_script:
            self.asm_visualizer_script = os.path.normpath(os.path.abspath(expand_path(self.asm_visualizer_script)))

        # 设置输出目录|Set output directory
        if not self.output_dir:
            self.output_dir = os.path.dirname(os.path.abspath(self.asm_file))
        else:
            self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))

        # 设置前缀|Set prefix
        if not self.prefix:
            self.prefix = os.path.splitext(os.path.basename(self.asm_file))[0]

        # 创建输出目录|Create output directory
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # 设置默认日志文件(落99_logs目录)|Set default log file (in 99_logs)
        if not self.log_file:
            self.log_file = os.path.join(self.output_dir, "99_logs", f"{self.prefix}_haphic.log")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.asm_file):
            errors.append(f"基因组组装文件不存在|Assembly file does not exist: {self.asm_file}")

        if not os.path.exists(self.hic_file):
            errors.append(f"Hi-C文件不存在|Hi-C file does not exist: {self.hic_file}")

        # 检查文件格式|Check file formats
        if not self.asm_file.endswith(('.fa', '.fasta', '.FA', '.FASTA')):
            errors.append("基因组文件必须是FASTA格式|Assembly file must be in FASTA format")

        if self.hic_file_type == "bam":
            if not self.hic_file.endswith(('.bam', '.BAM')):
                errors.append("Hi-C文件必须是BAM格式|Hi-C file must be in BAM format when --hic-file-type is bam")
        elif self.hic_file_type == "fastq":
            if os.path.isfile(self.hic_file) and not self.hic_file.endswith(('.fastq.gz', '.fq.gz', '.fastq', '.fq')):
                errors.append("Hi-C文件必须是FASTQ格式|Hi-C file must be in FASTQ format when --hic-file-type is fastq")

        # 检查HapHiC工具|Check HapHiC tool
        if self.haphic_bin != "haphic" and not os.path.exists(self.haphic_bin):
            errors.append(f"HapHiC可执行文件不存在|HapHiC executable does not exist: {self.haphic_bin}")

        # 检查BWA工具 (如果需要)|Check BWA tool (if needed)
        if self.hic_file_type == "fastq":
            if self.bwa_bin != "bwa" and not os.path.exists(self.bwa_bin):
                errors.append(f"BWA可执行文件不存在|BWA executable does not exist: {self.bwa_bin}")
            if self.samtools_bin != "samtools" and not os.path.exists(self.samtools_bin):
                errors.append(f"Samtools可执行文件不存在|Samtools executable does not exist: {self.samtools_bin}")

        # 检查数值参数|Check numeric parameters
        if self.nchrs <= 0:
            errors.append("染色体数必须大于0|Number of chromosomes must be greater than 0")

        if self.mapq_threshold < 0:
            errors.append("MAPQ阈值不能为负数|MAPQ threshold cannot be negative")

        if self.edit_distance < 0:
            errors.append("编辑距离不能为负数|Edit distance cannot be negative")

        if self.min_inflation <= 0:
            errors.append("最小膨胀参数必须大于0|Min inflation must be greater than 0")

        if self.max_inflation <= self.min_inflation:
            errors.append("最大膨胀参数必须大于最小膨胀参数|Max inflation must be greater than min inflation")

        if self.threads < 1:
            errors.append("线程数必须大于等于1|Number of threads must be >= 1")

        if self.processes < 1:
            errors.append("进程数必须大于等于1|Number of processes must be >= 1")

        # 检查高级选项|Check advanced options
        if self.quick_view and self.nchrs < 1:
            # Quick view mode doesn't care about nchrs, but needs at least 1 for parameter validation
            pass

        if self.remove_allelic_links and self.remove_allelic_links <= 0:
            errors.append("等位基因链接移除参数必须大于0|Remove allelic links parameter must be greater than 0")

        if self.phasing_weight < 0 or self.phasing_weight > 1:
            errors.append("分相权重必须在0-1之间|Phasing weight must be between 0 and 1")

        # 检查RE位点参数|Check restriction site parameters
        if not self.re_sites:
            errors.append("限制性内切酶位点不能为空|Restriction enzyme site cannot be empty")

        if self.bin_size <= 0:
            errors.append("装箱大小必须大于0|Bin size must be greater than 0")

        if self.min_len <= 0:
            errors.append("最小长度必须大于0|Min length must be greater than 0")

        # 检查Juicebox工具(如果启用)|Check Juicebox tools (if enabled)
        # 仅对路径(含/)做存在性检查,命令名(无/)交给conda run处理
        # Only check existence for paths (with /); command names handled by conda run
        if self.generate_juicebox:
            if '/' in self.matlock_bin and not os.path.exists(self.matlock_bin):
                errors.append(f"matlock可执行文件不存在|matlock executable not found: {self.matlock_bin}")
            if '/' in self.agp2assembly_script and not os.path.exists(self.agp2assembly_script):
                errors.append(f"agp2assembly脚本不存在|agp2assembly script not found: {self.agp2assembly_script}")
            if '/' in self.asm_visualizer_script and not os.path.exists(self.asm_visualizer_script):
                errors.append(f"asm-visualizer脚本不存在|asm-visualizer script not found: {self.asm_visualizer_script}")

        # 检查输出目录可写性|Check output directory writability
        try:
            test_file = os.path.join(self.output_dir, '.haphic_test')
            with open(test_file, 'w') as f:
                f.write('test')
            os.remove(test_file)
        except (PermissionError, OSError) as e:
            errors.append(f"输出目录不可写|Output directory not writable: {e}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_restriction_sites(self):
        """获取限制性内切酶位点列表|Get restriction enzyme sites list"""
        return [site.strip() for site in self.re_sites.split(',') if site.strip()]

    def get_output_files(self):
        """获取输出文件路径|Get output file paths"""
        base_path = os.path.join(self.output_dir, self.prefix)
        build_dir = resolve_legacy_path(self.output_dir, "04_build")

        return {
            "corrected_asm": f"{base_path}_corrected.asm.fa",
            "scaffolds_agp": os.path.join(build_dir, f"{self.prefix}.agp"),
            "scaffolds_raw_agp": os.path.join(build_dir, f"{self.prefix}.raw.agp"),
            "scaffolds_fasta": os.path.join(build_dir, f"{self.prefix}.fa"),
            "juicebox_script": os.path.join(build_dir, f"{self.prefix}_juicebox.sh"),
            "contact_map": f"{base_path}_contact_map.pdf",
            "hic_file": os.path.join(resolve_legacy_path(self.output_dir, "06_juicebox"), f"{self.prefix}.hic"),
            "assembly_file": os.path.join(resolve_legacy_path(self.output_dir, "06_juicebox"), f"{self.prefix}.assembly"),
            "log_file": self.log_file
        }

    def get_summary(self):
        """获取配置摘要|Get configuration summary"""
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

    def get_software_info(self) -> dict:
        """
        获取软件版本信息|Get software version information

        直调工具完整路径探测版本(不走conda run):haphic是python脚本但--version不需重依赖;
        执行期仍由build_conda_command包装(脚本运行需env激活)。两套各司其职。
        |Probe versions by calling tool full path directly (no conda run): haphic is a
        python script but --version needs no heavy deps; execution still goes through
        build_conda_command (script run needs env activation). Two mechanisms, separate roles.

        Returns:
            软件信息字典|Software information dictionary
        """
        import subprocess

        info = {
            'pipeline': {
                'name': 'biopytools haphic',
                'version': '1.0.0'
            },
            'tools': {},
            'parameters': {
                'nchrs': self.nchrs,
                'threads': self.threads,
                'processes': self.processes,
                'RE': self.re_sites,
                'correct_nrounds': self.correct_nrounds,
                'hic_file_type': self.hic_file_type
            }
        }

        # 检测各工具版本|Detect tool versions
        # haphic/samtools/samblaster 用 --version;bwa/matlock 无版本子命令(落unknown)
        # haphic/samtools/samblaster use --version; bwa/matlock have no version subcommand (unknown)
        tools_to_check = {
            'haphic': self.haphic_bin,
            'bwa': self.bwa_bin,
            'samtools': self.samtools_bin,
            'matlock': self.matlock_bin,
            'samblaster': self.samblaster_bin
        }

        for tool_name, tool_path in tools_to_check.items():
            try:
                result = subprocess.run(
                    [tool_path, '--version'] if tool_name in ['haphic', 'samtools', 'samblaster'] else [tool_path],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                version = result.stdout.strip() or result.stderr.strip()
                info['tools'][tool_name] = {
                    'version': version.split('\n')[0] if version else 'unknown',
                    'path': tool_path
                }
            except Exception:
                info['tools'][tool_name] = {
                    'version': 'unknown',
                    'path': tool_path
                }

        return info