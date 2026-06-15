"""
YaHS配置管理模块|YaHS Configuration Management Module
"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List
import os
import shutil


@dataclass
class YaHSConfig:
    """
    YaHS流程配置类|YaHS Pipeline Configuration Class

    Attributes:
        ref_fa: 参考基因组FASTA文件|Reference genome FASTA file
        hic_r1: Hi-C R1测序文件|Hi-C R1 sequencing file
        hic_r2: Hi-C R2测序文件|Hi-C R2 sequencing file
        output_dir: 输出目录|Output directory
        threads: 线程数|Number of threads
    """

    # 必需参数|Required parameters
    ref_fa: str
    hic_r1: str
    hic_r2: str

    # 输出配置|Output configuration
    output_dir: str = './yahs_output'

    # 资源配置|Resource configuration
    threads: int = 12
    java_ram: str = '32G'  # 默认值，会被自动调整|Default, will be auto-adjusted
    sam_ram: str = '4G'  # 默认值，会被自动调整|Default, will be auto-adjusted

    # YaHS 核心参数|YaHS core parameters
    enzyme_seq: str = 'GATC'
    min_len: int = 10000
    min_mapq: int = 30
    no_contig_ec: bool = False
    no_scaffold_ec: bool = False

    # YaHS 可选参数|YaHS optional parameters
    resolutions: Optional[str] = None  # 默认由YaHS自动设置
    rounds_per_resolution: int = 1
    telo_motif: Optional[str] = None

    # 工具路径（使用完整路径作为默认值）|Tool paths (use full paths as defaults)
    yahs_bin: str = field(default_factory=lambda: '~/miniforge3/envs/yahs_v.1.2.2/bin/yahs')
    juicer_bin: str = field(default_factory=lambda: '~/miniforge3/envs/yahs_v.1.2.2/bin/juicer')
    juicer_jar: Optional[str] = field(default_factory=lambda: '~/software/juicer/scripts/juicer_tools.jar')
    bwa_bin: str = field(default_factory=lambda: '~/miniforge3/envs/Population_genetics/bin/bwa')
    samtools_bin: str = field(default_factory=lambda: '~/miniforge3/envs/GATK_v.4.6.2.0/bin/samtools')
    java_cmd: str = field(default_factory=lambda: 'java')

    # 执行控制|Execution control
    step: Optional[str] = None  # None = 运行所有步骤|None = run all steps
    force_rerun: bool = False
    keep_temp: bool = False

    def __post_init__(self):
        """
        初始化后处理|Post-initialization processing

        展开~路径，创建输出目录，设置工具路径
        """
        from ..common.paths import expand_path

        # 展开用户输入路径|Expand user input paths
        self.ref_fa = expand_path(self.ref_fa)
        self.hic_r1 = expand_path(self.hic_r1)
        self.hic_r2 = expand_path(self.hic_r2)
        self.output_dir = expand_path(self.output_dir)

        # 展开工具路径（所有工具路径都需要展开）|Expand all tool paths
        self.yahs_bin = expand_path(self.yahs_bin)
        self.juicer_bin = expand_path(self.juicer_bin)
        if self.juicer_jar:
            self.juicer_jar = expand_path(self.juicer_jar)
        self.bwa_bin = expand_path(self.bwa_bin)
        self.samtools_bin = expand_path(self.samtools_bin)

        # 创建输出目录结构|Create output directory structure
        self.output_path = Path(self.output_dir)
        self._create_directory_structure()

        # 设置各步骤输出路径|Set step output paths
        self._setup_step_paths()

        # 设置资源参数|Set resource parameters
        self._setup_resources()

    def _create_directory_structure(self):
        """创建输出目录结构|Create output directory structure"""
        subdirs = [
            '00_pipeline_info',
            '01_indexing',
            '02_mapping',
            '03_scaffolding',
            '04_hic_standard',
            '05_jbat',
            '06_assessment',
            '99_logs'
        ]

        for subdir in subdirs:
            (self.output_path / subdir).mkdir(parents=True, exist_ok=True)

    def _setup_step_paths(self):
        """设置各步骤输出路径|Setup step output paths"""
        ref_name = Path(self.ref_fa).stem

        # 步骤1: 索引
        self.step1_dir = self.output_path / '01_indexing'
        self.ref_index = {
            'bwt': self.step1_dir / f'{ref_name}.bwt',
            'pac': self.step1_dir / f'{ref_name}.pac',
            'ann': self.step1_dir / f'{ref_name}.ann',
            'amb': self.step1_dir / f'{ref_name}.amb',
            'sa': self.step1_dir / f'{ref_name}.sa',
            'fai': self.step1_dir / f'{ref_name}.fai'
        }

        # 步骤2: 比对
        self.step2_dir = self.output_path / '02_mapping'
        self.final_bam = self.step2_dir / 'aligned_sorted_dedup.bam'
        self.bam_index = self.step2_dir / 'aligned_sorted_dedup.bam.bai'

        # 步骤3: YaHS挂载
        self.step3_dir = self.output_path / '03_scaffolding'
        self.yahs_out_prefix = self.step3_dir / 'yahs_out'
        self.final_scaffold_fa = self.step3_dir / 'yahs_out_scaffolds_final.fa'
        self.final_scaffold_agp = self.step3_dir / 'yahs_out_scaffolds_final.agp'
        self.yahs_bin_file = self.step3_dir / 'yahs_out.bin'

        # 步骤4: 标准Hi-C热图
        self.step4_dir = self.output_path / '04_hic_standard'
        self.final_hic = self.step4_dir / 'yahs_out_final.hic'

        # 步骤5: JBAT文件
        self.step5_dir = self.output_path / '05_jbat'
        self.jbat_hic = self.step5_dir / 'out_JBAT.hic'
        self.jbat_assembly = self.step5_dir / 'out_JBAT.assembly'
        self.jbat_liftover_agp = self.step5_dir / 'out_JBAT.liftover.agp'

        # 步骤6: 质量评估
        self.step6_dir = self.output_path / '06_assessment'
        self.assembly_metrics = self.step6_dir / 'assembly_metrics.txt'

        # 日志目录|Log directory
        self.log_dir = self.output_path / '99_logs'
        self.pipeline_log = self.log_dir / 'yahs_pipeline.log'

        # 将所有步骤目录解析为绝对路径|Resolve all step directories to absolute paths
        # 这确保所有基于它们的文件路径也是绝对路径|This ensures all file paths based on them are absolute
        self.step1_dir = self.step1_dir.resolve()
        self.step2_dir = self.step2_dir.resolve()
        self.step3_dir = self.step3_dir.resolve()
        self.step4_dir = self.step4_dir.resolve()
        self.step5_dir = self.step5_dir.resolve()
        self.step6_dir = self.step6_dir.resolve()
        self.log_dir = self.log_dir.resolve()

        # 重新定义基于step_dir的文件路径（现在使用绝对路径）|Redefine paths based on step_dir (now absolute)
        ref_name = Path(self.ref_fa).stem
        self.ref_index = {
            'bwt': self.step1_dir / f'{ref_name}.bwt',
            'pac': self.step1_dir / f'{ref_name}.pac',
            'ann': self.step1_dir / f'{ref_name}.ann',
            'amb': self.step1_dir / f'{ref_name}.amb',
            'sa': self.step1_dir / f'{ref_name}.sa',
            'fai': self.step1_dir / f'{ref_name}.fai'
        }
        self.final_bam = self.step2_dir / 'aligned_sorted_dedup.bam'
        self.bam_index = self.step2_dir / 'aligned_sorted_dedup.bam.bai'
        self.yahs_out_prefix = self.step3_dir / 'yahs_out'
        self.final_scaffold_fa = self.step3_dir / 'yahs_out_scaffolds_final.fa'
        self.final_scaffold_agp = self.step3_dir / 'yahs_out_scaffolds_final.agp'
        self.yahs_bin_file = self.step3_dir / 'yahs_out.bin'
        self.final_hic = self.step4_dir / 'yahs_out_final.hic'
        self.jbat_hic = self.step5_dir / 'out_JBAT.hic'
        self.jbat_assembly = self.step5_dir / 'out_JBAT.assembly'
        self.jbat_liftover_agp = self.step5_dir / 'out_JBAT.liftover.agp'
        self.assembly_metrics = self.step6_dir / 'assembly_metrics.txt'
        self.pipeline_log = self.log_dir / 'yahs_pipeline.log'

    def _setup_resources(self):
        """设置计算资源参数|Setup compute resource parameters"""
        import shutil

        # 自动检测可用资源|Auto-detect available resources
        total_mem_gb = self._get_total_memory_gb()
        cpu_cores = self._get_cpu_cores()

        # 调整默认资源分配|Adjust default resource allocation
        if self.threads == 12:  # 使用默认值|Using default value
            self.threads = max(2, cpu_cores - 2)

        if self.java_ram == '32G':  # 使用默认值|Using default value
            self.java_ram = f'{int(total_mem_gb * 0.6)}G'

        if self.sam_ram == '4G':  # 使用默认值|Using default value
            # 自动调整内存，但至少300G或总内存的50%，取较大值
            # Auto-adjust memory, but at least 300G or 50% of total, whichever is larger
            auto_ram = int(total_mem_gb * 0.5)
            self.sam_ram = f'{max(auto_ram, 300)}G'

        # 检测工具可用性|Check tool availability
        self._check_tool_availability()

    def _get_total_memory_gb(self) -> int:
        """获取总内存(GB)|Get total memory in GB"""
        try:
            with open('/proc/meminfo', 'r') as f:
                for line in f:
                    if line.startswith('MemTotal:'):
                        kb = int(line.split()[1])
                        return kb // (1024 * 1024)
        except Exception:
            pass
        return 64  # 默认64GB|Default 64GB

    def _get_cpu_cores(self) -> int:
        """获取CPU核心数|Get CPU core count"""
        try:
            return os.cpu_count() or 12
        except Exception:
            return 12

    def _check_tool_availability(self):
        """检查必需工具是否可用|Check required tools availability"""
        missing_tools = []

        # 检查常用工具（优先使用系统PATH）|Check common tools (prefer system PATH)
        common_tools = {
            'bwa': self.bwa_bin,
            'samtools': self.samtools_bin,
        }

        for tool_name, tool_path in common_tools.items():
            if shutil.which(tool_path) is None:
                missing_tools.append(f'{tool_name} ({tool_path})')

        # YaHS 和 juicer 在实际使用时检查，不在初始化时检查
        # YaHS and juicer are checked during actual use, not during initialization

        # 检查 Java（如果需要使用juicer_tools.jar）|Check Java (if using juicer_tools.jar)
        if self.juicer_jar:
            if shutil.which(self.java_cmd) is None:
                missing_tools.append(f'java ({self.java_cmd})')
            if not Path(self.juicer_jar).exists():
                missing_tools.append(f'juicer_tools.jar ({self.juicer_jar})')

        if missing_tools:
            raise ValueError(
                f"以下工具未找到|The following tools are not found:\n" +
                "\n".join(f"  - {tool}" for tool in missing_tools) +
                "\n\n请确保工具已安装或在PATH中可用|Please ensure tools are installed or available in PATH"
            )

    def validate(self) -> None:
        """
        验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 如果配置无效|If configuration is invalid
        """
        errors = []

        # 检查输入文件|Check input files
        if not Path(self.ref_fa).exists():
            errors.append(f"参考基因组文件不存在|Reference genome file not found: {self.ref_fa}")

        if not Path(self.hic_r1).exists():
            errors.append(f"Hi-C R1文件不存在|Hi-C R1 file not found: {self.hic_r1}")

        if not Path(self.hic_r2).exists():
            errors.append(f"Hi-C R2文件不存在|Hi-C R2 file not found: {self.hic_r2}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if self.min_len < 0:
            errors.append(f"最小长度不能为负数|Minimum length cannot be negative: {self.min_len}")

        if self.min_mapq < 0 or self.min_mapq > 255:
            errors.append(f"MAPQ必须在0-255之间|MAPQ must be between 0-255: {self.min_mapq}")

        # 检查酶切位点序列|Check enzyme sequence
        if not self.enzyme_seq or not all(c in 'ATCG' for c in self.enzyme_seq.upper()):
            errors.append(f"酶切位点序列无效|Invalid enzyme sequence: {self.enzyme_seq}")

        # 检查步骤参数|Check step parameter
        if self.step is not None:
            valid_steps = ['1', '2', '3', '4', '5', '6']
            if self.step not in valid_steps:
                errors.append(f"步骤参数必须为{valid_steps}之一|Step must be one of {valid_steps}: {self.step}")

        if errors:
            raise ValueError("\n".join(errors))

    def get_yahs_params(self) -> List[str]:
        """
        构建YaHS命令参数|Build YaHS command parameters

        Returns:
            YaHS参数列表|YaHS parameter list
        """
        params = []

        if self.enzyme_seq:
            params.extend(['-e', self.enzyme_seq])

        params.extend(['-q', str(self.min_mapq)])
        params.extend(['-l', str(self.min_len)])

        if self.no_contig_ec:
            params.append('--no-contig-ec')

        if self.no_scaffold_ec:
            params.append('--no-scaffold-ec')

        if self.resolutions:
            params.extend(['-r', self.resolutions])

        if self.rounds_per_resolution != 1:
            params.extend(['-R', str(self.rounds_per_resolution)])

        if self.telo_motif:
            params.extend(['--telo-motif', self.telo_motif])

        return params

    def get_software_info(self) -> dict:
        """
        获取软件版本信息|Get software version information

        Returns:
            软件信息字典|Software information dictionary
        """
        import subprocess

        info = {
            'pipeline': {
                'name': 'biopytools yahs',
                'version': '1.0.0'
            },
            'tools': {},
            'parameters': {
                'enzyme_seq': self.enzyme_seq,
                'min_len': self.min_len,
                'min_mapq': self.min_mapq,
                'no_contig_ec': self.no_contig_ec,
                'threads': self.threads
            }
        }

        # 检测各工具版本|Detect tool versions
        tools_to_check = {
            'yahs': self.yahs_bin,
            'juicer': self.juicer_bin,
            'bwa': self.bwa_bin,
            'samtools': self.samtools_bin
        }

        for tool_name, tool_path in tools_to_check.items():
            try:
                result = subprocess.run(
                    [tool_path, '--version'] if tool_name in ['bwa', 'samtools'] else [tool_path],
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

        # Java版本
        try:
            result = subprocess.run(
                [self.java_cmd, '-version'],
                capture_output=True,
                text=True,
                timeout=5
            )
            version = result.stderr.strip() or result.stdout.strip()
            info['tools']['java'] = {
                'version': version.split('\n')[0] if version else 'unknown',
                'path': self.java_cmd
            }
        except Exception:
            pass

        return info
