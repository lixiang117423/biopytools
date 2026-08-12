"""
CentIER着丝粒鉴定配置管理模块|CentIER Centromere Identification Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional
from ..common.paths import expand_path, get_tool_path


@dataclass
class CentIERConfig:
    """CentIER着丝粒鉴定配置类|CentIER Centromere Identification Configuration Class"""

    # 必需文件|Required files
    genome_fasta: str

    # 路径配置|Path configuration
    centier_path: str = field(
        default_factory=lambda: get_tool_path(
            'centier',
            '~/software/CentIER/CentIER-main',
            'CENTIER_PATH'
        )
    )
    output_dir: str = './centier_output'

    # 可选文件|Optional files
    gff_annotation: Optional[str] = None

    # Hi-C数据文件(可选)|Hi-C data files (optional)
    matrix1: Optional[str] = None
    matrix2: Optional[str] = None
    bed1: Optional[str] = None
    bed2: Optional[str] = None

    # Hi-C FASTQ 自动模式(可选,提供即启用)|Hi-C FASTQ auto mode (optional, enables auto mode)
    fastq_r1: Optional[str] = None
    fastq_r2: Optional[str] = None
    genome_id: Optional[str] = None          # bowtie2 索引命名,未给则从 genome 推导|for bowtie2 index naming
    restriction_enzyme: str = 'MboI'          # 限制性内切酶|Restriction enzyme
    bowtie2_idx: Optional[str] = None         # 无则 HiC-Pro 自动建|auto-built if None
    bin_sizes: str = '100000 20000'           # centier 需要的两个分辨率|two resolutions centier needs
    max_memory_gb: int = 200                  # HiC-Pro SORT_RAM(GB)|HiC-Pro memory
    force_hicpro: bool = False                # 强制重跑 HiC-Pro|Force rerun HiC-Pro
    hic_matrix_type: str = 'raw'              # raw|iced,选 HiC-Pro 产物的子目录|matrix subdir to use
    strict_chrname: bool = False              # True 时 ChrN 预检失败即中止|Abort if chr naming not ChrN

    # 分析参数|Analysis parameters
    threads: int = 12
    kmer_size: int = 21
    center_tolerance: int = 15
    step_len: int = 10000
    mul_cents: bool = False
    mingap: int = 2
    signal_threshold: float = 0.7

    # 步骤控制|Step control
    step: Optional[int] = None  # 1-6, None表示运行全部步骤|None means run all steps

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 展开路径|Expand paths
        self.centier_path = expand_path(self.centier_path)
        self.genome_fasta = expand_path(self.genome_fasta)
        self.output_dir = expand_path(self.output_dir)

        # 展开可选路径|Expand optional paths
        if self.gff_annotation:
            self.gff_annotation = expand_path(self.gff_annotation)
        if self.matrix1:
            self.matrix1 = expand_path(self.matrix1)
        if self.matrix2:
            self.matrix2 = expand_path(self.matrix2)
        if self.bed1:
            self.bed1 = expand_path(self.bed1)
        if self.bed2:
            self.bed2 = expand_path(self.bed2)

        # 展开 Hi-C FASTQ 模式路径|Expand Hi-C FASTQ mode paths
        if self.fastq_r1:
            self.fastq_r1 = expand_path(self.fastq_r1)
        if self.fastq_r2:
            self.fastq_r2 = expand_path(self.fastq_r2)
        if self.bowtie2_idx:
            self.bowtie2_idx = expand_path(self.bowtie2_idx)

        # 创建输出目录|Create output directory
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # Hi-C 模式:推导 genome_id 和 sample_name|Hi-C mode: derive genome_id and sample_name
        if self.fastq_r1:
            if not self.genome_id:
                self.genome_id = self._extract_genome_id()
            self.sample_name = self._extract_sample_name()
            # 创建 02_centier 子目录|Create 02_centier subdir
            (self.output_path / '02_centier').mkdir(parents=True, exist_ok=True)

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.genome_fasta):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome_fasta}")

        # 检查CentIER路径|Check CentIER path
        if not os.path.exists(self.centier_path):
            errors.append(f"CentIER路径不存在|CentIER path does not exist: {self.centier_path}")

        # 检查可选文件|Check optional files
        if self.gff_annotation and not os.path.exists(self.gff_annotation):
            errors.append(f"GFF注释文件不存在|GFF annotation file not found: {self.gff_annotation}")

        # Hi-C FASTQ 自动模式校验|Hi-C FASTQ auto mode validation
        if self.fastq_r1 or self.fastq_r2:
            if not (self.fastq_r1 and self.fastq_r2):
                errors.append("Hi-C FASTQ 模式需要 R1 和 R2 成对提供|"
                              "Hi-C FASTQ mode requires both R1 and R2")
            else:
                if not os.path.exists(self.fastq_r1):
                    errors.append(f"R1 文件不存在|R1 file not found: {self.fastq_r1}")
                if not os.path.exists(self.fastq_r2):
                    errors.append(f"R2 文件不存在|R2 file not found: {self.fastq_r2}")
            bins = self.bin_sizes.split()
            if '100000' not in bins:
                errors.append("bin_sizes 必须包含 100000|bin_sizes must include 100000")
            if '20000' not in bins:
                errors.append("bin_sizes 必须包含 20000|bin_sizes must include 20000")
            if self.hic_matrix_type not in ('raw', 'iced'):
                errors.append(f"hic_matrix_type 必须为 raw 或 iced|"
                              f"hic_matrix_type must be raw or iced: {self.hic_matrix_type}")

        # 检查Hi-C数据完整性|Check Hi-C data completeness
        hic_files = [self.matrix1, self.matrix2, self.bed1, self.bed2]
        if any(hic_files) and not all(hic_files):
            errors.append("Hi-C分析需要所有4个文件(medtrx1, matrix2, bed1, bed2)|"
                         "Hi-C analysis requires all 4 files (matrix1, matrix2, bed1, bed2)")

        if any(hic_files):
            for i, f in enumerate(hic_files):
                if f and not os.path.exists(f):
                    errors.append(f"Hi-C文件不存在|Hi-C file not found: {f}")

        # 检查步骤参数|Check step parameter
        if self.step is not None and self.step not in [1, 2, 3, 4, 5, 6]:
            errors.append(f"无效的步骤编号|Invalid step number: {self.step} (应为1-6|should be 1-6)")

        # 检查参数范围|Check parameter ranges
        if self.kmer_size <= 0:
            errors.append("kmer_size必须为正数|kmer_size must be positive")

        if self.step_len <= 0:
            errors.append("step_len必须为正数|step_len must be positive")

        if not (0 <= self.signal_threshold <= 1):
            errors.append("signal_threshold必须在0-1之间|signal_threshold must be between 0 and 1")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_centier_script_path(self) -> str:
        """获取centIER.py脚本路径|Get centIER.py script path"""
        return os.path.join(self.centier_path, 'centIER.py')

    def get_bin_path(self) -> str:
        """获取bin目录路径|Get bin directory path"""
        return os.path.join(self.centier_path, 'bin')

    def _extract_genome_id(self) -> str:
        """从基因组文件名推导 genome_id|Derive genome_id from genome filename"""
        p = Path(self.genome_fasta)
        stem = p.stem
        if stem.endswith('.fa'):
            stem = stem[:-3]
        if stem in ('genome', 'assembly'):
            stem = p.parent.name
        return stem

    def _extract_sample_name(self) -> str:
        """从 R1 文件名推导 sample_name(与 HiC-Pro 同算法)|Derive sample_name matching HiC-Pro's algorithm"""
        name = Path(self.fastq_r1).stem  # 去 .gz|remove .gz
        if name.endswith('.fastq') or name.endswith('.fq'):
            name = Path(name).stem
        name = (name.replace('_R1', '').replace('_R2', '')
                    .replace('_1', '').replace('_2', '')
                    .replace('.clean', ''))
        return name

    def get_centier_output_dir(self) -> Path:
        """CentIER 结果目录(Hi-C 模式 02_centier/,手动模式根)|CentIER output dir"""
        if self.fastq_r1:
            return self.output_path / '02_centier'
        return self.output_path
