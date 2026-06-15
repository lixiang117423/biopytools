"""
Samplot Configuration
Samplot配置类
"""

import os
from typing import List, Optional

from common.paths import expand_path, get_tool_path


class SamplotPlotConfig:
    """Samplot plot配置类|Samplot plot configuration class"""

    def __init__(
        self,
        bams: List[str],
        chrom: str,
        start: int,
        end: int,
        sv_type: Optional[str] = None,
        output_file: Optional[str] = None,
        output_dir: str = ".",
        reference: Optional[str] = None,
        max_depth: int = 1,
        window: Optional[int] = None,
        z: int = 4,
        plot_height: Optional[int] = None,
        plot_width: int = 8,
        dpi: int = 300,
        long_read: int = 1000,
        coverage_only: bool = False,
        same_yaxis_scales: bool = False,
        titles: Optional[List[str]] = None,
        samplot_path: str = "~/miniforge3/envs/samplot_v.1.3.0/bin/samplot",
    ):
        """
        初始化配置|Initialize configuration

        Args:
            bams: BAM/CRAM文件路径列表|BAM/CRAM file paths
            chrom: 染色体名称|Chromosome name
            start: 起始位置|Start position
            end: 结束位置|End position
            sv_type: SV类型(DEL/DUP/INV/BND)|SV type
            output_file: 输出文件名|Output file name
            output_dir: 输出目录|Output directory
            reference: 参考基因组(CRAM必须)|Reference genome (required for CRAM)
            max_depth: 最大正常pair数|Max normal pairs to plot
            window: 窗口大小|Window size
            z: 标准差倍数|Number of stdevs from mean
            plot_height: 图高|Plot height
            plot_width: 图宽|Plot width
            dpi: DPI|Dots per inch
            long_read: 长读长最小长度|Min length for long-read
            coverage_only: 仅显示覆盖度|Show only coverage
            same_yaxis_scales: 统一Y轴|Use same Y-axis scales
            titles: 样本标题列表|Sample title list
            samplot_path: samplot可执行文件路径|samplot binary path
        """
        self.bams = bams
        self.chrom = chrom
        self.start = start
        self.end = end
        self.sv_type = sv_type
        self.output_file = output_file
        self.output_dir = output_dir
        self.reference = reference
        self.max_depth = max_depth
        self.window = window
        self.z = z
        self.plot_height = plot_height
        self.plot_width = plot_width
        self.dpi = dpi
        self.long_read = long_read
        self.coverage_only = coverage_only
        self.same_yaxis_scales = same_yaxis_scales
        self.titles = titles
        self.samplot_path = samplot_path

    def validate(self):
        """
        验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时|When configuration parameters are invalid
        """
        # 展开路径|Expand paths
        self.samplot_path = expand_path(self.samplot_path)
        self.output_dir = expand_path(self.output_dir)

        if self.output_file:
            self.output_file = expand_path(self.output_file)
        if self.reference:
            self.reference = expand_path(self.reference)

        # 展开BAM路径|Expand BAM paths
        self.bams = [expand_path(b) for b in self.bams]

        if self.titles:
            self.titles = [expand_path(t) if os.path.exists(expand_path(t)) else t for t in self.titles]

        # 检查BAM文件|Check BAM files
        for bam in self.bams:
            if not os.path.exists(bam):
                raise ValueError(f"BAM文件不存在|BAM file not found: {bam}")

        # 检查参考基因组|Check reference genome
        if self.reference and not os.path.exists(self.reference):
            raise ValueError(f"参考基因组不存在|Reference genome not found: {self.reference}")

        # 检查位置参数|Check position parameters
        if self.start < 0:
            raise ValueError(f"起始位置不能为负数|Start position cannot be negative: {self.start}")
        if self.end <= self.start:
            raise ValueError(f"结束位置必须大于起始位置|End must be greater than start: {self.start}-{self.end}")

        # 创建输出目录|Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

        return True

    def __repr__(self):
        return (
            f"SamplotPlotConfig(\n"
            f"  bams={self.bams!r},\n"
            f"  chrom={self.chrom!r},\n"
            f"  start={self.start},\n"
            f"  end={self.end},\n"
            f"  sv_type={self.sv_type!r},\n"
            f"  output_dir={self.output_dir!r}\n"
            f")"
        )


class SamplotVcfConfig:
    """Samplot vcf配置类|Samplot vcf configuration class"""

    def __init__(
        self,
        bams: List[str],
        vcf: str,
        output_dir: str = "samplot-out",
        output_type: str = "png",
        threads: int = 1,
        downsample: int = 1,
        min_bp: int = 20,
        max_mb: Optional[int] = None,
        sample_ids: Optional[List[str]] = None,
        plot_all: bool = False,
        min_call_rate: Optional[float] = None,
        max_hets: Optional[int] = None,
        min_entries: int = 6,
        max_entries: int = 10,
        gff3: Optional[str] = None,
        samplot_path: str = "~/miniforge3/envs/samplot_v.1.3.0/bin/samplot",
    ):
        """
        初始化配置|Initialize configuration

        Args:
            bams: BAM/CRAM文件路径列表|BAM/CRAM file paths
            vcf: VCF文件路径|VCF file path
            output_dir: 输出目录|Output directory
            output_type: 输出格式(png/pdf/eps/jpg)|Output format
            threads: 线程数|Number of threads
            downsample: 下采样数|Downsample count
            min_bp: 最小SV长度|Min SV length in bp
            max_mb: 最大SV长度(MB)|Max SV length in MB
            sample_ids: 样本ID列表|Sample ID list
            plot_all: 绘制所有样本|Plot all samples
            min_call_rate: 最小call rate|Min call rate
            max_hets: 最大杂合数|Max heterozygotes
            min_entries: 最小样本数|Min entries to plot
            max_entries: 最大样本数|Max entries to plot
            gff3: GFF3注释文件|GFF3 annotation file
            samplot_path: samplot可执行文件路径|samplot binary path
        """
        self.bams = bams
        self.vcf = vcf
        self.output_dir = output_dir
        self.output_type = output_type
        self.threads = threads
        self.downsample = downsample
        self.min_bp = min_bp
        self.max_mb = max_mb
        self.sample_ids = sample_ids
        self.plot_all = plot_all
        self.min_call_rate = min_call_rate
        self.max_hets = max_hets
        self.min_entries = min_entries
        self.max_entries = max_entries
        self.gff3 = gff3
        self.samplot_path = samplot_path

    def validate(self):
        """
        验证配置参数|Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时|When configuration parameters are invalid
        """
        # 展开路径|Expand paths
        self.samplot_path = expand_path(self.samplot_path)
        self.output_dir = expand_path(self.output_dir)
        self.vcf = expand_path(self.vcf)
        self.bams = [expand_path(b) for b in self.bams]

        if self.gff3:
            self.gff3 = expand_path(self.gff3)

        # 检查VCF文件|Check VCF file
        if not self.vcf or not os.path.exists(self.vcf):
            raise ValueError(f"VCF文件不存在|VCF file not found: {self.vcf}")

        # 检查BAM文件|Check BAM files
        for bam in self.bams:
            if not os.path.exists(bam):
                raise ValueError(f"BAM文件不存在|BAM file not found: {bam}")

        # 检查GFF3文件|Check GFF3 file
        if self.gff3 and not os.path.exists(self.gff3):
            raise ValueError(f"GFF3文件不存在|GFF3 file not found: {self.gff3}")

        # 检查输出格式|Check output type
        valid_types = ["png", "pdf", "eps", "jpg"]
        if self.output_type not in valid_types:
            raise ValueError(f"无效的输出格式|Invalid output type: {self.output_type}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            raise ValueError(f"线程数必须为正数|Threads must be positive: {self.threads}")
        if self.min_bp < 0:
            raise ValueError(f"最小SV长度不能为负数|Min bp cannot be negative: {self.min_bp}")

        # 创建输出目录|Create output directory
        os.makedirs(self.output_dir, exist_ok=True)

        return True

    def __repr__(self):
        return (
            f"SamplotVcfConfig(\n"
            f"  bams={self.bams!r},\n"
            f"  vcf={self.vcf!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  output_type={self.output_type!r},\n"
            f"  threads={self.threads}\n"
            f")"
        )
