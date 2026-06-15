"""
PlotSR配置管理模块|PlotSR Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import List, Optional


@dataclass
class PlotSRConfig:
    """PlotSR配置类|PlotSR Configuration Class"""

    # 输入输出参数|Input and output parameters
    genomes: List[str]  # 基因组FASTA文件列表|List of genome FASTA files
    output_dir: str  # 输出目录|Output directory

    # 基因组参数|Genome parameters
    names: Optional[List[str]] = None  # 基因组名称列表|List of genome names

    # 性能参数|Performance parameters
    threads: int = 12  # 线程数|Number of threads

    # minimap2参数|minimap2 parameters
    minimap2_preset: str = "asm5"  # minimap2预设参数|minimap2 preset (asm5/asm10/asm20)

    # SyRI参数|SyRI parameters
    syri_filter: bool = True  # 是否过滤SyRI结果|Whether to filter SyRI output

    # PlotSR数据过滤参数|PlotSR data filtering parameters
    min_sr_size: int = 10000  # 最小结构变异大小|Minimum structural variant size
    nosyn: bool = False  # 不绘制同源区域|Do not plot syntenic regions
    noinv: bool = False  # 不绘制倒位|Do not plot inversions
    notr: bool = False  # 不绘制易位|Do not plot translocations
    nodup: bool = False  # 不绘制重复|Do not plot duplications

    # PlotSR可视化参数|PlotSR visualization parameters
    font_size: int = 6  # 字体大小|Font size
    dpi: int = 300  # 图片DPI|Image DPI
    space_ratio: float = 0.7  # 同源染色体间距(0.1-0.75)|Space for homologous chromosomes
    vertical: bool = False  # 垂直排列染色体|Plot vertical chromosomes
    itx: bool = False  # 染色体间交互模式|Inter-chromosomal plotting mode

    # 输出格式|Output format
    output_format: str = "pdf"  # 输出格式|Output format (pdf/png/svg)

    # 流程控制参数|Pipeline control parameters
    skip_existing: bool = True  # 跳过已完成的步骤|Skip completed steps

    # 染色体过滤参数|Chromosome filtering parameters
    chromosomes: Optional[List[str]] = None  # 指定要显示的染色体（数字或名称）|Specific chromosomes to display (numbers or names)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 规范化路径|Normalize paths
        self.genomes = [os.path.abspath(g) for g in self.genomes]
        self.output_dir = os.path.abspath(self.output_dir)

        # 创建输出目录|Create output directories
        Path(self.output_dir).mkdir(parents=True, exist_ok=True)
        for subdir in ['alignment', 'syri', 'plotsr']:
            Path(self.output_dir, subdir).mkdir(parents=True, exist_ok=True)

        # 自动生成基因组名称|Auto-generate genome names
        if self.names is None:
            self.names = self._extract_names_from_paths()

        # 验证基因组数量|Validate genome count
        if len(self.genomes) < 2:
            raise ValueError(f"至少需要2个基因组|At least 2 genomes required: {len(self.genomes)}")

    def _extract_names_from_paths(self) -> List[str]:
        """
        从文件路径提取基因组名称|Extract genome names from file paths

        Returns:
            list: 基因组名称列表|List of genome names
        """
        import re
        names = []
        for genome_path in self.genomes:
            basename = os.path.basename(genome_path)
            # 移除扩展名|Remove extension
            name = re.sub(r'\.(fa|fasta|fna)(\.gz)?$', '', basename)
            names.append(name)
        return names

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查基因组文件|Check genome files
        for i, genome in enumerate(self.genomes):
            if not os.path.exists(genome):
                errors.append(f"基因组文件不存在|Genome file not found [{i}]: {genome}")
            if not os.path.isfile(genome):
                errors.append(f"基因组路径不是文件|Genome path is not a file [{i}]: {genome}")

        # 检查参数范围|Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Thread count must be positive: {self.threads}")

        if not (0.1 <= self.space_ratio <= 0.75):
            errors.append(f"染色体间距必须在0.1-0.75之间|Space ratio must be 0.1-0.75: {self.space_ratio}")

        if self.min_sr_size < 0:
            errors.append(f"最小结构变异大小不能为负数|Min SR size cannot be negative: {self.min_sr_size}")

        if self.output_format not in ['pdf', 'png', 'svg']:
            errors.append(f"输出格式不支持|Unsupported output format: {self.output_format}")

        # 检查基因组名称数量|Check genome names count
        if self.names and len(self.names) != len(self.genomes):
            errors.append(
                f"基因组名称数量与基因组数量不匹配|"
                f"Genome names count mismatch: {len(self.names)} != {len(self.genomes)}"
            )

        if errors:
            raise ValueError("\n".join(errors))

        return True
