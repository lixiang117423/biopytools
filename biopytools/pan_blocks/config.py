"""泛基因组Block构建配置模块|Pan-Genome Block Construction Configuration Module"""

from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List, Dict
import os

from ..common.paths import get_tool_path, expand_path


@dataclass
class PanBlocksConfig:
    """泛基因组Block构建配置类|Pan-Genome Block Construction Configuration Class"""

    genome_list: str
    output_dir: str = "./pan_blocks_output"

    threads: int = 12
    parallel_alignments: int = 4
    min_alignment_length: int = 10000

    genome_order_file: Optional[str] = None
    chromosome: Optional[str] = None

    nucmer_path: str = field(default_factory=lambda: get_tool_path(
        'nucmer', '~/miniforge3/envs/pan-blocks/bin/nucmer', 'NUCMER_PATH'
    ))
    delta_filter_path: str = field(default_factory=lambda: get_tool_path(
        'delta-filter', '~/miniforge3/envs/pan-blocks/bin/delta-filter', 'DELTA_FILTER_PATH'
    ))
    show_coords_path: str = field(default_factory=lambda: get_tool_path(
        'show-coords', '~/miniforge3/envs/pan-blocks/bin/show-coords', 'SHOW_COORDS_PATH'
    ))
    bedtools_path: str = field(default_factory=lambda: get_tool_path(
        'bedtools', '~/.local/bin/bedtools', 'BEDTOOLS_PATH'
    ))
    minimap2_path: str = field(default_factory=lambda: get_tool_path(
        'minimap2', '~/miniforge3/envs/pan-blocks/bin/minimap2', 'MINIMAP2_PATH'
    ))

    plot_width: int = 20
    plot_height: int = 10
    plot_format: str = "svg"

    def __post_init__(self):
        self.output_dir = expand_path(self.output_dir)
        self.genome_list = expand_path(self.genome_list)

        self.nucmer_path = expand_path(self.nucmer_path)
        self.delta_filter_path = expand_path(self.delta_filter_path)
        self.show_coords_path = expand_path(self.show_coords_path)
        self.bedtools_path = expand_path(self.bedtools_path)
        self.minimap2_path = expand_path(self.minimap2_path)

        if self.genome_order_file:
            self.genome_order_file = expand_path(self.genome_order_file)

        self.output_path = Path(self.output_dir)
        self.coords_dir = self.output_path / "01_coords"
        self.blocks_dir = self.output_path / "02_blocks"
        self.plots_dir = self.output_path / "03_plots"
        self.logs_dir = self.output_path / "99_logs"
        self.info_dir = self.output_path / "00_pipeline_info"

        for d in [self.coords_dir, self.blocks_dir, self.plots_dir, self.logs_dir, self.info_dir]:
            d.mkdir(parents=True, exist_ok=True)

        self.genomes: List[Dict[str, str]] = self._parse_genome_list()
        self.genome_order_list: List[str] = self._parse_genome_order()
        self.chrlen_data: Dict[str, Dict[str, int]] = self._compute_chrlen()

    def _parse_genome_list(self) -> List[Dict[str, str]]:
        """解析基因组列表文件|Parse genome list file"""
        genomes = []
        with open(self.genome_list) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t')
                if len(parts) < 2:
                    continue
                name, path = parts[0], parts[1]
                path = expand_path(path)
                genomes.append({'name': name, 'path': path})
        return genomes

    def _parse_genome_order(self) -> List[str]:
        """解析基因组优先级顺序|Parse genome priority order"""
        if self.genome_order_file:
            order = []
            with open(self.genome_order_file) as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        order.append(line)
            return order
        return [g['name'] for g in self.genomes]

    def _compute_chrlen(self) -> Dict[str, Dict[str, int]]:
        """计算各基因组染色体长度|Compute chromosome lengths for each genome"""
        from .utils import parse_fasta_lengths
        data = {}
        for genome in self.genomes:
            name = genome['name']
            path = genome['path']
            if os.path.exists(path):
                lengths = parse_fasta_lengths(path)
                data[name] = lengths
            else:
                raise FileNotFoundError(f"基因组文件不存在|Genome file not found: {path}")
        return data

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.genome_list):
            errors.append(f"基因组列表文件不存在|Genome list file not found: {self.genome_list}")

        if len(self.genomes) < 2:
            errors.append(f"需要至少2个基因组|Need at least 2 genomes, got {len(self.genomes)}")

        for genome in self.genomes:
            if not os.path.exists(genome['path']):
                errors.append(f"基因组文件不存在|Genome file not found: {genome['path']}")

        for name, path in [('nucmer', self.nucmer_path), ('delta-filter', self.delta_filter_path),
                           ('show-coords', self.show_coords_path), ('bedtools', self.bedtools_path)]:
            if not os.path.isfile(path):
                errors.append(f"{name} 不存在|{name} not found: {path}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")
        if self.parallel_alignments <= 0:
            errors.append(f"并行比对数必须为正整数|Parallel alignments must be positive: {self.parallel_alignments}")

        if errors:
            raise ValueError("\n".join(errors))
        return True

    def get_target_chromosomes(self) -> List[str]:
        """获取目标染色体列表|Get target chromosome list"""
        if self.chromosome:
            return [self.chromosome]
        first_genome = self.genome_order_list[0]
        return sorted(self.chrlen_data[first_genome].keys())
