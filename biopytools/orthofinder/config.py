"""
OrthoFinder泛基因组分析配置管理模块|OrthoFinder Pangenome Analysis Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


@dataclass
class PangenomeConfig:
    """泛基因组分析配置类|Pangenome Analysis Configuration Class"""

    # 输入输出|Input/Output
    input_dir: str
    output_dir: str = './orthofinder_pangenome_output'
    project_name: Optional[str] = None

    # OrthoFinder基础参数|OrthoFinder basic parameters
    threads: int = 12
    search_program: str = 'blast'  # blast, diamond, diamond_ultra_sens, mmseqs
    mcl_inflation: float = 1.2
    sequence_type: str = 'protein'  # protein, dna

    # 泛基因组分类参数|Pangenome classification parameters
    softcore_missing_threshold: int = 1
    dispensable_missing_threshold: int = 1

    # 稀释分析参数|Rarefaction analysis parameters
    enable_rarefaction: bool = True
    rarefaction_iterations: int = 100

    # 单拷贝基因分析参数|Single copy gene analysis parameters
    enable_single_copy_analysis: bool = True
    extract_sequences: bool = True
    single_copy_output_format: str = 'both'  # 'by_orthogroup', 'by_genome', 'both'

    # 分析模式|Analysis mode
    basic_analysis_only: bool = True
    generate_trees: bool = False
    msa_program: str = 'mafft'
    tree_program: str = 'fasttree'

    # 断点续跑参数|Resume parameters
    resume_from_existing: bool = True
    skip_orthofinder: bool = False
    force_overwrite: bool = False

    # 可视化参数|Visualization parameters
    generate_plots: bool = True
    plot_format: str = 'png'
    figure_dpi: int = 300

    # 工具路径|Tool paths
    orthofinder_path: str = field(
        default_factory=lambda: get_tool_path('orthofinder', '/share/org/YZWL/yzwl_lixg/miniforge3/envs/orthofinder_v.3.1.5/bin/orthofinder', 'ORTHOFINDER_PATH')
    )

    # 内部属性|Internal attributes
    base_name: str = 'pangenome_analysis'

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.input_dir = expand_path(self.input_dir)
        self.output_dir = expand_path(self.output_dir)
        # orthofinder_path是命令名，不使用expand_path避免被转为cwd下的绝对路径
        if self.orthofinder_path and os.path.sep in self.orthofinder_path:
            self.orthofinder_path = expand_path(self.orthofinder_path)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        if not self.project_name:
            self.project_name = f"pangenome_{Path(self.input_dir).name}"

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.input_dir):
            errors.append(f"输入目录不存在|Input directory does not exist: {self.input_dir}")

        input_path = Path(self.input_dir)
        fasta_files = list(input_path.glob('*.fa*'))
        if not fasta_files:
            errors.append(f"输入目录中未找到FASTA文件|No FASTA files found in input directory: {self.input_dir}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        if not 0.1 <= self.mcl_inflation <= 10.0:
            errors.append(f"MCL inflation参数超出合理范围|MCL inflation parameter out of range: {self.mcl_inflation}")

        if self.softcore_missing_threshold < 0:
            errors.append(f"Softcore缺失阈值必须为非负整数|Softcore missing threshold must be non-negative: {self.softcore_missing_threshold}")

        if self.dispensable_missing_threshold < 0:
            errors.append(f"Dispensable缺失阈值必须为非负整数|Dispensable missing threshold must be non-negative: {self.dispensable_missing_threshold}")

        if self.rarefaction_iterations <= 0:
            errors.append(f"稀释分析迭代次数必须为正整数|Rarefaction iterations must be positive: {self.rarefaction_iterations}")

        valid_search_programs = ['blast', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl']
        if self.search_program not in valid_search_programs:
            errors.append(f"不支持的搜索程序|Unsupported search program: {self.search_program}")

        if errors:
            raise ValueError("\n".join(errors))

        return True
