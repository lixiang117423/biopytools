"""
PICRUSt2配置管理模块|PICRUSt2 Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List

from ..common.paths import get_tool_path, expand_path

VALID_TRAITS_SPLIT = {'EC', 'KO', 'GO', 'PFAM', 'BIGG', 'CAZY', 'GENE_NAMES'}
VALID_TRAITS_SINGLE = {'COG', 'EC', 'KO', 'PFAM', 'TIGRFAM'}
VALID_PLACEMENT_TOOLS = {'epa-ng', 'sepp'}
VALID_HSP_METHODS = {'mp', 'emp_prob', 'pic', 'scp', 'subtree_average'}
VALID_PIPELINES = {'auto', 'split', 'single'}


@dataclass
class Picrust2Config:
    """PICRUSt2配置类|PICRUSt2 Configuration Class"""

    # ===== 必需参数|Required parameters =====
    study_fasta: str
    input_table: str
    output_dir: str

    # ===== 流程参数|Pipeline parameters =====
    threads: int = 12
    max_nsti: float = 2.0
    stratified: bool = False
    in_traits: str = "EC,KO"
    placement_tool: str = "epa-ng"
    hsp_method: str = "mp"
    edge_exponent: float = 0.5
    min_align: float = 0.8
    min_reads: int = 1
    min_samples: int = 1

    # ===== 流程控制|Pipeline control =====
    pipeline: str = "auto"
    skip_pathways: bool = False
    coverage: bool = False
    skip_norm: bool = False
    remove_intermediate: bool = False
    verbose: bool = False
    per_sequence_contrib: bool = False
    skip_minpath: bool = False
    no_gap_fill: bool = False

    # ===== 工具路径|Tool paths =====
    picrust2_path: str = field(
        default_factory=lambda: get_tool_path(
            'picrust2',
            '~/miniforge3/envs/picrust_v.2.6.3/bin/picrust2_pipeline.py',
            'PICRUST2_PATH'
        )
    )
    picrust2_single_path: str = field(
        default_factory=lambda: get_tool_path(
            'picrust2_single',
            '~/miniforge3/envs/picrust_v.2.6.3/bin/picrust2_pipeline_singleRef.py',
            'PICRUST2_SINGLE_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.study_fasta = expand_path(self.study_fasta)
        self.input_table = expand_path(self.input_table)
        self.output_dir = expand_path(self.output_dir)
        self.picrust2_path = expand_path(self.picrust2_path)
        self.picrust2_single_path = expand_path(self.picrust2_single_path)

        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 标准子目录|Standard subdirectories
        self.info_dir = os.path.join(self.output_dir, "00_pipeline_info")
        self.placement_dir = os.path.join(self.output_dir, "01_placement")
        self.hsp_dir = os.path.join(self.output_dir, "02_hsp")
        self.metagenome_dir = os.path.join(self.output_dir, "03_metagenome")
        self.pathway_dir = os.path.join(self.output_dir, "04_pathway")
        self.log_dir = os.path.join(self.output_dir, "99_logs")

        for dir_path in [self.info_dir, self.placement_dir,
                         self.hsp_dir, self.metagenome_dir,
                         self.pathway_dir, self.log_dir]:
            Path(dir_path).mkdir(parents=True, exist_ok=True)

        # PICRUSt2 工作目录（wrapper在此运行原始pipeline，之后重组输出）|PICRUSt2 work dir
        self.work_dir = os.path.join(self.output_dir, "work")

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        if not os.path.exists(self.study_fasta):
            errors.append(f"代表序列文件不存在|Study FASTA not found: {self.study_fasta}")

        if not os.path.exists(self.input_table):
            errors.append(f"特征表文件不存在|Input table not found: {self.input_table}")

        if self.threads <= 0:
            errors.append(f"线程数必须为正整数|Thread count must be positive: {self.threads}")

        if self.max_nsti < 0:
            errors.append(f"NSTI阈值不能为负|Max NSTI must be non-negative: {self.max_nsti}")

        if self.pipeline not in VALID_PIPELINES:
            errors.append(f"无效的流程类型|Invalid pipeline type: {self.pipeline} (可选|options: {VALID_PIPELINES})")

        if self.placement_tool not in VALID_PLACEMENT_TOOLS:
            errors.append(f"无效的放置工具|Invalid placement tool: {self.placement_tool}")

        if self.hsp_method not in VALID_HSP_METHODS:
            errors.append(f"无效的HSP方法|Invalid HSP method: {self.hsp_method}")

        # 验证功能数据库参数|Validate trait parameters
        requested_traits = [t.strip().upper() for t in self.in_traits.split(',')]
        valid_traits = self._get_valid_traits()
        for trait in requested_traits:
            if trait not in valid_traits:
                errors.append(f"不支持的功能数据库|Unsupported trait: {trait} (可选|options: {valid_traits})")

        # 验证流程脚本可用|Validate pipeline script exists
        pipeline_script = self._resolve_pipeline_script()
        if not pipeline_script:
            errors.append("PICRUSt2 pipeline脚本未找到|PICRUSt2 pipeline script not found")
        elif not os.path.exists(pipeline_script):
            errors.append(f"PICRUSt2脚本不存在|PICRUSt2 script not found: {pipeline_script}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def _get_valid_traits(self) -> set:
        """获取当前流程支持的功能数据库|Get valid traits for current pipeline"""
        if self.pipeline == 'single':
            return VALID_TRAITS_SINGLE
        return VALID_TRAITS_SPLIT

    def _resolve_pipeline_script(self) -> Optional[str]:
        """解析使用的pipeline脚本路径|Resolve which pipeline script to use"""
        if self.pipeline == 'split':
            return self.picrust2_path
        elif self.pipeline == 'single':
            return self.picrust2_single_path
        else:
            # auto: 优先双域|auto: prefer split (dual-domain)
            if os.path.exists(self.picrust2_path):
                return self.picrust2_path
            elif os.path.exists(self.picrust2_single_path):
                return self.picrust2_single_path
            return None
