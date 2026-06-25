"""
HiCanu组装配置管理模块|HiCanu Assembly Configuration Management Module
"""

import os
from ..common.paths import expand_path
from .utils import build_conda_command
import glob
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional, List


@dataclass
class HiCanuConfig:
    """HiCanu组装配置类|HiCanu Assembly Configuration Class"""

    # 必需参数|Required parameters
    reads_file: str
    genome_size: str
    prefix: str

    # 路径配置|Path configuration
    canu_path: str = '~/miniforge3/envs/canu_v.2.3/bin/canu'
    base_dir: str = './hicanu_output'  # 基础输出目录|Base output directory

    # 组装参数|Assembly parameters
    min_read_length: int = 1000
    min_overlap_length: int = 500
    corrected_error_rate: Optional[float] = None
    raw_error_rate: Optional[float] = None
    max_input_coverage: Optional[int] = None

    # 计算资源|Computing resources
    threads: int = 12
    memory: str = '100G'
    use_grid: bool = False
    grid_options: Optional[str] = None

    # 执行控制|Execution control
    stage: str = 'assemble'  # haplotype, correct, trim, assemble, trim-assemble
    dry_run: bool = False
    keep_intermediate: bool = False
    resume: bool = True  # 断点续传（默认启用）|Resume from previous run (enabled by default)

    # 内部属性|Internal attributes
    work_dir: str = None
    raw_dir: str = None
    fasta_dir: str = None
    log_dir: str = None
    stat_dir: str = None
    output_dir: str = None  # canu的输出目录（raw_dir）|Canu output directory (raw_dir)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 关键：展开所有包含~的路径|CRITICAL: Expand all paths with ~
        self.canu_path = expand_path(self.canu_path)

        # 设置基础目录|Setup base directory
        self.base_dir = os.path.normpath(os.path.abspath(self.base_dir))
        self.work_dir = os.path.join(self.base_dir, self.prefix)

        # 定义子目录|Define subdirectories
        self.raw_dir = os.path.join(self.work_dir, "01.raw_output")
        self.fasta_dir = os.path.join(self.work_dir, "02.fasta")
        self.log_dir = os.path.join(self.work_dir, "03.logs")
        self.stat_dir = os.path.join(self.work_dir, "04.statistics")

        # canu的输出目录设置为raw_dir|Set canu output directory to raw_dir
        self.output_dir = self.raw_dir
        self.output_path = Path(self.output_dir)

        # 创建目录结构|Create directory structure
        if not self.dry_run:
            for dir_path in [self.work_dir, self.raw_dir, self.fasta_dir, self.log_dir, self.stat_dir]:
                Path(dir_path).mkdir(parents=True, exist_ok=True)

        # 标准化路径|Normalize paths
        self.reads_file = os.path.normpath(os.path.abspath(self.reads_file))

        # 验证 genome_size 格式|Validate genome_size format
        self._validate_genome_size()

        # 验证 stage 参数|Validate stage parameter
        valid_stages = ['haplotype', 'correct', 'trim', 'assemble', 'trim-assemble']
        if self.stage not in valid_stages:
            raise ValueError(
                f"无效的stage参数|Invalid stage parameter: {self.stage}. "
                f"必须是|Must be one of: {', '.join(valid_stages)}"
            )

    def _validate_genome_size(self):
        """验证基因组大小格式|Validate genome size format"""
        genome_size_lower = self.genome_size.lower()
        if not any(genome_size_lower.endswith(suffix) for suffix in ['g', 'm', 'k']):
            try:
                # 如果没有后缀，尝试转换为数字
                float(self.genome_size)
            except ValueError:
                raise ValueError(
                    f"无效的genomeSize格式|Invalid genomeSize format: {self.genome_size}. "
                    f"应该是数字加上后缀，如|Should be a number with suffix, e.g., '120m', '1g', '5000k'"
                )

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 检查必需文件|Check required files
        if not os.path.exists(self.reads_file):
            errors.append(f"reads文件不存在|reads file does not exist: {self.reads_file}")

        # 检查Canu路径|Check Canu path
        if not os.path.exists(self.canu_path):
            errors.append(f"Canu路径不存在|Canu path does not exist: {self.canu_path}")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_canu_command(self) -> List[str]:
        """构建Canu命令|Build Canu command"""
        # 构建参数列表（不包含canu_path本身）|Build argument list (excluding canu_path itself)
        args = [
            '-d', self.output_dir,
            '-p', self.prefix,
            f'genomeSize={self.genome_size}',
            f'-{self.stage}',
        ]

        # 添加reads文件类型|Add reads file type
        args.extend(['-pacbio-hifi', self.reads_file])

        # 添加可选参数|Add optional parameters
        if self.min_read_length != 1000:
            args.append(f'minReadLength={self.min_read_length}')

        if self.min_overlap_length != 500:
            args.append(f'minOverlapLength={self.min_overlap_length}')

        if self.corrected_error_rate is not None:
            args.append(f'correctedErrorRate={self.corrected_error_rate}')

        if self.raw_error_rate is not None:
            args.append(f'rawErrorRate={self.raw_error_rate}')

        if self.max_input_coverage is not None:
            args.append(f'maxInputCoverage={self.max_input_coverage}')

        # Grid选项|Grid options
        if self.use_grid:
            args.append('useGrid=true')
            if self.grid_options:
                args.append(f'gridOptions={self.grid_options}')
        else:
            args.append('useGrid=false')

        # 关键：使用build_conda_command包装（符合开发规范第13章）|CRITICAL: Use build_conda_command wrapper (per dev guide section 13)
        # 传递完整路径给build_conda_command，它会自动检测conda环境并包装|Pass full path to build_conda_command, it will auto-detect conda env and wrap
        return build_conda_command(self.canu_path, args)

    def get_output_files(self) -> dict:
        """获取预期的输出文件|Get expected output files"""
        return {
            'contigs': os.path.join(self.output_dir, f'{self.prefix}.contigs.fasta'),
            'unitigs': os.path.join(self.output_dir, f'{self.prefix}.unitigs.fasta'),
            'report': os.path.join(self.output_dir, f'{self.prefix}.report'),
            'log': os.path.join(self.output_dir, f'{self.prefix}.log'),
        }

    def get_completed_steps(self) -> dict:
        """
        检查各步骤完成状态|Check completion status of each step

        Returns:
            dict: 各步骤完成状态|Completion status of each step
        """
        steps = {
            'assembly': False,      # Canu组装完成|Canu assembly completed
            'fasta_copied': False,  # FASTA文件已复制|FASTA files copied
            'reads_mapped': False,  # Contig-reads映射已生成|Contig-reads mapping generated
        }

        output_files = self.get_output_files()

        # 检查组装是否完成|Check if assembly is completed
        if os.path.exists(output_files['contigs']):
            file_size = os.path.getsize(output_files['contigs'])
            if file_size > 0:
                steps['assembly'] = True

        # 检查FASTA是否已复制|Check if FASTA files are copied
        fasta_files = glob.glob(os.path.join(self.fasta_dir, "*.fasta"))
        if fasta_files:
            steps['fasta_copied'] = True

        # 检查contig-reads映射是否已生成|Check if contig-reads mapping is generated
        mapping_file = os.path.join(self.fasta_dir, f'{self.prefix}.contig_reads.tsv')
        if os.path.exists(mapping_file):
            file_size = os.path.getsize(mapping_file)
            if file_size > 0:
                steps['reads_mapped'] = True

        return steps
