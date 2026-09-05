"""NCBI datasets 下载配置模块|NCBI datasets download configuration module"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from biopytools.common.paths import expand_path, get_tool_path

# 合法筛选值|Valid filter values
VALID_ASSEMBLY_LEVELS = ('complete', 'chromosome', 'scaffold', 'contig')
ASSEMBLY_SOURCE_MAP = {'refseq': 'RefSeq', 'genbank': 'GenBank'}
# 清单 TSV 列|Manifest TSV columns
SUMMARY_FIELDS = ('accession', 'organism', 'assembly_level', 'assembly_status')
# 额外下载内容(默认只下 genome)|Extra include options (default: genome only)
INCLUDE_OPTIONS = {
    'include_gff3': 'gff3',
    'include_protein': 'protein',
    'include_cds': 'cds',
    'include_seq_report': 'seq-report',
}


@dataclass
class NCBIDatasetsConfig:
    """NCBI datasets 批量下载配置类|NCBI datasets batch download configuration class

    输入 taxon 编号, 下载该 taxon 下所有 assembly 的基因组序列(可加注释/蛋白等)
    |Input a taxon ID and download genome sequences of all assemblies under it
    """

    # 必需参数|Required parameters
    taxon: int

    # 输出设置|Output settings
    output_dir: str = './output'

    # 筛选参数|Filter parameters
    assembly_source: Optional[str] = None   # refseq / genbank
    assembly_level: Optional[str] = None    # 逗号分隔|comma-separated
    reference: bool = False                 # 只要参考基因组|reference genomes only
    annotated: bool = False                 # 只要有注释的|annotated only

    # 下载内容(默认只下 genome 序列)|Download content (genome only by default)
    include_gff3: bool = False
    include_protein: bool = False
    include_cds: bool = False
    include_seq_report: bool = False

    # 运行模式|Run mode
    dry_run: bool = False                   # 只查清单不下载|summary only
    organize: bool = True                   # 下载后整理到 02_organized|organize after download

    # 工具路径|Tool path
    datasets_path: Optional[str] = None     # 显式指定|explicit override

    log_level: str = 'INFO'

    # 内部属性|Internal attributes (computed in __post_init__)
    base_name: str = field(default='', init=False)
    output_path: Path = field(default=None, init=False)
    info_dir: Path = field(default=None, init=False)
    download_dir: Path = field(default=None, init=False)
    logs_dir: Path = field(default=None, init=False)
    organized_dir: Path = field(default=None, init=False)
    manifest_file: Path = field(default=None, init=False)
    zip_file: Path = field(default=None, init=False)
    dataset_dir: Path = field(default=None, init=False)
    log_file: Path = field(default=None, init=False)
    versions_file: Path = field(default=None, init=False)
    index_file: Path = field(default=None, init=False)
    resolved_datasets_path: str = field(default='', init=False)

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # taxon 归一化为正整数|Normalize taxon to positive int
        try:
            self.taxon = int(self.taxon)
        except (TypeError, ValueError):
            raise ValueError(
                f"taxon 必须是正整数|taxon must be a positive integer: {self.taxon}")

        # 路径展开(严禁硬编码绝对路径)|Expand paths (no hardcoded absolute paths)
        self.output_dir = expand_path(self.output_dir)
        self.output_path = Path(self.output_dir)

        # by-step 输出目录(§12.2)|by-step output layout
        self.info_dir = self.output_path / '00_pipeline_info'
        self.download_dir = self.output_path / '01_download'
        self.logs_dir = self.output_path / '99_logs'
        for d in (self.output_path, self.info_dir, self.download_dir, self.logs_dir):
            d.mkdir(parents=True, exist_ok=True)

        # 整理步骤目录(§12.2 条件性产物)|organized step dirs
        self.organized_dir = self.output_path / '02_organized'
        for d in (self.organized_dir,):
            d.mkdir(parents=True, exist_ok=True)

        self.base_name = str(self.taxon)
        self.manifest_file = self.info_dir / f'{self.base_name}_assemblies.tsv'
        self.zip_file = self.download_dir / f'{self.base_name}_genomes.zip'
        self.dataset_dir = self.download_dir / f'{self.base_name}.ncbi_dataset'
        self.log_file = self.logs_dir / f'{self.base_name}_ncbi_datasets.log'
        self.versions_file = self.info_dir / 'software_versions.yml'
        self.index_file = self.organized_dir / 'files.tsv'

        # 工具路径优先级: 显式参数 > 环境变量/配置文件 > 默认 ~/bin/datasets
        # |Tool path priority: explicit > env var / config file > default ~/bin/datasets
        if self.datasets_path:
            self.resolved_datasets_path = expand_path(self.datasets_path)
        else:
            self.resolved_datasets_path = get_tool_path(
                'datasets', '~/bin/datasets', 'DATASETS_PATH')

        # 筛选参数归一化|Normalize filters
        if self.assembly_source:
            self.assembly_source = self.assembly_source.strip().lower()
        if self.assembly_level:
            self.assembly_level = self.assembly_level.strip().lower()

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        if self.taxon <= 0:
            errors.append(
                f"taxon 必须是正整数|taxon must be a positive integer: {self.taxon}")
        if self.assembly_source not in (None, 'refseq', 'genbank'):
            errors.append(
                f"assembly-source 只能为 refseq 或 genbank|assembly-source must be "
                f"refseq or genbank: {self.assembly_source}")
        if self.assembly_level:
            for level in self._parsed_assembly_levels():
                if level not in VALID_ASSEMBLY_LEVELS:
                    errors.append(
                        f"无效的 assembly-level: {level}(可选|valid: "
                        f"complete, chromosome, scaffold, contig)")
        if errors:
            raise ValueError('\n'.join(errors))
        return True

    def _parsed_assembly_levels(self):
        """解析逗号分隔的 assembly-level|Parse comma-separated assembly levels"""
        return [x.strip() for x in self.assembly_level.split(',') if x.strip()]
