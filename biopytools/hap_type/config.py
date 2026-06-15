"""
单倍型分析配置管理模块|Haplotype Analysis Configuration Management Module
"""

import re
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Tuple


@dataclass
class HapTypeConfig:
    """单倍型分析配置类|Haplotype Analysis Configuration Class"""

    # 必需参数|Required parameters
    vcf_file: str = ""
    region: str = ""
    output: str = ""

    # geneHapR兼容参数|geneHapR-compatible parameters
    hetero_remove: bool = False
    na_drop: bool = True
    hap_prefix: str = "H"
    pad: int = 3

    # 区间格式正则|Interval pattern regex
    _interval_re = re.compile(r'^([^:]+):\s*(\d+)\s*[-]\s*(\d+)$')

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        self.vcf_path = Path(self.vcf_file) if self.vcf_file else None

        if self.region:
            match = self._interval_re.match(self.region.strip())
            if match:
                self.chrom = match.group(1)
                self.start = int(match.group(2))
                self.end = int(match.group(3))
                self.bed_path = None
            else:
                self.chrom = None
                self.start = None
                self.end = None
                self.bed_path = Path(self.region)
        else:
            self.chrom = None
            self.start = None
            self.end = None
            self.bed_path = None

        if self.output:
            self.output_path = Path(self.output)
        else:
            self.output_path = None

        self._bed_records: List[Tuple[str, int, int, str]] = field(init=False, default_factory=list)

    @property
    def is_bed_mode(self) -> bool:
        """是否为BED文件模式|Whether BED file mode is active"""
        return self.bed_path is not None

    def parse_bed(self) -> List[Tuple[str, int, int, str]]:
        """解析BED文件，自动识别BED3/4/6格式|Parse BED file, auto-detect BED3/4/6 format

        Returns:
            list: [(chrom, start, end, name), ...]  start/end为0-based转为1-based
        """
        if not self.bed_path:
            return []

        if not self.bed_path.exists():
            raise FileNotFoundError(f"BED文件不存在|BED file not found: {self.bed_path}")

        records = []
        with open(self.bed_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith("#"):
                    continue
                cols = line.split("\t")
                if len(cols) < 3:
                    continue

                chrom = cols[0]
                start = int(cols[1]) + 1  # BED 0-based → 1-based
                end = int(cols[2])         # BED end是开区间，直接用

                if len(cols) >= 4 and cols[3]:
                    name = cols[3]
                else:
                    name = ""

                records.append((chrom, start, end, name))

        self._bed_records = records
        return records

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []
        if not self.vcf_file:
            errors.append("VCF文件不能为空|VCF file cannot be empty")
        elif not self.vcf_path.exists():
            errors.append(f"VCF文件不存在|VCF file does not exist: {self.vcf_file}")
        if not self.region:
            errors.append("必须指定 -r/--region|Must specify -r/--region")
        if self.is_bed_mode and not self.bed_path.exists():
            errors.append(f"BED文件不存在|BED file not found: {self.region}")
        if self.region and not self.is_bed_mode and not self._interval_re.match(self.region.strip()):
            errors.append(f"区间格式不正确|Incorrect interval format: {self.region}")
        if errors:
            raise ValueError("\n".join(errors))
        return True

    @staticmethod
    def _region_prefix(chrom: str, start: int, end: int, name: str = "") -> str:
        """生成区间前缀: Chr_Start_End[_name]|Generate region prefix"""
        prefix = f"{chrom}_{start}_{end}"
        if name:
            prefix = f"{prefix}_{name}"
        return prefix

    def get_output_dir(self) -> Path:
        """获取输出目录|Get output directory"""
        if self.output_path:
            return self.output_path
        return Path.cwd()

    def get_output_files(self, chrom: str, start: int, end: int, name: str = "") -> dict:
        """获取输出文件路径|Get output file paths

        Args:
            chrom: 染色体|Chromosome
            start: 起始位置(1-based)|Start position
            end: 终止位置(1-based)|End position
            name: 区间名称(可选)|Interval name (optional)
        """
        region_name = self._region_prefix(chrom, start, end, name)
        output_dir = self.get_output_dir()
        base = str(output_dir / region_name)
        return {
            'hap_result': f"{base}.hapResult.txt",
            'hap_result_xlsx': f"{base}.hapResult.xlsx",
            'hap_summary': f"{base}.hapSummary.txt",
            'hap_summary_xlsx': f"{base}.hapSummary.xlsx",
            'sample_hap': f"{base}.sampleHap.txt",
            'sample_hap_xlsx': f"{base}.sampleHap.xlsx",
            'log': f"{base}.log",
        }
