"""vcf2splitstree 配置|vcf2splitstree configuration"""
from dataclasses import dataclass, field
from pathlib import Path

from ..common.paths import expand_path


@dataclass
class Vcf2SplitstreeConfig:
    """VCF→SplitsTree6 距离矩阵转换配置|VCF to SplitsTree6 distance matrix config

    只做一件事:读 VCF 变异文件,输出 SplitsTree6 GUI 可直接打开的
    p-distance 距离矩阵 CSV(SplitsTree6 打开 CSV 会自动跑 NeighborNet)。
    |Convert VCF variants into a p-distance CSV that SplitsTree6 GUI can open
    directly (it auto-runs NeighborNet on CSV distances).
    """

    # 必需|Required
    input: str                    # VCF 变异文件(.vcf / .vcf.gz)|VCF file

    # 输出|Output
    output_dir: str = './vcf2splitstree_output'

    log_level: str = 'INFO'

    # 内部属性|Internal attributes
    output_csv: Path = field(default=None, init=False)
    log_file: Path = field(default=None, init=False)

    def __post_init__(self):
        """展开路径、建目录|Expand paths, make dirs"""
        self.input = expand_path(self.input)
        self.output_dir = expand_path(self.output_dir)
        out = Path(self.output_dir)
        logs_dir = out / '99_logs'
        logs_dir.mkdir(parents=True, exist_ok=True)
        # 输出文件名与输入同名(去 .vcf/.gz)|same stem as input (strip .vcf/.gz)
        name = Path(self.input).name
        for suffix in ('.vcf.gz', '.vcf', '.VCF'):
            if name.lower().endswith(suffix.lower()):
                name = name[: -len(suffix)]
                break
        self.output_csv = out / f'{name}.distances.csv'
        self.log_file = logs_dir / 'vcf2splitstree.log'

    def validate(self):
        """校验配置|Validate configuration"""
        errors = []
        if not Path(self.input).exists():
            errors.append(f"输入不存在|Input not found: {self.input}")
        low = self.input.lower()
        if not (low.endswith('.vcf') or low.endswith('.vcf.gz')):
            errors.append(f"输入必须是 VCF 文件(.vcf/.vcf.gz)|Input must be a "
                          f"VCF file: {self.input}")
        if errors:
            raise ValueError('\n'.join(errors))
        return True
