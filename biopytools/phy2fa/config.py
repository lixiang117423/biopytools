"""phy2fa 配置|phy2fa configuration"""
from dataclasses import dataclass, field
from pathlib import Path

from ..common.paths import expand_path


@dataclass
class Phy2FaConfig:
    """Phylip→FASTA 转换配置|Phylip to FASTA conversion configuration"""

    # 必需|Required
    input: str                    # Phylip 文件(.phy/.phylip,可 .gz)|Phylip file

    # 输出|Output
    output_dir: str = './phy2fa_output'

    # 可选|Optional
    line_width: int = 60          # FASTA 序列换行宽度(0=单行)|FASTA line wrap width
    log_level: str = 'INFO'

    # 内部属性|Internal attributes
    output_fasta: Path = field(default=None, init=False)
    log_file: Path = field(default=None, init=False)

    def __post_init__(self):
        """展开路径、建目录、推导输出文件名|Expand paths, make dirs, derive output"""
        self.input = expand_path(self.input)
        self.output_dir = expand_path(self.output_dir)
        out = Path(self.output_dir)
        logs_dir = out / '99_logs'
        logs_dir.mkdir(parents=True, exist_ok=True)

        name = Path(self.input).name
        for suffix in ('.phy.gz', '.phylip.gz', '.phy', '.phylip'):
            if name.lower().endswith(suffix.lower()):
                name = name[: -len(suffix)]
                break
        self.output_fasta = out / f'{name}.fa'
        self.log_file = logs_dir / 'phy2fa.log'

    def validate(self):
        """校验配置(错误收集后一次抛)|Validate; raise ValueError listing all errors"""
        errors = []
        if not Path(self.input).exists():
            errors.append(f"输入不存在|Input not found: {self.input}")
        low = self.input.lower()
        if not (low.endswith('.phy') or low.endswith('.phylip') or
                low.endswith('.phy.gz') or low.endswith('.phylip.gz')):
            errors.append(f"输入必须是 Phylip 文件(.phy/.phylip/.gz)|Input must be "
                          f"a Phylip file: {self.input}")
        if self.line_width < 0:
            errors.append(f"换行宽度不能为负|Line width must be >= 0: {self.line_width}")
        if errors:
            raise ValueError('\n'.join(errors))
        return True
