"""check_reads 配置|check_reads configuration"""
from dataclasses import dataclass, field
from pathlib import Path
from typing import List

from ..common.paths import expand_path


@dataclass
class CheckReadsConfig:
    """fastq 完整性检查配置|FASTQ integrity check configuration

    校验输入目录的 fastq:gz 完整性、0 字节、R1-R2 配对完整性
    |Check fastq gz integrity, empty files, and R1-R2 pairing completeness
    """

    # 必需|Required
    input_dir: str               # fastq 目录(逗号分隔多个)|dir(s), comma-separated

    # 输出|Output
    output_dir: str = "./check_reads_output"

    # 运行|Runtime
    threads: int = 12
    log_level: str = "INFO"

    # 内部属性|Internal attributes
    report_file: Path = field(default=None, init=False)
    log_file: Path = field(default=None, init=False)

    def __post_init__(self):
        """展开路径、建目录|Expand paths, make dirs"""
        self.output_dir = expand_path(self.output_dir)
        out = Path(self.output_dir)
        logs_dir = out / "99_logs"
        logs_dir.mkdir(parents=True, exist_ok=True)
        self.report_file = out / "check_reads_report.txt"
        self.log_file = logs_dir / "check_reads.log"

    @property
    def input_dirs(self) -> List[str]:
        """展开逗号分隔的输入目录|Expand comma-separated input dirs"""
        return [expand_path(p.strip())
                for p in self.input_dir.split(",") if p.strip()]

    def validate(self):
        """校验配置|Validate configuration"""
        errors = []
        if not self.input_dirs:
            errors.append("输入目录不能为空|Input dir is required")
        for d in self.input_dirs:
            if not Path(d).is_dir():
                errors.append(f"输入目录不存在|Input dir not found: {d}")
        if self.threads <= 0:
            errors.append(f"线程数必须为正数|Threads must be positive: {self.threads}")
        if errors:
            raise ValueError("\n".join(errors))
        return True
