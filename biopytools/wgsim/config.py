"""
Wgsim配置类|Wgsim Configuration Class
"""

import os
from common.paths import get_tool_path, expand_path


class WgsimConfig:
    """Wgsim配置类|Wgsim Configuration Class"""

    def __init__(
        self,
        input_path: str,
        output_dir: str,
        num_reads: int = 50000000,
        read_length: int = 150,
        seed: int = 0,
        error_rate: float = 0.020,
        mutation_rate: float = 0.001,
        outer_distance: int = 500,
        inner_distance: int = 0,
        wgsim_path: str = None,
        extensions: list = None,
    ):
        """
        初始化配置|Initialize configuration

        Args:
            input_path: 输入文件或目录路径|Input file or directory path
            output_dir: 输出目录|Output directory
            num_reads: 模拟reads数量|Number of reads to simulate
            read_length: reads长度|Read length
            seed: 随机种子|Random seed
            error_rate: 测序错误率|Sequencing error rate
            mutation_rate: 突变率|Mutation rate
            outer_distance: 外部距离|Outer distance
            inner_distance: 内部距离|Inner distance
            wgsim_path: wgsim可执行文件路径|wgsim binary path
            extensions: 输入文件扩展名列表|Input file extension list
        """
        self.input_path = input_path
        self.output_dir = output_dir
        self.num_reads = num_reads
        self.read_length = read_length
        self.seed = seed
        self.error_rate = error_rate
        self.mutation_rate = mutation_rate
        self.outer_distance = outer_distance
        self.inner_distance = inner_distance
        self.wgsim_path = wgsim_path
        self.extensions = extensions or ['.fna', '.fa', '.fasta']

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        self.input_path = expand_path(self.input_path)
        self.output_dir = expand_path(self.output_dir)

        if not self.input_path:
            errors.append("输入路径不能为空|Input path cannot be empty")
        elif not os.path.exists(self.input_path):
            errors.append(f"输入路径不存在|Input path not found: {self.input_path}")

        if self.wgsim_path:
            self.wgsim_path = expand_path(self.wgsim_path)
        else:
            self.wgsim_path = get_tool_path(
                'wgsim',
                '~/miniforge3/envs/GATK_v.4.6.2.0/bin/wgsim',
                'WGSIM_PATH'
            )

        if not os.path.exists(self.wgsim_path):
            errors.append(f"wgsim不存在|wgsim not found: {self.wgsim_path}")

        if self.num_reads <= 0:
            errors.append(f"reads数量必须为正数|Number of reads must be positive: {self.num_reads}")

        if self.read_length <= 0:
            errors.append(f"reads长度必须为正数|Read length must be positive: {self.read_length}")

        if errors:
            raise ValueError("\n".join(errors))

        os.makedirs(self.output_dir, exist_ok=True)
        return True

    def __repr__(self):
        return (
            f"WgsimConfig(\n"
            f"  input_path={self.input_path!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  num_reads={self.num_reads!r},\n"
            f"  read_length={self.read_length!r},\n"
            f"  wgsim_path={self.wgsim_path!r}\n"
            f")"
        )
