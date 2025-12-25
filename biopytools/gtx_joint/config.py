"""
GTX Joint Calling配置类 | GTX Joint Calling Configuration Class
"""

import os
import re


class GTXJointConfig:
    """GTX Joint Calling配置类 | GTX Joint Calling Configuration Class"""

    def __init__(
        self,
        gtx_exec: str = None,
        reference: str = None,
        gvcf_dir: str = None,
        output_dir: str = None,
        threads: int = 88,
        tmp_dir: str = "./tmp",
        script_name: str = "run_gtx_joint.sh",
        faketime: str = "2020-10-20 00:00:00",
        chr_pattern: str = None,
        window_size: int = None
    ):
        """
        初始化配置 | Initialize configuration

        Args:
            gtx_exec: GTX可执行文件路径 | GTX executable path
            reference: 参考基因组文件路径 | Reference genome file path
            gvcf_dir: GVCF文件所在目录 | GVCF files directory
            output_dir: 输出结果目录 | Output directory
            threads: 线程数 | Number of threads
            tmp_dir: 临时目录 | Temporary directory
            script_name: 输出脚本文件名 | Output script filename
            faketime: faketime时间 | faketime string
            chr_pattern: 染色体过滤正则 | Chromosome filter pattern
            window_size: 区间大小(bp) | Window size in bp
        """
        self.gtx_exec = gtx_exec or "/share/org/YZWL/yzwl_lixg/software/gtx/bin/gtx"
        self.reference = reference
        self.gvcf_dir = gvcf_dir
        self.output_dir = output_dir
        self.threads = threads
        self.tmp_dir = tmp_dir
        self.script_name = script_name
        self.faketime = faketime
        self.chr_pattern = chr_pattern
        self.window_size = window_size

        # 运行时设置的变量 | Runtime variables
        self.gvcf_files = []
        self.chromosomes = []
        self.use_faketime = False

    def validate(self):
        """
        验证配置参数 | Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时 | When configuration parameters are invalid
        """
        # 验证GTX可执行文件
        if not os.path.isfile(self.gtx_exec):
            raise ValueError(f"GTX可执行文件不存在 | GTX executable not found: {self.gtx_exec}")
        if not os.access(self.gtx_exec, os.X_OK):
            raise ValueError(f"GTX文件不可执行 | GTX file is not executable: {self.gtx_exec}")

        # 验证参考基因组
        if not self.reference:
            raise ValueError("参考基因组文件路径不能为空 | Reference genome path cannot be empty")
        if not os.path.isfile(self.reference):
            raise ValueError(f"参考基因组文件不存在 | Reference genome not found: {self.reference}")
        # 检查索引文件
        if not os.path.isfile(f"{self.reference}.fai"):
            raise ValueError(
                f"参考基因组索引不存在 | Reference genome index not found: {self.reference}.fai\n"
                f"请运行 | Please run: samtools faidx {self.reference}"
            )

        # 验证GVCF目录
        if not self.gvcf_dir:
            raise ValueError("GVCF目录路径不能为空 | GVCF directory path cannot be empty")
        if not os.path.isdir(self.gvcf_dir):
            raise ValueError(f"GVCF目录不存在 | GVCF directory not found: {self.gvcf_dir}")

        # 验证输出目录
        if not self.output_dir:
            raise ValueError("输出目录路径不能为空 | Output directory path cannot be empty")

        # 验证染色体模式
        if self.chr_pattern:
            try:
                re.compile(self.chr_pattern)
            except re.error as e:
                raise ValueError(f"染色体过滤正则表达式无效 | Invalid chromosome pattern: {e}")

        # 验证窗口大小
        if self.window_size is not None and self.window_size <= 0:
            raise ValueError(f"窗口大小必须大于0 | Window size must be positive: {self.window_size}")

    def __repr__(self):
        """配置的字符串表示 | String representation of configuration"""
        return (
            f"GTXJointConfig(\n"
            f"  gtx_exec={self.gtx_exec!r},\n"
            f"  reference={self.reference!r},\n"
            f"  gvcf_dir={self.gvcf_dir!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  threads={self.threads},\n"
            f"  tmp_dir={self.tmp_dir!r},\n"
            f"  script_name={self.script_name!r},\n"
            f"  chr_pattern={self.chr_pattern!r},\n"
            f"  window_size={self.window_size}\n"
            f")"
        )
