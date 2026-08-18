"""
FASTQ GC过滤核心逻辑模块|FASTQ GC Filter Core Logic Module
基于seqkit+awk的高性能实现|High-performance implementation using seqkit+awk
"""

import subprocess
import tempfile
import os
from pathlib import Path
from typing import Tuple

from ..common.conda_runner import build_conda_command, run_pipeline


class FastqFilter:
    """FASTQ文件过滤器|FASTQ File Filter"""

    def __init__(self, config, logger):
        """
        初始化过滤器|Initialize filter

        Args:
            config: FastqGcFilterConfig配置对象|FastqGcFilterConfig object
            logger: 日志器对象|Logger object
        """
        self.config = config
        self.logger = logger

    @staticmethod
    def format_number(num: int) -> str:
        """
        格式化数字|Format number

        Args:
            num: 整数|Integer

        Returns:
            格式化后的字符串|Formatted string
        """
        if num >= 1_000_000:
            return f"{num / 1_000_000:.2f}M"
        elif num >= 1_000:
            return f"{num / 1_000:.2f}K"
        return str(num)

    def _build_awk_filter(self) -> str:
        """
        构建awk过滤表达式|Build awk filter expression

        seqkit fx2tab -n -l -g 输出格式: read_name<TAB>length<TAB>GC%
        $1: read名称|read name
        $2: 序列长度|sequence length
        $3: GC含量百分比|GC content percentage

        Returns:
            awk过滤字符串|awk filter string
        """
        conditions = []

        # GC含量过滤条件($3是GC含量)|GC content filter condition ($3 is GC content)
        gc_condition = f"$3 >= {self.config.min_gc} && $3 <= {self.config.max_gc}"
        conditions.append(gc_condition)

        # 序列长度过滤条件($2是序列长度)|Sequence length filter condition ($2 is sequence length)
        length_condition = f"$2 >= {self.config.min_length}"
        conditions.append(length_condition)

        if self.config.max_length is not None:
            max_length_condition = f"$2 <= {self.config.max_length}"
            conditions.append(max_length_condition)

        # 组合所有条件|Combine all conditions
        filter_expr = " && ".join(conditions)

        # 打印read名称(seqkit输出不包含@符号,直接打印$1)|Print read name (seqkit output has no @, print $1 directly)
        return f"{filter_expr} {{print $1}}"

    def filter_fastq(self) -> Tuple[int, int]:
        """
        根据GC含量和序列长度筛选FASTQ文件(使用seqkit+awk)|Filter FASTQ by GC content and sequence length (using seqkit+awk)

        Returns:
            tuple: (总reads数|total reads, 通过筛选数|passed reads)
        """
        self.logger.info(f"开始过滤FASTQ文件|Starting FASTQ filtering")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"输出文件|Output file: {self.config.output_file}")
        self.logger.info(f"GC含量范围|GC content range: {self.config.min_gc}% - {self.config.max_gc}%")
        self.logger.info(f"序列长度范围|Sequence length range: {self.config.min_length} - {self.config.max_length or 'unlimited'}")

        # 创建临时文件存放read名称列表|Create temp file for read name list
        tmp_root = os.path.join(os.path.dirname(self.config.output_file), 'tmp')
        os.makedirs(tmp_root, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode='w', suffix='.list', delete=False, dir=tmp_root) as tmp_list:
            tmp_list_path = tmp_list.name

        try:
            # 步骤1: 使用seqkit提取read名称、GC含量和长度|Step 1: Extract read names, GC content and length using seqkit
            self.logger.info(f"步骤1: 提取GC含量和序列长度|Step 1: Extracting GC content and sequence length")
            self.logger.info(f"步骤2: 过滤符合条件的reads|Step 2: Filtering qualified reads")

            awk_filter = self._build_awk_filter()

            # 执行seqkit + awk管道,生成read名称列表|Execute seqkit + awk pipeline to generate read name list
            # 方案B: seqkit经build_conda_command解析为域环境二进制, awk保持裸名, 严禁 conda run | conda run
            # |Solution B: seqkit resolved to domain env binary, awk stays bare; conda run in pipes is forbidden
            success, err = run_pipeline(
                [
                    [self.config.seqkit_path, 'fx2tab', '-n', '-l', '-g', self.config.input_file],
                    ['awk', awk_filter],
                ],
                self.logger,
                "seqkit提取GC含量+awk过滤|seqkit GC extraction + awk filtering",
                output_file=tmp_list_path,
            )
            if not success:
                self.logger.error(f"seqkit|awk管道执行失败|seqkit|awk pipeline failed: {err}")
                return 0, 0

            # 统计通过筛选的read数量|Count passed reads
            with open(tmp_list_path, 'r') as f:
                passed_reads = sum(1 for _ in f)

            self.logger.info(f"步骤3: 提取序列(使用多线程加速)|Step 3: Extracting sequences (multi-threaded)")

            # 步骤3: 使用seqkit grep根据read名称列表提取序列|Step 3: Use seqkit grep to extract sequences by read name list
            cmd_extract = build_conda_command(self.config.seqkit_path, [
                'grep',
                '-f', tmp_list_path,  # 从文件读取read名称列表|Read read name list from file
                '-n',  # 只匹配序列名称行|Only match sequence name lines
                '-j', '24',  # 使用24线程加速|Use 24 threads for acceleration
                self.config.input_file,
            ])

            # 如果输出文件是gzip格式,添加-o参数|If output is gzip format, add -o parameter
            if self.config.output_file.endswith('.gz'):
                cmd_extract.extend(['-o', self.config.output_file])
                self.logger.info(f"命令|Command: {' '.join(cmd_extract)}")
                result = subprocess.run(cmd_extract, capture_output=True, text=True)
            else:
                # 直接输出到文件|Output directly to file
                self.logger.info(f"命令|Command: {' '.join(cmd_extract)}")
                with open(self.config.output_file, 'w') as f_out:
                    result = subprocess.run(cmd_extract, stdout=f_out, stderr=subprocess.PIPE, text=True)

            if result.returncode != 0:
                self.logger.error(f"seqkit grep执行失败|seqkit grep execution failed: {result.stderr}")
                return 0, 0

            # 获取总reads数(使用seqkit stats)|Get total reads count (using seqkit stats)
            cmd_stats = build_conda_command(self.config.seqkit_path, ['stats', self.config.input_file])
            self.logger.info(f"命令|Command: {' '.join(cmd_stats)}")
            result_stats = subprocess.run(cmd_stats, capture_output=True, text=True, check=True)

            # 解析输出获取total reads|Parse output to get total reads
            total_reads = 0
            for line in result_stats.stdout.split('\n'):
                if 'total_seq' in line:
                    # seqkit stats输出格式: "total_seq: 123456"
                    total_reads = int(line.split(':')[1].strip())
                    break

            return total_reads, passed_reads

        finally:
            # 删除临时文件|Delete temp file
            if os.path.exists(tmp_list_path):
                os.unlink(tmp_list_path)
