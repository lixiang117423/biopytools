"""
甲烷循环分析配置模块 | Methane Cycle Analysis Configuration Module
"""

import os
from typing import Optional


class MCycConfig:
    """甲烷循环分析配置类 | Methane Cycle Analysis Configuration Class"""

    def __init__(
        self,
        input_list: str,
        output_dir: Optional[str] = None,
        mcyc_base: Optional[str] = None,
        staging_dir_name: str = "mcyc_staging_fastq",
        thread_count: int = 4,
        skip_diamond: bool = False,
        keep_temp_files: bool = False
    ):
        """
        初始化配置 | Initialize configuration

        Args:
            input_list: 输入文件列表路径 | Input file list path
            output_dir: 输出目录 | Output directory
            mcyc_base: MCycDB数据库基础路径 | MCycDB database base path
            staging_dir_name: 临时数据目录名称 | Temporary data directory name
            thread_count: 线程数 | Thread count
            skip_diamond: 是否跳过diamond比对 | Whether to skip diamond alignment
            keep_temp_files: 是否保留临时文件 | Whether to keep temporary files
        """
        self.input_list = input_list
        self.output_dir = output_dir if output_dir else os.getcwd()
        self.staging_dir_name = staging_dir_name
        self.thread_count = thread_count
        self.skip_diamond = skip_diamond
        self.keep_temp_files = keep_temp_files

        # MCycDB相关路径 | MCycDB related paths
        if mcyc_base is None:
            self.mcyc_base = "/share/org/YZWL/yzwl_lixg/software/MCycDB/MCycDB-main"
        else:
            self.mcyc_base = mcyc_base

        # 数据库文件路径 | Database file paths
        self.fasta_file = os.path.join(self.mcyc_base, "MCycDB_2021.faa")
        self.map_file = os.path.join(self.mcyc_base, "id2gene.map")
        self.dmnd_file = os.path.join(self.mcyc_base, "MCycDB_2021.dmnd")
        self.diamond_db_base = os.path.join(self.mcyc_base, "MCycDB_2021")

        # 工作目录和输出文件路径 | Working directory and output file paths
        self.work_dir = os.getcwd()
        self.staging_dir = os.path.join(self.work_dir, self.staging_dir_name)
        self.raw_output = os.path.join(self.work_dir, "MCyc_Raw_Temp.txt")
        self.sample_info_path = os.path.join(self.work_dir, "sample_info.txt")

        # 最终输出文件 | Final output files
        self.raw_counts_output = os.path.join(self.output_dir, "Matrix_1_RawCounts.txt")
        self.tpm_output = os.path.join(self.output_dir, "Matrix_2_TPM.txt")
        self.clr_output = os.path.join(self.output_dir, "Matrix_3_CLR.txt")

        # Perl脚本路径 | Perl script path
        self.perl_script_src = os.path.join(self.mcyc_base, "MCycDB_FunctionProfiler.PL")
        self.perl_script_dst = os.path.join(self.work_dir, "MCycDB_FunctionProfiler.PL")

        # 确保输出目录存在 | Ensure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        if not self.input_list:
            raise ValueError("输入文件列表路径不能为空 | Input file list path cannot be empty")

        if not os.path.exists(self.input_list):
            raise FileNotFoundError(f"输入文件列表不存在 | Input file list not found: {self.input_list}")

        # 检查MCycDB相关文件 | Check MCycDB related files
        required_files = [self.fasta_file, self.map_file]
        for file_path in required_files:
            if not os.path.exists(file_path):
                raise FileNotFoundError(f"MCycDB文件不存在 | MCycDB file not found: {file_path}")

        if self.thread_count < 1:
            raise ValueError("线程数必须大于0 | Thread count must be greater than 0")

    def __str__(self) -> str:
        """配置的字符串表示 | String representation of configuration"""
        return f"""MCycConfig:
    输入文件列表 | Input List: {self.input_list}
    输出目录 | Output Directory: {self.output_dir}
    工作目录 | Working Directory: {self.work_dir}
    临时数据目录 | Staging Directory: {self.staging_dir}
    MCycDB基础路径 | MCycDB Base: {self.mcyc_base}
    线程数 | Thread Count: {self.thread_count}
    跳过Diamond | Skip Diamond: {self.skip_diamond}
    保留临时文件 | Keep Temp Files: {self.keep_temp_files}
"""