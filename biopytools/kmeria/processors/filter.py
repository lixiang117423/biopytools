"""k-mer矩阵过滤处理器|K-mer Matrix Filter Processor"""

import os
from glob import glob
from ..utils import CommandRunner, format_number


class FilterProcessor:
    """k-mer矩阵过滤处理器|K-mer Matrix Filter Processor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """运行k-mer过滤|Run k-mer filtering"""
        # 检查输出是否已存在|Check if output already exists
        existing_files = glob(os.path.join(self.config.output_dir, 'filtered_*.bin'))
        existing_files += glob(os.path.join(self.config.output_dir, 'filtered_*.txt'))

        if not self.config.force and existing_files:
            total_size = sum(os.path.getsize(f) for f in existing_files)
            self.logger.info(f"跳过矩阵过滤|Skip matrix filtering (输出文件已存在|output files exist: {len(existing_files)} files, {format_number(total_size)})")
            return True

        self.logger.info("开始过滤k-mer矩阵|Starting k-mer matrix filtering")

        kmeria_bin = os.path.join(self.config.kmeria_path, 'bin', 'kmeria')

        cmd = [
            kmeria_bin,
            'filter',
            '-i', self.config.input_dir,
            '-o', self.config.output_dir,
            '-d', self.config.depth_file,
            '-t', str(self.config.threads),
            '-c', str(self.config.max_abund),
            '-s', str(self.config.missing_ratio),
            '-p', str(self.config.ploidy)
        ]

        success, _ = self.cmd_runner.run_command(
            cmd,
            description="过滤k-mer矩阵|Filtering k-mer matrix"
        )

        if success:
            self.logger.info("k-mer矩阵过滤完成|K-mer matrix filtering completed")

            # 统计输出文件|Statistics output files
            filtered_files = glob(os.path.join(self.config.output_dir, 'filtered_*.*'))
            self.logger.info(f"生成|Generated {len(filtered_files)} 个过滤文件|filtered files")

        return success
