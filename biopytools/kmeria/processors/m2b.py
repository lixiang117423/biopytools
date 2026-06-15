"""矩阵转BIMBAM格式处理器|Matrix to BIMBAM Processor"""

import os
from glob import glob
from ..utils import CommandRunner, format_number


class M2bProcessor:
    """矩阵转BIMBAM格式处理器|Matrix to BIMBAM Processor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """运行矩阵转换|Run matrix conversion"""
        # 检查输出是否已存在|Check if output already exists
        bimbam_files = glob(os.path.join(self.config.output_dir, '*.bimbam.gz'))
        bimbam_files += glob(os.path.join(self.config.output_dir, '*.bimbam'))

        # 统计输入文件数量|Count input files
        input_files = glob(os.path.join(self.config.input_dir, 'filtered_*.*'))
        expected_count = len(input_files)

        if not self.config.force and bimbam_files:
            # 如果BIMBAM文件数量合理（至少是输入文件数量的一半），认为已完成
            if len(bimbam_files) >= expected_count * 0.5:
                total_size = sum(os.path.getsize(f) for f in bimbam_files)
                self.logger.info(f"跳过格式转换|Skip format conversion (输出文件已存在|output files exist: {len(bimbam_files)} files, {format_number(total_size)})")
                return True

        self.logger.info("开始转换k-mer矩阵为BIMBAM格式|Starting k-mer matrix to BIMBAM conversion")

        kmeria_bin = os.path.join(self.config.kmeria_path, 'bin', 'kmeria')

        cmd = [
            kmeria_bin,
            'm2b',
            '--in', self.config.input_dir,
            '--out', self.config.output_dir,
            '--threads', str(self.config.threads)
        ]

        # 添加归一化参数|Add normalization parameters
        if not self.config.normalize:
            cmd.append('--no-normalize')
        if self.config.quantile_norm:
            cmd.append('--quantile-norm')
            cmd.extend(['--lower-quantile', str(self.config.lower_quantile)])
            cmd.extend(['--upper-quantile', str(self.config.upper_quantile)])

        # 添加压缩参数|Add compression parameters
        if not self.config.compress:
            cmd.append('--no-compress')
        else:
            cmd.extend(['--level', str(self.config.compression_level)])
            cmd.extend(['--bgzf-threads', str(self.config.bgzf_threads)])

        success, _ = self.cmd_runner.run_command(
            cmd,
            description="转换k-mer矩阵|Converting k-mer matrix"
        )

        if success:
            self.logger.info("k-mer矩阵转换完成|K-mer matrix conversion completed")

            # 统计输出文件|Statistics output files
            bimbam_files = glob(os.path.join(self.config.output_dir, '*.bimbam*'))
            self.logger.info(f"生成|Generated {len(bimbam_files)} 个BIMBAM文件|BIMBAM files")

        return success
