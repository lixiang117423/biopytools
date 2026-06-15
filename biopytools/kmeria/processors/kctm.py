"""k-mer矩阵构建处理器|K-mer Matrix Construction Processor"""

import os
from pathlib import Path
from ..utils import CommandRunner, read_sample_list, format_number


class KctmProcessor:
    """k-mer矩阵构建处理器|K-mer Matrix Construction Processor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """运行k-mer矩阵构建|Run k-mer matrix construction"""
        # 检查输出是否已存在|Check if output already exists
        output_prefix = os.path.join(self.config.output_dir, 'kmer_matrix')
        from glob import glob
        existing_matrix_files = glob(f"{output_prefix}.*.bin")

        if not self.config.force and existing_matrix_files:
            total_size = sum(os.path.getsize(f) for f in existing_matrix_files)
            self.logger.info(f"跳过矩阵构建|Skip matrix construction (输出文件已存在|output files exist: {len(existing_matrix_files)} files, {format_number(total_size)})")
            return True

        self.logger.info("开始构建k-mer矩阵|Starting k-mer matrix construction")

        # 查找k-mer文件|Find k-mer files
        kmer_files = self._find_kmer_files()

        if not kmer_files:
            self.logger.error("未找到k-mer文件|No k-mer files found")
            return False

        self.logger.info(f"找到|Found {len(kmer_files)} 个k-mer文件|k-mer files")

        # 创建文件列表|Create file list
        list_file = os.path.join(self.config.output_dir, 'kmer_files.txt')
        with open(list_file, 'w') as f:
            for kf in kmer_files:
                f.write(f"{kf}\n")

        # 构建kctm命令|Build kctm command
        # 注意kmatrix.cpp中 -t 为文本模式，-j 为线程数，-b 为k-mer批次大小，-s 为每块条目数
        # Note: in kmatrix.cpp, -t means text mode, -j for threads, -b for k-mer batch, -s for entries per block
        kmeria_bin = os.path.join(self.config.kmeria_path, 'bin', 'kmeria')

        cmd = [
            kmeria_bin,
            'kctm',
            '-i', list_file,
            '-o', output_prefix,
            '-j', str(self.config.threads),
            '-b', str(self.config.kmer_batch_size),
            '-s', str(self.config.block_size)
        ]

        # 执行命令|Execute command
        success, _ = self.cmd_runner.run_command(
            cmd,
            description="构建k-mer矩阵|Building k-mer matrix"
        )

        if success:
            self.logger.info("k-mer矩阵构建完成|K-mer matrix construction completed")

            # 统计输出文件|Statistics output files
            matrix_files = list(Path(self.config.output_dir).glob('kmer_matrix.*'))
            self.logger.info(f"生成|Generated {len(matrix_files)} 个矩阵文件|matrix files")
            for mf in matrix_files:
                size = os.path.getsize(mf)
                self.logger.debug(f"  {mf.name} ({format_number(size)} bytes)")

        return success

    def _find_kmer_files(self):
        """查找k-mer文件|Find k-mer files"""
        from glob import glob

        pattern = os.path.join(self.config.input_dir, '*.bin')
        files = glob(pattern)

        self.logger.debug(f"在目录中搜索|Searching in directory: {self.config.input_dir}")
        self.logger.debug(f"文件模式|File pattern: {pattern}")

        return sorted(files)
