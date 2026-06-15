"""k-mer计数处理器|K-mer Count Processor"""

import os
from pathlib import Path
from typing import List, Dict
from ..utils import CommandRunner, find_fastq_files, format_number


class CountProcessor:
    """k-mer计数处理器|K-mer Count Processor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """
        运行k-mer计数|Run k-mer counting

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始k-mer计数|Starting k-mer counting")

        # 读取样本列表|Read sample list
        from ..utils import read_sample_list
        samples = read_sample_list(self.config.samples_file)
        self.logger.info(f"样本数量|Sample count: {len(samples)}")

        # 查找FASTQ文件|Find FASTQ files
        sample_files = find_fastq_files(
            self.config.fastq_dir,
            samples,
            self.config.file_pattern
        )

        if len(sample_files) == 0:
            self.logger.error("未找到任何FASTQ文件|No FASTQ files found")
            return False

        self.logger.info(f"找到|Found {len(sample_files)} 个样本的FASTQ文件|samples with FASTQ files")

        # 批量处理|Batch processing
        success_count = 0
        failed_samples = []

        batch_size = self.config.batch_size
        total_batches = (len(samples) + batch_size - 1) // batch_size

        for batch_idx in range(total_batches):
            start_idx = batch_idx * batch_size
            end_idx = min((batch_idx + 1) * batch_size, len(samples))
            batch_samples = samples[start_idx:end_idx]

            self.logger.info(f"处理批次|Processing batch {batch_idx + 1}/{total_batches}: "
                           f"{len(batch_samples)} 个样本|samples")

            for sample in batch_samples:
                if sample not in sample_files:
                    self.logger.warning(f"跳过样本{sample}：未找到FASTQ文件|Skip sample {sample}: FASTQ files not found")
                    failed_samples.append(sample)
                    continue

                if self._count_sample(sample, sample_files[sample]):
                    success_count += 1
                else:
                    failed_samples.append(sample)

        # 统计结果|Statistics
        self.logger.info(f"k-mer计数完成|K-mer counting completed")
        self.logger.info(f"成功|Success: {success_count}/{len(samples)}")

        if failed_samples:
            self.logger.warning(f"失败样本|Failed samples ({len(failed_samples)}): {', '.join(failed_samples[:10])}"
                               if len(failed_samples) <= 10 else
                               f"失败样本|Failed samples ({len(failed_samples)}): {', '.join(failed_samples[:10])}...")

        return success_count > 0

    def _count_sample(self, sample: str, fastq_files: List[str]) -> bool:
        """
        计数单个样本|Count single sample

        Args:
            sample: 样本名|Sample name
            fastq_files: FASTQ文件列表|FASTQ file list

        Returns:
            是否成功|Whether successful
        """
        output_file = os.path.join(
            self.config.output_dir,
            f"{sample}_k{self.config.kmer_size}.bin"
        )

        # 检查输出文件是否已存在|Check if output file already exists
        if not self.config.force and os.path.exists(output_file):
            file_size = os.path.getsize(output_file)
            if file_size > 0:
                self.logger.info(f"跳过样本|Skip sample: {sample} (输出文件已存在|output file exists: {output_file}, {format_number(file_size)})")
                return True

        # 构建kmeria count命令|Build kmeria count command
        kmeria_bin = os.path.join(self.config.kmeria_path, 'bin', 'kmeria')

        cmd = [
            kmeria_bin,
            'count',
            '-k', str(self.config.kmer_size),
            '-t', str(self.config.threads),
            '-o', output_file
        ]

        # 添加可选参数|Add optional parameters
        if self.config.count_separate_strands:
            cmd.append('-C')

        if self.config.text_output:
            cmd.append('-T')

        # 添加输入文件|Add input files
        cmd.extend(fastq_files)

        # 执行命令|Execute command
        success, _ = self.cmd_runner.run_command(
            cmd,
            description=f"计数样本|Counting sample: {sample}"
        )

        if success:
            self.logger.debug(f"样本{sample}计数完成|Sample {sample} counting completed")

            # 统计输出文件大小|Statistics output file size
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file)
                self.logger.debug(f"输出文件|Output file: {output_file} ({format_number(file_size)} bytes)")

        return success
