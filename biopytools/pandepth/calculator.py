"""
PanDepth覆盖度计算核心逻辑|PanDepth Coverage Calculation Core Logic
"""

import os
from .utils import CommandRunner, BAMFileFinder


class PanDepthCalculator:
    """PanDepth覆盖度计算器|PanDepth Coverage Calculator"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.bam_finder = BAMFileFinder(logger)

    def _resolve_bam_path(self, bam_file: str) -> str:
        """解析BAM文件路径，处理软链接|Resolve BAM file path, handle symlinks

        Args:
            bam_file: BAM文件路径|BAM file path

        Returns:
            str: 解析后的绝对路径|Resolved absolute path
        """
        # 检查是否是软链接|Check if it's a symlink
        if os.path.islink(bam_file):
            self.logger.debug(f"检测到软链接|Detected symlink: {bam_file}")

            # 获取真实路径|Get real path
            real_path = os.path.realpath(bam_file)

            if not os.path.exists(real_path):
                self.logger.error(f"软链接目标文件不存在|Symlink target does not exist")
                self.logger.error(f"  软链接|Symlink: {bam_file}")
                self.logger.error(f"  指向|Points to: {real_path}")
                return None

            self.logger.debug(f"软链接已解析|Symlink resolved to: {real_path}")
            return real_path

        return bam_file

    def build_pandepth_command(self, bam_file: str, output_prefix: str) -> str:
        """构建PanDepth命令|Build PanDepth command

        Args:
            bam_file: BAM文件路径|BAM file path
            output_prefix: 输出文件前缀|Output file prefix

        Returns:
            str: 完整的PanDepth命令|Complete PanDepth command
        """
        cmd_parts = [self.config.pandepth_path]

        # 输入输出选项|Input/Output options
        cmd_parts.extend(['-i', bam_file])
        cmd_parts.extend(['-o', output_prefix])

        # 目标区域选项|Target region options
        if self.config.gff_file:
            cmd_parts.extend(['-g', self.config.gff_file])
            cmd_parts.extend(['-f', self.config.feature_type])

        if self.config.bed_file:
            cmd_parts.extend(['-b', self.config.bed_file])

        if self.config.window_size:
            cmd_parts.extend(['-w', str(self.config.window_size)])

        # 过滤选项|Filter options
        cmd_parts.extend(['-q', str(self.config.min_mapq)])
        cmd_parts.extend(['-d', str(self.config.min_depth)])
        cmd_parts.extend(['-x', str(self.config.exclude_flag)])

        # 其他选项|Other options
        cmd_parts.extend(['-t', str(self.config.threads)])

        if self.config.reference:
            cmd_parts.extend(['-r', self.config.reference])

        if self.config.enable_gc:
            cmd_parts.append('-c')

        if self.config.output_all_sites:
            cmd_parts.append('-a')

        return ' '.join(cmd_parts)

    def calculate_coverage_single_bam(self, bam_file: str, sample_name: str = None) -> bool:
        """计算单个BAM文件的覆盖度|Calculate coverage for single BAM file

        Args:
            bam_file: BAM文件路径|BAM file path
            sample_name: 样本名称(用于输出)|Sample name (for output)

        Returns:
            bool: 是否成功|Whether successful
        """
        # 解析BAM文件路径（处理软链接）|Resolve BAM file path (handle symlinks)
        resolved_bam = self._resolve_bam_path(bam_file)
        if resolved_bam is None:
            return False

        if not sample_name:
            # 从BAM文件名提取样本名|Extract sample name from BAM filename
            sample_name = os.path.splitext(os.path.basename(bam_file))[0]

        self.logger.info(f"处理样本|Processing sample: {sample_name}")
        self.logger.info(f"BAM文件|BAM file: {resolved_bam}")

        # 构建输出前缀|Build output prefix
        output_prefix = os.path.join(self.config.output_dir, sample_name)

        # 构建命令|Build command
        cmd = self.build_pandepth_command(resolved_bam, output_prefix)

        # 执行命令|Execute command
        success = self.cmd_runner.run(
            cmd,
            f"PanDepth覆盖度计算|PanDepth coverage calculation for {sample_name}"
        )

        if success:
            self.logger.info(f"样本处理成功|Sample processed successfully: {sample_name}")
            # 列出输出文件|List output files
            self._list_output_files(output_prefix)
        else:
            self.logger.error(f"样本处理失败|Sample processing failed: {sample_name}")

        return success

    def calculate_coverage_batch(self, bam_files: list) -> dict:
        """批量计算多个BAM文件的覆盖度|Calculate coverage for multiple BAM files in batch

        Args:
            bam_files: BAM文件路径列表|List of BAM file paths

        Returns:
            dict: 处理结果统计|Processing result statistics
        """
        results = {
            'total': len(bam_files),
            'success': 0,
            'failed': 0,
            'failed_samples': []
        }

        self.logger.info(f"开始批量处理 {len(bam_files)} 个样本|Starting batch processing of {len(bam_files)} samples")

        for i, bam_file in enumerate(bam_files, 1):
            sample_name = os.path.splitext(os.path.basename(bam_file))[0]

            self.logger.info(f"\n进度|Progress: [{i}/{len(bam_files)}]")

            success = self.calculate_coverage_single_bam(bam_file, sample_name)

            if success:
                results['success'] += 1
            else:
                results['failed'] += 1
                results['failed_samples'].append(sample_name)

        return results

    def _list_output_files(self, output_prefix: str):
        """列出生成的输出文件|List generated output files

        Args:
            output_prefix: 输出文件前缀|Output file prefix
        """
        self.logger.info(f"输出文件|Output files:")

        # 可能的输出文件|Possible output files
        possible_extensions = [
            '.chr.stat.gz',  # 染色体统计|Chromosome statistics
            '.gene.stat.gz',  # 基因统计|Gene statistics
            '.bed.stat.gz',  # BED区域统计|BED region statistics
            '.win.stat.gz'   # 窗口统计|Window statistics
        ]

        for ext in possible_extensions:
            output_file = output_prefix + ext
            if os.path.exists(output_file):
                file_size = os.path.getsize(output_file)
                self.logger.info(f"  - {output_file} ({file_size} bytes)")

    def run(self):
        """运行覆盖度计算|Run coverage calculation

        Returns:
            dict: 处理结果统计|Processing result statistics
        """
        # 查找BAM文件|Find BAM files
        if self.config.is_batch_mode:
            bam_files = self.bam_finder.find_bam_files(self.config.input_path)

            if not bam_files:
                self.logger.error("未找到任何BAM文件|No BAM files found")
                return {'total': 0, 'success': 0, 'failed': 0, 'failed_samples': []}

            # 批量处理|Batch processing
            results = self.calculate_coverage_batch(bam_files)

            # 输出统计结果|Output statistics
            self.logger.info("\n" + "=" * 60)
            self.logger.info("批量处理统计结果|Batch Processing Statistics")
            self.logger.info("=" * 60)
            self.logger.info(f"总样本数|Total samples: {results['total']}")
            self.logger.info(f"成功处理|Successfully processed: {results['success']}")
            self.logger.info(f"处理失败|Failed: {results['failed']}")

            if results['failed_samples']:
                self.logger.warning(f"失败样本|Failed samples: {', '.join(results['failed_samples'])}")

            return results

        else:
            # 单文件处理|Single file processing
            success = self.calculate_coverage_single_bam(self.config.input_path)

            return {
                'total': 1,
                'success': 1 if success else 0,
                'failed': 0 if success else 1,
                'failed_samples': [] if success else [os.path.basename(self.config.input_path)]
            }
