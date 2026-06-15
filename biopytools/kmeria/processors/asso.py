"""k-mer关联分析处理器|K-mer Association Analysis Processor"""

import os
import time
from glob import glob
from ..utils import CommandRunner, format_number


class AssoProcessor:
    """k-mer关联分析处理器|K-mer Association Analysis Processor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run(self) -> bool:
        """运行k-mer关联分析|Run k-mer association analysis"""
        self.logger.info("开始k-mer关联分析|Starting k-mer association analysis")

        # 检查输出是否已存在|Check if output already exists
        # 包括多种可能的输出文件格式|Include multiple possible output file formats
        assoc_files = glob(os.path.join(self.config.output_dir, '*.assoc.txt'))
        assoc_files += glob(os.path.join(self.config.output_dir, '*.assoc'))
        assoc_files += glob(os.path.join(self.config.output_dir, '*_assoc*'))
        assoc_files += glob(os.path.join(self.config.output_dir, '*.ps'))
        assoc_files += glob(os.path.join(self.config.output_dir, '*.reml'))

        if not self.config.force and assoc_files:
            total_size = sum(os.path.getsize(f) for f in assoc_files)
            self.logger.info(f"跳过关联分析|Skip association analysis (输出文件已存在|output files exist: {len(assoc_files)} files, {format_number(total_size)})")
            return True

        # 清理空的 filtered_kmer_files.bimbam.gz 文件
        # Remove empty filtered_kmer_files.bimbam.gz file
        # (文件小于100字节视为空|Files smaller than 100 bytes are considered empty)
        empty_bimbam = os.path.join(self.config.input_dir, 'filtered_kmer_files.bimbam.gz')
        if os.path.exists(empty_bimbam):
            file_size = os.path.getsize(empty_bimbam)
            if file_size < 100:
                self.logger.info(f"删除空的 BIMBAM 文件|Removing empty BIMBAM file: {empty_bimbam} ({file_size} bytes)")
                os.remove(empty_bimbam)

        kmeria_bin = os.path.join(self.config.kmeria_path, 'bin', 'kmeria')

        cmd = [
            kmeria_bin,
            'asso',
            '--tool', 'bimbamAsso' if self.config.use_bimbam_tools else 'gemma',
            '-i', self.config.input_dir,
            '-p', self.config.pheno_file,
            '-o', self.config.output_dir,
            '-n', str(self.config.pheno_col),
            '-t', str(self.config.threads)
        ]

        # bimbamAsso 工具需要额外参数|bimbamAsso tool needs extra parameters
        if self.config.use_bimbam_tools:
            # 指定输入文件为 gzip 压缩格式|Specify input files are gzipped
            cmd.append('--bimbam-gzip')

            # 如果没有提供 kinship 文件且需要 kinship，则自动生成|Generate kinship if needed and not provided
            if not self.config.kinship_file and self.config.use_kinship:
                cmd.append('--gen-kinship')

            # 生成样本文件|Generate sample file
            sample_file = self._generate_sample_file()
            if sample_file:
                cmd.extend(['-s', sample_file])

        # 添加可选文件|Add optional files
        if self.config.covar_file:
            cmd.extend(['-c', self.config.covar_file])
        if self.config.kinship_file:
            cmd.extend(['-k', self.config.kinship_file])

        # 添加精度参数|Add precision parameters (仅 bimbamAsso | bimbamAsso only)
        if self.config.use_bimbam_tools:
            cmd.extend(['--kin-precision', str(self.config.kinship_precision)])
            cmd.extend(['--out-precision', str(self.config.output_precision)])
            cmd.extend(['--kin-method', str(self.config.kinship_method)])

        # 分析方法|Analysis method (-m)
        if self.config.analysis_method:
            cmd.extend(['-m', self.config.analysis_method])

        # 关闭kinship（如果指定）|Disable kinship if requested
        if not self.config.use_kinship:
            cmd.append('--no-kinship')

        # 禁用GLS|Disable GLS
        if self.config.disable_gls:
            cmd.append('--disable-gls')

        # 输出特征值/特征向量|Output eigenvalues/eigenvectors
        if self.config.write_eigen:
            cmd.append('--write-eigen')

        # 质量控制|Quality control
        if self.config.minor_allele_freq is not None:
            cmd.extend(['--maf', str(self.config.minor_allele_freq)])
        if self.config.missing_threshold is not None:
            cmd.extend(['--miss', str(self.config.missing_threshold)])

        # 分块分析|Block analysis
        if self.config.start_marker is not None:
            cmd.extend(['--start-marker', str(self.config.start_marker)])
        if self.config.end_marker is not None:
            cmd.extend(['--end-marker', str(self.config.end_marker)])

        # 其他选项|Other options
        if self.config.generate_plots:
            cmd.append('--generate-plots')
        if self.config.compress_output:
            cmd.append('--compress')
        if self.config.verbose:
            cmd.append('--verbose')
        if self.config.dry_run:
            cmd.append('--dry-run')
        if self.config.no_validate:
            cmd.append('--no-validate')
        if self.config.no_check_deps:
            cmd.append('--no-check-deps')
        if self.config.no_cleanup:
            cmd.append('--no-cleanup')

        # 记录命令开始时间，用于事后精准识别本次产生的输出文件
        # Record start time to precisely identify files produced by this run
        cmd_start_time = time.time()

        success, _ = self.cmd_runner.run_command(
            cmd,
            description="执行k-mer关联分析|Executing k-mer association analysis"
        )

        if success:
            self.logger.info("k-mer关联分析完成|K-mer association analysis completed")

            # 修复软件bug：bimbamAsso可能将文件输出到工作目录而非指定目录
            # Fix software bug: bimbamAsso may output files to working directory instead of specified directory
            # 仅搬移本次命令开始后新建的文件，避免误移用户已有同名文件
            # Only move files created during this run to avoid clobbering user's pre-existing files
            import shutil
            cwd = os.getcwd()
            moved = 0

            for ext in ('ps', 'reml', 'log'):
                for f in glob(os.path.join(cwd, f'*.{ext}')):
                    try:
                        mtime = os.path.getmtime(f)
                    except OSError:
                        continue
                    # 仅搬移本次运行产生的文件（容差1秒）|Only files from this run (1s tolerance)
                    if mtime >= cmd_start_time - 1:
                        dest = os.path.join(self.config.output_dir, os.path.basename(f))
                        if not os.path.exists(dest):
                            shutil.move(f, dest)
                            moved += 1

            if moved:
                self.logger.info(f"已将|Moved {moved} 个工作目录残留文件移动到输出目录|files from cwd to output dir")

            # 统计输出文件|Statistics output files
            # 包括所有可能的输出格式|Include all possible output formats
            assoc_files = glob(os.path.join(self.config.output_dir, '*.ps'))
            assoc_files += glob(os.path.join(self.config.output_dir, '*.reml'))
            assoc_files += glob(os.path.join(self.config.output_dir, '*assoc*'))

            if assoc_files:
                total_size = sum(os.path.getsize(f) for f in assoc_files)
                self.logger.info(f"生成|Generated {len(assoc_files)} 个关联结果文件|association result files ({format_number(total_size)})")
            else:
                self.logger.warning("未找到关联结果文件|No association result files found in output directory")

        return success

    def _generate_sample_file(self):
        """从表型文件生成样本文件|Generate sample file from phenotype file"""
        import tempfile

        try:
            # 读取表型文件，提取样本名|Read phenotype file and extract sample names
            samples = []
            with open(self.config.pheno_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('#'):
                        parts = line.split()
                        if parts:
                            samples.append(parts[0])

            if not samples:
                self.logger.warning("无法从表型文件提取样本名|Cannot extract sample names from phenotype file")
                return None

            # 创建样本文件|Create sample file
            # bimbamAsso 要求格式：每行两个相同的样本名（用空格分隔）
            # Format required by bimbamAsso: two identical sample names per line (space-separated)
            sample_file = os.path.join(self.config.output_dir, 'samples.txt')
            with open(sample_file, 'w') as f:
                for sample in samples:
                    f.write(f"{sample} {sample}\n")

            self.logger.info(f"生成样本文件|Generated sample file: {sample_file} ({len(samples)} samples)")
            return sample_file

        except Exception as e:
            self.logger.error(f"生成样本文件失败|Failed to generate sample file: {e}")
            return None
