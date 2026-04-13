"""
FASTP处理核心模块|FASTP Processing Core Module
"""

import os
import shutil
import subprocess
from pathlib import Path
from .utils import CommandRunner


class FastpCore:
    """FASTP核心处理器|FASTP Core Processor"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        """
        初始化FASTP核心处理器|Initialize FASTP core processor

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def create_output_directories(self):
        """创建输出目录|Create output directories"""
        self.config.output_path.mkdir(parents=True, exist_ok=True)
        self.config.report_path.mkdir(parents=True, exist_ok=True)

        # 如果启用pair，创建临时目录和unpaired目录|Create temp and unpaired directories if pair is enabled
        if self.config.enable_pair:
            self.config.temp_fastp_path = self.config.output_path / "temp_after_fastp"
            self.config.temp_fastp_path.mkdir(parents=True, exist_ok=True)

            self.config.unpaired_path = self.config.output_path / "unpaired"
            self.config.unpaired_path.mkdir(parents=True, exist_ok=True)

            self.logger.info(f"创建临时目录|Created temp directories:")
            self.logger.info(f"  - temp_after_fastp: {self.config.temp_fastp_path}")
            self.logger.info(f"创建unpaired目录|Created unpaired directory: {self.config.unpaired_path}")

        self.logger.info(f"创建清洁数据输出目录|Created clean data output directory: {self.config.output_path}")
        self.logger.info(f"创建报告文件输出目录|Created report output directory: {self.config.report_path}")

    def validate_fastp(self) -> bool:
        """
        验证fastp可执行性|Validate fastp executable

        Returns:
            验证是否成功|Whether validation succeeded
        """
        if not self.cmd_runner.check_executable(self.config.fastp_path):
            self.logger.error(
                f"Fastp未找到或无执行权限|"
                f"Fastp not found or no execute permission: {self.config.fastp_path}"
            )
            return False

        self.logger.info(f"Fastp可执行文件验证成功|Fastp executable validated: {self.config.fastp_path}")
        return True

    def process_sample(self, sample_name: str, read1_file: Path, read2_file: Path = None) -> bool:
        """
        处理单个样本|Process single sample

        流程|Flow:
        - 启用pair: fastp -> seqkit pair -> 最终输出
        - 不启用pair: fastp -> 最终输出

        Args:
            sample_name: 样本名称|Sample name
            read1_file: Read1文件路径|Read1 file path
            read2_file: Read2文件路径（单末端模式为None）|Read2 file path (None for single-end mode)

        Returns:
            处理是否成功|Whether processing succeeded
        """
        # 单末端模式|Single-end mode
        if self.config.single_end:
            return self._process_single_end(sample_name, read1_file)

        # 双末端模式|Paired-end mode
        return self._process_paired_end(sample_name, read1_file, read2_file)

    def _process_single_end(self, sample_name: str, read1_file: Path) -> bool:
        """处理单末端样本|Process single-end sample"""
        # 构建输出文件路径|Construct output file paths
        final_read1 = self.config.output_path / f"{sample_name}.clean.fq.gz"
        html_report = self.config.report_path / f"{sample_name}.html"
        json_report = self.config.report_path / f"{sample_name}.json"

        # 检查输出文件是否已存在|Check if output files already exist
        if final_read1.exists():
            if not self.config.force:
                self.logger.info(f"样本 {sample_name} 的输出文件已存在，跳过处理|Output files exist for sample {sample_name}, skipping")
                return True
            else:
                self.logger.info(f"强制模式：覆盖已存在的文件|Force mode: overwriting existing files for sample {sample_name}")

        # 打印处理信息|Print processing information
        self.logger.info("---")
        self.logger.info(f"正在处理样本（单末端模式）|Processing sample (single-end mode): {sample_name}")
        self.logger.info(f"  Input (输入): {read1_file}")
        self.logger.info(f"  Output (输出): {final_read1}")
        self.logger.info(f"  HTML report (报告): {html_report}")
        self.logger.info("---")

        # 构建fastp命令（单末端）|Build fastp command (single-end)
        cmd = [
            self.config.fastp_path,
            "-i", str(read1_file),
            "-o", str(final_read1),
            "-h", str(html_report),
            "-j", str(json_report),
            "-w", str(self.config.threads),
            "-q", str(self.config.quality_threshold),
            "-l", str(self.config.min_length),
            "-u", str(self.config.unqualified_percent),
            "-n", str(self.config.n_base_limit)
        ]

        # 检查是否为dry-run模式|Check if dry-run mode
        if self.config.dry_run:
            self.logger.info(f"[DRY-RUN] 将执行命令|Would execute command for: {sample_name}")
            return True

        # 执行fastp命令|Execute fastp command
        fastp_success = self.cmd_runner.run(cmd, f"FASTP质控处理|FASTP quality control -> {sample_name}")

        if fastp_success:
            self.logger.info(f"样本 {sample_name} 处理完成|Sample {sample_name} processing completed")
            return True
        else:
            self.logger.error(f"样本 {sample_name} FASTP处理失败|Sample {sample_name} FASTP processing failed")
            return False

    def _process_paired_end(self, sample_name: str, read1_file: Path, read2_file: Path) -> bool:
        """
        处理双末端样本|Process paired-end sample

        流程|Flow:
        - 启用pair: fastp -> seqkit pair -> 最终输出
        - 不启用pair: fastp -> 最终输出
        """
        html_report = self.config.report_path / f"{sample_name}.html"
        json_report = self.config.report_path / f"{sample_name}.json"

        # 检查最终输出文件是否已存在|Check if final output files already exist
        final_read1 = self.config.output_path / f"{sample_name}_1.clean.fq.gz"
        final_read2 = self.config.output_path / f"{sample_name}_2.clean.fq.gz"
        if final_read1.exists() and final_read2.exists():
            if not self.config.force:
                self.logger.info(f"样本 {sample_name} 的输出文件已存在，跳过处理|Output files exist for sample {sample_name}, skipping")
                return True
            else:
                self.logger.info(f"强制模式：覆盖已存在的文件|Force mode: overwriting existing files for sample {sample_name}")

        if self.config.enable_pair:
            # fastp -> seqkit pair -> 最终输出
            return self._process_with_pair(sample_name, read1_file, read2_file, html_report, json_report)
        else:
            # 直接fastp|Direct fastp
            return self._process_fastp_only(sample_name, read1_file, read2_file, final_read1, final_read2, html_report, json_report)

    def _process_with_pair(self, sample_name: str, read1_file: Path, read2_file: Path,
                           html_report: Path, json_report: Path) -> bool:
        """
        使用fastp + seqkit pair处理样本|Process sample with fastp + seqkit pair

        步骤|Steps:
        1. fastp质控 -> temp_after_fastp/
        2. seqkit pair配对修复 -> 最终输出 | seqkit pair -> final output
        """
        # 步骤1: fastp质控|Step 1: fastp QC
        self.logger.info("=" * 60)
        self.logger.info(f"步骤1: FASTP质控|Step 1: FASTP QC -> {sample_name}")
        self.logger.info("=" * 60)

        fastp_read1 = self.config.temp_fastp_path / f"{sample_name}_1.fastp.fq.gz"
        fastp_read2 = self.config.temp_fastp_path / f"{sample_name}_2.fastp.fq.gz"

        if not self._run_fastp(read1_file, read2_file, fastp_read1, fastp_read2,
                              html_report, json_report, sample_name):
            return False

        # 步骤2: seqkit pair配对修复|Step 2: seqkit pair
        self.logger.info("=" * 60)
        self.logger.info(f"步骤2: SEQKIT PAIR配对修复|Step 2: SEQKIT PAIR -> {sample_name}")
        self.logger.info("=" * 60)

        final_read1 = self.config.output_path / f"{sample_name}_1.clean.fq.gz"
        final_read2 = self.config.output_path / f"{sample_name}_2.clean.fq.gz"

        if not self._run_seqkit_pair(fastp_read1, fastp_read2, final_read1, final_read2, sample_name):
            return False

        # 清理fastp临时文件|Clean up fastp temp files
        fastp_read1.unlink(missing_ok=True)
        fastp_read2.unlink(missing_ok=True)
        self.logger.debug(f"已删除fastp临时文件|Deleted fastp temp files for {sample_name}")

        self.logger.info(f"样本 {sample_name} 处理完成|Sample {sample_name} processing completed")
        return True

    def _process_fastp_only(self, sample_name: str, read1_file: Path, read2_file: Path,
                           final_read1: Path, final_read2: Path,
                           html_report: Path, json_report: Path) -> bool:
        """只使用fastp处理样本|Process sample with fastp only"""
        self.logger.info("=" * 60)
        self.logger.info(f"FASTP质控|FASTP Quality Control -> {sample_name}")
        self.logger.info("=" * 60)

        return self._run_fastp(read1_file, read2_file, final_read1, final_read2,
                             html_report, json_report, sample_name)

    def _run_fastp(self, read1: Path, read2: Path, out1: Path, out2: Path,
                  html_report: Path, json_report: Path, sample_name: str) -> bool:
        """运行fastp命令|Run fastp command"""
        self.logger.info(f"  Read 1 (输入|input): {read1}")
        self.logger.info(f"  Read 2 (输入|input): {read2}")
        self.logger.info(f"  Read 1 (输出|output): {out1}")
        self.logger.info(f"  Read 2 (输出|output): {out2}")
        self.logger.info(f"  HTML report (报告): {html_report}")
        self.logger.info("---")

        cmd = [
            self.config.fastp_path,
            "-i", str(read1),
            "-I", str(read2),
            "-o", str(out1),
            "-O", str(out2),
            "-h", str(html_report),
            "-j", str(json_report),
            "-w", str(self.config.threads),
            "-q", str(self.config.quality_threshold),
            "-l", str(self.config.min_length),
            "-u", str(self.config.unqualified_percent),
            "-n", str(self.config.n_base_limit)
        ]

        if self.config.dry_run:
            self.logger.info(f"[DRY-RUN] 将执行命令|Would execute command for: {sample_name}")
            return True

        success = self.cmd_runner.run(cmd, f"FASTP质控处理|FASTP quality control -> {sample_name}")

        if success:
            self.logger.info(f"FASTP处理成功|FASTP processing succeeded for {sample_name}")
        else:
            self.logger.error(f"FASTP处理失败|FASTP processing failed for {sample_name}")

        return success

    def _detect_id_regexp(self, read1: Path, read2: Path) -> str:
        """
        自动检测R1/R2文件的read name配对后缀格式，生成seqkit pair的--id-regexp
        Auto-detect read name pairing suffix format from R1/R2 files

        原理：从两个文件中各读取前N条read name，找到公共前缀，比对差异后缀，
        据此推断配对标识方式（/1 /2, .1 .2, _1 _2, Casava空格格式等）。
        Read first N read names from both files, find the common prefix,
        compare the differing suffixes, and infer the pairing format.

        支持格式|Supported formats:
        - /1 /2 (BGI/MGISEQ, old Illumina)
        - .1 .2 (SRA/NCBI fastq-dump)
        - _1 _2 (underscore separated)
        - space + 1:... / 2:... (Illumina Casava 1.8+)
        - f r (old 454/Sanger)
        - 无后缀 (modern bcl2fastq/DRAGEN, same names in both files)

        Args:
            read1: R1文件路径|R1 file path
            read2: R2文件路径|R2 file path

        Returns:
            --id-regexp参数字符串，None表示无需特殊处理
            --id-regexp string, None means no special handling needed
        """
        import gzip
        import re

        def first_names(path, n=20):
            names = []
            with gzip.open(str(path), 'rt') as f:
                for i, line in enumerate(f):
                    if i % 4 == 0:
                        names.append(line[1:].rstrip())
                        if len(names) >= n:
                            break
            return names

        r1_names = first_names(read1)
        r2_names = first_names(read2)

        if not r1_names or not r2_names:
            return None

        # 找到所有read name对的最长公共前缀|Find longest common prefix across all pairs
        prefix_len = len(r1_names[0])
        for n1, n2 in zip(r1_names, r2_names):
            i = 0
            while i < min(len(n1), len(n2)) and n1[i] == n2[i]:
                i += 1
            prefix_len = min(prefix_len, i)

        # R1/R2的read name完全相同（如bcl2fastq输出）|Names identical (e.g., bcl2fastq)
        if all(len(n) <= prefix_len for n in r1_names):
            self.logger.info("  配对格式检测|Pair format detected: 无后缀(no suffix), 使用默认regex")
            return None

        # 提取差异后缀并验证一致性|Extract differing suffixes and verify consistency
        r1_suffix = r1_names[0][prefix_len:]
        r2_suffix = r2_names[0][prefix_len:]

        for n1, n2 in zip(r1_names, r2_names):
            if n1[prefix_len:] != r1_suffix or n2[prefix_len:] != r2_suffix:
                return None

        # 在后缀中定位第一个不同字符|Find first differing character in suffix
        for i, (c1, c2) in enumerate(zip(r1_suffix, r2_suffix)):
            if c1 != c2:
                sep = r1_suffix[:i]  # 分隔符|separator
                # 如果分隔符为空，向前回退一个字符作为分隔符
                # If separator is empty, include the preceding character as separator
                if not sep and prefix_len > 0:
                    prefix_len -= 1
                    sep = r1_names[0][prefix_len]
                rest = r1_suffix[i + 1:]  # 1/2之后的内容|content after 1/2

                # 根据分隔符和后续内容构建正则|Build regex based on separator and rest
                if sep == " ":
                    # Casava 1.8+: "EAS139:xxx 1:N:18:1"
                    regexp = r"^(\S+)\s+[12]"
                    desc = f"Casava空格格式(space + 1/2)"
                elif rest:
                    # 后缀还有更多内容（如1:N:18:1），不加$锚定
                    # Suffix has more content, no $ anchor
                    regexp = f"^(\\S+?){re.escape(sep)}[12]"
                    desc = f"分隔符'{sep}' + 1/2 + 更多内容"
                else:
                    # 后缀仅是分隔符+1/2（如/1, _1, .1），用$锚定避免多匹配
                    # Suffix is just separator + 1/2, use $ anchor
                    regexp = f"^(\\S+?){re.escape(sep)}[12]$"
                    desc = f"分隔符'{sep}' + 1/2"

                self.logger.info(f"  配对格式检测|Pair format detected: {desc}, regex={regexp}")
                return regexp

        return None

    def _run_seqkit_pair(self, read1: Path, read2: Path,
                         final_read1: Path, final_read2: Path,
                         sample_name: str) -> bool:
        """
        运行seqkit pair命令|Run seqkit pair command

        seqkit pair将fastp输出中未配对的reads分离，确保最终R1/R2严格配对。
        seqkit pair separates unpaired reads from fastp output to ensure strict R1/R2 pairing.

        Args:
            read1: fastp输出的R1文件|fastp output R1 file
            read2: fastp输出的R2文件|fastp output R2 file
            final_read1: 最终配对后的R1输出路径|Final paired R1 output path
            final_read2: 最终配对后的R2输出路径|Final paired R2 output path
            sample_name: 样本名称|Sample name

        Returns:
            执行是否成功|Whether execution succeeded
        """
        self.logger.info(f"  Input R1 (输入): {read1}")
        self.logger.info(f"  Input R2 (输入): {read2}")
        self.logger.info(f"  Final R1 (最终输出): {final_read1}")
        self.logger.info(f"  Final R2 (最终输出): {final_read2}")
        self.logger.info("---")

        # 使用临时输出目录，因为seqkit pair -O会保持输入文件名不变
        # Use temp output dir because seqkit pair -O keeps input file names unchanged
        pair_output_dir = self.config.output_path / "temp_after_pair"
        pair_output_dir.mkdir(parents=True, exist_ok=True)

        # 自动检测R1/R2配对后缀格式|Auto-detect R1/R2 pairing suffix format
        id_regexp = self._detect_id_regexp(read1, read2)

        cmd = [
            self.config.seqkit_path, "pair",
            "-1", str(read1),
            "-2", str(read2),
            "-O", str(pair_output_dir),
            "-u",  # 保存未配对reads|save unpaired reads
            "-f",  # 覆盖已有输出|force overwrite
            "-j", str(self.config.threads),
        ]
        if id_regexp:
            cmd.extend(["--id-regexp", id_regexp])

        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        if self.config.dry_run:
            self.logger.info(f"[DRY-RUN] 将执行命令|Would execute seqkit pair for: {sample_name}")
            return True

        try:
            subprocess.run(cmd, check=True)
        except subprocess.CalledProcessError as e:
            self.logger.error(f"SEQKIT PAIR执行失败|SEQKIT PAIR execution failed with return code {e.returncode}")
            # 清理临时目录|Clean up temp directory
            shutil.rmtree(pair_output_dir, ignore_errors=True)
            return False
        except FileNotFoundError:
            self.logger.error(f"seqkit未找到|seqkit not found: {self.config.seqkit_path}")
            shutil.rmtree(pair_output_dir, ignore_errors=True)
            return False
        except Exception as e:
            self.logger.error(f"SEQKIT PAIR处理出错|Error during seqkit pair processing: {e}")
            shutil.rmtree(pair_output_dir, ignore_errors=True)
            return False

        # seqkit pair -O 保持输入文件名不变，在pair_output_dir中查找配对输出
        # seqkit pair -O keeps input file names unchanged, find paired output in pair_output_dir
        paired_r1 = pair_output_dir / read1.name
        paired_r2 = pair_output_dir / read2.name

        if not paired_r1.exists() or not paired_r2.exists():
            self.logger.error(
                f"seqkit pair输出文件未找到|seqkit pair output files not found: "
                f"R1={paired_r1} (exists={paired_r1.exists()}), R2={paired_r2} (exists={paired_r2.exists()})"
            )
            shutil.rmtree(pair_output_dir, ignore_errors=True)
            return False

        # 将未配对reads移动到unpaired目录|Move unpaired reads to unpaired directory
        for f in pair_output_dir.iterdir():
            if f.is_file() and "unpaired" in f.name:
                target = self.config.unpaired_path / f.name
                shutil.move(str(f), str(target))
                self.logger.info(f"未配对reads已保存|Unpaired reads saved: {target.name}")

        # 移动配对文件到最终位置|Move paired files to final location
        try:
            shutil.move(str(paired_r1), str(final_read1))
            shutil.move(str(paired_r2), str(final_read2))
        except Exception as e:
            self.logger.error(f"移动文件失败|Failed to move files: {e}")
            shutil.rmtree(pair_output_dir, ignore_errors=True)
            return False

        # 清理临时目录|Clean up temp directory
        shutil.rmtree(pair_output_dir, ignore_errors=True)

        self.logger.info(f"SEQKIT PAIR配对修复完成|SEQKIT PAIR completed for {sample_name}")
        return True
