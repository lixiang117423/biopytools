"""
HapHiC序列比对模块 | HapHiC Sequence Alignment Module
使用官方推荐的BWA比对流程
"""

import os
import subprocess
import glob
import logging
from pathlib import Path
from typing import List, Tuple, Optional


class BWAAligner:
    """BWA比对器 | BWA Aligner - 使用官方推荐的流程"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def run_alignment(self, is_realignment: bool = False) -> str:
        """
        运行BWA比对流程 | Run BWA alignment pipeline using official method

        Args:
            is_realignment: 是否是纠错后的重新比对 | Whether this is realignment after correction

        Returns:
            str: 最终过滤后的BAM文件路径 | Final filtered BAM file path
        """
        if self.config.hic_file_type != "fastq":
            return self.config.bam_file

        self.logger.info("🔗 开始BWA比对 | Starting BWA alignment (official method)")

        # 创建比对输出目录 | Create alignment output directory
        if is_realignment:
            alignment_dir = os.path.join(self.config.output_dir, "00.mapping_corrected")
            self.logger.info("📂 使用纠错后重新比对目录 | Using corrected realignment directory")
        else:
            alignment_dir = os.path.join(self.config.output_dir, "00.mapping")
        os.makedirs(alignment_dir, exist_ok=True)

        # 获取FASTQ文件对 | Get FASTQ file pairs
        if hasattr(self.config, 'hic2_file') and self.config.hic2_file:
            # 使用提供的两个文件 | Use the two provided files
            fastq_pairs = [(self.config.hic_file, self.config.hic2_file)]
            self.logger.info(f"📄 使用提供的FASTQ文件对 | Using provided FASTQ pair: {os.path.basename(self.config.hic_file)}")
        else:
            # 查找FASTQ文件 | Find FASTQ files
            fastq_pairs = self._find_fastq_pairs()

        if not fastq_pairs:
            raise ValueError("❌ 未找到FASTQ文件对 | No FASTQ pairs found")

        self.logger.info(f"📄 找到 {len(fastq_pairs)} 对FASTQ文件 | Found {len(fastq_pairs)} FASTQ pairs")

        # 在当前工作目录创建软链接 | Create symbolic links in current directory
        original_cwd = os.getcwd()
        os.chdir(alignment_dir)

        try:
            # 创建基因组文件的软链接 | Create symlink for assembly file
            asm_link = "assembly.fa"
            if os.path.exists(asm_link):
                os.remove(asm_link)
            os.symlink(self.config.asm_file, asm_link)
            self.logger.info(f"🔗 创建基因组软链接 | Created assembly symlink: {asm_link}")

            # 创建FASTQ文件的软链接 | Create symlinks for FASTQ files
            fastq_links = []
            for i, (read1, read2) in enumerate(fastq_pairs):
                read1_link = f"read1_{i+1}.fastq.gz"
                read2_link = f"read2_{i+1}.fastq.gz"

                # 创建软链接 | Create symlinks
                if os.path.exists(read1_link):
                    os.remove(read1_link)
                if os.path.exists(read2_link):
                    os.remove(read2_link)

                os.symlink(read1, read1_link)
                os.symlink(read2, read2_link)
                fastq_links.append((read1_link, read2_link))

                self.logger.info(f"🔗 创建FASTQ软链接 | Created FASTQ symlinks: {read1_link}, {read2_link}")

            # 执行BWA比对流程 | Execute BWA alignment pipeline
            final_bam = self._run_official_alignment_pipeline(asm_link, fastq_links, is_realignment)

            self.logger.info(f"✅ BWA比对完成 | BWA alignment completed: {final_bam}")
            return os.path.join(alignment_dir, final_bam)

        finally:
            # 返回原工作目录 | Return to original directory
            os.chdir(original_cwd)

    def _find_fastq_pairs(self) -> List[Tuple[str, str]]:
        """查找FASTQ文件对 | Find FASTQ file pairs"""
        fastq_pairs = []

        if os.path.isfile(self.config.hic_file):
            # 单个文件 | Single file
            if self.config.hic_file.endswith('.fastq.gz') or self.config.hic_file.endswith('.fq.gz'):
                base_name = self.config.hic_file
                read1 = base_name
                # 查找对应的read2文件
                if self.config.read2_pattern in base_name:
                    read2 = base_name.replace(self.config.read2_pattern, self.config.read1_pattern)
                else:
                    read2 = base_name.replace(self.config.read1_pattern, self.config.read2_pattern)

                if os.path.exists(read2):
                    fastq_pairs.append((read1, read2))

        elif os.path.isdir(self.config.hic_file):
            # 目录 | Directory
            # 查找所有read1文件 | Find all read1 files
            read1_pattern = self.config.hic_file + "/*" + self.config.read1_pattern
            read1_files = glob.glob(read1_pattern)

            for read1 in read1_files:
                # 查找对应的read2文件 | Find corresponding read2 file
                read2 = read1.replace(self.config.read1_pattern, self.config.read2_pattern)
                if os.path.exists(read2):
                    fastq_pairs.append((read1, read2))

        return fastq_pairs

    def _run_official_alignment_pipeline(self, asm_link: str, fastq_links: List[Tuple[str, str]], is_realignment: bool = False) -> str:
        """
        运行官方推荐的BWA比对流程 | Run official BWA alignment pipeline

        官方流程:
        1. bwa mem -5SP -t threads asm.fa read1.fq.gz read2.fq.gz | samblaster | samtools view -@ threads -S -h -b -F 3340 -o HiC.bam
        2. HapHiC/utils/filter_bam HiC.bam 1 --nm 3 --threads threads | samtools view -b -@ threads -o HiC.filtered.bam

        Args:
            is_realignment: 是否是纠错后的重新比对 | Whether this is realignment after correction
        """

        # (1) BWA索引 | BWA indexing
        self._ensure_bwa_index(asm_link)

        # (1.5) 检查并调整线程数以避免内存问题 | Check and adjust threads to avoid memory issues
        self._adjust_threads_for_memory()

        # (2) BWA比对 + samblaster + samtools view | BWA alignment + samblaster + samtools view
        self.logger.info(f"🔗 步骤1: BWA比对 + 去重 + 初步过滤 (使用 {self.config.threads} 线程) | Step 1: BWA alignment + deduplication + initial filtering (using {self.config.threads} threads)")

        # 根据是否是重新比对设置文件名 | Set filename based on whether this is realignment
        if is_realignment:
            raw_bam = "HiC.corrected.bam"
            self.logger.info("📝 创建纠错后的BAM文件 | Creating corrected BAM file")
        else:
            raw_bam = "HiC.bam"

        # 构建BWA比对命令 | Build BWA alignment command (使用官方参数)
        bwa_cmd = [
            self.config.bwa_bin, "mem",
            "-5SP",  # 5: mark shorter split hits as secondary, S: output mate pairs, P: output in sam format
            "-t", str(self.config.threads),
            asm_link
        ]

        # 添加所有FASTQ文件 | Add all FASTQ files
        for read1_link, read2_link in fastq_links:
            bwa_cmd.extend([read1_link, read2_link])

        try:
            if self.config.dry_run:
                self.logger.info(f"[DRY RUN] BWA alignment command: {' '.join(bwa_cmd)}")
                if self.config.use_samblaster:
                    self.logger.info(f"[DRY RUN] ... | {self.config.samblaster_bin} | samtools view -@ {self.config.threads} -S -h -b -F 3340 -o {raw_bam}")
                else:
                    self.logger.info(f"[DRY RUN] ... | samtools view -@ {self.config.threads} -S -h -b -F 3340 -o {raw_bam}")
            else:
                # 使用官方的管道命令格式 | Use official pipeline command format
                if self.config.use_samblaster:
                    # 构建完整的管道命令字符串 | Build complete pipeline command string
                    bwa_cmd_str = " ".join(bwa_cmd)
                    # 完全按照官方的格式 | Exactly following official format
                    full_cmd = f"{bwa_cmd_str} | {self.config.samblaster_bin} | {self.config.samtools_bin} view -@ {self.config.threads} -S -h -b -F 3340 -o {raw_bam}"

                    self.logger.info(f"🔄 执行BWA比对管道 | Executing BWA alignment pipeline")

                    result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True)

                    if result.returncode != 0:
                        error_msg = f"BWA比对管道失败 | BWA alignment pipeline failed:\n"
                        error_msg += f"输出 | Output: {result.stderr}\n" if result.stderr else ""
                        error_msg += f"命令 | Command: {full_cmd}"
                        self.logger.error(error_msg)
                        raise RuntimeError(error_msg)

                else:
                    # 不使用samblaster，直接用samtools | Don't use samblaster, use samtools directly
                    bwa_cmd_str = " ".join(bwa_cmd)
                    # 完全按照官方的格式 | Exactly following official format
                    full_cmd = f"{bwa_cmd_str} | {self.config.samtools_bin} view -@ {self.config.threads} -S -h -b -F 3340 -o {raw_bam}"

                    self.logger.info(f"🔄 执行BWA比对管道（无samblaster）| Executing BWA alignment pipeline (no samblaster)")

                    result = subprocess.run(full_cmd, shell=True, capture_output=True, text=True)

                    if result.returncode != 0:
                        error_msg = f"BWA比对管道失败 | BWA alignment pipeline failed:\n"
                        error_msg += f"输出 | Output: {result.stderr}\n" if result.stderr else ""
                        error_msg += f"命令 | Command: {full_cmd}"
                        self.logger.error(error_msg)
                        raise RuntimeError(error_msg)

                self.logger.info(f"✅ BWA比对完成，生成 {raw_bam}")

                # 检查文件大小 | Check file size
                if os.path.exists(raw_bam):
                    size_mb = os.path.getsize(raw_bam) / (1024 * 1024)
                    self.logger.info(f"📊 BAM文件大小: {size_mb:.2f} MB")

        except Exception as e:
            self.logger.error(f"❌ BWA比对失败 | BWA alignment failed: {e}")
            raise

        # (3) HapHiC filter_bam过滤 | HapHiC filter_bam filtering
        if self.config.use_haphic_filter:
            self.logger.info("🔍 步骤2: HapHiC过滤 | Step 2: HapHiC filtering")
            # 根据是否是重新比对设置文件名 | Set filename based on whether this is realignment
            if is_realignment:
                filtered_bam = "HiC.corrected.filtered.bam"
                self.logger.info("📝 创建纠错后的过滤BAM文件 | Creating corrected filtered BAM file")
            else:
                filtered_bam = "HiC.filtered.bam"

            if self.config.dry_run:
                self.logger.info(f"[DRY RUN] {self.config.haphic_filter_bam_bin} {raw_bam} 1 --nm 3 --threads {self.config.threads}")
                self.logger.info(f"[DRY RUN] ... | samtools view -b -@ {self.config.threads} -o {filtered_bam}")
            else:
                try:
                    # 构建HapHiC过滤命令 | Build HapHiC filter command
                    filter_cmd = [
                        self.config.haphic_filter_bam_bin,
                        raw_bam,  # 输入BAM文件
                        str(self.config.mapq_threshold),  # MAPQ阈值 (默认1)
                        "--nm", str(self.config.edit_distance),  # 编辑距离 (默认3)
                        "--threads", str(self.config.threads)
                    ]

                    # samtools view命令 | samtools view command
                    samtools_filter_cmd = [
                        self.config.samtools_bin, "view",
                        "-b", "-@", str(self.config.threads),
                        "-o", filtered_bam
                    ]

                    # 创建过滤管道 | Create filtering pipeline
                    filter_proc = subprocess.Popen(filter_cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
                    samtools_filter_proc = subprocess.Popen(samtools_filter_cmd, stdin=filter_proc.stdout, stderr=subprocess.PIPE, text=True)

                    # 关闭管道 | Close pipes
                    filter_proc.stdout.close()

                    # 等待完成 | Wait for completion
                    filter_stdout, filter_stderr = filter_proc.communicate()
                    samtools_filter_stdout, samtools_filter_stderr = samtools_filter_proc.communicate()

                    # 检查返回码 | Check return codes
                    if filter_proc.returncode != 0 or samtools_filter_proc.returncode != 0:
                        error_msg = f"HapHiC过滤失败 | HapHiC filtering failed:\n"
                        error_msg += f"Filter: {filter_stderr}\n" if filter_stderr else ""
                        error_msg += f"Samtools: {samtools_filter_stderr}" if samtools_filter_stderr else ""
                        self.logger.error(error_msg)
                        # 如果过滤失败，使用原始BAM文件 | If filtering fails, use original BAM
                        self.logger.warning("⚠️ 过滤失败，使用原始BAM文件 | Filtering failed, using original BAM file")
                        filtered_bam = raw_bam
                    else:
                        self.logger.info(f"✅ HapHiC过滤完成，生成 {filtered_bam}")

                        # 检查过滤后的文件大小 | Check filtered file size
                        if os.path.exists(filtered_bam):
                            size_mb = os.path.getsize(filtered_bam) / (1024 * 1024)
                            self.logger.info(f"📊 过滤后BAM文件大小: {size_mb:.2f} MB")

                except Exception as e:
                    self.logger.error(f"❌ HapHiC过滤失败: {e}")
                    self.logger.warning("⚠️ 使用原始BAM文件继续 | Continuing with original BAM file")
                    filtered_bam = raw_bam
        else:
            filtered_bam = raw_bam
            self.logger.info("⏭️ 跳过HapHiC过滤步骤 | Skipping HapHiC filtering step")

        # HapHiC不需要BAM索引，直接返回过滤后的文件
        # HapHiC does not need BAM index, directly return filtered file
        return filtered_bam

    def _adjust_threads_for_memory(self):
        """简单调整线程数以避免内存不足 | Simple adjustment of threads to avoid memory issues"""
        # BWA每个线程大约需要1-2GB内存，限制最大线程数
        max_safe_threads = 64  # 安全的最大线程数
        if self.config.threads > max_safe_threads:
            self.logger.warning(f"⚠️ 线程数过多 ({self.config.threads})，自动调整为 {max_safe_threads} 以避免内存不足 | Too many threads ({self.config.threads}), automatically adjusting to {max_safe_threads} to avoid memory issues")
            self.config.threads = max_safe_threads

    def _check_memory_usage(self):
        """检查系统内存使用情况并给出优化建议 | Check system memory usage and provide optimization suggestions"""
        try:
            import psutil
            memory = psutil.virtual_memory()
            available_gb = memory.available / (1024**3)

            self.logger.info(f"💾 系统内存状态 | System Memory Status:")
            self.logger.info(f"  总内存 | Total Memory: {memory.total / (1024**3):.1f} GB")
            self.logger.info(f"  可用内存 | Available Memory: {available_gb:.1f} GB ({memory.percent:.1f}% used)")

            # 建议线程数
            recommended_threads = max(1, min(self.config.threads, int(available_gb / 2)))
            if recommended_threads < self.config.threads:
                self.logger.warning(f"⚠️ 建议减少线程数从 {self.config.threads} 到 {recommended_threads} 以避免内存不足 | Suggest reducing threads from {self.config.threads} to {recommended_threads} to avoid memory issues")
                self.config.threads = recommended_threads

            if available_gb < 4:
                self.logger.warning("⚠️ 可用内存较少，建议关闭其他程序或增加系统内存 | Low available memory, suggest closing other programs or adding system memory")

        except ImportError:
            self.logger.info("💾 无法检查内存使用情况 (psutil未安装) | Cannot check memory usage (psutil not installed)")
        except Exception as e:
            self.logger.warning(f"⚠️ 检查内存使用时出错 | Error checking memory usage: {e}")

    def _ensure_bwa_index(self, asm_file: str):
        """确保BWA索引存在 | Ensure BWA index exists"""
        index_files = [
            asm_file + ".bwt",
            asm_file + ".pac",
            asm_file + ".ann",
            asm_file + ".amb",
            asm_file + ".sa"
        ]

        if all(os.path.exists(f) for f in index_files):
            self.logger.info(f"✅ BWA索引已存在 | BWA index already exists: {asm_file}")
            return

        self.logger.info(f"📝 创建BWA索引 | Creating BWA index: {asm_file}")
        cmd = [self.config.bwa_bin, "index", asm_file]

        if self.config.dry_run:
            self.logger.info(f"[DRY RUN] {' '.join(cmd)}")
        else:
            try:
                result = subprocess.run(cmd, capture_output=True, text=True, check=True)
                self.logger.info("✅ BWA索引创建完成 | BWA index created")
            except subprocess.CalledProcessError as e:
                self.logger.error(f"❌ BWA索引创建失败 | BWA index creation failed: {e.stderr}")
                raise

    def _create_bam_index(self, bam_file: str):
        """创建BAM索引 | Create BAM index"""
        # 检查两种可能的索引文件名
        index_file = bam_file + ".bai"
        index_file_alt = bam_file.rsplit('.', 1)[0] + ".bai"  # 例如 HiC.bam -> HiC.bai

        if os.path.exists(index_file) or os.path.exists(index_file_alt):
            existing_index = index_file if os.path.exists(index_file) else index_file_alt
            self.logger.info(f"✅ BAM索引已存在 | BAM index already exists: {existing_index}")
            return

        self.logger.info(f"📝 创建BAM索引 | Creating BAM index: {bam_file}")

        if self.config.dry_run:
            self.logger.info(f"[DRY RUN] {self.config.samtools_bin} index {bam_file}")
        else:
            try:
                result = subprocess.run([self.config.samtools_bin, "index", bam_file],
                                      capture_output=True, text=True, check=True)
                self.logger.info("✅ BAM索引创建完成 | BAM index created")
            except subprocess.CalledProcessError as e:
                self.logger.warning(f"⚠️ BAM索引创建失败 | BAM index creation failed: {e.stderr}")
                # 检查是否有其他错误信息，比如染色体不连续
                if "Chromosome blocks not continuous" in str(e.stderr):
                    self.logger.warning("⚠️ 检测到BAM文件染色体块不连续的问题，这可能是由于排序问题引起的")
                    self.logger.warning("⚠️ 尝试使用samtools sort重新排序BAM文件 | Attempting to re-sort BAM file with samtools sort")
                    try:
                        sorted_bam = bam_file.replace(".bam", ".sorted.bam")
                        subprocess.run([self.config.samtools_bin, "sort", "-@", str(self.config.threads), "-o", sorted_bam, bam_file],
                                     capture_output=True, text=True, check=True)
                        subprocess.run([self.config.samtools_bin, "index", sorted_bam],
                                     capture_output=True, text=True, check=True)
                        self.logger.info(f"✅ 重新排序并索引完成: {sorted_bam}")
                        return sorted_bam
                    except subprocess.CalledProcessError as sort_error:
                        self.logger.error(f"❌ BAM排序失败 | BAM sorting failed: {sort_error.stderr}")
                # 不是致命错误，可以继续 | Not fatal error, can continue