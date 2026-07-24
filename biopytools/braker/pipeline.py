"""
BRAKER3基因组注释核心流程模块|BRAKER3 Genome Annotation Core Pipeline Module
"""

import os
import subprocess
import fcntl
from pathlib import Path
from .utils import BrakerLogger, CommandRunner, check_file_exists, check_step_completed, format_number, fix_duplicate_gtf_transcript_ids


class BrakerPipeline:
    """BRAKER3基因组注释流程|BRAKER3 Genome Annotation Pipeline"""

    def __init__(self, config, logger=None):
        """
        初始化流程|Initialize pipeline

        Args:
            config: BrakerConfig配置对象|BrakerConfig object
            logger: 日志器|Logger object (optional)
        """
        from .config import BrakerConfig
        if not isinstance(config, BrakerConfig):
            raise TypeError("config必须是BrakerConfig实例|config must be BrakerConfig instance")

        self.config = config
        self.logger = logger

        # 创建日志器|Create logger
        if not self.logger:
            log_file = os.path.join(self.config.log_dir, "braker_pipeline.log")
            logger_manager = BrakerLogger(log_file)
            self.logger = logger_manager.get_logger()

        # 创建命令执行器|Create command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir)

        # 运行状态|Run status
        self.completed_steps = {}
        self._run_lock_fd = None  # 项目级运行锁fd,进程退出自动释放|Project run lock fd, auto-released on process exit
        self._check_completed_steps()

    def _check_completed_steps(self):
        """检查已完成步骤|Check completed steps"""
        # 检查重复序列屏蔽|Check repeat masking
        # 只检查屏蔽后的基因组文件作为完成标志
        # Only check masked genome file as completion marker
        masked_genome = os.path.join(self.config.repeat_dir, f"{Path(self.config.genome).name}.masked")
        if not os.path.exists(masked_genome):
            # 尝试其他可能的命名|Try other possible namings
            possible_names = [
                f"{Path(self.config.genome).name}.masked",
                f"{Path(self.config.genome).stem}.masked.fa",
                f"{Path(self.config.genome).stem}.masked"
            ]
            for name in possible_names:
                candidate = os.path.join(self.config.repeat_dir, name)
                if os.path.exists(candidate):
                    masked_genome = candidate
                    break
        repeat_done = os.path.exists(masked_genome)
        if repeat_done and self.logger:
            self.logger.debug(f"检测到重复序列屏蔽已完成|Detected completed repeat masking: {masked_genome}")
        self.completed_steps['repeat_masking'] = repeat_done

        # 检查三代转录本处理|Check long-read processing
        long_reads_done = check_step_completed(
            self.config.long_reads_dir,
            ["isoseq.sorted.bam"],
            self.logger
        )
        self.completed_steps['long_reads'] = long_reads_done

        # 检查二代RNA-seq处理|Check short-read processing
        # 检查排序后的BAM文件和索引|Check sorted BAM file and index
        short_reads_done = check_step_completed(
            self.config.short_reads_dir,
            ["rnaseq.sorted.bam", "rnaseq.sorted.bam.bai"],
            self.logger
        )
        self.completed_steps['short_reads'] = short_reads_done

        # 检查BRAKER3运行|Check BRAKER3 run
        braker_done = check_step_completed(
            self.config.braker_dir,
            ["braker.gtf", "braker.gff3"],
            self.logger
        )
        self.completed_steps['braker'] = braker_done

        # 检查BRAKER3运行|Check BRAKER3 run

    def _acquire_run_lock(self):
        """
        获取项目级排他锁,防止同一项目并发运行互相踩踏|Acquire project-level exclusive lock to prevent concurrent-run collision

        设计原因|Design reason:
            braker_safe_dir 的名字是 (braker_dir+genome) 的 md5 哈希,确定性且会被复用
            (为支持断点续传 --useexisting)。同一项目若被并发执行两次,两个 braker.pl 进程会
            落入同一个工作目录,在 braker.pl 内部 ln -s traingenes.gtf 等非原子操作上竞争,
            触发 "File exists" 致命错误(braker.pl:6107)。本锁让第二个并发任务在流程入口立即失败,
            而不是白跑数小时后在 BRAKER 内部互踩。flock 锁绑定到打开的文件描述符,进程退出
            (含被 kill)时由内核自动释放,不会产生死锁或残留锁文件问题。

            The braker_safe_dir name is a deterministic md5 hash and is reused across runs
            (to support --useexisting resume). Two concurrent runs of the same project land
            in the same working dir and race on non-atomic BRAKER-internal ops (e.g. the
            braker.pl `ln -s traingenes.gtf`), causing fatal "File exists" errors. This lock
            makes the second concurrent run fail fast at pipeline entry instead of colliding
            hours later. flock is bound to the open fd and auto-released by the kernel on
            process exit (including kill), so no deadlock or stale-lock risk.
        """
        # 锁文件放在 safe_dir 同级(safe_dir 可能被清理重建,锁文件不能随之删除)
        # Lock file is a sibling of safe_dir (safe_dir may be rmtree'd and rebuilt;
        # the lock file must survive that)
        lock_file = self.config.braker_safe_dir + ".lock"
        fd = open(lock_file, "w")
        try:
            # LOCK_NB: 非阻塞,锁被占用立即抛 BlockingIOError 而非挂起等待
            # LOCK_NB: non-blocking; raises BlockingIOError immediately if held
            fcntl.flock(fd, fcntl.LOCK_EX | fcntl.LOCK_NB)
        except BlockingIOError:
            fd.close()
            raise RuntimeError(
                f"该项目正在被另一个进程运行,请等待其完成后再重试|"
                f"This project is already running in another process; wait for it to finish and retry. "
                f"Lock: {lock_file}"
            )
        except OSError:
            fd.close()
            raise  # 其它系统级错误向上抛|Other OS-level errors propagate
        self._run_lock_fd = fd  # 持有fd保持锁,进程退出时自动释放|Hold fd to keep lock, auto-released on exit
        self.logger.info(f"已获取项目运行锁(防止并发冲突)|Acquired project run lock: {lock_file}")

    def run_pipeline(self):
        """运行完整流程|Run complete pipeline"""
        self.logger.info("=" * 80)
        self.logger.info("开始BRAKER3基因组注释流程|Starting BRAKER3 Genome Annotation Pipeline")
        self.logger.info("=" * 80)

        try:
            # 流程入口加项目级锁:防止同一项目并发运行导致 braker_safe_dir 内部互踩
            # Acquire project-level lock at entry: prevent concurrent runs from colliding
            # inside the shared braker_safe_dir
            self._acquire_run_lock()

            # 步骤1: 重复序列屏蔽|Step 1: Repeat masking
            if not self.config.skip_repeat:
                masked_genome = self._step1_repeat_masking()
            else:
                self.logger.info("跳过重复序列屏蔽步骤|Skipping repeat masking step")
                masked_genome = self.config.genome

            # 步骤2: 三代转录本处理|Step 2: Long-read processing
            lr_bam = None
            if self.config.isoseq and not self.config.skip_long_reads:
                lr_bam = self._step2_long_reads(masked_genome)
            elif self.config.isoseq:
                self.logger.info("跳过三代转录本处理步骤|Skipping long-read processing step")

            # 步骤3: 二代RNA-seq处理|Step 3: Short-read processing
            sr_bam = None
            if self.config.rnaseq_dirs and not self.config.skip_short_reads:
                sr_bam = self._step3_short_reads(masked_genome)
            elif self.config.rnaseq_dirs:
                self.logger.info("跳过二代RNA-seq处理步骤|Skipping short-read processing step")

            # 步骤3.5: 证据驱动还原被误mask的真基因区(方案2)|Step 3.5: Rescue
            bam_files = [b for b in [sr_bam, lr_bam] if b]
            if not self.config.skip_rescue:
                from .repeat_refine import rescue_masked_regions
                refined = rescue_masked_regions(
                    masked_genome, self.config.prot_seq, bam_files,
                    self.config.repeat_dir, self.config, self.cmd_runner, self.logger)
                if refined:
                    self.logger.info(f"使用还原后的基因组|Using refined genome: {refined}")
                    masked_genome = refined
            else:
                self.logger.info("跳过证据还原步骤|Skipping rescue step")

            # 步骤4: 运行BRAKER3|Step 4: Run BRAKER3
            braker_results = self._step4_braker_annotation(
                masked_genome, sr_bam, lr_bam
            )

            self.logger.info("=" * 80)
            self.logger.info("BRAKER3基因组注释流程完成|BRAKER3 Genome Annotation Pipeline Completed")
            self.logger.info(f"最终注释文件|Final annotation: {braker_results['gtf']}")
            self.logger.info("=" * 80)

            return braker_results['gtf']

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {e}")
            raise

    def _step1_repeat_masking(self):
        """
        步骤1: 重复序列屏蔽|Step 1: Repeat masking

        Returns:
            str: 屏蔽后的基因组文件路径|Path to masked genome file
        """
        self.logger.info("-" * 80)
        self.logger.info("步骤1: 重复序列屏蔽|Step 1: Repeat Masking")
        self.logger.info("-" * 80)

        if self.completed_steps['repeat_masking']:
            self.logger.info("检测到重复序列屏蔽已完成，跳过|Repeat masking already completed, skipping")
            # 查找屏蔽后的基因组文件|Find masked genome file
            possible_names = [
                f"{Path(self.config.genome).name}.masked",
                f"{Path(self.config.genome).stem}.masked.fa",
                f"{Path(self.config.genome).stem}.masked"
            ]
            for name in possible_names:
                masked_genome = os.path.join(self.config.repeat_dir, name)
                if os.path.exists(masked_genome):
                    self.logger.info(f"使用已有的屏蔽基因组|Using existing masked genome: {masked_genome}")
                    return masked_genome
            # 如果没找到，列出目录中的文件供调试|If not found, list directory for debugging
            self.logger.warning(f"未找到屏蔽后的基因组文件，重新运行|Masked genome not found, re-running")

        # 1.1 运行RepeatModeler|Run RepeatModeler
        self.logger.info("步骤1.1: 运行RepeatModeler构建重复序列库|Step 1.1: Running RepeatModeler")

        # 使用基因组绝对路径，因为命令在output_dir下执行|Use absolute genome path since command runs in output_dir
        genome_abs = os.path.abspath(self.config.genome)
        database_name = f"{self.config.species}_repeats"
        build_cmd = f"{self.config.build_database_bin} -name {database_name} {genome_abs}"

        if not self.cmd_runner.run_command(build_cmd, "构建RepeatModeler数据库|Build RepeatModeler database"):
            raise RuntimeError("RepeatModeler数据库构建失败|RepeatModeler database construction failed")

        repeatmodeler_cmd = f"{self.config.repeatmodeler_bin} -database {database_name} -threads {self.config.threads}"

        if not self.cmd_runner.run_command(repeatmodeler_cmd, "运行RepeatModeler|Run RepeatModeler"):
            raise RuntimeError("RepeatModeler运行失败|RepeatModeler execution failed")

        # 1.2 运行RepeatMasker|Run RepeatMasker
        self.logger.info("步骤1.2: 运行RepeatMasker屏蔽重复序列|Step 1.2: Running RepeatMasker")

        # 查找 RepeatModeler 输出的 consensi.fa.classified 文件
        # Find consensi.fa.classified file output by RepeatModeler
        # 支持动态目录名（如 RM_进程ID.时间戳）
        # Support dynamic directory names (e.g., RM_processid.timestamp)
        # RepeatModeler输出到主输出目录，而非子目录
        # RepeatModeler outputs to main output dir, not subdirectory
        consensi_files = list(Path(self.config.output_dir).glob("RM_*/consensi.fa.classified"))

        # 如果没找到，尝试递归搜索
        # If not found, try recursive search
        if not consensi_files:
            self.logger.debug("未在 RM_*/ 目录中找到，尝试递归搜索|Not found in RM_*/ directories, trying recursive search")
            consensi_files = list(Path(self.config.output_dir).rglob("*/consensi.fa.classified"))

        if not consensi_files:
            raise RuntimeError(f"找不到 RepeatModeler 输出的 consensi.fa.classified 文件|Cannot find consensi.fa.classified from RepeatModeler in {self.config.output_dir}")

        consensi_fa = os.path.abspath(str(consensi_files[0]))
        self.logger.info(f"使用重复序列库|Using repeat library: {consensi_fa}")

        # repeat 库过滤:剔除被误判的基因家族(方案1)|Filter repeat library (scheme 1)
        lib_for_masking = consensi_fa
        if not self.config.skip_repeat_filter:
            from .repeat_refine import filter_repeat_library
            filtered = filter_repeat_library(
                consensi_fa, self.config.repeat_dir, self.config,
                self.cmd_runner, self.logger)
            if filtered:
                lib_for_masking = filtered
                self.logger.info(f"RepeatMasker 使用过滤后的库|Using filtered library: {filtered}")
            else:
                self.logger.warning("repeat 库过滤未产出结果,回退原库|Filter produced nothing, fallback to original library")
        else:
            self.logger.info("跳过 repeat 库过滤|Skipping repeat library filtering")

        mask_option = "-xsmall" if self.config.soft_masking else ""

        repeat_dir_abs = os.path.abspath(self.config.repeat_dir)
        repeatmasker_cmd = (
            f"{self.config.repeatmasker_bin} -lib {lib_for_masking} "
            f"{mask_option} -dir {repeat_dir_abs} -pa {self.config.threads} "
            f"{genome_abs}"
        )

        if not self.cmd_runner.run_command(repeatmasker_cmd, "运行RepeatMasker|Run RepeatMasker"):
            raise RuntimeError("RepeatMasker运行失败|RepeatMasker execution failed")

        # 查找屏蔽后的基因组|Find masked genome
        masked_genome = os.path.join(self.config.repeat_dir, f"{Path(self.config.genome).name}.masked")
        if not os.path.exists(masked_genome):
            # 尝试其他可能的命名|Try other possible namings
            possible_names = [
                f"{Path(self.config.genome).name}.masked",
                f"{Path(self.config.genome).stem}.masked.fa",
                f"{Path(self.config.genome).stem}.masked"
            ]
            for name in possible_names:
                candidate = os.path.join(self.config.repeat_dir, name)
                if os.path.exists(candidate):
                    masked_genome = candidate
                    break
            else:
                raise RuntimeError(f"找不到屏蔽后的基因组文件|Cannot find masked genome file in {self.config.repeat_dir}")

        self.logger.info(f"重复序列屏蔽完成|Repeat masking completed: {masked_genome}")
        return masked_genome

    def _step2_long_reads(self, genome):
        """
        步骤2: 三代转录本处理|Step 2: Long-read processing

        Args:
            genome: 基因组文件|Genome file

        Returns:
            str: sorted BAM文件路径|Path to sorted BAM file
        """
        self.logger.info("-" * 80)
        self.logger.info("步骤2: 三代全长转录本处理|Step 2: Long-read Transcript Processing")
        self.logger.info("-" * 80)

        if self.completed_steps['long_reads']:
            self.logger.info("检测到三代转录本处理已完成，跳过|Long-read processing already completed, skipping")
            bam_file = os.path.join(self.config.long_reads_dir, "isoseq.sorted.bam")
            if os.path.exists(bam_file):
                return bam_file

        # 确保路径为绝对路径，因为命令在output_dir下执行|Ensure absolute paths since command runs in output_dir
        genome = os.path.abspath(genome)
        isoseq = os.path.abspath(self.config.isoseq)

        # 使用minimap2比对|Align with minimap2
        sam_out = os.path.abspath(os.path.join(self.config.long_reads_dir, "isoseq.sam"))
        minimap_cmd = (
            f"{self.config.minimap2_bin} -ax splice -uf -C5 "
            f"-t {self.config.threads} {genome} {isoseq} > {sam_out}"
        )

        if not self.cmd_runner.run_command(minimap_cmd, "三代转录本比对|Long-read alignment"):
            raise RuntimeError("三代转录本比对失败|Long-read alignment failed")

        # 将SAM转为sorted BAM|Convert SAM to sorted BAM
        bam_out = os.path.abspath(os.path.join(self.config.long_reads_dir, "isoseq.sorted.bam"))
        view_cmd = f"samtools view -bS {sam_out} > {bam_out.replace('.sorted.bam', '.bam')}"

        if not self.cmd_runner.run_command(view_cmd, "SAM转BAM|Convert SAM to BAM"):
            raise RuntimeError("SAM转BAM失败|SAM to BAM conversion failed")

        unsorted_bam = bam_out.replace('.sorted.bam', '.bam')
        sort_cmd = f"samtools sort -@ {self.config.threads} -o {bam_out} {unsorted_bam}"

        if not self.cmd_runner.run_command(sort_cmd, "排序三代BAM|Sort long-read BAM"):
            raise RuntimeError("三代BAM排序失败|Long-read BAM sorting failed")

        index_cmd = f"samtools index {bam_out}"

        if not self.cmd_runner.run_command(index_cmd, "建立三代BAM索引|Index long-read BAM"):
            raise RuntimeError("三代BAM索引建立失败|Long-read BAM index creation failed")

        # 清理中间文件|Clean up intermediate files
        os.remove(unsorted_bam)
        os.remove(sam_out)

        self.logger.info(f"三代转录本处理完成|Long-read processing completed: {bam_out}")

        return bam_out

    def _step3_short_reads(self, genome):
        """
        步骤3: 二代RNA-seq处理|Step 3: Short-read processing

        Args:
            genome: 基因组文件|Genome file

        Returns:
            str: BAM文件路径|Path to BAM file
        """
        self.logger.info("-" * 80)
        self.logger.info("步骤3: 二代RNA-seq数据处理|Step 3: Short-read RNA-seq Processing")
        self.logger.info("-" * 80)

        if self.completed_steps['short_reads']:
            self.logger.info("检测到二代RNA-seq处理已完成，跳过|Short-read processing already completed, skipping")
            # 使用排序后的BAM文件|Use sorted BAM file
            bam_file = os.path.join(self.config.short_reads_dir, "rnaseq.sorted.bam")
            if os.path.exists(bam_file):
                self.logger.info(f"使用已有的RNA-seq比对文件|Using existing RNA-seq alignment: {bam_file}")
                return bam_file

        # 导入文件识别函数|Import file finder functions
        from .utils import find_rnaseq_files_in_directory

        # 确保路径为绝对路径，因为命令在output_dir下执行|Ensure absolute paths since command runs in output_dir
        genome = os.path.abspath(genome)

        # 构建HISAT2索引|Build HISAT2 index
        index_prefix = os.path.abspath(os.path.join(self.config.short_reads_dir, self.config.species))
        hisat_build_cmd = f"{self.config.hisat2_build_bin} {genome} {index_prefix}"

        if not self.cmd_runner.run_command(hisat_build_cmd, "构建HISAT2索引|Build HISAT2 index"):
            raise RuntimeError("HISAT2索引构建失败|HISAT2 index construction failed")

        # 处理所有RNA-seq目录|Process all RNA-seq directories
        all_paired_files = []
        for rnaseq_dir in self.config.rnaseq_dirs:
            self.logger.info(f"处理RNA-seq目录|Processing RNA-seq directory: {rnaseq_dir}")
            paired_files = find_rnaseq_files_in_directory(
                rnaseq_dir,
                read1_pattern=self.config.read1_pattern,
                read2_pattern=self.config.read2_pattern,
                logger=self.logger
            )
            all_paired_files.extend(paired_files)

        if not all_paired_files:
            raise RuntimeError(f"找不到RNA-seq文件|Cannot find RNA-seq files in any directory")

        self.logger.info(f"共找到 {len(all_paired_files)} 对RNA-seq文件|Found {len(all_paired_files)} pairs of RNA-seq files")

        # 合并所有R1和R2文件列表|Merge all R1 and R2 file lists
        r1_files = [os.path.abspath(f[0]) for f in all_paired_files]
        r2_files = [os.path.abspath(f[1]) for f in all_paired_files]

        # 当有多对FASTQ文件时，先合并为单文件再传给HISAT2
        # HISAT2多文件输入时内部在/tmp创建命名管道(mkfifo)，不尊重TMPDIR，
        # 某些计算节点/tmp不支持创建FIFO导致失败
        # 当只有一对文件时跳过合并
        # When multiple FASTQ pairs, merge into single files before passing to HISAT2
        # HISAT2 creates named pipes in /tmp for multi-file input, ignores TMPDIR,
        # some compute nodes don't support FIFO creation in /tmp
        # Skip merging when there's only one pair
        if len(r1_files) > 1:
            merged_r1 = os.path.abspath(os.path.join(self.config.short_reads_dir, "merged_R1.fq.gz"))
            merged_r2 = os.path.abspath(os.path.join(self.config.short_reads_dir, "merged_R2.fq.gz"))

            if os.path.exists(merged_r1) and os.path.exists(merged_r2):
                self.logger.info(
                    f"使用已有的合并文件|Using existing merged files: "
                    f"{os.path.basename(merged_r1)}, {os.path.basename(merged_r2)}"
                )
                r1_input = merged_r1
                r2_input = merged_r2
            else:
                self.logger.info(
                    f"合并 {len(r1_files)} 对FASTQ文件|Concatenating {len(r1_files)} FASTQ pairs"
                )
                # 合并R1|Merge R1
                cat_r1_cmd = f"cat {' '.join(r1_files)} > {merged_r1}"
                if not self.cmd_runner.run_command(cat_r1_cmd, "合并R1文件|Merge R1 files"):
                    raise RuntimeError("R1文件合并失败|Failed to merge R1 files")
                # 合并R2|Merge R2
                cat_r2_cmd = f"cat {' '.join(r2_files)} > {merged_r2}"
                if not self.cmd_runner.run_command(cat_r2_cmd, "合并R2文件|Merge R2 files"):
                    raise RuntimeError("R2文件合并失败|Failed to merge R2 files")
                r1_input = merged_r1
                r2_input = merged_r2
        else:
            r1_input = r1_files[0]
            r2_input = r2_files[0]

        # HISAT2比对（支持多个样本）
        # 拆分为两步执行，避免管道缓冲区溢出导致samtools解析错误
        # HISAT2 alignment (support multiple samples)
        # Split into two steps to avoid pipe buffer overflow causing samtools parse errors
        sam_file = os.path.abspath(os.path.join(self.config.short_reads_dir, "rnaseq.sam"))

        hisat_align_cmd = (
            f"{self.config.hisat2_bin} -x {index_prefix} "
            f"-1 {r1_input} -2 {r2_input} "
            f"--rna-strandness RF -p {self.config.threads} -S {sam_file}"
        )

        if not self.cmd_runner.run_command(hisat_align_cmd, "HISAT2比对|HISAT2 alignment"):
            raise RuntimeError("HISAT2比对失败|HISAT2 alignment failed")

        # SAM转BAM并排序|Convert SAM to sorted BAM
        sorted_bam = os.path.abspath(os.path.join(self.config.short_reads_dir, "rnaseq.sorted.bam"))
        sort_cmd = f"samtools sort -@ {self.config.threads} -o {sorted_bam} {sam_file}"

        if not self.cmd_runner.run_command(sort_cmd, "SAM转排序BAM|Sort SAM to BAM"):
            raise RuntimeError("BAM文件排序失败|BAM file sorting failed")

        # 建立索引|Create index
        index_cmd = f"samtools index {sorted_bam}"

        if not self.cmd_runner.run_command(index_cmd, "建立BAM索引|Create BAM index"):
            raise RuntimeError("BAM索引建立失败|BAM index creation failed")

        # 清理中间SAM文件|Clean up intermediate SAM file
        if os.path.exists(sam_file):
            os.remove(sam_file)

        self.logger.info(f"二代RNA-seq处理完成|Short-read processing completed: {sorted_bam}")

        return sorted_bam

    def _step4_braker_annotation(self, genome, bam_file, lr_bam):
        """
        步骤4: 运行BRAKER3注释|Step 4: Run BRAKER3 annotation

        Args:
            genome: 基因组文件|Genome file
            bam_file: 二代RNA-seq BAM文件|Short-read RNA-seq BAM file
            lr_bam: 三代转录本BAM文件|Long-read BAM file

        Returns:
            dict: BRAKER结果路径|BRAKER result paths
        """
        self.logger.info("-" * 80)
        self.logger.info("步骤4: BRAKER3基因组注释|Step 4: BRAKER3 Genome Annotation")
        self.logger.info("-" * 80)

        if self.completed_steps['braker']:
            self.logger.info("检测到BRAKER3注释已完成，跳过|BRAKER3 annotation already completed, skipping")
            return {
                'gtf': os.path.join(self.config.braker_dir, "braker.gtf"),
                'aa': os.path.join(self.config.braker_dir, "braker.aa"),
                'hints': os.path.join(self.config.braker_dir, "hintsfile.gff")
            }

        # 构建BRAKER命令|Build BRAKER command
        braker_cmd_parts = []

        # Singularity执行|Singularity execution
        if self.config.use_singularity:
            # 设置 Singularity 临时目录为工作目录，避免 /tmp 空间不足
            # Set Singularity temporary directory to work dir to avoid /tmp space issues
            # 注意：Singularity临时目录必须使用绝对路径，否则会报错
            # Note: Singularity temp directory MUST use absolute path, otherwise will fail
            # 获取原工作目录并拼接绝对路径|Get original working directory and join for absolute path
            original_cwd = getattr(self.config, 'working_dir', None)
            if original_cwd and not os.path.isabs(self.config.braker_dir):
                # braker_dir是相对路径，需要拼接原工作目录得到绝对路径
                # braker_dir is relative, need to join with original working directory for absolute path
                singularity_tmpdir = os.path.join(original_cwd, self.config.braker_dir, "singularity_tmp")
            else:
                # braker_dir已经是绝对路径，直接使用
                # braker_dir is already absolute, use directly
                singularity_tmpdir = os.path.join(self.config.braker_dir, "singularity_tmp")

            singularity_tmpdir = os.path.abspath(singularity_tmpdir)
            os.makedirs(singularity_tmpdir, exist_ok=True)
            singularity_env = f"SINGULARITY_TMPDIR={singularity_tmpdir} SINGULARITY_CACHEDIR={singularity_tmpdir}"
            # 创建 AUGUSTUS 配置目录用于挂载|Create AUGUSTUS config directory for mounting
            # 需要挂载整个 config 目录，因为 species 目录依赖其他配置文件
            # Need to mount entire config directory because species depends on other config files
            # AUGUSTUS配置目录也需要绝对路径用于bind mount
            # AUGUSTUS config directory also needs absolute path for bind mount
            if original_cwd and not os.path.isabs(self.config.braker_dir):
                augustus_config_dir = os.path.join(original_cwd, self.config.braker_dir, "augustus_config")
            else:
                augustus_config_dir = os.path.join(self.config.braker_dir, "augustus_config")

            augustus_config_dir = os.path.abspath(augustus_config_dir)
            os.makedirs(augustus_config_dir, exist_ok=True)

            # 首先从容器复制完整的配置结构（如果尚未复制）
            # First copy complete config structure from container (if not yet copied)
            # 检查关键目录是否存在来判断是否需要复制
            # Check if key directories exist to determine if copy is needed
            key_dirs = ['cgp', 'extrinsic', 'model', 'parameters', 'profile', 'species']
            needs_copy = False
            for key_dir in key_dirs:
                if not os.path.exists(os.path.join(augustus_config_dir, key_dir)):
                    needs_copy = True
                    break

            if needs_copy:
                self.logger.info(f"从容器复制 AUGUSTUS 配置|Copying AUGUSTUS config from container to {augustus_config_dir}")
                # 使用 singularity exec 从容器复制配置文件
                # Use singularity exec to copy config files from container
                # 使用 bash -c 来正确处理路径和通配符
                # Use bash -c to properly handle paths and wildcards
                copy_cmd = f"{self.config.singularity_bin} exec {self.config.singularity_image} bash -c 'cp -r /opt/Augustus/config/. {augustus_config_dir}/'"
                result = subprocess.run(copy_cmd, shell=True, capture_output=True, text=True)
                if result.returncode != 0:
                    self.logger.error(f"复制配置失败|Failed to copy config: {result.stderr}")
                    raise RuntimeError(f"无法从容器复制 AUGUSTUS 配置|Cannot copy AUGUSTUS config from container")

            self.logger.debug(f"使用 AUGUSTUS 配置目录|Using AUGUSTUS config directory: {augustus_config_dir}")

            # 添加 singularity 执行命令，绑定挂载整个 AUGUSTUS 配置目录
            # Add singularity exec command with bind mount for entire AUGUSTUS config directory
            braker_cmd_parts.append(f"env {singularity_env} {self.config.singularity_bin} exec")
            braker_cmd_parts.append(f"--bind {augustus_config_dir}:/opt/Augustus/config")

            # 如果使用了相对路径（说明路径包含中文），需要挂载原工作目录
            # If relative paths are used (means path contains Chinese), need to mount original working directory
            original_cwd = getattr(self.config, 'working_dir', None)

            # 收集所有需要挂载的目录
            # Collect all directories that need to be mounted
            dirs_to_mount = set()

            if original_cwd:
                dirs_to_mount.add(original_cwd)

            # 添加符号链接的目标目录（~/tmp也需要挂载以让容器能访问）
            # Add symlink target directories (~/tmp also needs to be mounted for container access)
            user_tmp = os.path.join(os.path.expanduser("~"), "tmp")
            dirs_to_mount.add(user_tmp)

            # 挂载所有必要的目录
            # Mount all necessary directories
            for mount_dir in dirs_to_mount:
                self.logger.info(f"挂载目录到容器|Mounting directory to container: {mount_dir}")
                braker_cmd_parts.append(f"--bind {mount_dir}:{mount_dir}")

            braker_cmd_parts.append(self.config.singularity_image)

        braker_cmd_parts.append(f"perl {self.config.braker_in_image}")

        # 转换所有输入路径为安全绝对路径（通过符号链接）
        # Convert all input paths to safe absolute paths (via symlinks)
        from .config import get_safe_absolute_path

        # genome路径需要安全绝对路径
        # genome path needs safe absolute path
        genome_safe = get_safe_absolute_path(genome, getattr(self.config, 'working_dir', None))
        braker_cmd_parts.append(f"--genome={genome_safe}")

        braker_cmd_parts.append(f"--species={self.config.species}")
        braker_cmd_parts.append(f"--threads={self.config.threads}")

        # 添加证据数据|Add evidence data
        # 添加二代RNA-seq BAM文件|Add short-read RNA-seq BAM
        if bam_file:
            bam_safe = get_safe_absolute_path(bam_file, getattr(self.config, 'working_dir', None))
            braker_cmd_parts.append(f"--bam={bam_safe}")

        # 添加三代转录本比对文件|Add long-read alignment file
        # BRAKER3 接受额外的 BAM 文件作为证据
        # BRAKER3 accepts additional BAM files as evidence
        if lr_bam:
            lr_bam_safe = get_safe_absolute_path(lr_bam, getattr(self.config, 'working_dir', None))
            braker_cmd_parts.append(f"--bam={lr_bam_safe}")

        # 添加蛋白质序列|Add protein sequences
        if self.config.prot_seq:
            prot_safe = get_safe_absolute_path(self.config.prot_seq, getattr(self.config, 'working_dir', None))
            braker_cmd_parts.append(f"--prot_seq={prot_safe}")

        # 添加参数|Add parameters
        if self.config.use_fungus:
            braker_cmd_parts.append("--fungus")

        if self.config.soft_masking:
            braker_cmd_parts.append("--softmasking")

        # 输出GFF3格式文件|Output GFF3 format file
        braker_cmd_parts.append("--gff3")

        # --useexisting策略：
        # 1. braker.gtf存在 → 上次成功完成，使用--useexisting跳过所有步骤
        # 2. braker.gtf不存在但GeneMark-ET anchored GTF存在 → 上次TSEBRA失败，
        #    修复anchored GTF中的重复transcript_id后使用--useexisting，让BRAKER3
        #    跳过GeneMark-ET训练，仅重新执行TSEBRA合并
        # 3. braker.gtf和anchored GTF都不存在 → 从头运行
        # --useexisting strategy:
        # 1. braker.gtf exists → previous run completed, use --useexisting to skip all
        # 2. braker.gtf missing but anchored GTF exists → previous TSEBRA failed,
        #    fix duplicate transcript_ids in anchored GTF, use --useexisting to skip
        #    GeneMark-ET training, only re-run TSEBRA merge
        # 3. braker.gtf and anchored GTF both missing → run from scratch
        braker_gtf = os.path.join(self.config.braker_safe_dir, "braker.gtf")
        anchored_gtf = os.path.join(
            self.config.braker_safe_dir, "GeneMark-ET", "genemark.f.multi_anchored.gtf"
        )
        if os.path.exists(braker_gtf):
            # safe_dir 有 braker.gtf 但本步仍执行 → 04_braker_dir/braker.gtf 不存在(用户删了目标目录重跑)
            # → 上次完成残留,genome 可能已变,清理 safe_dir 全新运行,避免 --useexisting 用旧状态
            # safe_dir has braker.gtf but target missing => stale residue, fresh run
            import shutil
            shutil.rmtree(self.config.braker_safe_dir)
            os.makedirs(self.config.braker_safe_dir, exist_ok=True)
            self.logger.info("清理 braker_safe_dir 上次完成残留,全新运行|Cleaned stale safe_dir, fresh run")
        elif os.path.exists(anchored_gtf):
            # 上次TSEBRA阶段失败，GeneMark-ET已完成但合并失败
            # Previous run failed at TSEBRA, GeneMark-ET completed but merge failed
            # 先修复anchored GTF中的重复transcript_id
            # First fix duplicate transcript_ids in anchored GTF
            self.logger.info(
                "检测到上次运行在TSEBRA阶段失败|Detected previous run failed at TSEBRA stage"
            )
            import shutil
            # 清理旧的AUGUSTUS species目录，避免BRAKER3报错
            # Clean old AUGUSTUS species dir to avoid BRAKER3 error
            if self.config.use_singularity:
                species_dir = os.path.join(augustus_config_dir, "species", self.config.species)
                if os.path.exists(species_dir):
                    shutil.rmtree(species_dir)
                    self.logger.info(f"清理旧species目录|Cleaned stale species dir: {species_dir}")
            braker_cmd_parts.append("--useexisting")
        else:
            # 从头运行，无需清理
            # Run from scratch, no cleanup needed
            pass

        # 使用安全的工作目录路径（不含非ASCII字符）
        # Use safe working directory path (no non-ASCII characters)
        braker_cmd_parts.append(f"--workingdir={self.config.braker_safe_dir}")

        braker_cmd = " ".join(braker_cmd_parts)

        # 在运行BRAKER3之前，检查并修复GeneMark-ET GTF中的重复transcript_id
        # GeneMark-ET可能产生重复transcript_id导致TSEBRA合并失败
        # Before running BRAKER3, check and fix duplicate transcript_ids in GeneMark-ET GTF
        # GeneMark-ET may produce duplicate transcript_ids causing TSEBRA merge failure
        genemark_gtf_files = [
            os.path.join(self.config.braker_safe_dir, "GeneMark-ET", "genemark.f.multi_anchored.gtf"),
            os.path.join(self.config.braker_safe_dir, "GeneMark-ET", "genemark.f.single_anchored.gtf"),
        ]
        for gtf in genemark_gtf_files:
            if os.path.exists(gtf) and os.path.getsize(gtf) > 0:
                fix_duplicate_gtf_transcript_ids(gtf, self.logger)

        # 清理上次运行遗留的符号链接，避免BRAKER3内部ln -s失败
        # Clean stale symlinks from previous run to avoid BRAKER3 ln -s failure
        stale_symlinks = ["traingenes.gtf"]
        for link_name in stale_symlinks:
            link_path = os.path.join(self.config.braker_safe_dir, link_name)
            if os.path.islink(link_path):
                os.remove(link_path)
                self.logger.info(f"清理旧符号链接|Cleaned stale symlink: {link_name}")

        # 不在原工作目录执行，让BRAKER3使用安全路径
        # Don't execute in original working directory, let BRAKER3 use safe path
        # 这样BRAKER3内部生成的所有路径都是安全路径
        # This way all paths generated inside BRAKER3 are safe paths
        needs_copyback = getattr(self.config, 'braker_needs_copyback', False)
        if needs_copyback:
            self.logger.info(f"使用安全工作目录（临时目录）|Using safe working directory (temp): {self.config.braker_safe_dir}")
            self.logger.info(f"运行完成后将结果复制到|Will copy results to: {self.config.braker_dir}")
        else:
            self.logger.info(f"使用工作目录|Using working directory: {self.config.braker_safe_dir}")

        if not self.cmd_runner.run_singularity_command(braker_cmd, "运行BRAKER3|Run BRAKER3"):
            raise RuntimeError("BRAKER3运行失败|BRAKER3 execution failed")

        # 如果使用了临时目录，需要将结果复制回原始目录
        # If used temporary directory, need to copy results back to original directory
        if needs_copyback:
            self.logger.info("复制BRAKER3结果到原始输出目录|Copying BRAKER3 results to original output directory")
            import shutil

            # 定义需要复制的关键结果文件和目录|Define key result files and dirs to copy
            key_files = [
                "braker.gtf",
                "braker.gff3",
                "braker.aa",
                "braker.codingseq",
                "braker.log",
                "hintsfile.gff",
            ]

            safe_dir = self.config.braker_safe_dir
            target_dir = self.config.braker_dir

            # 复制关键文件（逐个复制，容错处理）|Copy key files one by one with error handling
            for filename in key_files:
                src = os.path.join(safe_dir, filename)
                dst = os.path.join(target_dir, filename)
                if os.path.exists(src):
                    try:
                        shutil.copy2(src, dst)
                        self.logger.debug(f"复制|Copied: {filename}")
                    except Exception as e:
                        self.logger.warning(f"复制文件失败|Failed to copy file {filename}: {e}")
                else:
                    self.logger.warning(f"临时目录中缺少文件|Missing file in temp dir: {filename}")

            # 复制关键子目录（Augustus、GeneMark-ETP），逐项容错
            # Copy key subdirectories (Augustus, GeneMark-ETP) with per-item error handling
            key_subdirs = ["Augustus", "GeneMark-ETP"]
            for subdir in key_subdirs:
                src_subdir = os.path.join(safe_dir, subdir)
                dst_subdir = os.path.join(target_dir, subdir)
                if not os.path.exists(src_subdir):
                    continue
                try:
                    if os.path.exists(dst_subdir):
                        shutil.rmtree(dst_subdir)
                    shutil.copytree(src_subdir, dst_subdir)
                    self.logger.debug(f"复制目录|Copied dir: {subdir}/")
                except Exception as e:
                    self.logger.warning(f"复制目录失败|Failed to copy dir {subdir}/: {e}")
                    # 尝试逐项复制以跳过无效符号链接|Try per-item copy to skip broken symlinks
                    try:
                        os.makedirs(dst_subdir, exist_ok=True)
                        for root, dirs, files in os.walk(src_subdir):
                            rel_root = os.path.relpath(root, src_subdir)
                            dst_root = os.path.join(dst_subdir, rel_root)
                            for d in dirs:
                                src_d = os.path.join(root, d)
                                dst_d = os.path.join(dst_root, d)
                                if os.path.exists(src_d) and not os.path.islink(src_d):
                                    os.makedirs(dst_d, exist_ok=True)
                            for f in files:
                                src_f = os.path.join(root, f)
                                dst_f = os.path.join(dst_root, f)
                                if os.path.exists(src_f) and not os.path.islink(src_f):
                                    try:
                                        shutil.copy2(src_f, dst_f)
                                    except Exception:
                                        pass
                        self.logger.info(f"部分复制目录完成|Partially copied dir: {subdir}/")
                    except Exception as e2:
                        self.logger.error(f"目录复制完全失败|Dir copy completely failed: {subdir}/: {e2}")

            self.logger.info(f"结果复制完成|Results copied to: {target_dir}")

        # 检查输出文件|Check output files
        gtf_file = os.path.join(self.config.braker_dir, "braker.gtf")
        gff3_file = os.path.join(self.config.braker_dir, "braker.gff3")
        aa_file = os.path.join(self.config.braker_dir, "braker.aa")
        hints_file = os.path.join(self.config.braker_dir, "hintsfile.gff")

        if not os.path.exists(gtf_file):
            # 如果原始目录没有文件，检查临时目录
            # If no files in original directory, check temp directory
            if needs_copyback:
                gtf_file_temp = os.path.join(self.config.braker_safe_dir, "braker.gtf")
                if os.path.exists(gtf_file_temp):
                    self.logger.warning(f"使用临时目录的结果文件|Using result file from temp directory: {gtf_file_temp}")
                    gtf_file = gtf_file_temp
                    gff3_file = os.path.join(self.config.braker_safe_dir, "braker.gff3")
                    aa_file = os.path.join(self.config.braker_safe_dir, "braker.aa")
                    hints_file = os.path.join(self.config.braker_safe_dir, "hintsfile.gff")
                else:
                    raise RuntimeError(f"BRAKER3输出文件不存在|BRAKER3 output file not found")
            else:
                raise RuntimeError(f"BRAKER3输出文件不存在|BRAKER3 output file not found: {gtf_file}")

        self.logger.info(f"BRAKER3注释完成|BRAKER3 annotation completed: {gtf_file}")
        self.logger.info(f"GFF3格式文件|GFF3 format file: {gff3_file}")

        return {
            'gtf': gtf_file,
            'gff3': gff3_file,
            'aa': aa_file,
            'hints': hints_file
        }
