"""
YaHS流程步骤实现|YaHS Pipeline Steps Implementation

实现Hi-C scaffolding的各个步骤
Implements individual steps for Hi-C scaffolding pipeline
"""

import os
from pathlib import Path
from typing import Optional
from .utils import CommandRunner, check_file_exists


class YaHSSteps:
    """
    YaHS流程步骤执行器|YaHS Pipeline Steps Executor

    封装Hi-C scaffolding的6个主要步骤
    Encapsulates the 6 main steps of Hi-C scaffolding
    """

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        """
        初始化步骤执行器|Initialize steps executor

        Args:
            config: 配置对象|Configuration object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def step_1_indexing(self) -> bool:
        """
        步骤1: 构建基因组索引|Step 1: Build genome index

        为参考基因组构建BWA和SAMtools索引
        Build BWA and SAMtools indices for reference genome

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 1/6: 构建基因组索引|Step 1/6: Building genome index")
        self.logger.info("=" * 60)

        step_dir = self.config.step1_dir
        step_dir.mkdir(parents=True, exist_ok=True)

        # 在输出目录中创建参考基因组的副本|Create reference genome copy in output directory
        ref_fa = str(Path(self.config.ref_fa).resolve())  # 原始参考基因组路径|Original reference path
        ref_name = Path(ref_fa).name  # 使用完整文件名（包含.fa）|Use full filename (including .fa)
        local_ref = str(step_dir / ref_name)  # 输出目录中的参考基因组|Local reference in output dir

        # 创建本地参考基因组（硬链接或副本）|Create local reference (hardlink or copy)
        # 注意：不能使用软链接，因为BWA会解析软链接并在原始位置查找索引文件
        # Note: Cannot use symlink because BWA resolves it and looks for index in original location
        local_ref_path = Path(local_ref)
        if not local_ref_path.exists():
            import shutil
            try:
                # 尝试创建硬链接（节省空间，且BWA可以正确找到索引）|Try hardlink (saves space, BWA can find index)
                local_ref_path.hardlink_to(ref_fa)
                self.logger.info(f"创建参考基因组硬链接|Created hardlink: {ref_name}")
            except (OSError, AttributeError):
                # 如果硬链接失败（AttributeError是Python<3.10没有hardlink_to方法），复制文件
                # If hardlink fails (AttributeError for Python<3.10 without hardlink_to), copy file
                shutil.copy2(ref_fa, local_ref)
                self.logger.info(f"复制参考基因组|Copied reference: {ref_name}")
        else:
            self.logger.info(f"本地参考基因组已存在|Local reference exists: {ref_name}")

        # 检查BWA索引是否已存在|Check if BWA index already exists
        if check_file_exists(str(step_dir / f'{ref_name}.bwt')):
            self.logger.info("BWA索引已存在，跳过|BWA index exists, skipping")
        else:
            # 在本地副本上构建BWA索引|Build BWA index on local copy
            self.logger.info(f"构建BWA索引|Building BWA index: {ref_name}")

            # 使用绝对路径，避免conda run工作目录问题
            # Use absolute path to avoid conda run working directory issues
            local_ref_abs = str(Path(local_ref).resolve())

            cmd = [
                self.config.bwa_bin,
                'index',
                local_ref_abs  # 使用绝对路径|Use absolute path
            ]

            success, _, stderr, _ = self.cmd_runner.run_command(
                cmd,
                description="BWA索引|BWA indexing"
            )

            if not success:
                self.logger.error(f"BWA索引构建失败|BWA index failed")
                return False

            self.logger.info("BWA索引构建完成|BWA index completed")

        # 检查SAMtools索引|Check SAMtools index
        if check_file_exists(str(step_dir / f'{ref_name}.fai')):
            self.logger.info("SAMtools索引已存在，跳过|SAMtools index exists, skipping")
        else:
            self.logger.info(f"构建SAMtools索引|Building SAMtools index")

            # 使用绝对路径|Use absolute path
            local_ref_abs = str(Path(local_ref).resolve())

            cmd = [self.config.samtools_bin, 'faidx', local_ref_abs]  # 使用绝对路径|Use absolute path

            success, stdout, stderr, _ = self.cmd_runner.run_command(
                cmd,
                description="SAMtools索引|SAMtools indexing"
            )

            if not success:
                self.logger.error(f"SAMtools索引构建失败|SAMtools index failed")
                return False

            self.logger.info("SAMtools索引构建完成|SAMtools index completed")

        self.logger.info("步骤 1 完成|Step 1 completed")
        return True

    def step_2_mapping(self) -> bool:
        """
        步骤2: Hi-C数据比对|Step 2: Hi-C data mapping

        使用BWA MEM进行比对，然后排序和去重
        Align using BWA MEM, then sort and mark duplicates

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 2/6: Hi-C数据比对|Step 2/6: Hi-C data mapping")
        self.logger.info("=" * 60)

        # 检查是否已完成
        if check_file_exists(str(self.config.final_bam), 1000000):
            self.logger.info("BAM文件已存在，跳过|BAM file exists, skipping")
            return True

        # 使用绝对路径|Use absolute paths
        ref_fa_orig = str(Path(self.config.ref_fa).resolve())
        hic_r1 = str(Path(self.config.hic_r1).resolve())
        hic_r2 = str(Path(self.config.hic_r2).resolve())
        step2_dir = self.config.step2_dir.resolve()  # 解析为绝对路径|Resolve to absolute path

        # 使用步骤1中创建的本地参考基因组副本|Use local reference copy from step 1
        ref_name = Path(ref_fa_orig).name
        ref_fa = str(Path(self.config.step1_dir / ref_name).resolve())  # 使用绝对路径|Use absolute path

        # 创建临时目录|Create temp directories
        tmp_nsort = step2_dir / 'tmp_nsort'
        tmp_sort = step2_dir / 'tmp_sort'
        tmp_nsort.mkdir(exist_ok=True)
        tmp_sort.mkdir(exist_ok=True)

        # 所有文件路径都使用绝对路径|All file paths use absolute paths
        aligned_bam = step2_dir / 'aligned.bam'
        aligned_nsorted = step2_dir / 'aligned_nsorted.bam'
        aligned_fixmate = step2_dir / 'aligned_fixmate.bam'
        aligned_sorted = step2_dir / 'aligned_sorted.bam'

        try:
            # 2.1 BWA MEM比对|BWA MEM alignment
            self.logger.info("[1/6] BWA MEM比对|BWA MEM alignment")

            # 先输出到SAM（文本格式），避免conda run | conda run管道问题
            # First output to SAM (text format) to avoid conda run | conda run pipe issue
            temp_sam = step2_dir / 'aligned.sam'

            # 检查是否需要运行BWA（检查最终产物）|Check if BWA is needed (check final output)
            # 如果BAM或后续产物已存在，说明BWA已完成|If BAM or later outputs exist, BWA is done
            if (check_file_exists(str(aligned_bam), 1000000) or
                check_file_exists(str(aligned_nsorted), 1000000) or
                check_file_exists(str(aligned_fixmate), 1000000) or
                check_file_exists(str(aligned_sorted), 1000000) or
                check_file_exists(str(self.config.final_bam), 1000000)):
                self.logger.info(f"BWA比对已完成（产物存在），跳过|BWA alignment completed (outputs exist), skipping")
            else:
                bwa_cmd = [
                    self.config.bwa_bin,
                    'mem',
                    '-5SP',
                    '-t', str(self.config.threads),
                    ref_fa,
                    hic_r1,
                    hic_r2
                ]

                # 记录命令
                from .utils import build_conda_command
                wrapped_bwa_cmd = build_conda_command(self.config.bwa_bin, bwa_cmd[1:])
                self.logger.info(f"执行|Executing: BWA MEM比对|BWA MEM alignment")
                self.logger.info(f"命令|Command: {' '.join(wrapped_bwa_cmd)} > {temp_sam}")

                # 直接运行BWA并输出到SAM文件
                import subprocess
                import os
                with open(temp_sam, 'wb') as f_out:
                    result = subprocess.run(
                        wrapped_bwa_cmd,
                        stdout=f_out,
                        stderr=subprocess.PIPE,
                        env=os.environ.copy()
                    )

                if result.returncode != 0 or not check_file_exists(str(temp_sam), 1000000):
                    self.logger.error("BWA比对失败|BWA alignment failed")
                    if result.stderr:
                        self.logger.error(f"错误|Error: {result.stderr.decode()[:500]}")
                    return False

            # 2.1.1 SAM转BAM|SAM to BAM conversion
            self.logger.info("[1.5/6] SAM转BAM|Converting SAM to BAM")

            # 检查BAM文件或后续产物是否已存在|Check if BAM or later outputs already exist
            # 如果aligned.bam或任何后续产物存在，说明SAM转BAM已完成|If aligned.bam or later outputs exist, SAM to BAM is done
            if (check_file_exists(str(aligned_bam), 1000000) or
                check_file_exists(str(aligned_nsorted), 1000000) or
                check_file_exists(str(aligned_fixmate), 1000000) or
                check_file_exists(str(aligned_sorted), 1000000) or
                check_file_exists(str(self.config.final_bam), 1000000)):
                self.logger.info(f"BAM或后续产物已存在，跳过SAM转BAM|BAM or later outputs exist, skipping SAM to BAM")
                # 删除临时SAM文件（如果存在）|Delete temporary SAM file (if exists)
                if temp_sam.exists():
                    temp_sam.unlink()
            elif not temp_sam.exists():
                # SAM文件不存在且后续产物也不存在，说明BWA没运行，需要报错
                # SAM doesn't exist and later outputs don't exist either, BWA was skipped but outputs missing
                self.logger.error(f"SAM和后续产物都不存在，无法继续|SAM and later outputs both missing, cannot continue")
                return False
            else:
                cmd = [
                    self.config.samtools_bin,
                    'view',
                    '-@', str(self.config.threads),
                    '-b',
                    '-h',
                    str(temp_sam),
                    '-o', str(aligned_bam)
                ]

                success, stdout, stderr, _ = self.cmd_runner.run_command(
                    cmd,
                    description="SAM转BAM|SAM to BAM conversion"
                )

                if not success or not check_file_exists(str(aligned_bam), 1000000):
                    self.logger.error("SAM转BAM失败|SAM to BAM conversion failed")
                    if stderr:
                        self.logger.error(f"错误|Error: {stderr[:500]}")
                    return False

                # 删除临时SAM文件|Delete temporary SAM file
                temp_sam.unlink()
                self.logger.info("SAM转BAM完成|SAM to BAM conversion completed")

            # 2.2 按名称排序|Sort by read name
            self.logger.info("[2/6] 按read name排序|Sorting by read name")

            # 检查是否已完成|Check if already done
            if (check_file_exists(str(aligned_nsorted), 1000000) or
                check_file_exists(str(aligned_fixmate), 1000000) or
                check_file_exists(str(aligned_sorted), 1000000) or
                check_file_exists(str(self.config.final_bam), 1000000)):
                self.logger.info(f"按名称排序产物已存在，跳过|Name sort outputs exist, skipping")
            else:
                # 计算每个线程内存，留有余地|Calculate per-thread memory with margin
                # samtools的-m参数是每个线程的内存限制|-m is per-thread memory limit
                # 公式：sam_ram / (threads * 1.1) 留10%余量|Formula: sam_ram / (threads * 1.1) for 10% margin
                total_mem_gb = int(self.config.sam_ram.rstrip("G"))
                per_thread_gb = int(total_mem_gb / (self.config.threads * 1.1))
                per_thread_mem = f'{per_thread_gb}G'

                cmd = [
                    self.config.samtools_bin,
                    'sort',
                    '-n',
                    '-@', str(self.config.threads),
                    '-m', per_thread_mem,
                    '-T', str(tmp_nsort / 'split'),
                    '-o', str(aligned_nsorted),
                    str(aligned_bam)
                ]

                success, _, stderr, _ = self.cmd_runner.run_command(
                    cmd,
                    description="按名称排序|Sort by name"
                )

                if not success:
                    self.logger.error("按名称排序失败|Sort by name failed")
                    return False

                # 删除中间文件|Delete intermediate file
                aligned_bam.unlink()
                self.logger.info("按名称排序完成|Sort by name completed")

            # 2.3 修复mate pair信息|Fix mate pair information
            self.logger.info("[3/6] 修复mate pair信息|Fixing mate pair")

            if (check_file_exists(str(aligned_fixmate), 1000000) or
                check_file_exists(str(aligned_sorted), 1000000) or
                check_file_exists(str(self.config.final_bam), 1000000)):
                self.logger.info(f"修复mate产物已存在，跳过|Fixmate outputs exist, skipping")
            else:
                cmd = [
                    self.config.samtools_bin,
                    'fixmate',
                    '-m',
                    '-@', str(self.config.threads),
                    str(aligned_nsorted),
                    str(aligned_fixmate)
                ]

                success, _, stderr, _ = self.cmd_runner.run_command(
                    cmd,
                    description="修复mate pair|Fix mate pairs"
                )

                if not success:
                    self.logger.error("修复mate pair失败|Fix mate pairs failed")
                    return False

                aligned_nsorted.unlink()
                self.logger.info("修复mate pair完成|Fix mate pairs completed")

            # 2.4 按坐标排序|Sort by coordinate
            self.logger.info("[4/6] 按坐标排序|Sorting by coordinate")

            if (check_file_exists(str(aligned_sorted), 1000000) or
                check_file_exists(str(self.config.final_bam), 1000000)):
                self.logger.info(f"按坐标排序产物已存在，跳过|Coordinate sort outputs exist, skipping")
            else:
                # 计算每个线程内存，留有余地|Calculate per-thread memory with margin
                total_mem_gb = int(self.config.sam_ram.rstrip("G"))
                per_thread_gb = int(total_mem_gb / (self.config.threads * 1.1))
                per_thread_mem = f'{per_thread_gb}G'

                cmd = [
                    self.config.samtools_bin,
                    'sort',
                    '-@', str(self.config.threads),
                    '-m', per_thread_mem,
                    '-T', str(tmp_sort / 'split'),
                    '-o', str(aligned_sorted),
                    str(aligned_fixmate)
                ]

                success, _, stderr, _ = self.cmd_runner.run_command(
                    cmd,
                    description="按坐标排序|Sort by coordinate"
                )

                if not success:
                    self.logger.error("按坐标排序失败|Sort by coordinate failed")
                    return False

                aligned_fixmate.unlink()
                self.logger.info("按坐标排序完成|Sort by coordinate completed")

            # 2.5 标记并移除PCR重复|Mark and remove PCR duplicates
            self.logger.info("[5/6] 标记并移除PCR重复|Marking and removing duplicates")

            if check_file_exists(str(self.config.final_bam), 1000000):
                self.logger.info(f"最终BAM已存在，跳过标记重复|Final BAM exists, skipping markdup: {self.config.final_bam.name}")
            else:
                cmd = [
                    self.config.samtools_bin,
                    'markdup',
                    '-r',
                    '-@', str(self.config.threads),
                    str(aligned_sorted),
                    str(self.config.final_bam)
                ]

                success, _, stderr, _ = self.cmd_runner.run_command(
                    cmd,
                    description="标记重复|Mark duplicates"
                )

                if not success:
                    self.logger.error("标记重复失败|Mark duplicates failed")
                    return False

                aligned_sorted.unlink()
                self.logger.info("标记重复完成|Mark duplicates completed")

            # 2.6 构建BAM索引|Build BAM index
            self.logger.info("[6/6] 构建BAM索引|Building BAM index")

            if check_file_exists(str(self.config.final_bam) + '.bai'):
                self.logger.info(f"BAM索引已存在，跳过|BAM index exists, skipping")
            else:
                cmd = [
                    self.config.samtools_bin,
                    'index',
                    '-@', str(self.config.threads),
                    str(self.config.final_bam)
                ]

                success, _, stderr, _ = self.cmd_runner.run_command(
                    cmd,
                    description="BAM索引|BAM indexing"
                )

                if not success:
                    self.logger.warning("BAM索引构建失败（非致命）|BAM index failed (non-fatal)")

            # 2.7 统计比对结果|Calculate alignment statistics
            self.logger.info("比对统计信息|Alignment statistics")

            cmd = [self.config.samtools_bin, 'flagstat', '-@', str(self.config.threads), str(self.config.final_bam)]

            success, stdout, stderr, _ = self.cmd_runner.run_command(
                cmd,
                description="统计比对|Calculate statistics",
                check=False
            )

            if success and stdout:
                self.logger.info(f"比对统计|Alignment stats:\n{stdout}")

            self.logger.info("步骤 2 完成|Step 2 completed")
            return True

        except Exception as e:
            self.logger.error(f"步骤2执行异常|Step 2 execution error: {e}")
            return False

        finally:
            # 清理临时目录|Clean up temp directories
            if not self.config.keep_temp:
                import shutil
                if tmp_nsort.exists():
                    shutil.rmtree(tmp_nsort, ignore_errors=True)
                if tmp_sort.exists():
                    shutil.rmtree(tmp_sort, ignore_errors=True)

    def step_3_scaffolding(self) -> bool:
        """
        步骤3: YaHS染色体挂载|Step 3: YaHS scaffolding

        使用YaHS进行Hi-C scaffolding
        Perform Hi-C scaffolding using YaHS

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 3/6: YaHS染色体挂载|Step 3/6: YaHS scaffolding")
        self.logger.info("=" * 60)

        # 检查是否已完成
        if check_file_exists(str(self.config.final_scaffold_fa)):
            self.logger.info("Scaffold文件已存在，跳过|Scaffold file exists, skipping")
            return True

        # 构建YaHS命令|Build YaHS command
        yahs_params = self.config.get_yahs_params()

        # 使用步骤1中创建的本地参考基因组副本|Use local reference copy from step 1
        ref_fa_orig = str(Path(self.config.ref_fa).resolve())
        ref_name = Path(ref_fa_orig).name
        local_ref = str(Path(self.config.step1_dir / ref_name).resolve())  # 使用绝对路径|Use absolute path

        cmd = [
            self.config.yahs_bin,
            '-o', str(self.config.yahs_out_prefix),
            local_ref,  # 使用本地副本|Use local copy
            str(self.config.final_bam)
        ] + yahs_params

        self.logger.info(f"运行YaHS|Running YaHS")
        self.logger.debug(f"YaHS参数|YaHS params: {' '.join(yahs_params)}")

        success, stdout, stderr, _ = self.cmd_runner.run_command(
            cmd,
            description="YaHS scaffolding",
            timeout=7200  # 2小时超时
        )

        if not success:
            self.logger.error("YaHS运行失败|YaHS execution failed")
            if stderr:
                self.logger.error(f"YaHS错误|YaHS error:\n{stderr}")
            return False

        # 检查输出文件|Check output files
        if not check_file_exists(str(self.config.final_scaffold_fa)):
            self.logger.error("YaHS未生成最终scaffold文件|YaHS did not generate final scaffold")
            return False

        self.logger.info(f"Scaffold文件生成成功|Scaffold file generated: {self.config.final_scaffold_fa}")
        self.logger.info("步骤 3 完成|Step 3 completed")
        return True

    def step_4_hic_standard(self) -> bool:
        """
        步骤4: 生成标准Hi-C热图|Step 4: Generate standard Hi-C heatmap

        生成标准.hic文件用于可视化
        Generate standard .hic file for visualization

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 4/6: 生成标准Hi-C热图|Step 4/6: Generate standard Hi-C heatmap")
        self.logger.info("=" * 60)

        if not self.config.juicer_jar:
            self.logger.warning("未指定juicer_tools.jar，跳过此步骤|juicer_tools.jar not specified, skipping")
            return True

        # 检查是否已完成
        if check_file_exists(str(self.config.final_hic), 100000):
            self.logger.info("标准.hic文件已存在，跳过|Standard .hic file exists, skipping")
            return True

        # 检查bin文件|Check bin file
        if not check_file_exists(str(self.config.yahs_bin_file)):
            self.logger.error("YaHS .bin文件未找到|YaHS .bin file not found")
            return False

        try:
            # 4.1 生成并排序比对链接|Generate and sort alignments
            self.logger.info("[1/3] 生成并排序比对链接|Generating and sorting alignments")

            # 使用juicer pre生成链接
            # 使用步骤1中创建的本地参考基因组副本|Use local reference copy from step 1
            ref_fa_orig = str(Path(self.config.ref_fa).resolve())
            ref_name = Path(ref_fa_orig).name
            local_ref = str(Path(self.config.step1_dir / ref_name).resolve())  # 使用绝对路径|Use absolute path

            juicer_pre_cmd = [
                self.config.juicer_bin,
                'pre',
                str(self.config.yahs_bin_file),
                str(self.config.final_scaffold_agp),
                local_ref + '.fai'  # 使用本地副本的.fai文件|Use local copy's .fai file
            ]

            success, stdout, stderr, _ = self.cmd_runner.run_command(
                juicer_pre_cmd,
                description="juicer pre",
                check=False
            )

            if not success or not stdout:
                self.logger.error("juicer pre失败|juicer pre failed")
                return False

            # 保存到临时文件并排序|Save to temp file and sort
            import tempfile
            with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt', dir=self.config.step4_dir) as f:
                f.write(stdout)
                temp_file = f.name

            sorted_file = self.config.step4_dir / 'alignments_sorted.txt'

            sort_cmd = [
                'sort',
                '-k2,2d',
                '-k6,6d',
                '-k3,3n',
                '-k7,7n',
                '--parallel', str(self.config.threads),
                '-S', self.config.sam_ram,
                '-T', str(self.config.step4_dir / 'tmp_sort'),
                '-o', str(sorted_file),
                temp_file
            ]

            success, _, stderr, _ = self.cmd_runner.run_command(
                sort_cmd,
                description="排序比对链接|Sort alignments"
            )

            # 删除临时文件|Delete temp file
            Path(temp_file).unlink(missing_ok=True)

            if not success:
                self.logger.error("排序失败|Sort failed")
                return False

            # 4.2 准备染色体大小文件|Prepare chromosome sizes file
            self.logger.info("[2/3] 准备染色体大小文件|Preparing chromosome sizes file")

            # 为scaffold FASTA建立索引
            cmd = [self.config.samtools_bin, 'faidx', str(self.config.final_scaffold_fa)]
            self.cmd_runner.run_command(cmd, description="索引scaffold|Index scaffold", check=False)

            # 提取染色体大小
            chrom_sizes = self.config.step4_dir / 'chrom.sizes'
            with open(str(Path(self.config.final_scaffold_fa).resolve()) + '.fai', 'r') as f_in, open(chrom_sizes, 'w') as f_out:
                for line in f_in:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        f_out.write(f"{fields[0]}\t{fields[1]}\n")

            # 4.3 生成.hic文件|Generate .hic file
            self.logger.info("[3/3] 运行juicer_tools|Running juicer_tools")

            cmd = [
                self.config.java_cmd,
                f'-Xmx{self.config.java_ram}',
                f'-Xms8G',
                '-jar', self.config.juicer_jar,
                'pre',
                str(sorted_file),
                str(self.config.final_hic) + '.part',
                str(chrom_sizes)
            ]

            success, stdout, stderr, _ = self.cmd_runner.run_command(
                cmd,
                description="juicer_tools pre",
                timeout=3600
            )

            if success and Path(str(self.config.final_hic) + '.part').exists():
                Path(str(self.config.final_hic) + '.part').rename(self.config.final_hic)
                self.logger.info("标准.hic文件生成成功|Standard .hic file generated")
            else:
                self.logger.error(".hic文件生成失败|.hic file generation failed")
                return False

            # 清理大文件|Clean up large files
            if not self.config.keep_temp:
                sorted_file.unlink(missing_ok=True)
                chrom_sizes.unlink(missing_ok=True)

            self.logger.info("步骤 4 完成|Step 4 completed")
            return True

        except Exception as e:
            self.logger.error(f"步骤4执行异常|Step 4 execution error: {e}")
            return False

    def step_5_jbat(self) -> bool:
        """
        步骤5: 生成JBAT文件|Step 5: Generate JBAT files

        生成用于Juicebox JBAT手动校正的文件
        Generate files for Juicebox JBAT manual curation

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 5/6: 生成JBAT文件|Step 5/6: Generate JBAT files")
        self.logger.info("=" * 60)

        if not self.config.juicer_jar:
            self.logger.warning("未指定juicer_tools.jar，跳过此步骤|juicer_tools.jar not specified, skipping")
            return True

        # 检查是否已完成
        if (check_file_exists(str(self.config.jbat_hic), 100000) and
            check_file_exists(str(self.config.jbat_assembly))):
            self.logger.info("JBAT文件已存在，跳过|JBAT files exist, skipping")
            return True

        try:
            # 5.1 生成JBAT文本和assembly文件|Generate JBAT text and assembly files
            self.logger.info("[1/4] 生成JBAT文件|Generating JBAT files")

            jbat_prefix = self.config.step5_dir / 'out_JBAT'

            # 使用步骤1中创建的本地参考基因组副本|Use local reference copy from step 1
            ref_fa_orig = str(Path(self.config.ref_fa).resolve())
            ref_name = Path(ref_fa_orig).name
            local_ref = str(Path(self.config.step1_dir / ref_name).resolve())  # 使用绝对路径|Use absolute path

            cmd = [
                self.config.juicer_bin,
                'pre',
                '-a',
                '-o', str(jbat_prefix),
                str(self.config.yahs_bin_file),
                str(self.config.final_scaffold_agp),
                local_ref + '.fai'  # 使用本地副本的.fai文件|Use local copy's .fai file
            ]

            success, stdout, stderr, _ = self.cmd_runner.run_command(
                cmd,
                description="juicer pre -a (JBAT)",
                timeout=600
            )

            if not success:
                self.logger.error("JBAT文件生成失败|JBAT file generation failed")
                return False

            # 检查输出文件|Check output files
            if not check_file_exists(str(self.config.jbat_assembly)):
                self.logger.error("JBAT assembly文件未找到|JBAT assembly file not found")
                return False

            self.logger.info("JBAT文件生成完成|JBAT files generated")

            # 5.2 提取assembly大小|Extract assembly size
            self.logger.info("[2/4] 提取assembly大小信息|Extracting assembly size")

            jbat_log = self.config.step5_dir / 'juicer_pre.log'

            # 运行juicer pre并保存日志|Run juicer pre and save log
            with open(jbat_log, 'w') as f_log:
                result = self.cmd_runner.run_command(
                    cmd,
                    description="juicer pre with log",
                    check=False
                )

                # 这里需要重新运行来获取日志，或者从上面的stderr中提取
                # 简化处理：使用默认assembly大小
                total_bp = sum(len(seq) for seq in self._parse_fasta(str(self.config.final_scaffold_fa)))
                chrom_sizes = self.config.step5_dir / 'jbat_chrom_sizes.txt'

                with open(chrom_sizes, 'w') as f:
                    f.write(f"assembly\t{total_bp}\n")

            # 5.3 排序JBAT文件|Sort JBAT file
            self.logger.info("[3/4] 排序JBAT文件|Sorting JBAT file")

            jbat_txt = self.config.step5_dir / 'out_JBAT.txt'
            sorted_file = self.config.step5_dir / 'out_JBAT_sorted.txt'

            if jbat_txt.exists():
                sort_cmd = [
                    'sort',
                    '-k2,2d',
                    '-k6,6d',
                    '-k3,3n',
                    '-k7,7n',
                    '--parallel', str(self.config.threads),
                    '-S', self.config.sam_ram,
                    '-T', str(self.config.step5_dir / 'tmp'),
                    '-o', str(sorted_file),
                    str(jbat_txt)
                ]

                success, _, stderr, _ = self.cmd_runner.run_command(
                    sort_cmd,
                    description="排序JBAT|Sort JBAT"
                )

                if not success:
                    self.logger.warning("排序失败，使用未排序文件|Sort failed, using unsorted file")
                    sorted_file = jbat_txt
            else:
                sorted_file = None

            # 5.4 生成JBAT .hic文件|Generate JBAT .hic file
            self.logger.info("[4/4] 生成JBAT .hic文件|Generating JBAT .hic file")

            if sorted_file and sorted_file.exists():
                cmd = [
                    self.config.java_cmd,
                    f'-Xmx{self.config.java_ram}',
                    f'-Xms8G',
                    '-jar', self.config.juicer_jar,
                    'pre',
                    str(sorted_file),
                    str(self.config.jbat_hic) + '.part',
                    str(chrom_sizes)
                ]

                success, stdout, stderr, _ = self.cmd_runner.run_command(
                    cmd,
                    description="juicer_tools pre (JBAT)",
                    timeout=3600
                )

                if success and Path(str(self.config.jbat_hic) + '.part').exists():
                    Path(str(self.config.jbat_hic) + '.part').rename(self.config.jbat_hic)
                    self.logger.info("JBAT .hic文件生成成功|JBAT .hic file generated")
                else:
                    self.logger.error("JBAT .hic文件生成失败|JBAT .hic file generation failed")
                    return False

            # 清理大文件|Clean up large files
            if not self.config.keep_temp:
                if sorted_file and sorted_file != jbat_txt:
                    sorted_file.unlink(missing_ok=True)
                jbat_txt.unlink(missing_ok=True)

            self.logger.info("步骤 5 完成|Step 5 completed")
            return True

        except Exception as e:
            self.logger.error(f"步骤5执行异常|Step 5 execution error: {e}")
            return False

    def step_6_assessment(self) -> bool:
        """
        步骤6: 组装质量评估|Step 6: Assembly quality assessment

        计算scaffold统计指标（N50, N90等）
        Calculate scaffold statistics (N50, N90, etc.)

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤 6/6: 组装质量评估|Step 6/6: Assembly quality assessment")
        self.logger.info("=" * 60)

        try:
            scaffold_fa = self.config.final_scaffold_fa

            if not check_file_exists(str(scaffold_fa)):
                self.logger.error("Scaffold文件未找到|Scaffold file not found")
                return False

            # 提取序列长度|Extract sequence lengths
            self.logger.info("计算统计指标|Calculating statistics")

            lengths = []
            current_seq = []

            with open(scaffold_fa, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        if current_seq:
                            lengths.append(len(''.join(current_seq)))
                        current_seq = []
                    else:
                        current_seq.append(line.strip())

                if current_seq:
                    lengths.append(len(''.join(current_seq)))

            if not lengths:
                self.logger.error("无法解析FASTA文件|Failed to parse FASTA file")
                return False

            # 排序|Sort
            lengths.sort(reverse=True)

            # 计算统计指标|Calculate statistics
            total_length = sum(lengths)
            n_scaffolds = len(lengths)
            max_scaffold = lengths[0] if lengths else 0

            # 计算N50, N90|Calculate N50, N90
            n50, n90 = self._calculate_nx(lengths, total_length, [0.5, 0.9])

            # 计算L50, L90|Calculate L50, L90
            l50, l90 = self._calculate_lx(lengths, total_length, [0.5, 0.9])

            # 输出报告|Output report
            self.logger.info("组装质量统计报告|Assembly Quality Statistics")
            self.logger.info("=" * 60)
            self.logger.info(f"Scaffold总数|Total scaffolds: {n_scaffolds}")
            self.logger.info(f"总长度|Total length: {self._format_bp(total_length)} bp")
            self.logger.info(f"最长Scaffold|Longest scaffold: {self._format_bp(max_scaffold)} bp")
            self.logger.info(f"N50: {self._format_bp(n50)} bp")
            self.logger.info(f"N90: {self._format_bp(n90)} bp")
            self.logger.info(f"L50: {l50}")
            self.logger.info(f"L90: {l90}")
            self.logger.info("=" * 60)

            # 保存到文件|Save to file
            with open(self.config.assembly_metrics, 'w') as f:
                f.write("Assembly Quality Statistics\n")
                f.write("=" * 60 + "\n")
                f.write(f"Total scaffolds: {n_scaffolds}\n")
                f.write(f"Total length: {total_length} bp\n")
                f.write(f"Longest scaffold: {max_scaffold} bp\n")
                f.write(f"N50: {n50} bp\n")
                f.write(f"N90: {n90} bp\n")
                f.write(f"L50: {l50}\n")
                f.write(f"L90: {l90}\n")

            self.logger.info("步骤 6 完成|Step 6 completed")
            return True

        except Exception as e:
            self.logger.error(f"步骤6执行异常|Step 6 execution error: {e}")
            return False

    def _calculate_nx(self, lengths: list, total: int, thresholds: list) -> list:
        """计算Nx值|Calculate Nx values"""
        results = []
        for threshold in thresholds:
            target = total * threshold
            cum = 0
            nx = 0
            for length in lengths:
                cum += length
                if cum >= target:
                    nx = length
                    break
            results.append(nx)
        return results

    def _calculate_lx(self, lengths: list, total: int, thresholds: list) -> list:
        """计算Lx值|Calculate Lx values"""
        results = []
        for threshold in thresholds:
            target = total * threshold
            cum = 0
            lx = 0
            for i, length in enumerate(lengths):
                cum += length
                if cum >= target:
                    lx = i + 1
                    break
            results.append(lx)
        return results

    def _format_bp(self, bp: int) -> str:
        """格式化碱基对数量|Format base pairs"""
        if bp >= 1_000_000_000:
            return f"{bp / 1_000_000_000:.2f}G"
        elif bp >= 1_000_000:
            return f"{bp / 1_000_000:.2f}M"
        elif bp >= 1_000:
            return f"{bp / 1_000:.2f}K"
        return str(bp)

    def _parse_fasta(self, fasta_file: str):
        """解析FASTA文件生成器|Parse FASTA file generator"""
        with open(fasta_file, 'r') as f:
            seq = []
            for line in f:
                if line.startswith('>'):
                    if seq:
                        yield ''.join(seq)
                    seq = []
                else:
                    seq.append(line.strip())
            if seq:
                yield ''.join(seq)
