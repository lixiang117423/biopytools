"""
三代转录组比对核心模块|Long RNA-seq Alignment Core Module
"""

import re
import subprocess
from pathlib import Path
from typing import Dict, Tuple


class LongRNASeqAligner:
    """三代转录组比对器|Long RNA-seq Aligner"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def align(self) -> bool:
        """
        执行三代转录组比对流程|Execute long RNA-seq alignment pipeline

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始三代转录组比对分析|Starting long RNA-seq alignment")
        self.logger.info("=" * 60)
        self.logger.info(f"样本名称|Sample name: {self.config.sample_name}")
        self.logger.info(f"参考基因组|Reference genome: {self.config.ref_genome}")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"输入类型|Input type: {self.config.input_type}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"最大intron长度|Max intron length: {self.config.max_intron}")
        self.logger.info(f"最小MAPQ|Min MAPQ: {self.config.min_mapq}")

        try:
            # 步骤1: 准备FASTQ文件（如果是BAM则需要转换，如果是FASTQ则直接使用）|Step 1: Prepare FASTQ file
            self.logger.info("步骤1: 准备FASTQ文件|Step 1: Preparing FASTQ file")
            fastq_files = self._prepare_fastq()
            if not fastq_files:
                return False

            # 步骤2: minimap2比对|Step 2: Align with minimap2
            self.logger.info("步骤2: minimap2比对|Step 2: Aligning with minimap2")
            sorted_bam = self._align_with_minimap2(fastq_files)
            if not sorted_bam:
                return False

            # 步骤3: 建立索引|Step 3: Build index
            self.logger.info("步骤3: 建立索引|Step 3: Building index")
            if not self._build_index(sorted_bam):
                return False

            # 步骤4: 生成统计信息|Step 4: Generate statistics
            self.logger.info("步骤4: 生成统计信息|Step 4: Generating statistics")
            if not self._generate_statistics(sorted_bam):
                return False

            # 步骤5: 质量控制|Step 5: Quality control
            self.logger.info("步骤5: 质量控制|Step 5: Quality control")
            self._quality_check(sorted_bam)

            self.logger.info("=" * 60)
            self.logger.info("三代转录组比对完成|Long RNA-seq alignment completed")
            self.logger.info("=" * 60)
            self.logger.info(f"输出BAM|Output BAM: {sorted_bam}")
            self.logger.info(f"统计文件|Statistics: {self.config.stats_dir}")

            return True

        except Exception as e:
            self.logger.error(f"比对失败|Alignment failed: {e}")
            import traceback
            traceback.print_exc()
            return False

    def _prepare_fastq(self):
        """
        准备FASTQ文件|Prepare FASTQ file

        Returns:
            单个FASTQ路径或元组(r1, r2)|Single FASTQ path or tuple (r1, r2)
        """
        from .utils import run_command

        # 如果是FASTQ文件，直接使用|If FASTQ file, use directly
        if self.config.input_type in ['fastq_single', 'fastq_paired']:
            self.logger.info(f"检测到FASTQ输入文件|Detected FASTQ input file: {self.config.input_file}")

            if self.config.input_type == 'fastq_paired':
                # 查找R2文件|Find R2 file
                r2_file = self.config._find_r2_file()
                if r2_file:
                    self.logger.info(f"双端测序|Paired-end sequencing:")
                    self.logger.info(f"  R1|Read 1: {self.config.input_file}")
                    self.logger.info(f"  R2|Read 2: {r2_file}")
                    return (self.config.input_file, r2_file)
                else:
                    self.logger.error("未找到配对的R2文件|R2 file not found")
                    return None
            else:
                self.logger.info(f"单端测序|Single-end sequencing: {self.config.input_file}")
                return self.config.input_file

        # 如果是BAM文件，转换为FASTQ|If BAM file, convert to FASTQ
        elif self.config.input_type == 'bam':
            self.logger.info("检测到BAM输入文件，正在转换为FASTQ|Detected BAM input file, converting to FASTQ")
            return self._bam_to_fastq()

        return None

    def _bam_to_fastq(self):
        """
        将BAM转换为FASTQ|Convert BAM to FASTQ

        Returns:
            Path: FASTQ文件路径|FASTQ file path
        """
        from .utils import run_command

        fastq_file = self.config.tmp_dir / f"{self.config.sample_name}.fastq"

        self.logger.info(f"转换BAM到FASTQ|Converting BAM to FASTQ: {fastq_file}")

        cmd = [
            self.config.samtools_path,
            "fastq",
            "-@", str(self.config.threads),
            self.config.input_file
        ]

        try:
            result = run_command(cmd, self.logger, check=True, capture_output=True)

            with open(fastq_file, 'w') as f:
                f.write(result.stdout)

            # 统计reads数量|Count reads
            read_count = 0
            with open(fastq_file, 'r') as f:
                for i, line in enumerate(f):
                    if i % 4 == 0:  # FASTQ每4行一条reads
                        read_count += 1

            self.logger.info(f"总reads数|Total reads: {read_count}")

            return fastq_file

        except Exception as e:
            self.logger.error(f"BAM转FASTQ失败|Failed to convert BAM to FASTQ: {e}")
            return None

    def _align_with_minimap2(self, fastq_files) -> Path:
        """
        使用minimap2比对|Align with minimap2

        Args:
            fastq_files: FASTQ文件（路径或元组(r1, r2)）|FASTQ file (path or tuple (r1, r2))

        Returns:
            Path: 排序后的BAM文件路径|Sorted BAM file path
        """
        from .utils import run_command

        # 分步骤文件
        sam_file = self.config.tmp_dir / f"{self.config.sample_name}.sam"
        bam_file = self.config.tmp_dir / f"{self.config.sample_name}.bam"
        sorted_bam = self.config.result_dir / f"{self.config.sample_name}.sorted.bam"
        tmp_prefix = self.config.tmp_dir / self.config.sample_name

        self.logger.info(f"比对参数|Alignment parameters: {' '.join(self.config.get_minimap2_params())}")

        try:
            # 步骤1: 运行minimap2，输出SAM文件|Step 1: Run minimap2, output SAM file
            self.logger.info("步骤1/3: 运行minimap2生成SAM文件|Step 1/3: Running minimap2 to generate SAM file")

            minimap2_cmd = [
                self.config.minimap2_path
            ] + self.config.get_minimap2_params() + [
                self.config.ref_genome
            ]

            # 根据输入类型添加FASTQ文件|Add FASTQ files based on input type
            if isinstance(fastq_files, tuple):
                # 双端FASTQ|Paired-end FASTQ
                r1, r2 = fastq_files
                minimap2_cmd.extend([str(r1), str(r2)])
                self.logger.info(f"双端比对模式|Paired-end alignment mode")
            else:
                # 单端FASTQ|Single-end FASTQ
                minimap2_cmd.append(str(fastq_files))
                self.logger.info(f"单端比对模式|Single-end alignment mode")

            self.logger.info(f"运行命令|Running command: {' '.join(minimap2_cmd)} > {sam_file}")

            # 使用shell重定向输出到SAM文件|Use shell redirection to output to SAM file
            with open(sam_file, 'w') as f:
                result = subprocess.run(
                    minimap2_cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True
                )

            if not sam_file.exists():
                self.logger.error(f"SAM文件未生成|SAM file not generated: {sam_file}")
                return None

            self.logger.info(f"SAM文件生成完成|SAM file generated: {sam_file}")

            # 步骤2: 转换SAM为BAM|Step 2: Convert SAM to BAM
            self.logger.info("步骤2/3: 转换SAM为BAM|Step 2/3: Converting SAM to BAM")

            samtools_view_cmd = [
                self.config.samtools_path,
                "view",
                "-b",
                str(sam_file)
            ]

            self.logger.info(f"运行命令|Running command: {' '.join(samtools_view_cmd)} > {bam_file}")

            # 使用shell重定向输出到BAM文件|Use shell redirection to output to BAM file
            with open(bam_file, 'wb') as f:
                result = subprocess.run(
                    samtools_view_cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    check=True
                )

            if not bam_file.exists():
                self.logger.error(f"BAM文件未生成|BAM file not generated: {bam_file}")
                return None

            self.logger.info(f"BAM文件生成完成|BAM file generated: {bam_file}")

            # 删除SAM文件以节省空间|Delete SAM file to save space
            sam_file.unlink()
            self.logger.info(f"已删除临时SAM文件|Temporary SAM file deleted: {sam_file}")

            # 步骤3: 排序BAM文件|Step 3: Sort BAM file
            self.logger.info("步骤3/3: 排序BAM文件|Step 3/3: Sorting BAM file")

            samtools_sort_cmd = [
                self.config.samtools_path,
                "sort",
                "-@", str(self.config.threads),
                "-o", str(sorted_bam),
                str(bam_file)
            ]

            self.logger.info(f"运行命令|Running command: {' '.join(samtools_sort_cmd)}")
            result = subprocess.run(
                samtools_sort_cmd,
                stderr=subprocess.PIPE,
                check=True
            )

            if not sorted_bam.exists():
                self.logger.error(f"排序BAM文件未生成|Sorted BAM file not generated: {sorted_bam}")
                return None

            self.logger.info(f"BAM文件排序完成|BAM file sorted: {sorted_bam}")

            # 删除未排序的BAM文件|Delete unsorted BAM file
            bam_file.unlink()
            self.logger.info(f"已删除未排序BAM文件|Unsorted BAM file deleted: {bam_file}")

            self.logger.info(f"比对和排序完成|Alignment and sorting completed: {sorted_bam}")

            return sorted_bam

        except Exception as e:
            self.logger.error(f"比对失败|Alignment failed: {e}")
            import traceback
            traceback.print_exc()
            return None

    def _build_index(self, bam_file: Path) -> bool:
        """
        建立BAM索引|Build BAM index

        Args:
            bam_file: BAM文件路径|BAM file path

        Returns:
            bool: 是否成功|Whether successful
        """
        from .utils import run_command

        self.logger.info(f"建立索引|Building index: {bam_file}")

        cmd = [
            self.config.samtools_path,
            "index",
            "-@", str(self.config.threads),
            str(bam_file)
        ]

        try:
            run_command(cmd, self.logger, check=True)
            self.logger.info("索引创建完成|Index created")
            return True

        except Exception as e:
            self.logger.error(f"建立索引失败|Failed to build index: {e}")
            return False

    def _generate_statistics(self, bam_file: Path) -> bool:
        """
        生成比对统计信息|Generate alignment statistics

        Args:
            bam_file: BAM文件路径|BAM file path

        Returns:
            bool: 是否成功|Whether successful
        """
        from .utils import run_command

        stats_file = self.config.stats_dir / f"{self.config.sample_name}_stats.txt"
        detail_stats = self.config.stats_dir / f"{self.config.sample_name}_detail_stats.txt"

        try:
            # 基础统计|Basic statistics
            self.logger.info("生成基础统计|Generating basic statistics")
            cmd = [
                self.config.samtools_path,
                "flagstat",
                "-@", str(self.config.threads),
                str(bam_file)
            ]
            result = run_command(cmd, self.logger, check=True, capture_output=True)

            with open(stats_file, 'w') as f:
                f.write(result.stdout)

            # 详细统计|Detailed statistics
            self.logger.info("生成详细统计|Generating detailed statistics")
            cmd = [
                self.config.samtools_path,
                "stats",
                "-@", str(self.config.threads),
                str(bam_file)
            ]
            result = run_command(cmd, self.logger, check=True, capture_output=True)

            with open(detail_stats, 'w') as f:
                f.write(result.stdout)

            # 覆盖度统计|Coverage statistics
            self.logger.info("计算覆盖度统计|Calculating coverage statistics")
            cmd = [
                self.config.samtools_path,
                "depth",
                str(bam_file)
            ]
            result = run_command(cmd, self.logger, check=False, capture_output=True)

            # 计算平均深度|Calculate average depth
            lines = result.stdout.strip().split('\n')
            total_depth = 0
            count = 0
            for line in lines:
                parts = line.split('\t')
                if len(parts) >= 3:
                    total_depth += int(parts[2])
                    count += 1

            avg_depth = total_depth / count if count > 0 else 0

            with open(stats_file, 'a') as f:
                f.write(f"\n平均深度|Average depth: {avg_depth:.2f}\n")

            # 解析并显示统计结果|Parse and display statistics
            self._display_statistics(stats_file)

            return True

        except Exception as e:
            self.logger.error(f"生成统计信息失败|Failed to generate statistics: {e}")
            return False

    def _display_statistics(self, stats_file: Path):
        """
        显示统计信息|Display statistics

        Args:
            stats_file: 统计文件路径|Statistics file path
        """
        self.logger.info(f"比对统计结果|Alignment statistics results:")
        self.logger.info("-" * 60)

        with open(stats_file, 'r') as f:
            content = f.read()
            self.logger.info(content)

        # 提取关键指标|Extract key metrics
        with open(stats_file, 'r') as f:
            for line in f:
                if "mapped (" in line:
                    parts = line.split()
                    if len(parts) >= 5:
                        total_reads = parts[0]
                        mapped_reads = parts[0]
                        mapping_rate = parts[4]

                        self.logger.info(f"总reads数|Total reads: {total_reads}")
                        self.logger.info(f"比对上的reads|Mapped reads: {mapped_reads}")
                        self.logger.info(f"比对率|Mapping rate: {mapping_rate}")

        # 生成样本汇总|Generate sample summary
        summary_file = self.config.stats_dir / f"{self.config.sample_name}_summary.txt"
        with open(summary_file, 'w') as f:
            f.write(f"样本|Sample: {self.config.sample_name}\n")
            f.write(f"输入文件|Input file: {self.config.input_file}\n")
            f.write(f"输出文件|Output file: {self.config.result_dir / f'{self.config.sample_name}.sorted.bam'}\n")

            with open(stats_file, 'r') as stats:
                for line in stats:
                    if "mapped (" in line:
                        parts = line.split()
                        if len(parts) >= 5:
                            f.write(f"总reads数|Total reads: {parts[0]}\n")
                            f.write(f"比对reads数|Mapped reads: {parts[0]}\n")
                            f.write(f"比对率|Mapping rate: {parts[4]}\n")
                            break

    def _quality_check(self, bam_file: Path):
        """
        质量控制检查|Quality control check

        Args:
            bam_file: BAM文件路径|BAM file path
        """
        from .utils import run_command

        self.logger.info("执行质量检查|Performing quality check")

        # 检查文件完整性|Check file integrity
        cmd = [
            self.config.samtools_path,
            "quickcheck",
            str(bam_file)
        ]

        try:
            result = run_command(cmd, self.logger, check=False, capture_output=True)

            if result.returncode != 0:
                self.logger.error(f"BAM文件可能损坏|BAM file may be corrupted: {bam_file}")
            else:
                self.logger.info("BAM文件完整性检查通过|BAM file integrity check passed")

        except Exception as e:
            self.logger.error(f"质量检查失败|Quality check failed: {e}")

        # 检查比对率|Check mapping rate
        stats_file = self.config.stats_dir / f"{self.config.sample_name}_stats.txt"
        if stats_file.exists():
            with open(stats_file, 'r') as f:
                for line in f:
                    if "mapped (" in line:
                        parts = line.split()
                        if len(parts) >= 5:
                            mapping_rate_str = parts[4].strip('()%')
                            try:
                                mapping_rate = float(mapping_rate_str)
                                if mapping_rate < 70:
                                    self.logger.warning(
                                        f"比对率较低|Mapping rate is low ({mapping_rate}%)，请检查数据质量和参考基因组"
                                    )
                                else:
                                    self.logger.info(f"比对率正常|Mapping rate is normal: {mapping_rate}%")
                            except ValueError:
                                pass
                            break
