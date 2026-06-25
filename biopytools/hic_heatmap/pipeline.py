"""Hi-C热图分析流程|Hi-C heatmap analysis pipeline"""

import os
import subprocess
from pathlib import Path
from typing import List, Optional
from ..common.paths import get_tool_path


class HiCPipeline:
    """Hi-C热图分析流程类|Hi-C heatmap analysis pipeline class"""

    def __init__(self, config, logger):
        """初始化流程|Initialize pipeline

        Args:
            config: HiCConfig配置对象|HiCConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

        # 输出文件路径|Output file paths
        self.base_name = Path(self.config.genome).stem
        self.bam_file = self.config.bam_dir / f"{self.base_name}.sorted.bam"
        self.restriction_bed = self.config.matrix_dir / f"{self.base_name}_restriction_sites.bed"
        self.raw_matrix = self.config.matrix_dir / f"{self.base_name}_raw.h5"
        self.corrected_matrix = self.config.matrix_dir / f"{self.base_name}_corrected.h5"
        self.heatmap_file = self.config.plot_dir / f"{self.base_name}_heatmap.png"

    def run(self) -> bool:
        """运行完整流程|Run complete pipeline

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始Hi-C热图分析流程|Starting Hi-C heatmap analysis pipeline")
        self.logger.info("=" * 60)

        # 记录软件版本信息|Log software version information
        self._log_versions()

        # 步骤1：建立基因组索引|Step 1: Build genome index
        if not self._build_genome_index():
            return False

        # 步骤2：BWA比对|Step 2: BWA alignment
        if not self._run_bwa_alignment():
            return False

        # 步骤3：samtools排序|Step 3: samtools sort
        if not self._sort_bam():
            return False

        # 步骤4：生成限制性内切酶位点文件|Step 4: Generate restriction sites file
        if not self._find_restriction_sites():
            return False

        # 步骤5：生成Hi-C矩阵|Step 5: Build Hi-C matrix
        if not self._build_matrix():
            return False

        # 步骤6：矫正矩阵|Step 6: Correct matrix
        if not self._correct_matrix():
            return False

        # 步骤7：可视化|Step 7: Visualization
        if not self._plot_heatmap():
            return False

        self.logger.info("=" * 60)
        self.logger.info("Hi-C热图分析完成|Hi-C heatmap analysis completed")
        self.logger.info(f"结果保存在|Results saved to: {self.config.output_dir}")
        self.logger.info("=" * 60)

        return True

    def _log_versions(self):
        """记录软件版本信息|Log software version information"""
        self.logger.info("=" * 60)
        self.logger.info("软件版本信息|Software versions")
        self.logger.info("=" * 60)

        import subprocess

        # BWA版本|BWA version
        try:
            result = subprocess.run(
                [self.config.bwa_path],
                capture_output=True,
                text=True,
                timeout=10
            )
            version_line = [line for line in result.stderr.split('\n') if 'Version:' in line]
            if version_line:
                self.logger.info(f"BWA: {version_line[0].strip()}")
        except:
            self.logger.info("BWA: 版本检测失败|Version detection failed")

        # samtools版本|samtools version
        try:
            result = subprocess.run(
                [self.config.samtools_path, '--version'],
                capture_output=True,
                text=True,
                timeout=10
            )
            version_line = result.stdout.split('\n')[0]
            self.logger.info(f"samtools: {version_line}")
        except:
            self.logger.info("samtools: 版本检测失败|Version detection failed")

        # HiCExplorer版本|HiCExplorer version
        try:
            result = subprocess.run(
                [self.config.hicexplorer_path],
                capture_output=True,
                text=True,
                timeout=10
            )
            # HiCExplorer没有明确的版本信息，记录路径|HiCExplorer has no clear version, log path
            self.logger.info(f"HiCExplorer: {self.config.hicexplorer_path}")
        except:
            self.logger.info("HiCExplorer: 版本检测失败|Version detection failed")

        # 记录关键参数|Log key parameters
        self.logger.info("=" * 60)
        self.logger.info("分析参数|Analysis parameters")
        self.logger.info("=" * 60)
        self.logger.info(f"限制性内切酶|Restriction enzyme: {self.config.restriction_sequence}")
        self.logger.info(f"矩阵分辨率|Bin size: {self.config.bin_size}bp ({self.config.bin_size/1000}kb)")
        self.logger.info(f"矫正方法|Correction method: {self.config.correction_method}")
        self.logger.info(f"颜色方案|Color map: {self.config.color_map}")
        self.logger.info(f"DPI: {self.config.dpi}")
        self.logger.info(f"对数变换|Log1p: {self.config.log1p}")
        self.logger.info("=" * 60)

    def _log_command(self, cmd_list, description=""):
        """记录执行的命令|Log executed command

        Args:
            cmd_list: 命令列表|Command list
            description: 命令描述|Command description
        """
        # 将命令列表转换为可复制的shell命令|Convert command list to copyable shell command
        if isinstance(cmd_list, list):
            # 对于带空格的参数用引号包围|Quote arguments with spaces
            cmd_str = ' '.join(
                f'"{arg}"' if ' ' in arg else arg
                for arg in cmd_list
            )
        else:
            cmd_str = str(cmd_list)

        if description:
            self.logger.info(f"{description}命令|{description} command:")
        self.logger.info(f"  {cmd_str}")

    def _get_hic_tool_path(self, tool_name: str) -> str:
        """获取HiCExplorer工具的完整路径|Get full path to HiCExplorer tool

        Args:
            tool_name: 工具名称，如 'hicFindRestSite', 'hicBuildMatrix' 等|Tool name

        Returns:
            str: 工具的完整路径|Full path to the tool
        """
        return str(Path(self.config.hicexplorer_path).parent / tool_name)

    def _build_genome_index(self) -> bool:
        """建立基因组索引|Build genome index

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤1: 建立基因组索引|Step 1: Building genome index")

        # 检查索引是否已存在|Check if index exists
        # BWA索引文件包括：.amb, .ann, .bwt, .pac, .sa
        # 注意：索引文件名是 genome.fa.ext 而不是 genome.ext
        index_files = [
            Path(str(self.config.genome) + ext)
            for ext in ['.amb', '.ann', '.bwt', '.pac', '.sa']
        ]

        if self.config.skip_existing and all(f.exists() for f in index_files):
            self.logger.info("跳过|Skipping: 基因组索引已存在|Genome index exists")
            return True

        try:
            cmd = [self.config.bwa_path, 'index', self.config.genome]
            self._log_command(cmd, "BWA索引|BWA index")
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            self.logger.info("基因组索引已建立|Genome index built")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"建立基因组索引失败|Failed to build genome index: {e}")
            return False

    def _run_bwa_alignment(self) -> bool:
        """运行BWA比对|Run BWA alignment

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤2: BWA比对|Step 2: Running BWA alignment")

        # 输出SAM文件|Output SAM file
        sam_file = self.config.bam_dir / f"{self.base_name}.sam"

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.bam_file.exists() and self.bam_file.with_suffix('.bam.bai').exists():
            self.logger.info("跳过|Skipping: BAM文件已存在|BAM file exists")
            return True

        try:
            # BWA mem比对|BWA mem alignment
            self.logger.info(f"比对中|Aligning with BWA mem...")
            cmd_bwa = [
                self.config.bwa_path, 'mem',
                '-A', '1',
                '-B', '4',
                '-E', '50',
                '-L', '0',
                '-t', str(self.config.threads),
                self.config.genome,
                self.config.fastq_r1,
                self.config.fastq_r2
            ]

            self._log_command(cmd_bwa, "BWA mem比对|BWA mem alignment")

            with open(sam_file, 'w') as f_out:
                result = subprocess.run(
                    cmd_bwa,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=True
                )

            self.logger.info(f"BWA比对完成|BWA alignment completed: {sam_file}")

            # 转换为BAM|Convert to BAM
            self.logger.info("转换为BAM格式|Converting to BAM format...")
            cmd_view = [
                self.config.samtools_path, 'view',
                '-Sb',
                str(sam_file),
                '-@', str(self.config.threads)
            ]

            self._log_command(cmd_view, "samtools view转换|samtools view convert")

            bam_file = self.config.bam_dir / f"{self.base_name}.bam"

            with open(bam_file, 'wb') as f_out:
                subprocess.run(cmd_view, stdout=f_out, check=True)

            # 删除SAM文件|Delete SAM file
            sam_file.unlink()

            self.logger.info(f"BAM文件已生成|BAM file created: {bam_file}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"BWA比对失败|BWA alignment failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False

    def _sort_bam(self) -> bool:
        """排序BAM文件|Sort BAM file

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤3: 排序BAM文件|Step 3: Sorting BAM file")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.bam_file.exists() and self.bam_file.with_suffix('.bam.bai').exists():
            # 检查是否已过滤（通过标记文件）
            filtered_marker = self.config.bam_dir / f"{self.base_name}.filtered.done"
            if not filtered_marker.exists():
                # 需要重新过滤
                self.logger.info("检测到需要过滤BAM文件|Detected need to filter BAM file")
            else:
                self.logger.info("跳过|Skipping: 已排序的BAM文件存在|Sorted BAM exists")
                return True

        try:
            # 定义备份BAM文件路径|Define backup BAM file path
            backup_bam = self.config.bam_dir / f"{self.base_name}.before_filter.bam"

            # 如果BAM已存在但需要过滤，先备份
            if self.bam_file.exists():
                if not backup_bam.exists():
                    import shutil
                    shutil.copy2(self.bam_file, backup_bam)
                    if self.bam_file.with_suffix('.bam.bai').exists():
                        shutil.copy2(self.bam_file.with_suffix('.bam.bai'), backup_bam.with_suffix('.bam.bai'))

            unsorted_bam = self.config.bam_dir / f"{self.base_name}.bam"

            # 如果备份存在，直接从备份开始过滤
            if backup_bam.exists():
                self.logger.info("从备份文件过滤|Filtering from backup file")
                bam_to_filter = backup_bam
            else:
                self.logger.info("排序中|Sorting...")
                cmd_sort = [
                    self.config.samtools_path, 'sort',
                    '-@', str(self.config.threads),
                    '-o', str(self.bam_file),
                    str(unsorted_bam)
                ]

                self._log_command(cmd_sort, "samtools sort排序|samtools sort")

                subprocess.run(cmd_sort, check=True)

                # 删除未排序的BAM|Delete unsorted BAM
                if unsorted_bam != self.bam_file:
                    unsorted_bam.unlink()

                bam_to_filter = self.bam_file

            # 过滤掉有问题的reads（没有CIGAR信息的reads）
            # HiCExplorer无法处理cigartuples为None的reads
            self.logger.info("过滤有问题的reads|Filtering problematic reads...")
            filtered_bam = self.config.bam_dir / f"{self.base_name}.filtered.bam"
            cmd_filter = [
                self.config.samtools_path, 'view',
                '-h',  # 包含header
                '-F', '256',  # 过滤secondary alignments (可选)
                '-q', str(self.config.min_mapping_quality),  # 最低比对质量
                '-@', str(self.config.threads),
                str(bam_to_filter),
                '-o', str(filtered_bam)
            ]

            self._log_command(cmd_filter, "samtools view过滤|samtools view filter")

            subprocess.run(cmd_filter, check=True)
            # 替换原BAM文件
            import shutil
            shutil.move(str(filtered_bam), str(self.bam_file))
            self.logger.info(f"BAM文件已过滤|BAM file filtered: {self.bam_file}")

            # 删除备份文件
            if backup_bam.exists():
                backup_bam.unlink()
                backup_bai = backup_bam.with_suffix('.bam.bai')
                if backup_bai.exists():
                    backup_bai.unlink()

            # 创建过滤完成标记
            filtered_marker = self.config.bam_dir / f"{self.base_name}.filtered.done"
            filtered_marker.touch()

            # 建立索引|Build index
            self.logger.info("建立BAM索引|Building BAM index...")
            cmd_index = [self.config.samtools_path, 'index', str(self.bam_file)]
            self._log_command(cmd_index, "samtools index索引|samtools index")
            subprocess.run(cmd_index, check=True)

            self.logger.info(f"BAM文件已排序|BAM file sorted: {self.bam_file}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"BAM排序失败|BAM sorting failed: {e}")
            return False

    def _find_restriction_sites(self) -> bool:
        """生成限制性内切酶位点文件|Generate restriction sites file

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤4: 生成限制性内切酶位点文件|Step 4: Finding restriction sites")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.restriction_bed.exists() and self.restriction_bed.stat().st_size > 0:
            self.logger.info("跳过|Skipping: 限制性位点文件已存在|Restriction sites file exists")
            return True

        try:
            # 使用独立的hicFindRestSite工具|Use standalone hicFindRestSite tool
            cmd = [
                self._get_hic_tool_path('hicFindRestSite'),
                '--fasta', self.config.genome,
                '--searchPattern', self.config.restriction_sequence,
                '--outFile', str(self.restriction_bed)
            ]

            self._log_command(cmd, "hicFindRestSite查找限制性位点|hicFindRestSite")

            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            self.logger.info(f"限制性位点文件已生成|Restriction sites file created: {self.restriction_bed}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"生成限制性位点文件失败|Failed to find restriction sites: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False

    def _build_matrix(self) -> bool:
        """生成Hi-C矩阵|Build Hi-C matrix

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤5: 生成Hi-C矩阵|Step 5: Building Hi-C matrix")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.raw_matrix.exists() and self.raw_matrix.stat().st_size > 0:
            self.logger.info("跳过|Skipping: Hi-C矩阵已存在|Hi-C matrix exists")
            return True

        # 选择矩阵构建方法|Choose matrix building method
        if self.config.use_cooler:
            return self._build_matrix_with_cooler()
        else:
            return self._build_matrix_with_hicexplorer()

    def _build_matrix_with_hicexplorer(self) -> bool:
        """使用HiCExplorer生成Hi-C矩阵|Build Hi-C matrix with HiCExplorer

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("使用HiCExplorer hicBuildMatrix生成矩阵|Using HiCExplorer hicBuildMatrix")

        try:
            cmd = [
                self._get_hic_tool_path('hicBuildMatrix'),
                '--samFiles', str(self.bam_file), str(self.bam_file),
                '--outFileName', str(self.raw_matrix),
                '--QCfolder', str(self.config.qc_dir),
                '--restrictionCutFile', str(self.restriction_bed),
                '--restrictionSequence', self.config.restriction_sequence,
                '--danglingSequence', self.config.restriction_sequence,
                '--binSize', str(self.config.bin_size),
                '--threads', str(self.config.threads),
                '--minMappingQuality', str(self.config.min_mapping_quality),
                '--maxLibraryInsertSize', str(self.config.max_library_insert_size),
                '--inputBufferSize', str(self.config.input_buffer_size),
                '--skipDuplicationCheck',  # 跳过重复检查|Skip duplication check
            ]

            # 添加区域限制参数|Add region limit
            if self.config.region:
                cmd.extend(['--region', self.config.region])

            # 添加距离限制参数|Add distance limits
            if self.config.min_distance > 0:
                cmd.extend(['--minDistance', str(self.config.min_distance)])

            self._log_command(cmd, "hicBuildMatrix生成矩阵|hicBuildMatrix")

            self.logger.info(f"生成矩阵中|Building matrix (bin size: {self.config.bin_size}bp)...")
            result = subprocess.run(cmd, check=True, capture_output=True, text=True)

            self.logger.info(f"Hi-C矩阵已生成|Hi-C matrix created: {self.raw_matrix}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"生成Hi-C矩阵失败|Failed to build Hi-C matrix: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False

    def _build_matrix_with_cooler(self) -> bool:
        """使用cooler生成Hi-C矩阵（内存友好）|Build Hi-C matrix with cooler (memory efficient)

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("使用cooler生成矩阵（内存友好）|Using cooler to build matrix (memory efficient)")

        try:
            # 步骤1: 生成chrom.sizes文件|Step 1: Generate chrom.sizes file
            chrom_sizes_file = self.config.matrix_dir / f"{self.base_name}.chrom.sizes"
            if not chrom_sizes_file.exists():
                self.logger.info("生成chrom.sizes文件|Generating chrom.sizes file...")
                cmd_chrom = [
                    self.config.samtools_path, 'view', '-H', str(self.bam_file)
                ]
                result = subprocess.run(cmd_chrom, capture_output=True, text=True, check=True)

                # 解析SAM header提取染色体长度|Parse SAM header for chromosome lengths
                chrom_sizes = []
                for line in result.stdout.split('\n'):
                    if line.startswith('@SQ'):
                        parts = line.split('\t')
                        sn = None
                        ln = None
                        for part in parts:
                            if part.startswith('SN:'):
                                sn = part[3:]
                            elif part.startswith('LN:'):
                                ln = part[3:]
                        if sn and ln:
                            chrom_sizes.append(f"{sn}\t{ln}")

                with open(chrom_sizes_file, 'w') as f:
                    f.write('\n'.join(chrom_sizes))

                self.logger.info(f"chrom.sizes文件已生成|chrom.sizes file created: {chrom_sizes_file}")

            # 步骤2: 使用pairtools将BAM转换为pairs|Step 2: Convert BAM to pairs using pairtools
            pairs_file = self.config.matrix_dir / f"{self.base_name}.pairs.gz"
            if not pairs_file.exists():
                self.logger.info("将BAM转换为pairs格式|Converting BAM to pairs format...")

                # 使用pairtools parse|Use pairtools parse
                # 使用--walks-policy mask减少内存使用（mask复杂walks而不是报错）
                pairtools_path = get_tool_path('pairtools', '~/miniforge3/envs/pairtools_v.1.1.3/bin/pairtools', 'PAIRTOOLS_PATH')
                cmd_parse = [
                    pairtools_path, 'parse',
                    '-c', str(chrom_sizes_file),  # chrom.sizes文件路径
                    '-o', str(pairs_file),
                    '--min-mapq', str(self.config.min_mapping_quality),
                    '--walks-policy', 'mask',  # mask复杂walks减少内存
                    '--max-inter-align-gap', '0',
                    '--drop-seq',  # 不保存序列，减少内存
                    '--drop-sam',  # 不保存SAM tags，减少内存
                    str(self.bam_file)
                ]

                self._log_command(cmd_parse, "pairtools parse")
                subprocess.run(cmd_parse, check=True)
                self.logger.info(f"Pairs文件已生成|Pairs file created: {pairs_file}")

            # 步骤3: 使用cooler cload生成.cool文件|Step 3: Generate .cool file with cooler cload
            cool_file = self.config.matrix_dir / f"{self.base_name}.cool"
            if not cool_file.exists():
                self.logger.info("使用cooler生成.cool文件|Generating .cool file with cooler...")

                cmd_cload = [
                    self.config.cooler_path, 'cload', 'pairs',
                    '-c1', '2',  # chrom1在第2列
                    '-p1', '3',  # pos1在第3列
                    '-c2', '4',  # chrom2在第4列
                    '-p2', '5',  # pos2在第5列
                    '--chunksize', '1000000',  # 100万行一个chunk
                    f"{chrom_sizes_file}:{self.config.bin_size}",
                    str(pairs_file),
                    str(cool_file)
                ]

                self._log_command(cmd_cload, "cooler cload pairs")
                subprocess.run(cmd_cload, check=True)
                self.logger.info(f"Cool文件已生成|Cool file created: {cool_file}")

            # 步骤4: 转换.cool为.h5格式|Step 4: Convert .cool to .h5 format
            self.logger.info("转换.cool为.h5格式|Converting .cool to .h5 format...")

            cmd_convert = [
                self._get_hic_tool_path('hicConvertFormat'),
                '--matrices', str(cool_file),
                '--inputFormat', 'cool',
                '--outputFormat', 'h5',
                '--outFileName', str(self.raw_matrix)
            ]

            self._log_command(cmd_convert, "hicConvertFormat")
            subprocess.run(cmd_convert, check=True)

            self.logger.info(f"Hi-C矩阵已生成|Hi-C matrix created: {self.raw_matrix}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"使用cooler生成矩阵失败|Failed to build matrix with cooler: {e}")
            return False
        except Exception as e:
            self.logger.error(f"错误|Error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _correct_matrix(self) -> bool:
        """矫正Hi-C矩阵|Correct Hi-C matrix

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤6: 矫正Hi-C矩阵|Step 6: Correcting Hi-C matrix")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.corrected_matrix.exists() and self.corrected_matrix.stat().st_size > 0:
            self.logger.info("跳过|Skipping: 矫正矩阵已存在|Corrected matrix exists")
            return True

        try:
            # 先生成诊断图|Generate diagnostic plot first
            diagnostic_plot = self.config.qc_dir / f"{self.base_name}_correction_diagnostic.png"

            cmd_diagnostic = [
                self._get_hic_tool_path('hicCorrectMatrix'),
                'diagnostic_plot',
                '--matrix', str(self.raw_matrix),
                '--plotName', str(diagnostic_plot)
            ]

            self.logger.info("生成矫正诊断图|Generating correction diagnostic plot...")
            self._log_command(cmd_diagnostic, "hicCorrectMatrix诊断图|hicCorrectMatrix diagnostic_plot")
            subprocess.run(cmd_diagnostic, check=True, capture_output=True, text=True)

            # 执行矫正|Perform correction
            cmd_correct = [
                self._get_hic_tool_path('hicCorrectMatrix'),
                'correct',
                '--correctionMethod', self.config.correction_method,
                '--matrix', str(self.raw_matrix),
                '--outFileName', str(self.corrected_matrix)
            ]

            # 添加过滤阈值参数|Add filter threshold parameters (仅ICE使用|Only for ICE)
            if self.config.correction_method == 'ICE':
                if self.config.filter_threshold:
                    min_threshold, max_threshold = self.config.filter_threshold
                    cmd_correct.extend(['--filterThreshold', str(min_threshold), str(max_threshold)])
                    self.logger.info(f"使用过滤阈值|Using filter thresholds: min={min_threshold}, max={max_threshold}")

                # 添加inflation cutoff参数|Add inflation cutoff parameter (仅ICE使用|Only for ICE)
                if self.config.inflation_cutoff is not None:
                    cmd_correct.extend(['--inflationCutoff', str(self.config.inflation_cutoff)])
                    self.logger.info(f"使用inflation cutoff|Using inflation cutoff: {self.config.inflation_cutoff}")

                # 添加sequenced count cutoff参数|Add sequenced count cutoff parameter (仅ICE使用|Only for ICE)
                if self.config.sequenced_count_cutoff is not None:
                    cmd_correct.extend(['--sequencedCountCutoff', str(self.config.sequenced_count_cutoff)])
                    self.logger.info(f"使用sequenced count cutoff|Using sequenced count cutoff: {self.config.sequenced_count_cutoff}")

            self.logger.info(f"执行{self.config.correction_method}矫正|Performing {self.config.correction_method} correction...")
            self._log_command(cmd_correct, f"hicCorrectMatrix矫正 ({self.config.correction_method})|hicCorrectMatrix correct")
            subprocess.run(cmd_correct, check=True, capture_output=True, text=True)

            self.logger.info(f"矩阵已矫正|Matrix corrected: {self.corrected_matrix}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"矩阵矫正失败|Matrix correction failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False

    def _plot_heatmap(self) -> bool:
        """绘制Hi-C热图|Plot Hi-C heatmap

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("步骤7: 绘制Hi-C热图|Step 7: Plotting Hi-C heatmap")

        # 检查是否已完成|Check if already done
        if self.config.skip_existing and self.heatmap_file.exists() and self.heatmap_file.stat().st_size > 0:
            self.logger.info("跳过|Skipping: 热图已存在|Heatmap exists")
            return True

        try:
            cmd = [
                self._get_hic_tool_path('hicPlotMatrix'),
                '--matrix', str(self.corrected_matrix),
                '--outFileName', str(self.heatmap_file),
                '--title', f'Hi-C Heatmap - {self.base_name}',
                '--scoreName', 'Contact Frequency',
                '--colorMap', self.config.color_map,
                '--dpi', str(self.config.dpi)
            ]

            # 添加对数变换参数|Add log transformation
            if self.config.log1p:
                cmd.append('--log1p')

            self.logger.info("绘制热图中|Plotting heatmap...")
            self._log_command(cmd, "hicPlotMatrix绘制热图|hicPlotMatrix")
            subprocess.run(cmd, check=True, capture_output=True, text=True)

            self.logger.info(f"热图已生成|Heatmap created: {self.heatmap_file}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"绘制热图失败|Failed to plot heatmap: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False
