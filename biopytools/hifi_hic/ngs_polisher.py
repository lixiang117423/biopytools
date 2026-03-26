"""
NGS数据polish模块|NGS Data Polish Module

使用二代数据筛选高质量HiFi reads并重新组装
Use second-generation data to filter high-quality HiFi reads and reassemble
"""

import os
import subprocess
from pathlib import Path


class NGSPolisher:
    """NGS数据polish器|NGS Data Polisher"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def run_ngs_polish(self, completed_steps=None):
        """
        运行完整的NGS polish流程|Run complete NGS polish pipeline

        Args:
            completed_steps: 已完成步骤字典|Completed steps dictionary (for resume)

        Returns:
            bool: 是否成功|Whether successful
        """
        if completed_steps is None:
            completed_steps = {}

        try:
            self.logger.info("=" * 80)
            self.logger.info("开始NGS polish流程|Starting NGS Polish Pipeline")
            self.logger.info("=" * 80)

            # 1. BWA比对|BWA alignment
            if not completed_steps.get('bwa_alignment', False):
                if not self._run_bwa_alignment():
                    return False
            else:
                self.logger.info("检测到BWA比对已完成，跳过|BWA alignment already completed, skipping")

            # 2. Coverage filter筛选高质量contigs|Coverage filter for high-quality contigs
            if not completed_steps.get('coverage_filter', False):
                high_quality_list = self._run_coverage_filter()
                if not high_quality_list:
                    return False
            else:
                self.logger.info("检测到覆盖度过滤已完成，跳过|Coverage filter already completed, skipping")
                # 获取已生成的高质量contig列表|Get already generated high-quality contig list
                high_quality_list = os.path.join(self.config.ngs_polish_dir,
                                                 "02.coverage_filter",
                                                 f"{self.config.prefix}_high_quality.list")

            # 3. 提取高质量reads|Extract high-quality reads
            if not completed_steps.get('filtered_reads', False):
                filtered_reads = self._extract_high_quality_reads(high_quality_list)
                if not filtered_reads:
                    return False
            else:
                self.logger.info("检测到reads筛选已完成，跳过|Reads filtering already completed, skipping")
                filtered_reads = os.path.join(self.config.ngs_polish_dir,
                                             "03.filtered_reads",
                                             f"{self.config.prefix}_high_quality_reads.fq.gz")

            # 4. 重新组装|Reassembly
            if not completed_steps.get('reassembly', False):
                if not self._run_reassembly(filtered_reads):
                    return False
            else:
                self.logger.info("检测到重新组装已完成，跳过|Reassembly already completed, skipping")

            self.logger.info("=" * 80)
            self.logger.info("NGS polish流程完成|NGS Polish Pipeline Completed")
            self.logger.info("=" * 80)

            return True

        except Exception as e:
            self.logger.error(f"NGS polish流程失败|NGS polish pipeline failed: {str(e)}")
            return False

    def _run_bwa_alignment(self):
        """
        运行BWA比对|Run BWA alignment

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("[1/4] 运行BWA比对|Running BWA alignment")

        try:
            # 导入BWA模块|Import BWA module
            from biopytools.bwa.main import BWAAligner

            # 定义输出目录|Define output directory
            bwa_output_dir = os.path.join(self.config.ngs_polish_dir, "01.bwa_alignment")

            # 创建BWA比对器|Create BWA aligner
            bwa_aligner = BWAAligner(
                genome=os.path.join(self.config.fasta_dir, f"{self.config.prefix}.primary.fa"),
                input_dir=self.config.ngs_data,
                pattern=self.config.ngs_pattern,
                output_dir=bwa_output_dir,
                threads=self.config.threads
            )

            # 运行比对|Run alignment
            success = bwa_aligner.run_analysis()

            if success:
                self.logger.info(f"BWA比对完成|BWA alignment completed: {bwa_output_dir}")
                return True
            else:
                self.logger.error("BWA比对失败|BWA alignment failed")
                return False

        except ImportError as e:
            self.logger.error(f"无法导入BWA模块|Cannot import BWA module: {e}")
            return False
        except Exception as e:
            self.logger.error(f"BWA比对出错|BWA alignment error: {e}")
            return False

    def _run_coverage_filter(self):
        """
        运行覆盖度过滤|Run coverage filter

        Returns:
            str: 高质量contig列表文件路径|High-quality contig list file path
        """
        self.logger.info("[2/4] 运行覆盖度过滤|Running coverage filter")

        try:
            # 导入coverage-filter模块|Import coverage-filter module
            from biopytools.coverage_filter.main import CoverageFilter

            # 定义输出路径|Define output paths
            filter_output_dir = os.path.join(self.config.ngs_polish_dir, "02.coverage_filter")
            Path(filter_output_dir).mkdir(parents=True, exist_ok=True)

            # 获取BAM文件|Get BAM file
            bwa_output_dir = os.path.join(self.config.ngs_polish_dir, "01.bwa_alignment")
            bam_dir = os.path.join(bwa_output_dir, "bam")

            # 查找实际生成的BAM文件（以样本名命名）|Find actual BAM files (named by sample)
            import glob
            bam_files = glob.glob(os.path.join(bam_dir, "*.bam"))

            if not bam_files:
                self.logger.error(f"未找到任何BAM文件|No BAM files found in: {bam_dir}")
                return None

            # 如果有多个BAM文件，合并它们|If multiple BAM files, merge them
            if len(bam_files) == 1:
                bam_file = bam_files[0]
                self.logger.info(f"使用BAM文件|Using BAM file: {bam_file}")
            else:
                # 合并多个BAM文件|Merge multiple BAM files
                merged_bam = os.path.join(bam_dir, f"{self.config.prefix}.merged.bam")
                self.logger.info(f"发现{len(bam_files)}个BAM文件，正在合并|Found {len(bam_files)} BAM files, merging...")

                # 使用samtools merge合并|Merge using samtools merge
                bam_list = " ".join(bam_files)
                cmd = f"samtools merge -@ {self.config.threads} - {bam_list} | samtools sort -@ {self.config.threads} -o {merged_bam}"
                if not self._run_command(cmd):
                    self.logger.error("BAM文件合并失败|Failed to merge BAM files")
                    return None

                # 构建索引|Build index
                cmd = f"samtools index -@ {self.config.threads} {merged_bam}"
                if not self._run_command(cmd):
                    self.logger.error("合并BAM索引构建失败|Failed to index merged BAM")
                    return None

                bam_file = merged_bam
                self.logger.info(f"BAM文件合并完成|BAM files merged: {bam_file}")

            if not os.path.exists(bam_file):
                self.logger.error(f"BAM文件不存在|BAM file not found: {bam_file}")
                return None

            # 创建coverage filter|Create coverage filter
            coverage_filter = CoverageFilter(
                bam_file=bam_file,
                fasta_file=os.path.join(self.config.fasta_dir, f"{self.config.prefix}.primary.fa"),
                output_prefix=self.config.prefix,
                output_dir=filter_output_dir,
                threads=self.config.threads,
                high_coverage=self.config.high_cov,
                medium_cov_min=self.config.medium_cov_min
            )

            # 运行过滤|Run filter
            success = coverage_filter.run_filter()

            if success:
                high_quality_list = os.path.join(filter_output_dir,
                                                 f"{self.config.prefix}_high_quality.list")
                if os.path.exists(high_quality_list):
                    self.logger.info(f"覆盖度过滤完成|Coverage filter completed: {high_quality_list}")
                    return high_quality_list
                else:
                    self.logger.error("高质量contig列表文件不存在|High-quality contig list file not found")
                    return None
            else:
                self.logger.error("覆盖度过滤失败|Coverage filter failed")
                return None

        except ImportError as e:
            self.logger.error(f"无法导入coverage-filter模块|Cannot import coverage-filter module: {e}")
            return None
        except Exception as e:
            self.logger.error(f"覆盖度过滤出错|Coverage filter error: {e}")
            return None

    def _extract_high_quality_reads(self, high_quality_list):
        """
        提取高质量reads|Extract high-quality reads

        Args:
            high_quality_list: 高质量contig列表文件|High-quality contig list file

        Returns:
            str: 筛选后的reads文件路径|Filtered reads file path
        """
        self.logger.info("[3/4] 提取高质量HiFi reads|Extracting high-quality HiFi reads")

        try:
            # 定义输出路径|Define output paths
            filtered_reads_dir = os.path.join(self.config.ngs_polish_dir, "03.filtered_reads")
            Path(filtered_reads_dir).mkdir(parents=True, exist_ok=True)

            read_names_file = os.path.join(filtered_reads_dir, "high_quality_read_names.txt")
            filtered_reads = os.path.join(filtered_reads_dir,
                                          f"{self.config.prefix}_high_quality_reads.fq.gz")

            # 获取contig-reads映射文件|Get contig-reads mapping file
            # 注意：文件名格式是 {prefix}.p_ctg.contig_reads.tsv（去掉了.hic.和.bp.）
            # Note: filename format is {prefix}.p_ctg.contig_reads.tsv (removed .hic. and .bp.)
            # assembler.py会去掉.hic.前缀，所以这里统一使用不带.hic.的文件名
            # assembler.py removes .hic. prefix, so we use the filename without .hic. here
            mapping_file = os.path.join(self.config.fasta_dir,
                                       f"{self.config.prefix}.p_ctg.contig_reads.tsv")

            if not os.path.exists(mapping_file):
                self.logger.error(f"Contig-reads映射文件不存在|Contig-reads mapping file not found: {mapping_file}")
                return None

            # 从映射文件中提取read names|Extract read names from mapping file
            self.logger.info(f"从contig-reads映射中提取read names|Extracting read names from contig-reads mapping")
            cmd = f"grep -Fwf {high_quality_list} {mapping_file} | cut -f2 > {read_names_file}"
            if not self._run_command(cmd):
                return None

            # 统计read数量|Count reads
            read_count = subprocess.check_output(f"wc -l {read_names_file}", shell=True).decode().strip().split()[0]
            self.logger.info(f"找到|Found {read_count} 个高质量reads|high-quality reads")

            # 使用seqkit提取reads|Extract reads using seqkit
            self.logger.info(f"使用seqkit提取reads|Extracting reads using seqkit")
            # 使用多线程加速|Use multithreading to speed up
            # cmd = f"seqkit grep -n -f {read_names_file} {self.config.hifi_data} -j {self.config.threads} -o {filtered_reads}"
            cmd = f"seqkit grep -f {read_names_file} {self.config.hifi_data} -j {self.config.threads} -o {filtered_reads}"
            if not self._run_command(cmd):
                return None

            # 检查输出文件|Check output file
            if os.path.exists(filtered_reads) and os.path.getsize(filtered_reads) > 0:
                self.logger.info(f"高质量reads提取完成|High-quality reads extracted: {filtered_reads}")
                return filtered_reads
            else:
                self.logger.error("筛选后的reads文件为空|Filtered reads file is empty")
                return None

        except Exception as e:
            self.logger.error(f"提取高质量reads失败|Failed to extract high-quality reads: {e}")
            return None

    def _run_reassembly(self, filtered_reads):
        """
        重新组装|Reassembly

        Args:
            filtered_reads: 筛选后的reads文件|Filtered reads file

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("[4/4] 使用高质量reads重新组装|Reassembly with high-quality reads")

        try:
            # 导入assembler模块|Import assembler module
            from assembler import HifiasmAssembler

            # 定义输出目录|Define output directory
            reassembly_dir = os.path.join(self.config.ngs_polish_dir, "04.reassembly")
            Path(reassembly_dir).mkdir(parents=True, exist_ok=True)

            # 创建临时配置用于重新组装|Create temporary config for reassembly
            import copy
            reassembly_config = copy.copy(self.config)
            # 重要：更新hifi_data为筛选后的reads|Important: update hifi_data to filtered reads
            reassembly_config.hifi_data = filtered_reads
            reassembly_config.raw_dir = os.path.join(reassembly_dir, "01.raw_output")
            reassembly_config.fasta_dir = os.path.join(reassembly_dir, "02.fasta")
            Path(reassembly_config.raw_dir).mkdir(parents=True, exist_ok=True)
            Path(reassembly_config.fasta_dir).mkdir(parents=True, exist_ok=True)

            self.logger.info(f"使用筛选后的reads进行重新组装|Reassembly with filtered reads: {filtered_reads}")

            # 创建assembler|Create assembler
            assembler = HifiasmAssembler(reassembly_config, self.logger)

            # 运行组装|Run assembly
            if not assembler.run_assembly():
                self.logger.error("重新组装失败|Reassembly failed")
                return False

            # 转换GFA到FASTA|Convert GFA to FASTA
            fasta_results = assembler.convert_gfa_to_fasta()

            # 生成contig-reads映射|Generate contig-reads mapping
            assembler.generate_contig_reads_mapping()

            # 将最终结果复制到ngs_polish目录|Copy final results to ngs_polish directory
            self._copy_final_results(reassembly_dir)

            self.logger.info(f"重新组装完成|Reassembly completed")
            self.logger.info(f"最终结果保存在|Final results saved to: {reassembly_dir}")

            return True

        except Exception as e:
            self.logger.error(f"重新组装失败|Reassembly failed: {e}")
            return False

    def _copy_final_results(self, reassembly_dir):
        """
        复制最终结果到ngs_polish目录|Copy final results to ngs_polish directory

        Args:
            reassembly_dir: 重新组装目录|Reassembly directory
        """
        self.logger.info("统计最终组装结果|Summarizing final assembly results")

        # 所有文件保留在 04.reassembly/02.fasta/ 目录中
        # All files remain in 04.reassembly/02.fasta/ directory
        reassembly_fasta_dir = os.path.join(reassembly_dir, "02.fasta")
        fasta_files = [
            f"{self.config.prefix}.primary.fa",
            f"{self.config.prefix}.hap1.fa",
            f"{self.config.prefix}.hap2.fa",
            f"{self.config.prefix}.alternate.fa"
        ]

        existing_files = []
        for fasta_name in fasta_files:
            fasta_path = os.path.join(reassembly_fasta_dir, fasta_name)
            if os.path.exists(fasta_path):
                existing_files.append(fasta_name)
                # 获取文件大小|Get file size
                size_mb = os.path.getsize(fasta_path) / (1024 * 1024)
                self.logger.info(f"  ✓ {fasta_name} ({size_mb:.1f} MB)")

        self.logger.info(f"共生成|Total generated: {len(existing_files)}/{len(fasta_files)} 个FASTA文件")
        self.logger.info(f"文件位置|Files location: {reassembly_fasta_dir}")

    def _run_command(self, cmd):
        """
        执行命令|Run command

        Args:
            cmd: 命令字符串|Command string

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            self.logger.info(f"执行命令|Executing: {cmd}")
            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )
            if result.stdout:
                self.logger.info(result.stdout.strip())
            return True
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败|Command failed: {e}")
            if e.stderr:
                self.logger.error(f"错误信息|Error: {e.stderr}")
            return False
