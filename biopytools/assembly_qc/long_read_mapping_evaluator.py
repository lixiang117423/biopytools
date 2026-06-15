"""
三代数据Mapping评估模块|Long-read Mapping Evaluation Module

使用minimap2将三代测序数据（PacBio/Nanopore/HiFi）mapping到基因组
Use minimap2 to map third-generation sequencing data (PacBio/Nanopore/HiFi) to genome
"""

import os
import subprocess
import glob
import re
from pathlib import Path
from typing import Dict, Any, Optional, List


class LongReadMappingEvaluator:
    """三代数据Mapping评估器|Long-read Mapping Evaluator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.working_dir = Path(self.config.long_read_mapping_output_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)

        # 创建BAM输出子目录|Create BAM output subdirectory
        self.bam_dir = self.working_dir / "bam_files"
        self.bam_dir.mkdir(parents=True, exist_ok=True)

        # 使用直接工具路径（避免conda run的管道问题）|Use direct tool paths
        self.minimap2_cmd = [self.config.minimap2_path]
        self.samtools_cmd = [self.config.samtools_path]

    def evaluate(self) -> Optional[Dict[str, Any]]:
        """运行三代数据Mapping评估|Run long-read mapping evaluation"""
        if not self.config.enable_long_read_mapping:
            self.logger.info("跳过三代数据Mapping评估|Skipping long-read mapping evaluation")
            return None

        self.logger.info("开始三代数据Mapping评估|Starting long-read mapping evaluation")

        # 获取reads目录|Get reads directory
        reads_dir = self.config.long_reads
        if not reads_dir:
            self.logger.error("未提供三代数据Mapping所需的reads数据|No long-read data provided for mapping evaluation")
            return None

        self.logger.info("开始三代数据Mapping评估|Starting long-read mapping evaluation")

        # 发现样品|Discover samples
        samples = self._discover_samples(reads_dir)

        if not samples:
            self.logger.error(f"未找到有效的三代数据FASTQ文件|No valid long-read FASTQ files found in: {reads_dir}")
            return None

        self.logger.info(f"发现 {len(samples)} 个样品|Found {len(samples)} sample(s)")

        # 处理每个样品|Process each sample
        all_results = {}
        for sample_name, fastq_files in samples.items():
            self.logger.info(f"处理样品|Processing sample: {sample_name}")

            result = self._process_sample(sample_name, fastq_files)
            if result:
                all_results[sample_name] = result

        if not all_results:
            self.logger.error("所有样品处理失败|All samples failed")
            return None

        # 汇总结果|Summarize results
        return self._summarize_results(all_results)

    def _discover_samples(self, reads_dir: str) -> Dict[str, List[str]]:
        """
        发现三代数据样品|Discover long-read samples

        Args:
            reads_dir: reads目录|Reads directory

        Returns:
            样品名字典|Sample name dictionary: {sample_name: [fastq_files]}
        """
        self.logger.info(f"扫描三代数据reads目录|Scanning long-read reads directory: {reads_dir}")

        reads_path = Path(reads_dir)
        if not reads_path.exists():
            self.logger.error(f"reads目录不存在|Reads directory not found: {reads_dir}")
            return {}

        # 查找FASTQ文件|Find FASTQ files
        pattern = self.config.long_read_mapping_pattern
        fastq_files = list(reads_path.glob(pattern))

        if not fastq_files:
            self.logger.error(f"未找到匹配模式的FASTQ文件|No FASTQ files found matching pattern: {pattern}")
            return {}

        self.logger.info(f"找到 {len(fastq_files)} 个FASTQ文件|Found {len(fastq_files)} FASTQ files")

        # 按样品名分组|Group by sample name
        samples = {}
        for fastq_file in fastq_files:
            # 提取样品名|Extract sample name
            # 文件名格式：sample.fq.gz, sample.fastq.gz, sample_1.fastq.gz等
            # 去除扩展名和测序批次标识|Remove extensions and sequencing batch identifiers
            sample_name = self._extract_sample_name(fastq_file)

            if sample_name not in samples:
                samples[sample_name] = []
            samples[sample_name].append(str(fastq_file))

        self.logger.info(f"识别到 {len(samples)} 个独立样品|Identified {len(samples)} unique samples")
        for sample_name in sorted(samples.keys()):
            self.logger.debug(f"  样品|Sample: {sample_name}, 文件数|File count: {len(samples[sample_name])}")

        return samples

    def _extract_sample_name(self, fastq_file: Path) -> str:
        """
        从FASTQ文件名提取样品名|Extract sample name from FASTQ filename

        Args:
            fastq_file: FASTQ文件路径|FASTQ file path

        Returns:
            样品名|Sample name
        """
        # 去除扩展名|Remove extensions
        name = fastq_file.name
        for ext in ['.fq.gz', '.fastq.gz', '.fasta.gz', '.fa.gz', '.fq', '.fastq', '.fasta', '.fa']:
            if name.endswith(ext):
                name = name[:-len(ext)]
                break

        # 去除测序批次标识（如_1, _2等）|Remove sequencing batch identifiers (e.g., _1, _2)
        name = re.sub(r'_[\d]+$', '', name)

        return name

    def _process_sample(self, sample_name: str, fastq_files: List[str]) -> Optional[Dict[str, Any]]:
        """
        处理单个三代数据样品|Process single long-read sample

        Args:
            sample_name: 样品名|Sample name
            fastq_files: FASTQ文件列表|List of FASTQ files

        Returns:
            Mapping统计结果|Mapping statistics results
        """
        # 定义输出文件|Define output files
        sam_file = self.bam_dir / f"{sample_name}.sam"
        bam_file = self.bam_dir / f"{sample_name}.bam"
        sorted_bam = self.bam_dir / f"{sample_name}.sorted.bam"
        flagstat_file = self.bam_dir / f"{sample_name}.flagstat.txt"
        coverage_file = self.bam_dir / f"{sample_name}.coverage.txt"

        # 检查断点续传|Check resume
        if self.config.resume and sorted_bam.exists() and flagstat_file.exists() and coverage_file.exists():
            self.logger.info(f"跳过已完成的样品|Skipping completed sample: {sample_name}")
            return self._parse_flagstat_and_coverage(flagstat_file, coverage_file)

        # Step 1: Minimap2 mapping（直接输出BAM）
        # Minimap2 mapping (direct BAM output)
        if not self._run_minimap2(sample_name, fastq_files, bam_file):
            return None

        # Step 2: 生成基础统计|Generate basic statistics
        if not self._run_flagstat(bam_file, flagstat_file):
            return None

        # Step 3: 排序BAM（用于覆盖度计算）|Sort BAM (for coverage calculation)
        if not self._sort_bam(bam_file, sorted_bam):
            return None

        # Step 4: 构建索引（用于覆盖度计算）|Build index (for coverage calculation)
        if not self._build_bam_index(sorted_bam):
            return None

        # Step 5: 计算覆盖比例|Calculate coverage fraction
        if not self._calculate_coverage(sorted_bam, coverage_file):
            return None

        # Step 6: 清理中间文件|Cleanup intermediate files
        self._cleanup_intermediate_files(sam_file, bam_file)

        # 解析统计结果|Parse statistics results
        return self._parse_flagstat_and_coverage(flagstat_file, coverage_file)

    def _run_minimap2(
        self,
        sample_name: str,
        fastq_files: List[str],
        bam_file: Path
    ) -> bool:
        """
        运行Minimap2 mapping，直接输出BAM|Run Minimap2 mapping with direct BAM output

        Args:
            sample_name: 样品名|Sample name
            fastq_files: FASTQ文件列表|List of FASTQ files
            bam_file: 输出BAM文件|Output BAM file

        Returns:
            成功返回True，失败返回False|True on success, False on failure
        """
        # 根据数据类型选择preset|Select preset based on data type
        preset = self._get_preset()

        # Minimap2命令|Minimap2 command
        minimap2_args = [
            "-ax", preset,  # alignment mode + preset
            "-t", str(self.config.long_read_mapping_threads),
            self.config.genome,
        ] + fastq_files

        minimap2_cmd = self.minimap2_cmd + minimap2_args

        # samtools命令：从stdin读取SAM，输出BAM
        # samtools command: read SAM from stdin, output BAM
        samtools_args = [
            "view",
            "-@", str(self.config.long_read_mapping_threads),
            "-b",  # BAM格式|BAM format
            "-h",  # 在header中包含@PG|Include @PG in header
            "-o", str(bam_file),
            "-"  # 从stdin读取|Read from stdin
        ]
        samtools_cmd = self.samtools_cmd + samtools_args

        try:
            self.logger.info(f"三代数据Mapping（管道输出BAM）|Long-read mapping (piped to BAM):")
            self.logger.info(f"  样品|Sample: {sample_name}")
            self.logger.info(f"  数据类型|Data type: {self.config.long_read_type}")
            self.logger.info(f"  Minimap2命令|Minimap2 command: {' '.join(minimap2_cmd)}")
            self.logger.info(f"  SAMtools命令|SAMtools command: {' '.join(samtools_cmd)}")

            # 使用管道：minimap2 | samtools view -b
            # Use pipeline: minimap2 | samtools view -b
            process_minimap2 = subprocess.Popen(
                minimap2_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            process_samtools = subprocess.Popen(
                samtools_cmd,
                stdin=process_minimap2.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # 关闭minimap2的stdout，让samtools能接收到EOF
            # Close minimap2's stdout to allow samtools to receive EOF
            process_minimap2.stdout.close()

            # 等待两个进程完成|Wait for both processes to complete
            minimap2_stdout, minimap2_stderr = process_minimap2.communicate()
            samtools_stdout, samtools_stderr = process_samtools.communicate()

            minimap2_returncode = process_minimap2.returncode
            samtools_returncode = process_samtools.returncode

            if minimap2_returncode != 0:
                stderr_msg = minimap2_stderr if minimap2_stderr else ""
                self.logger.error(f"Minimap2 mapping失败|Minimap2 mapping failed: {stderr_msg}")
                return False

            if samtools_returncode != 0:
                stderr_msg = samtools_stderr if samtools_stderr else ""
                self.logger.error(f"SAMtools view失败|SAMtools view failed: {stderr_msg}")
                return False

            return True

        except Exception as e:
            self.logger.error(f"Minimap2运行异常|Unexpected error running Minimap2: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _get_preset(self) -> str:
        """
        根据数据类型获取minimap2 preset|Get minimap2 preset based on data type

        Returns:
            preset字符串|Preset string
        """
        data_type = self.config.long_read_type.lower()

        # 根据minimap2文档的preset选择|Preset selection based on minimap2 documentation
        preset_map = {
            "pacbio": "map-pb",      # PacBio CLR genomic reads
            "pb": "map-pb",
            "ont": "map-ont",        # Oxford Nanopore genomic reads
            "nanopore": "map-ont",
            "hifi": "map-hifi",      # PacBio HiFi/CCS genomic reads
            "ccs": "map-hifi",
        }

        return preset_map.get(data_type, "map-ont")  # 默认使用ont|Default to ont

    def _run_flagstat(self, bam_file: Path, flagstat_file: Path) -> bool:
        """运行samtools flagstat|Run samtools flagstat"""
        args = ["flagstat", str(bam_file)]
        cmd = self.samtools_cmd + args

        try:
            with open(flagstat_file, 'w') as f:
                result = subprocess.run(
                    cmd,
                    stdout=f,
                    stderr=subprocess.PIPE,
                    text=True,
                    check=False
                )

            if result.returncode != 0:
                self.logger.error(f"flagstat运行失败|flagstat failed: {result.stderr}")
                return False

            return True

        except Exception as e:
            self.logger.error(f"flagstat运行异常|Error running flagstat: {e}")
            return False

    def _sort_bam(self, input_bam: Path, output_bam: Path) -> bool:
        """排序BAM|Sort BAM"""
        args = [
            "sort",
            "-@", str(self.config.long_read_mapping_threads),
            "-o", str(output_bam),
            str(input_bam)
        ]
        cmd = self.samtools_cmd + args

        try:
            self.logger.info(f"排序BAM|Sorting BAM: {output_bam.name}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                self.logger.error(f"BAM排序失败|BAM sorting failed: {result.stderr}")
                return False
            self.logger.info(f"BAM排序完成|BAM sorting completed: {output_bam.name}")
            return True
        except Exception as e:
            self.logger.error(f"BAM排序异常|Error sorting BAM: {e}")
            return False

    def _build_bam_index(self, bam_file: Path) -> bool:
        """构建BAM索引|Build BAM index"""
        args = [
            "index",
            "-@", str(self.config.long_read_mapping_threads),
            str(bam_file)
        ]
        cmd = self.samtools_cmd + args

        try:
            self.logger.info(f"构建BAM索引|Building BAM index: {bam_file.name}")
            result = subprocess.run(cmd, capture_output=True, text=True, check=False)
            if result.returncode != 0:
                self.logger.error(f"BAM索引构建失败|BAM indexing failed: {result.stderr}")
                return False
            self.logger.info(f"BAM索引构建完成|BAM indexing completed: {bam_file.name}")
            return True
        except Exception as e:
            self.logger.error(f"BAM索引构建异常|Error building BAM index: {e}")
            return False

    def _calculate_coverage(self, bam_file: Path, coverage_file: Path) -> bool:
        """
        计算覆盖比例|Calculate coverage fraction
        使用samtools depth统计深度>0的位点比例
        Use samtools depth to calculate fraction of positions with depth > 0
        """
        args = ["depth", "-@", str(self.config.long_read_mapping_threads), str(bam_file)]
        cmd = self.samtools_cmd + args

        try:
            self.logger.info(f"计算覆盖比例|Calculating coverage fraction: {bam_file.name}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode != 0:
                self.logger.error(f"覆盖度计算失败|Coverage calculation failed: {result.stderr}")
                return False

            # 解析samtools depth输出|Parse samtools depth output
            # 格式：chr position depth
            total_positions = 0
            covered_positions = 0

            for line in result.stdout.split('\n'):
                if not line.strip():
                    continue
                parts = line.split('\t')
                if len(parts) >= 3:
                    total_positions += 1
                    depth = int(parts[2])
                    if depth > 0:
                        covered_positions += 1

            if total_positions == 0:
                self.logger.error(f"未找到覆盖度数据|No coverage data found")
                return False

            coverage_fraction = (covered_positions / total_positions) * 100

            # 保存结果|Save results
            with open(coverage_file, 'w') as f:
                f.write(f"{coverage_fraction:.2f}\n")

            self.logger.info(f"覆盖比例|Coverage fraction: {coverage_fraction:.2f}% ({covered_positions}/{total_positions} positions)")
            return True

        except Exception as e:
            self.logger.error(f"覆盖度计算异常|Error calculating coverage: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _parse_flagstat_and_coverage(self, flagstat_file: Path, coverage_file: Path) -> Optional[Dict[str, Any]]:
        """
        解析flagstat和coverage输出|Parse flagstat and coverage output
        同时获取mapping统计和覆盖比例|Get both mapping stats and coverage fraction
        """
        try:
            # 解析flagstat|Parse flagstat
            with open(flagstat_file, 'r') as f:
                lines = f.readlines()

            stats = {}

            # 解析flagstat输出|Parse flagstat output
            for line in lines:
                if 'in total' in line:
                    stats['total_reads'] = int(line.split()[0])
                elif 'mapped (' in line and 'primary' not in line:
                    stats['mapped_reads'] = int(line.split()[0])
                    percentage_part = line.split('(')[1].split('%')[0]
                    stats['mapping_rate'] = float(percentage_part)

            # 解析coverage|Parse coverage
            if coverage_file.exists():
                with open(coverage_file, 'r') as f:
                    coverage_line = f.read().strip()
                    if coverage_line:
                        stats['coverage_fraction'] = float(coverage_line)
                        self.logger.info(f"覆盖比例|Coverage fraction: {stats['coverage_fraction']:.2f}%")
                    else:
                        self.logger.warning("覆盖度文件为空|Coverage file is empty")
            else:
                self.logger.warning(f"覆盖度文件不存在|Coverage file not found: {coverage_file}")

            self.logger.info(f"Mapping统计|Mapping statistics:")
            self.logger.info(f"  总reads数|Total reads: {stats.get('total_reads', 'N/A')}")
            self.logger.info(f"  Mapped reads|Mapped reads: {stats.get('mapped_reads', 'N/A')}")
            self.logger.info(f"  Mapping rate|Mapping rate: {stats.get('mapping_rate', 'N/A'):.2f}%")

            return stats

        except Exception as e:
            self.logger.error(f"解析flagstat/coverage文件失败|Failed to parse flagstat/coverage files: {e}")
            return None

    def _cleanup_intermediate_files(self, sam_file: Path, raw_bam: Path):
        """清理中间文件|Cleanup intermediate files"""
        try:
            if sam_file.exists():
                sam_file.unlink()
            if raw_bam.exists():
                raw_bam.unlink()
        except Exception as e:
            self.logger.warning(f"清理中间文件失败|Failed to cleanup intermediate files: {e}")

    def _summarize_results(self, all_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """
        汇总所有样品的结果|Summarize results from all samples

        Args:
            all_results: 所有样品结果|All sample results

        Returns:
            汇总结果|Summary results
        """
        total_reads = sum(r.get('total_reads', 0) for r in all_results.values())
        total_mapped = sum(r.get('mapped_reads', 0) for r in all_results.values())

        overall_mapping_rate = (total_mapped / total_reads * 100) if total_reads > 0 else 0

        summary = {
            'total_reads': total_reads,
            'mapped_reads': total_mapped,
            'overall_mapping_rate': overall_mapping_rate,
            'samples': all_results
        }

        self.logger.info(f"总体Mapping统计|Overall mapping statistics:")
        self.logger.info(f"  总reads数|Total reads: {total_reads}")
        self.logger.info(f"  Mapped reads|Mapped reads: {total_mapped}")
        self.logger.info(f"  Overall mapping rate|Overall mapping rate: {overall_mapping_rate:.2f}%")

        return summary
