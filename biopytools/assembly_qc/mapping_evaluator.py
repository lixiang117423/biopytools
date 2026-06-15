"""
Mapping评估模块（完全独立实现）|Mapping Evaluation Module (Fully Independent Implementation)

直接调用bwa mem/samtools命令行工具
Directly calls bwa mem/samtools command-line tools
"""

import os
import subprocess
import glob
import re
from pathlib import Path
from typing import Dict, Any, Optional, List


class MappingEvaluator:
    """Mapping评估器|Mapping Evaluator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.working_dir = Path(self.config.mapping_output_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)

        # 创建BAM输出子目录|Create BAM output subdirectory
        self.bam_dir = self.working_dir / "bam_files"
        self.bam_dir.mkdir(parents=True, exist_ok=True)

        # 使用直接工具路径（避免conda run的管道问题）|Use direct tool paths (avoid conda run pipe issues)
        self.bwa_cmd = [self.config.bwa_path]
        self.samtools_cmd = [self.config.samtools_path]

    def evaluate(self) -> Optional[Dict[str, Any]]:
        """运行Mapping评估|Run mapping evaluation"""
        if not self.config.enable_mapping or self.config.skip_mapping:
            self.logger.info("跳过Mapping评估|Skipping mapping evaluation")
            return None

        reads_dir = self.config.ngs_reads
        if not reads_dir:
            self.logger.error("未提供Mapping评估所需的NGS reads数据|No NGS reads data provided for mapping evaluation")
            return None

        self.logger.info("开始Mapping评估|Starting mapping evaluation")

        # 发现样品|Discover samples
        samples = self._discover_samples(reads_dir)

        if not samples:
            self.logger.error(f"未找到有效的FASTQ文件对|No valid FASTQ pairs found in: {reads_dir}")
            return None

        self.logger.info(f"发现 {len(samples)} 个样品|Found {len(samples)} sample(s)")

        # 检查并构建基因组索引|Check and build genome index
        if not self._check_and_build_index():
            return None

        # 运行比对和统计|Run alignment and statistics
        try:
            results = {}
            total_reads = 0
            total_mapped = 0

            for sample_name, read1, read2 in samples:
                sample_result = self._process_sample(sample_name, read1, read2)
                if sample_result:
                    results[sample_name] = sample_result
                    total_reads += sample_result['total_reads']
                    total_mapped += sample_result['mapped_reads']

            # 汇总统计|Summary statistics
            if results:
                overall_mapping_rate = (total_mapped / total_reads * 100) if total_reads > 0 else 0.0

                # 计算平均覆盖度（所有样品coverage_fraction的平均值）|Calculate mean coverage (average of all samples' coverage_fraction)
                coverage_fractions = [r['coverage_fraction'] for r in results.values() if r.get('coverage_fraction', 0) > 0]
                mean_coverage = sum(coverage_fractions) / len(coverage_fractions) if coverage_fractions else 0.0

                summary = {
                    'total_samples': len(results),
                    'total_reads': total_reads,
                    'total_mapped': total_mapped,
                    'overall_mapping_rate': overall_mapping_rate,
                    'mean_coverage': mean_coverage,
                    'samples': results,
                }

                self.logger.info(f"Mapping统计|Mapping statistics:")
                self.logger.info(f"  样品数|Samples: {summary['total_samples']}")
                self.logger.info(f"  总reads数|Total reads: {total_reads:,}")
                self.logger.info(f"  总比对reads数|Total mapped reads: {total_mapped:,}")
                self.logger.info(f"  总比对率|Overall mapping rate: {overall_mapping_rate:.2f}%")

                return summary
            else:
                return None

        except Exception as e:
            self.logger.error(f"Mapping评估异常|Mapping evaluation error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _discover_samples(self, reads_dir: str) -> List[tuple]:
        """发现样品FASTQ文件对|Discover sample FASTQ pairs"""
        samples = []
        processed_reads = set()

        # 查找read1文件|Find read1 files
        pattern = self.config.mapping_pattern.replace("_1.", "*_1.")
        read1_files = glob.glob(os.path.join(reads_dir, pattern))

        for read1 in sorted(read1_files):
            # 推断read2文件|Infer read2 file
            read2 = read1.replace("_1.", "_2.")
            read2 = read2.replace("_1.clean", "_2.clean")

            if not os.path.exists(read2):
                self.logger.warning(f"未找到对应的read2文件|Read2 not found for: {read1}")
                continue

            # 提取样品名|Extract sample name
            sample_name = os.path.basename(read1).split("_1")[0]

            # 避免重复|Avoid duplicates
            if sample_name in processed_reads:
                continue
            processed_reads.add(sample_name)

            samples.append((sample_name, read1, read2))

        return samples

    def _check_and_build_index(self) -> bool:
        """检查并构建基因组索引|Check and build genome index"""
        self.logger.info("检查基因组索引|Checking genome index")

        # BWA索引文件后缀|BWA index file extensions
        index_extensions = ['.amb', '.ann', '.bwt', '.pac', '.sa']
        genome_path = Path(self.config.genome)

        # 检查所有索引文件是否存在|Check if all index files exist
        missing_indices = []
        for ext in index_extensions:
            index_file = genome_path.parent / f"{genome_path.name}{ext}"
            if not index_file.exists():
                missing_indices.append(ext)

        if missing_indices:
            self.logger.info(f"缺少索引文件|Missing index files: {', '.join(missing_indices)}")
            self.logger.info("开始构建BWA索引|Building BWA index")
            return self._build_index()
        else:
            self.logger.info("基因组索引已存在|Genome index already exists")
            return True

    def _build_index(self) -> bool:
        """构建BWA索引|Build BWA index"""
        args = ["index", self.config.genome]
        cmd = self.bwa_cmd + args

        try:
            self.logger.info(f"运行BWA索引|Running BWA index: {' '.join(cmd)}")
            result = subprocess.run(
                cmd,
                capture_output=True,
                text=True,
                check=False
            )

            if result.returncode != 0:
                self.logger.error(f"BWA索引构建失败|BWA index building failed: {result.stderr}")
                return False

            self.logger.info("BWA索引构建成功|BWA index built successfully")
            return True

        except Exception as e:
            self.logger.error(f"BWA索引构建异常|Error building BWA index: {e}")
            return False

    def _process_sample(self, sample_name: str, read1: str, read2: str) -> Optional[Dict[str, Any]]:
        """处理单个样品|Process single sample"""
        self.logger.info(f"处理样品|Processing sample: {sample_name}")

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

        # Step 1: BWA mem比对（直接输出BAM）|BWA mem alignment (direct BAM output)
        if not self._run_bwa_mem(sample_name, read1, read2, bam_file):
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

    def _run_bwa_mem(self, sample_name: str, read1: str, read2: str, bam_file: Path) -> bool:
        """运行BWA mem，直接输出BAM|Run BWA mem with direct BAM output"""
        # 添加read group信息|Add read group info
        rg = f"@RG\\tID:{sample_name}\\tSM:{sample_name}\\tPL:ILLUMINA"

        # BWA mem命令|BWA mem command
        bwa_args = [
            "mem",
            "-t", str(self.config.mapping_threads),
            "-R", rg,
            self.config.genome,
            read1,
            read2
        ]
        bwa_cmd = self.bwa_cmd + bwa_args

        # samtools命令：从stdin读取SAM，输出BAM
        # samtools command: read SAM from stdin, output BAM
        samtools_args = [
            "view",
            "-@", str(self.config.mapping_threads),
            "-b",  # BAM格式|BAM format
            "-o", str(bam_file),
            "-"  # 从stdin读取|Read from stdin
        ]
        samtools_cmd = self.samtools_cmd + samtools_args

        try:
            self.logger.info(f"比对reads到基因组（管道输出BAM）|Aligning reads to genome (piped to BAM):")
            self.logger.info(f"  样品|Sample: {sample_name}")
            self.logger.info(f"  BWA命令|BWA command: {' '.join(bwa_cmd)}")
            self.logger.info(f"  SAMtools命令|SAMtools command: {' '.join(samtools_cmd)}")

            # 使用管道：bwa mem | samtools view -b
            # Use pipeline: bwa mem | samtools view -b
            process_bwa = subprocess.Popen(
                bwa_cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            process_samtools = subprocess.Popen(
                samtools_cmd,
                stdin=process_bwa.stdout,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True
            )

            # 关闭bwa的stdout，让samtools能接收到EOF
            # Close bwa's stdout to allow samtools to receive EOF
            process_bwa.stdout.close()

            # 等待两个进程完成|Wait for both processes to complete
            bwa_stdout, bwa_stderr = process_bwa.communicate()
            samtools_stdout, samtools_stderr = process_samtools.communicate()

            bwa_returncode = process_bwa.returncode
            samtools_returncode = process_samtools.returncode

            if bwa_returncode != 0:
                stderr_msg = bwa_stderr if bwa_stderr else ""
                # 检测内存错误|Detect memory errors
                if "memory" in stderr_msg.lower() or "MemoryError" in stderr_msg or "CondaMemoryError" in stderr_msg:
                    self.logger.error(f"BWA mem内存不足|BWA mem out of memory for {sample_name}")
                    self.logger.error("建议|Suggestion:")
                    self.logger.error("  1. 减少线程数|Reduce thread count: use -t parameter")
                    self.logger.error("  2. 增加系统内存|Increase system memory")
                    self.logger.error("  3. 使用分批比对|Use batch alignment for large datasets")
                    self.logger.error(f"当前使用线程数|Current threads: {self.config.mapping_threads}")
                else:
                    self.logger.error(f"BWA mem比对失败|BWA mem alignment failed for {sample_name}")
                    if stderr_msg:
                        self.logger.error(f"BWA STDERR: {stderr_msg}")
                return False

            if samtools_returncode != 0:
                stderr_msg = samtools_stderr if samtools_stderr else ""
                self.logger.error(f"SAMtools view失败|SAMtools view failed: {stderr_msg}")
                return False

            return True

        except Exception as e:
            self.logger.error(f"BWA mem运行异常|Unexpected error running BWA mem: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

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
            "-@", str(self.config.mapping_threads),
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
            "-@", str(self.config.mapping_threads),
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
        args = ["depth", "-@", str(self.config.mapping_threads), str(bam_file)]
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

            stats = {
                'total_reads': 0,
                'mapped_reads': 0,
                'mapping_rate': 0.0,
                'coverage_fraction': 0.0,
            }

            # 解析flagstat|Parse flagstat
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

