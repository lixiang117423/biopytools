"""
FASTQ配对检查与排序混入类|FASTQ Pair Checking and Sorting Mixin Class

提供配对检查、排序和过滤功能
Provides pair checking, sorting and filtering functionality
"""

import os
import subprocess
import shutil
from typing import Dict, Tuple, List
from pathlib import Path


class PairFixMixin:
    """
    FASTQ配对检查与排序混入类|FASTQ Pair Checking and Sorting Mixin Class

    提供配对检查、排序和过滤功能
    Provides pair checking, sorting and filtering functionality
    """

    def __init__(self):
        """初始化混入类|Initialize mixin"""
        # 这些属性将由使用此mixin的类提供
        # These attributes will be provided by classes using this mixin
        self.config = None
        self.logger = None
        self.cmd_runner = None

    def setup_pair_fix(self, config, logger, cmd_runner):
        """
        设置pair-fixing所需的依赖|Setup dependencies for pair-fixing

        Args:
            config: 配置对象|Configuration object
            logger: 日志对象|Logger object
            cmd_runner: 命令执行器|Command runner object
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def find_pairs(self) -> Dict[str, Tuple[str, str]]:
        """
        查找配对文件|Find paired files

        Returns:
            字典: {sample_name: (r1_path, r2_path)}
            Dict: {sample_name: (r1_path, r2_path)}
        """
        input_dir = self.config.input_dir
        suffix1 = self.config.read1_suffix
        suffix2 = self.config.read2_suffix

        self.logger.info(f"扫描输入目录: {input_dir}|Scanning input directory")

        r1_files = {}
        r2_files = {}

        # 查找R1文件|Find R1 files
        for file in os.listdir(input_dir):
            if file.endswith(suffix1):
                sample_name = file.replace(suffix1, '')
                r1_path = os.path.join(input_dir, file)
                r1_files[sample_name] = r1_path

        # 查找R2文件|Find R2 files
        for file in os.listdir(input_dir):
            if file.endswith(suffix2):
                sample_name = file.replace(suffix2, '')
                r2_path = os.path.join(input_dir, file)
                r2_files[sample_name] = r2_path

        # 找出完整配对|Find complete pairs
        pairs = {}
        r1_only = set(r1_files.keys()) - set(r2_files.keys())
        r2_only = set(r2_files.keys()) - set(r1_files.keys())

        for sample_name in r1_files.keys():
            if sample_name in r2_files:
                pairs[sample_name] = (r1_files[sample_name], r2_files[sample_name])

        # 报告统计信息|Report statistics
        self.logger.info(f"找到 {len(pairs)} 对完整文件|Found {len(pairs)} complete pairs")

        if r1_only:
            self.logger.warning(f"发现 {len(r1_only)} 个只有R1的样本|Found {len(r1_only)} R1-only samples")
            if len(r1_only) <= 5:
                self.logger.warning(f"  样本列表|Sample list: {list(r1_only)}")
            else:
                self.logger.warning(f"  部分样本|Partial samples: {list(r1_only)[:5]}...")

        if r2_only:
            self.logger.warning(f"发现 {len(r2_only)} 个只有R2的样本|Found {len(r2_only)} R2-only samples")
            if len(r2_only) <= 5:
                self.logger.warning(f"  样本列表|Sample list: {list(r2_only)}")
            else:
                self.logger.warning(f"  部分样本|Partial samples: {list(r2_only)[:5]}...")

        return pairs

    def run_pair_fixing(self) -> bool:
        """
        运行配对检查和排序流程|Run pair checking and sorting pipeline

        步骤|Steps:
        1. 配对检查（抽样前1000条reads）|Pair checking (sample first 1000 reads)
        2. 如果配对不对应，进行排序和过滤|If not paired, sort and filter

        Returns:
            是否成功|Success status
        """
        self.logger.info("=" * 60)
        self.logger.info("步骤1: 配对检查与排序|Step 1: Pair Check and Sort")
        self.logger.info("=" * 60)

        # 显示配置信息|Show configuration
        self.logger.info("配置信息|Configuration:")
        self.logger.info(f"  输入目录|Input directory: {self.config.input_dir}")
        sorted_dir = self.config.output_path / "fixed" / "sorted"
        sorted_dir.mkdir(parents=True, exist_ok=True)
        self.logger.info(f"  输出目录|Output directory: {sorted_dir}")
        self.logger.info(f"  工具|Tool: seqkit")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info("-" * 60)

        # 查找配对文件|Find paired files
        pairs = self.find_pairs()

        if not pairs:
            self.logger.error("未找到任何配对文件|No paired files found")
            return False

        # 对每个样本进行检查和处理|Check and process each sample
        total_samples = len(pairs)
        self.logger.info(f"开始处理 {total_samples} 个样本|Start processing {total_samples} samples")
        self.logger.info("-" * 60)

        success_count = 0
        failed_samples = []
        already_ok_count = 0
        need_sort_count = 0

        for idx, (sample_name, (r1_file, r2_file)) in enumerate(pairs.items(), 1):
            self.logger.info(f"[{idx}/{total_samples}] 检查样本 {sample_name}|Checking sample {sample_name}")

            # 先检查配对是否已经对应|Check if pairs are already matched
            is_ok, check_result = self.check_paired_files(r1_file, r2_file)

            if is_ok:
                # 配对已经对应，直接复制到sorted目录|Pairs already matched, copy directly
                self.logger.info(f"  配对已对应，跳过排序|Pairs already matched, skipping sort")
                shutil.copy(r1_file, sorted_dir / Path(r1_file).name)
                shutil.copy(r2_file, sorted_dir / Path(r2_file).name)
                already_ok_count += 1
                success_count += 1
            else:
                # 需要排序|Need to sort
                self.logger.info(f"  配对不对应，需要排序|Pairs not matched, sorting required")
                self.logger.info(f"  检查结果|Check result: {check_result}")
                need_sort_count += 1

                success = self.process_sample_with_seqkit_sort(sample_name, r1_file, r2_file, sorted_dir)

                if success:
                    success_count += 1
                else:
                    failed_samples.append(sample_name)

        # 处理结果汇总|Processing summary
        self.logger.info("=" * 60)
        self.logger.info("配对检查与排序结果汇总|Pair Check and Sort Summary:")
        self.logger.info(f"  总样本数|Total samples: {total_samples}")
        self.logger.info(f"  配对已对应无需排序|Already matched: {already_ok_count}/{total_samples}")
        self.logger.info(f"  需要排序|Required sorting: {need_sort_count}/{total_samples}")
        self.logger.info(f"  成功|Success: {success_count}/{total_samples}")
        self.logger.info(f"  失败|Failed: {len(failed_samples)}/{total_samples}")

        if failed_samples:
            self.logger.warning(f"失败样本列表|Failed samples: {', '.join(failed_samples)}")

        self.logger.info("=" * 60)

        # 更新输入目录为sorted目录，供fastp使用|Update input_dir to sorted directory for fastp
        self.config.input_dir = str(sorted_dir)
        self.logger.info(f"配对文件已准备好，供fastp使用|Paired files ready for fastp: {self.config.input_dir}")
        self.logger.info("")

        return success_count == total_samples

    def check_paired_files(self, r1_file: str, r2_file: str) -> tuple:
        """
        检查R1和R2文件的配对是否已经对应（检查所有reads）
        Check if R1 and R2 files are already properly paired (check all reads)

        Args:
            r1_file: R1文件路径|R1 file path
            r2_file: R2文件路径|R2 file path

        Returns:
            (是否配对对应, 检查结果信息)|(Is paired matched, Check result message)
        """
        try:
            self.logger.debug(f"  检查所有reads的配对对应关系...|Checking all reads pairing...")

            # 使用seqkit提取所有read name|Extract all read names using seqkit
            temp_dir = Path(self.config.output_path) / "temp"
            temp_dir.mkdir(parents=True, exist_ok=True)

            names_r1 = temp_dir / "r1_names.txt"
            names_r2 = temp_dir / "r2_names.txt"

            # 提取R1的所有read name|Extract all R1 read names
            if not self.cmd_runner.run(["seqkit", "seq", "-n", "-i", r1_file, "-o", str(names_r1)], "提取R1 read names|Extract R1 read names"):
                return False, "提取R1 read names失败|Failed to extract R1 read names"

            # 提取R2的所有read name|Extract all R2 read names
            if not self.cmd_runner.run(["seqkit", "seq", "-n", "-i", r2_file, "-o", str(names_r2)], "提取R2 read names|Extract R2 read names"):
                return False, "提取R2 read names失败|Failed to extract R2 read names"

            # 读取所有read name|Read all read names
            with open(names_r1, 'r') as f:
                r1_names = [line.strip() for line in f if line.strip()]
            with open(names_r2, 'r') as f:
                r2_names = [line.strip() for line in f if line.strip()]

            # 清理临时文件|Clean temporary files
            names_r1.unlink()
            names_r2.unlink()

            # 检查数量|Check counts
            if len(r1_names) != len(r2_names):
                return False, f"reads数量不一致|Reads count mismatch: R1={len(r1_names)}, R2={len(r2_names)}"

            # 检查所有对应位置的read name是否相同|Check if all read names match at same positions
            mismatch_count = 0
            total_count = len(r1_names)

            for i in range(total_count):
                # read name格式通常为: @READNAME/1 或 @READNAME/2
                # 去掉 /1 或 /2 后缀进行比较
                name1 = r1_names[i].rsplit('/', 1)[0] if '/' in r1_names[i] else r1_names[i]
                name2 = r2_names[i].rsplit('/', 1)[0] if '/' in r2_names[i] else r2_names[i]

                if name1 != name2:
                    mismatch_count += 1
                    # 只显示前5个不匹配的例子
                    if mismatch_count <= 5:
                        self.logger.debug(f"    位置{i+1}不匹配|Position {i+1} mismatch:")
                        self.logger.debug(f"      R1: {r1_names[i]}")
                        self.logger.debug(f"      R2: {r2_names[i]}")

                    # 如果已经发现10个不匹配，可以提前结束
                    if mismatch_count >= 10:
                        self.logger.debug(f"    ... (已发现{mismatch_count}个不匹配，停止检查)|... (found {mismatch_count} mismatches, stopping)")
                        break

            # 判断结果|Determine result
            if mismatch_count == 0:
                return True, f"配对完全对应|Fully matched ({total_count} reads)"
            else:
                mismatch_rate = (mismatch_count / total_count) * 100
                return False, f"发现{mismatch_count}个顺序不匹配|Found {mismatch_count} order mismatches out of {total_count} ({mismatch_rate:.2f}%)"

        except Exception as e:
            self.logger.warning(f"  配对检查失败|Pair check failed: {e}")
            return False, f"检查失败|Check failed: {e}"

    def process_sample_with_seqkit_sort(self, sample_name: str, r1_file: str, r2_file: str, output_dir: Path) -> bool:
        """
        使用seqkit处理单个样本（按read name排序+配对过滤）
        Process single sample with seqkit (sort by read name + pair filtering)

        Args:
            sample_name: 样本名称|Sample name
            r1_file: R1文件路径|R1 file path
            r2_file: R2文件路径|R2 file path
            output_dir: 输出目录|Output directory

        Returns:
            是否成功|Success status
        """
        # 创建临时目录|Create temporary directory
        temp_dir = output_dir / "temp"
        temp_dir.mkdir(parents=True, exist_ok=True)

        # 定义输出文件|Define output files
        output_r1 = output_dir / f"{sample_name}{self.config.read1_suffix}"
        output_r2 = output_dir / f"{sample_name}{self.config.read2_suffix}"

        # 定义临时文件前缀|Define temporary file prefix
        temp_prefix = temp_dir / sample_name

        try:
            # 步骤1: 按read name排序|Step 1: Sort by read name
            self.logger.info(f"  [1/3] 按read name排序...|Sorting by read name...")

            sorted_r1 = f"{temp_prefix}_sorted_1.fq.gz"
            sorted_r2 = f"{temp_prefix}_sorted_2.fq.gz"

            # R1排序|Sort R1
            cmd_sort_r1 = f"zcat {r1_file} | paste - - - - | sort -k1,1 -S 10G --parallel={self.config.threads} | tr '\\t' '\\n' | gzip > {sorted_r1}"
            if not self.cmd_runner.run_shell(cmd_sort_r1, "排序R1|Sort R1"):
                return False

            # R2排序|Sort R2
            cmd_sort_r2 = f"zcat {r2_file} | paste - - - - | sort -k1,1 -S 10G --parallel={self.config.threads} | tr '\\t' '\\n' | gzip > {sorted_r2}"
            if not self.cmd_runner.run_shell(cmd_sort_r2, "排序R2|Sort R2"):
                return False

            # 步骤2: 统计排序后的reads数|Step 2: Count reads after sorting
            r1_count = self._count_reads(sorted_r1)
            r2_count = self._count_reads(sorted_r2)
            self.logger.info(f"  排序后 R1|After sorting R1: {r1_count} reads")
            self.logger.info(f"  排序后 R2|After sorting R2: {r2_count} reads")

            # 步骤3: 提取共有read name|Step 3: Extract common read names
            self.logger.info(f"  [2/3] 提取R1/R2共有的read name...|Extracting common read names...")

            names_r1 = f"{temp_prefix}_names_R1.txt"
            names_r2 = f"{temp_prefix}_names_R2.txt"

            # 使用seqkit提取read name|Extract read names using seqkit
            if not self.cmd_runner.run(["seqkit", "seq", "-n", "-i", sorted_r1, "-o", names_r1], "提取R1 read names|Extract R1 read names"):
                return False
            if not self.cmd_runner.run(["seqkit", "seq", "-n", "-i", sorted_r2, "-o", names_r2], "提取R2 read names|Extract R2 read names"):
                return False

            # 排序name文件|Sort name files
            subprocess.run(f"sort {names_r1} -o {names_r1}", shell=True, check=True)
            subprocess.run(f"sort {names_r2} -o {names_r2}", shell=True, check=True)

            # 计算共有reads|Count common reads
            common_names = f"{temp_prefix}_common_names.txt"
            subprocess.run(f"comm -12 {names_r1} {names_r2} > {common_names}", shell=True, check=True)

            with open(common_names, 'r') as f:
                common_count = len(f.readlines())

            only_r1 = subprocess.run(f"comm -23 {names_r1} {names_r2} | wc -l", shell=True, capture_output=True, text=True).stdout.strip()
            only_r2 = subprocess.run(f"comm -13 {names_r1} {names_r2} | wc -l", shell=True, capture_output=True, text=True).stdout.strip()

            self.logger.info(f"    共有reads|Common reads: {common_count}")
            self.logger.info(f"    仅R1有|R1 only: {only_r1}")
            self.logger.info(f"    仅R2有|R2 only: {only_r2}")

            # 计算共有比例|Calculate common ratio
            ratio = (common_count / (common_count + int(only_r1))) * 100 if (common_count + int(only_r1)) > 0 else 0
            self.logger.info(f"    共有比例|Common ratio: {ratio:.1f}%")

            # 共有比例低于50%时警告|Warn if common ratio < 50%
            if ratio < 50:
                self.logger.warning("  警告: 共有read比例低于50%，数据可能存在根本性错误！|Warning: Common read ratio < 50%, data may have fundamental errors!")
                self.logger.warning("  建议联系测序公司确认原始数据|Recommend contacting sequencing company to verify data")
                # 自动跳过此样本|Skip this sample automatically
                self.logger.warning(f"  自动跳过样本 {sample_name}|Automatically skipping sample {sample_name}")
                return False

            # 步骤4: 过滤只保留配对的reads|Step 4: Filter to keep only paired reads
            self.logger.info(f"  [3/3] 过滤保留配对reads...|Filtering to keep paired reads...")

            # 使用seqkit grep过滤|Filter using seqkit grep
            if not self.cmd_runner.run(["seqkit", "grep", "-f", common_names, sorted_r1, "-o", str(output_r1)], "过滤R1|Filter R1"):
                return False
            if not self.cmd_runner.run(["seqkit", "grep", "-f", common_names, sorted_r2, "-o", str(output_r2)], "过滤R2|Filter R2"):
                return False

            # 步骤5: 最终验证|Step 5: Final validation
            final_r1_count = self._count_reads(str(output_r1))
            final_r2_count = self._count_reads(str(output_r2))

            self.logger.info(f"  最终结果|Final results:")
            self.logger.info(f"    R1: {final_r1_count} reads -> {output_r1.name}")
            self.logger.info(f"    R2: {final_r2_count} reads -> {output_r2.name}")

            if final_r1_count == final_r2_count:
                self.logger.info(f"  配对验证通过|Pair validation passed")
                success = True
            else:
                self.logger.error(f"  配对验证失败，R1/R2数量仍不一致|Pair validation failed, R1/R2 counts still mismatch")
                success = False

        except Exception as e:
            self.logger.error(f"  处理样本{sample_name}时出错|Error processing sample {sample_name}: {e}")
            success = False

        finally:
            # 清理临时文件|Clean up temporary files
            self.logger.debug(f"  清理临时文件...|Cleaning temporary files...")
            if temp_dir.exists():
                shutil.rmtree(temp_dir)

        return success

    def _count_reads(self, fastq_file: str) -> int:
        """
        统计FASTQ文件中的reads数|Count reads in FASTQ file

        Args:
            fastq_file: FASTQ文件路径|FASTQ file path

        Returns:
            reads数量|Number of reads
        """
        try:
            result = subprocess.run(
                f"zcat {fastq_file} | wc -l",
                shell=True,
                capture_output=True,
                text=True,
                check=True
            )
            line_count = int(result.stdout.strip())
            return line_count // 4
        except Exception as e:
            self.logger.warning(f"  统计reads失败|Failed to count reads: {e}")
            return 0
