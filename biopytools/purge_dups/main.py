"""
Purge_Dups主程序模块|Purge_Dups Main Module
"""

import argparse
import sys
import os
from pathlib import Path
from .config import PurgeDupsConfig
from .utils import PurgeDupsLogger, CommandRunner, format_number


class PurgeDupsRunner:
    """Purge_Dups去冗余主类|Main Purge_Dups Deduplication Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = PurgeDupsConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = PurgeDupsLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 设置工具路径|Set tool paths
        self.purge_dups_bin = os.path.join(self.config.purge_dups_path, "bin", "purge_dups")
        self.split_fa_bin = os.path.join(self.config.purge_dups_path, "bin", "split_fa")
        self.get_seqs_bin = os.path.join(self.config.purge_dups_path, "bin", "get_seqs")
        self.calcuts_bin = os.path.join(self.config.purge_dups_path, "bin", "calcuts")
        self.pbcstat_bin = os.path.join(self.config.purge_dups_path, "bin", "pbcstat")
        self.ngscstat_bin = os.path.join(self.config.purge_dups_path, "bin", "ngscstat")
        self.minimap2_bin = os.path.join(self.config.purge_dups_path, "bin", "minimap2")

        # 获取基因组基础名称|Get genome base name
        # 去掉所有可能的扩展名（.fa, .fasta, .primary.fa等）|Remove all possible extensions
        input_path = Path(self.config.input)
        genome_base = input_path.stem
        # 如果还有.primary或其他常见后缀，继续去掉|If there are .primary or other common suffixes, continue removing
        if genome_base.endswith('.primary'):
            genome_base = genome_base[:-8]  # 去掉'.primary'|Remove '.primary'
        elif genome_base.endswith('.asm'):
            genome_base = genome_base[:-4]  # 去掉'.asm'|Remove '.asm'
        self.genome_base = genome_base

    def step1_calculate_coverage(self):
        """步骤1: 计算测序深度和统计信息|Step 1: Calculate coverage and statistics"""
        self.logger.info("开始步骤1: 计算测序深度|Starting step 1: Calculate coverage")

        # 选择reads类型|Select reads type
        if self.config.read_type in ['pacbio', 'hifi']:
            stat_cmd = self.pbcstat_bin
            minimap_opt = '-xmap-pb' if self.config.read_type == 'pacbio' else '-xmap-hifi'
        else:  # illumina
            stat_cmd = self.ngscstat_bin

        # 生成PAF文件|Generate PAF files
        paf_files = []
        for reads_file in self.config.reads_files:
            reads_base = Path(reads_file).stem
            paf_file = self.config.coverage_dir / f"{reads_base}.paf.gz"
            paf_file_abs = str(paf_file.absolute())  # 转换为绝对路径

            if self.config.read_type in ['pacbio', 'hifi']:
                cmd = f"{self.minimap2_bin} -t {self.config.threads} {minimap_opt} " \
                      f"-I 50G {self.config.input} {reads_file} | gzip -c > {paf_file_abs}"
            else:  # illumina需要先用bwa比对|illumina needs bwa alignment first
                self.logger.error("Illumina reads需要先用bwa比对生成BAM文件|"
                                  "Illumina reads require bwa alignment first")
                return False

            if not self.cmd_runner.run(cmd, f"比对reads文件|Aligning reads: {reads_base}"):
                return False
            paf_files.append(paf_file_abs)  # 使用绝对路径|Use absolute path

        # 计算测序深度统计|Calculate coverage statistics
        base_cov_file = self.config.coverage_dir / f"{self.genome_base}.base.cov"
        stat_file = self.config.coverage_dir / f"{self.genome_base}.stat"

        # 转换为绝对路径用于命令|Convert to absolute paths for command
        coverage_dir_abs = str(self.config.coverage_dir.absolute())
        paf_list = " ".join(paf_files)
        cmd = f"{stat_cmd} -O {coverage_dir_abs} {paf_list}"
        if not self.cmd_runner.run(cmd, "计算测序深度统计|Calculate coverage statistics"):
            return False

        # 检查输出文件|Check output files
        expected_base_cov = self.config.coverage_dir / "PB.base.cov"
        expected_stat = self.config.coverage_dir / "PB.stat"

        if expected_base_cov.exists():
            os.rename(expected_base_cov, base_cov_file)
        if expected_stat.exists():
            os.rename(expected_stat, stat_file)

        self.logger.info(f"步骤1完成|Step 1 completed: {base_cov_file}, {stat_file}")
        return True

    def step2_calculate_cutoffs(self):
        """步骤2: 计算覆盖度阈值|Step 2: Calculate coverage cutoffs"""
        self.logger.info("开始步骤2: 计算覆盖度阈值|Starting step 2: Calculate cutoffs")

        stat_file = self.config.coverage_dir / f"{self.genome_base}.stat"
        cutoffs_file = self.config.coverage_dir / "cutoffs"

        # 转换为绝对路径|Convert to absolute paths
        stat_file_abs = str(stat_file.absolute())
        cutoffs_abs = str(cutoffs_file.absolute())

        if not stat_file.exists():
            self.logger.error(f"找不到统计文件|Statistics file not found: {stat_file_abs}")
            return False

        # 检查是否手动指定阈值|Check if manual cutoffs specified
        if self.config.manual_cutoffs:
            if os.path.exists(self.config.manual_cutoffs):
                import shutil
                shutil.copy(self.config.manual_cutoffs, cutoffs_abs)
                self.logger.info(f"使用手动指定的阈值文件|Using manual cutoffs file: {self.config.manual_cutoffs}")
            else:
                self.logger.error(f"手动阈值文件不存在|Manual cutoffs file not found: {self.config.manual_cutoffs}")
                return False
        else:
            cmd = f"{self.calcuts_bin} -f {self.config.min_depth_fraction} {stat_file_abs} > {cutoffs_abs}"
            if not self.cmd_runner.run(cmd, "计算覆盖度阈值|Calculate coverage cutoffs"):
                return False

        self.logger.info(f"步骤2完成|Step 2 completed: {cutoffs_abs}")
        return True

    def step3_split_and_align(self):
        """步骤3: 分割基因组并做自比对|Step 3: Split genome and self-alignment"""
        self.logger.info("开始步骤3: 分割基因组并自比对|Starting step 3: Split genome and self-alignment")

        split_fa = self.config.split_aln_dir / f"{self.genome_base}.split"
        split_fa_abs = str(split_fa.absolute())  # 转换为绝对路径
        split_opt = "-n" if self.config.split_by_n else ""

        # 分割基因组|Split genome
        cmd = f"{self.split_fa_bin} {split_opt} {self.config.input} > {split_fa_abs}"
        if not self.cmd_runner.run(cmd, "分割基因组|Split genome"):
            return False

        # 自比对|Self-alignment
        paf_file = self.config.split_aln_dir / f"{self.genome_base}.split.self.paf.gz"
        paf_file_abs = str(paf_file.absolute())  # 转换为绝对路径
        cmd = f"{self.minimap2_bin} -xasm5 -DP -t {self.config.threads} " \
              f"{split_fa_abs} {split_fa_abs} | gzip -c > {paf_file_abs}"
        if not self.cmd_runner.run(cmd, "基因组自比对|Genome self-alignment"):
            return False

        self.logger.info(f"步骤3完成|Step 3 completed: {split_fa_abs}, {paf_file_abs}")
        return True

    def step4_purge_dups(self):
        """步骤4: 识别并去除冗余序列|Step 4: Identify and purge duplications"""
        self.logger.info("开始步骤4: 识别并去除冗余序列|Starting step 4: Identify and purge duplications")

        base_cov_file = self.config.coverage_dir / f"{self.genome_base}.base.cov"
        cutoffs_file = self.config.coverage_dir / "cutoffs"
        paf_file = self.config.split_aln_dir / f"{self.genome_base}.split.self.paf.gz"
        dups_bed = self.config.purge_dups_dir / "dups.bed"

        # 检查输入文件|Check input files
        for f in [base_cov_file, cutoffs_file, paf_file]:
            if not f.exists():
                self.logger.error(f"找不到输入文件|Input file not found: {f}")
                return False

        # 转换为绝对路径|Convert to absolute paths
        base_cov_abs = str(base_cov_file.absolute())
        cutoffs_abs = str(cutoffs_file.absolute())
        paf_abs = str(paf_file.absolute())
        dups_bed_abs = str(dups_bed.absolute())
        log_file = str((self.config.purge_dups_dir / 'purge_dups.log').absolute())

        # 构建purge_dups命令|Build purge_dups command
        cmd = f"{self.purge_dups_bin} "
        cmd += f"-c {base_cov_abs} "
        cmd += f"-T {cutoffs_abs} "
        cmd += f"-f {self.config.min_fraction} "
        cmd += f"-a {self.config.min_alignment_score} "
        cmd += f"-b {self.config.min_match_score} "
        cmd += f"-m {self.config.min_match_bases} "
        cmd += f"-M {self.config.max_gap_1st} "
        cmd += f"-G {self.config.max_gap_2nd} "
        cmd += f"-l {self.config.min_chain_score} "
        cmd += f"-E {self.config.max_contig_end_ext} "

        if self.config.two_round_chaining:
            cmd += "-2 "

        cmd += f"{paf_abs} > {dups_bed_abs} 2> {log_file}"

        if not self.cmd_runner.run(cmd, "识别冗余序列|Identify duplications"):
            return False

        self.logger.info(f"步骤4完成|Step 4 completed: {dups_bed_abs}")
        return True

    def step5_get_sequences(self):
        """步骤5: 获取纯化后的序列|Step 5: Get purged sequences"""
        self.logger.info("开始步骤5: 获取纯化后的序列|Starting step 5: Get purged sequences")

        dups_bed = self.config.purge_dups_dir / "dups.bed"
        output_prefix = self.config.seqs_dir / f"{self.genome_base}_purged"

        # 转换为绝对路径|Convert to absolute paths
        dups_bed_abs = str(dups_bed.absolute())
        output_prefix_abs = str(output_prefix.absolute())

        # 构建get_seqs命令|Build get_seqs command
        cmd = f"{self.get_seqs_bin} "
        cmd += f"-g {self.config.max_dup_gap} "
        cmd += f"-l {self.config.min_primary_length} "
        cmd += f"-m {self.config.min_length_ratio} "

        if self.config.ends_only:
            cmd += "-e "
        if self.config.split_contigs:
            cmd += "-s "
        if self.config.keep_high_coverage:
            cmd += "-c "
        if not self.config.add_prefix:
            cmd += "-a "

        cmd += f"-p {output_prefix_abs} "
        cmd += f"{dups_bed_abs} {self.config.input}"

        if not self.cmd_runner.run(cmd, "获取纯化后的序列|Get purged sequences"):
            return False

        # 检查输出文件|Check output files
        purged_fa = f"{output_prefix_abs}.purge.fa"
        haplotig_fa = f"{output_prefix_abs}.haplotigs.fa"

        if not os.path.exists(purged_fa):
            self.logger.warning(f"未找到纯化后的主序列文件|Purged primary file not found: {purged_fa}")
        else:
            self.logger.info(f"纯化后的主序列|Purged primary contigs: {purged_fa}")

        if not os.path.exists(haplotig_fa):
            self.logger.warning(f"未找到单倍型序列文件|Haplotigs file not found: {haplotig_fa}")
        else:
            self.logger.info(f"单倍型序列|Haplotigs: {haplotig_fa}")

        self.logger.info(f"步骤5完成|Step 5 completed")
        return True

    def run_single_step(self, step_num: int):
        """运行单个步骤|Run single step"""
        step_functions = {
            1: (self.step1_calculate_coverage, "计算测序深度|Calculate coverage"),
            2: (self.step2_calculate_cutoffs, "计算覆盖度阈值|Calculate cutoffs"),
            3: (self.step3_split_and_align, "分割基因组并自比对|Split genome and self-alignment"),
            4: (self.step4_purge_dups, "识别并去除冗余序列|Identify and purge duplications"),
            5: (self.step5_get_sequences, "获取纯化后的序列|Get purged sequences")
        }

        if step_num not in step_functions:
            self.logger.error(f"无效的步骤编号|Invalid step number: {step_num}")
            return False

        step_func, step_name = step_functions[step_num]
        self.logger.info(f"执行步骤{step_num}|Executing step {step_num}: {step_name}")

        success = step_func()
        if success:
            self.logger.info(f"步骤{step_num}完成|Step {step_num} completed: {step_name}")
        else:
            self.logger.error(f"步骤{step_num}失败|Step {step_num} failed: {step_name}")

        return success

    def run_full_pipeline(self):
        """运行完整的Purge_Dups流程|Run complete Purge_Dups pipeline"""
        self.logger.info("开始Purge_Dups去冗余流程|Starting Purge_Dups deduplication pipeline")

        steps = [
            (self.step1_calculate_coverage, "计算测序深度|Calculate coverage"),
            (self.step2_calculate_cutoffs, "计算覆盖度阈值|Calculate cutoffs"),
            (self.step3_split_and_align, "分割基因组并自比对|Split genome and self-alignment"),
            (self.step4_purge_dups, "识别并去除冗余序列|Identify and purge duplications"),
            (self.step5_get_sequences, "获取纯化后的序列|Get purged sequences")
        ]

        for i, (step_func, step_name) in enumerate(steps, 1):
            self.logger.info(f"执行步骤{i}|Executing step {i}: {step_name}")

            if not step_func():
                self.logger.error(f"步骤{i}失败|Step {i} failed: {step_name}")
                return False

            self.logger.info(f"步骤{i}完成|Step {i} completed: {step_name}")

        self.logger.info("Purge_Dups去冗余流程全部完成|Purge_Dups deduplication pipeline completed!")
        return True

    def run_analysis(self):
        """运行分析|Run analysis"""
        try:
            if self.config.step:
                success = self.run_single_step(self.config.step)
                if not success:
                    sys.exit(1)
            else:
                success = self.run_full_pipeline()
                if not success:
                    sys.exit(1)

        except Exception as e:
            self.logger.error(f"程序执行出错|Program execution error: {str(e)}")
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='Purge_Dups基因组去冗余工具(基于测序深度去除单倍型和重叠序列)|'
                    'Purge_Dups Genome Deduplication Tool (Remove haplotigs and overlaps based on read depth)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='基因组组装文件(FASTA格式)|Genome assembly file in FASTA format')
    parser.add_argument('-r', '--reads', required=True,
                       help='测序文件(PacBio/HiFi/illumina)|Sequencing reads file')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output-dir',
                       default='./purge_dups_output',
                       help='输出目录|Output directory')
    parser.add_argument('-t', '--threads',
                       type=int, default=12,
                       help='线程数|Number of threads')

    # reads类型|Reads type
    parser.add_argument('--read-type',
                       choices=['pacbio', 'hifi', 'illumina'],
                       default='hifi',
                       help='测序数据类型|Sequencing data type')

    # purge_dups参数|purge_dups parameters
    parser.add_argument('--min-fraction',
                       type=float, default=0.8,
                       help='最小比例阈值|Minimum fraction threshold')
    parser.add_argument('--two-round-chaining',
                       action='store_true', default=True,
                       help='启用两轮链式匹配|Enable two-round chaining')
    parser.add_argument('--no-two-round-chaining',
                       action='store_false', dest='two_round_chaining',
                       help='禁用两轮链式匹配|Disable two-round chaining')

    # get_seqs参数|get_seqs parameters
    parser.add_argument('--ends-only',
                       action='store_true', default=True,
                       help='只去除contig末端的冗余|Only remove duplications at contig ends')
    parser.add_argument('--no-ends-only',
                       action='store_false', dest='ends_only',
                       help='也去除contig中间的冗余|Also remove duplications in contig middle')
    parser.add_argument('--min-primary-length',
                       type=int, default=10000,
                       help='最小主contig长度|Minimum primary contig length')

    # 步骤控制|Step control
    parser.add_argument('-s', '--step', type=int, choices=[1, 2, 3, 4, 5],
                       help='只运行指定步骤|Run only specified step '
                            '(1: 计算深度|coverage, 2: 计算阈值|cutoffs, '
                            '3: 分割比对|split&align, 4: 去冗余|purge, 5: 获取序列|get seqs)')

    # 其他参数|Other parameters
    parser.add_argument('--split-by-n',
                       action='store_true',
                       help='split_fa按N分割|split_fa split by N')
    parser.add_argument('--manual-cutoffs',
                       help='手动指定阈值文件|Manual cutoffs file')

    args = parser.parse_args()

    # 创建运行器并运行|Create runner and run
    runner = PurgeDupsRunner(
        input=args.input,
        reads=args.reads,
        output_dir=args.output_dir,
        threads=args.threads,
        read_type=args.read_type,
        min_fraction=args.min_fraction,
        two_round_chaining=args.two_round_chaining,
        ends_only=args.ends_only,
        min_primary_length=args.min_primary_length,
        step=args.step,
        split_by_n=args.split_by_n,
        manual_cutoffs=args.manual_cutoffs
    )

    runner.run_analysis()


if __name__ == "__main__":
    main()
