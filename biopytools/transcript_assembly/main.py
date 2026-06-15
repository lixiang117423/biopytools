"""
转录本从头组装主程序模块|Transcript De Novo Assembly Main Module
"""

import argparse
import glob
import sys
import os
import time
from .config import TranscriptAssemblyConfig
from .utils import TranscriptAssemblyLogger, CommandRunner
from .assembly import (
    HISAT2Indexer, HISAT2Aligner, SamtoolsProcessor,
    StringTieAssembler, StringTieMerger, GFFReadExtractor, SampleParser
)


class TranscriptAssembler:
    """转录本从头组装主类|Main Transcript De Novo Assembly Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = TranscriptAssemblyConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = TranscriptAssemblyLogger(
            self.config.output_path,
            verbose=self.config.verbose,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各个处理器|Initialize processors
        self.sample_parser = SampleParser(self.logger)
        self.hisat2_indexer = HISAT2Indexer(self.config, self.logger, self.cmd_runner)
        self.hisat2_aligner = HISAT2Aligner(self.config, self.logger, self.cmd_runner)
        self.samtools_processor = SamtoolsProcessor(self.config, self.logger, self.cmd_runner)
        self.stringtie_assembler = StringTieAssembler(self.config, self.logger, self.cmd_runner)
        self.stringtie_merger = StringTieMerger(self.config, self.logger, self.cmd_runner)
        self.gffread_extractor = GFFReadExtractor(self.config, self.logger, self.cmd_runner)

    def _parse_samples(self):
        """解析输入样本|Parse input samples"""
        samples = self.sample_parser.parse_input_samples(
            self.config.input_dir, self.config.fastq_pattern
        )

        if not samples:
            self.logger.error("未找到有效的样本文件|No valid sample files found")
            sys.exit(1)

        self.logger.info(f"找到 {len(samples)} 个样本|Found {len(samples)} samples:")
        for sample in samples:
            self.logger.info(f"  - {sample['name']}: {sample['fastq1']}, {sample['fastq2']}")

        self.config.samples = samples
        return samples

    def step1_build_index(self) -> str:
        """步骤1: 构建HISAT2基因组索引|Step 1: Build HISAT2 genome index"""
        self.logger_manager.step("步骤1: 构建HISAT2基因组索引|Step 1: Building HISAT2 genome index")
        return self.hisat2_indexer.build_index()

    def step2_align_samples(self, index_prefix: str):
        """步骤2: HISAT2比对|Step 2: HISAT2 alignment"""
        self.logger_manager.step("步骤2: HISAT2比对（--dta）|Step 2: HISAT2 alignment (--dta)")

        samples = self.config.samples
        output_dir = self.config.output_dir

        for sample_info in samples:
            sample_name = sample_info["name"]
            output_sam = os.path.join(output_dir, "02_hisat2_align", f"{sample_name}.sam")

            self.logger.info(f"比对样本|Aligning sample: {sample_name}")

            if not self.hisat2_aligner.align(
                index_prefix, sample_info["fastq1"], sample_info["fastq2"], output_sam
            ):
                self.logger.warning(f"样本比对失败，跳过|Sample alignment failed, skipping: {sample_name}")

    def step3_sort_bam(self):
        """步骤3: SAM转排序BAM|Step 3: Convert SAM to sorted BAM"""
        self.logger_manager.step("步骤3: SAM转排序BAM|Step 3: Convert SAM to sorted BAM")

        samples = self.config.samples
        output_dir = self.config.output_dir

        for sample_info in samples:
            sample_name = sample_info["name"]
            input_sam = os.path.join(output_dir, "02_hisat2_align", f"{sample_name}.sam")
            output_bam = os.path.join(output_dir, "03_bam_sort", f"{sample_name}.sorted.bam")

            self.logger.info(f"排序样本|Sorting sample: {sample_name}")

            if not self.samtools_processor.sort_and_index(input_sam, output_bam):
                self.logger.warning(f"样本排序失败，跳过|Sample sorting failed, skipping: {sample_name}")

    def step4_assemble_per_sample(self):
        """步骤4: StringTie逐样本从头组装|Step 4: StringTie per-sample de novo assembly"""
        self.logger_manager.step("步骤4: StringTie逐样本组装|Step 4: StringTie per-sample de novo assembly")

        samples = self.config.samples
        output_dir = self.config.output_dir

        for sample_info in samples:
            sample_name = sample_info["name"]
            bam_file = os.path.join(output_dir, "03_bam_sort", f"{sample_name}.sorted.bam")
            output_gtf = os.path.join(output_dir, "04_stringtie", f"{sample_name}.gtf")

            self.logger.info(f"组装样本|Assembling sample: {sample_name}")

            if not self.stringtie_assembler.assemble(bam_file, output_gtf):
                self.logger.warning(f"样本组装失败，跳过|Sample assembly failed, skipping: {sample_name}")

    def step5_merge_gtf(self):
        """步骤5: StringTie合并GTF|Step 5: StringTie merge GTFs"""
        self.logger_manager.step("步骤5: StringTie合并GTF|Step 5: StringTie merge GTFs")

        output_dir = self.config.output_dir
        gtf_dir = os.path.join(output_dir, "04_stringtie")

        # 收集所有GTF文件|Collect all GTF files
        gtf_files = sorted(glob.glob(os.path.join(gtf_dir, "*.gtf")))

        if not gtf_files:
            self.logger.error("未找到样本GTF文件，无法合并|No sample GTF files found, cannot merge")
            return False

        self.logger.info(f"找到 {len(gtf_files)} 个样本GTF文件|Found {len(gtf_files)} sample GTF files")

        # 创建GTF列表文件|Create GTF list file
        gtf_list_file = os.path.join(output_dir, "04_stringtie", "gtf_list.txt")
        with open(gtf_list_file, 'w') as f:
            for gtf_file in gtf_files:
                f.write(gtf_file + '\n')

        self.logger.info(f"GTF列表文件|GTF list file: {gtf_list_file}")

        # stringtie --merge
        output_merged_gtf = os.path.join(output_dir, "05_merge", "merged.gtf")

        return self.stringtie_merger.merge(
            gtf_list_file, self.config.genome_file, output_merged_gtf
        )

    def step6_extract_transcripts(self):
        """步骤6: gffread提取转录本序列|Step 6: gffread extract transcript sequences"""
        self.logger_manager.step("步骤6: 提取转录本序列|Step 6: Extract transcript sequences")

        output_dir = self.config.output_dir
        merged_gtf = os.path.join(output_dir, "05_merge", "merged.gtf")

        if not os.path.exists(merged_gtf):
            self.logger.error(f"合并GTF文件不存在|Merged GTF file not found: {merged_gtf}")
            return False

        output_fa = os.path.join(output_dir, "06_transcripts", "transcripts.fa")

        return self.gffread_extractor.extract_transcripts(
            merged_gtf, self.config.genome_file, output_fa
        )

    def run_analysis(self):
        """运行完整的转录本从头组装流程|Run complete transcript de novo assembly pipeline"""
        start_time = time.time()

        try:
            self.logger.info("=" * 60)
            self.logger.info("转录本从头组装流程开始|Transcript de novo assembly pipeline started")
            self.logger.info("=" * 60)
            self.logger.info(f"基因组文件|Genome file: {self.config.genome_file}")
            self.logger.info(f"输入目录|Input directory: {self.config.input_dir}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")
            if self.config.fastq_pattern:
                self.logger.info(f"文件模式|File pattern: {self.config.fastq_pattern}")

            # 解析样本|Parse samples
            self.logger_manager.step("解析输入样本|Parsing input samples")
            self._parse_samples()

            if self.config.step:
                # 运行指定步骤|Run specified step
                step_map = {
                    1: ("构建HISAT2索引|Build HISAT2 index", self._run_step1),
                    2: ("HISAT2比对|HISAT2 alignment", self._run_step2),
                    3: ("SAM转BAM|Convert SAM to BAM", self._run_step3),
                    4: ("StringTie组装|StringTie assembly", self._run_step4),
                    5: ("合并GTF|Merge GTFs", self._run_step5),
                    6: ("提取转录本|Extract transcripts", self._run_step6),
                }
                step_name, step_func = step_map[self.config.step]
                self.logger.info(f"执行步骤{self.config.step}|Executing step {self.config.step}: {step_name}")
                success = step_func()
                if not success and self.config.step != 6:
                    self.logger.error(f"步骤{self.config.step}执行失败|Step {self.config.step} failed")
                    sys.exit(1)
            else:
                # 运行完整流程|Run full pipeline
                index_prefix = self._run_step1()
                if index_prefix is None:
                    sys.exit(1)

                self._run_step2(index_prefix)
                self._run_step3()
                success = self._run_step4()
                if not success:
                    self.logger.error("组装步骤失败，无法继续|Assembly step failed, cannot continue")
                    sys.exit(1)

                success = self._run_step5()
                if not success:
                    self.logger.error("合并GTF步骤失败|Merge GTF step failed")
                    sys.exit(1)

                success = self._run_step6()
                if not success:
                    self.logger.error("提取转录本步骤失败|Extract transcript step failed")
                    sys.exit(1)

            # 输出总结信息|Output summary
            elapsed_time = time.time() - start_time
            self.logger.info("=" * 60)
            self.logger.info("流程总结|Pipeline Summary")
            self.logger.info("=" * 60)
            self.logger.info(f"总运行时间|Total runtime: {elapsed_time:.2f}秒|seconds")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

            # 列出关键输出文件|List key output files
            transcript_fa = os.path.join(self.config.output_dir, "06_transcripts", "transcripts.fa")
            merged_gtf = os.path.join(self.config.output_dir, "05_merge", "merged.gtf")
            if os.path.exists(transcript_fa):
                self.logger.info(f"转录本序列|Transcript sequences: {transcript_fa}")
            if os.path.exists(merged_gtf):
                self.logger.info(f"合并GTF|Merged GTF: {merged_gtf}")

            self.logger.info("流程成功完成|Pipeline completed successfully")
            self.logger.info("=" * 60)

        except KeyboardInterrupt:
            self.logger.warning("操作被用户中断|Operation interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"流程执行出错|Pipeline execution error: {e}", exc_info=True)
            sys.exit(1)

    def _run_step1(self):
        """运行步骤1|Run step 1"""
        try:
            return self.step1_build_index()
        except RuntimeError as e:
            self.logger.error(str(e))
            return None

    def _run_step2(self, index_prefix=None):
        """运行步骤2|Run step 2"""
        if index_prefix is None:
            index_prefix = os.path.join(
                self.config.output_dir, "01_hisat2_index",
                os.path.splitext(os.path.basename(self.config.genome_file))[0]
            )
        self.step2_align_samples(index_prefix)

    def _run_step3(self):
        """运行步骤3|Run step 3"""
        self.step3_sort_bam()

    def _run_step4(self):
        """运行步骤4|Run step 4"""
        self.step4_assemble_per_sample()
        return True

    def _run_step5(self):
        """运行步骤5|Run step 5"""
        return self.step5_merge_gtf()

    def _run_step6(self):
        """运行步骤6|Run step 6"""
        return self.step6_extract_transcripts()


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="转录本从头组装流程：HISAT2 + StringTie|Transcript de novo assembly pipeline: HISAT2 + StringTie",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -g genome.fasta -i ./clean_data -o ./transcript_output
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument("-g", "--genome", required=True,
                          help="基因组FASTA文件路径|Genome FASTA file path")
    required.add_argument("-i", "--input", required=True,
                          help="输入FASTQ文件目录|Input FASTQ file directory")
    required.add_argument("-o", "--output", required=True,
                          help="输出目录|Output directory")

    # 可选参数|Optional parameters
    optional = parser.add_argument_group('可选参数|Optional parameters')
    optional.add_argument("-p", "--pattern", default="*_1.clean.fq.gz",
                          help='FASTQ文件命名模式（*为样本名占位符）|FASTQ file naming pattern (* is sample name placeholder)')
    optional.add_argument("-t", "--threads", type=int, default=12,
                          help="线程数|Number of threads")
    optional.add_argument("--sample-timeout", type=int, default=43200,
                          help="单个样本处理超时时间（秒）|Sample processing timeout in seconds")

    # 步骤控制|Step control
    step_group = parser.add_argument_group('步骤控制|Step control')
    step_group.add_argument("-s", "--step", type=int, choices=[1, 2, 3, 4, 5, 6],
                            help='运行指定步骤|Run only specified step '
                                 '(1: 索引|index, 2: 比对|alignment, 3: 排序|sort, '
                                 '4: 组装|assembly, 5: 合并|merge, 6: 提取|extraction)')

    # 日志选项|Logging options
    logging_group = parser.add_argument_group('日志选项|Logging options')
    logging_group.add_argument('-v', '--verbose', action='count', default=0,
                               help='增加输出详细程度|Increase output verbosity')
    logging_group.add_argument('--quiet', action='store_true',
                               help='静默模式，仅输出错误信息|Quiet mode, only output errors')
    logging_group.add_argument('--log-level',
                               choices=['DEBUG', 'INFO', 'WARNING', 'ERROR', 'CRITICAL'],
                               default='INFO', help='日志级别|Log level')

    # 高级选项|Advanced options
    advanced = parser.add_argument_group('高级选项|Advanced options')
    advanced.add_argument('--force', action='store_true',
                          help='强制重新处理已完成的步骤|Force re-process completed steps')

    args = parser.parse_args()

    try:
        assembler = TranscriptAssembler(
            genome_file=args.genome,
            input_dir=args.input,
            output_dir=args.output,
            threads=args.threads,
            fastq_pattern=args.pattern,
            sample_timeout=args.sample_timeout,
            step=args.step,
            verbose=(args.verbose > 0),
            quiet=args.quiet,
            log_level=args.log_level,
            force=args.force
        )

        assembler.run_analysis()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("操作被用户中断|Operation interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
