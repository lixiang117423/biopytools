"""
转录组验证注释主程序模块|Transcriptome Validation Main Pipeline Module
"""

import os
import sys
import argparse
import time
from typing import List, Optional

from .config import RnaseqValConfig
from .utils import RnaseqValLogger, CommandRunner, SampleParser, build_conda_command, format_number
from .align_2nd import HISAT2Indexer, HISAT2Aligner
from .align_3rd import Minimap2Aligner
from .assemble_2nd import StringTieAssembler
from .assemble_3rd import FlairAssembler, IsoQuantRunner
from .compare import GFFCompareRunner
from .correct import (
    read_gffcompare_annotated,
    classify_transcripts,
    CorrectionExporter,
)
from .report import SummaryReporter


class RnaseqValPipeline:
    """转录组验证注释 Pipeline|Transcriptome Validation Pipeline"""

    def __init__(self, **kwargs):
        """初始化 Pipeline|Initialize pipeline

        Args:
            **kwargs: 传入 RnaseqValConfig 的参数|Parameters for RnaseqValConfig
        """
        # 初始化配置|Initialize config
        self.config = RnaseqValConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = RnaseqValLogger(
            self.config.log_dir,
            verbose=self.config.verbose,
            quiet=self.config.quiet,
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化样本解析器|Initialize sample parser
        self.sample_parser = SampleParser(self.logger)

        # 解析步骤列表|Parse step list
        self.steps = self.config.get_steps_list()

        # 内部状态|Internal state
        self.sr_bam_files: List[str] = []
        self.lr_bam_files: List[str] = []
        self.sr_merged_gtf: Optional[str] = None
        self.lr_gtf: Optional[str] = None
        self.gffcmp_prefix: Optional[str] = None

    def run(self):
        """运行 Pipeline|Run the pipeline"""
        start_time = time.time()

        try:
            self._print_header()
            self._detect_samples()
            self._execute_steps(start_time)

        except KeyboardInterrupt:
            self.logger.warning("用户中断|Interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"Pipeline 异常终止|Pipeline terminated unexpectedly: {e}", exc_info=True)
            sys.exit(1)

    def _print_header(self):
        """打印 Pipeline 头部信息|Print pipeline header"""
        # 创建 pipeline_info 目录并写入软件版本信息|Create pipeline_info dir and write software versions
        pipeline_info_dir = os.path.join(self.config.output_dir, "00_pipeline_info")
        os.makedirs(pipeline_info_dir, exist_ok=True)

        versions_file = os.path.join(pipeline_info_dir, "software_versions.yml")
        with open(versions_file, "w") as f:
            f.write("# RNA-seq Validation Pipeline - Software Versions\n")
            f.write(f"# Generated: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write("tools:\n")
            f.write("  python: \"3.x\"\n")
            f.write("  conda_env: \"%s\"\n" % self.config.conda_env)
            f.write("# Tool versions will be populated during pipeline execution\n")

        self.logger.info("=" * 60)
        self.logger.info("转录组验证注释 Pipeline|Transcriptome Validation Pipeline")
        self.logger.info("=" * 60)
        self.logger.info(f"基因组|Genome: {self.config.genome_fa}")
        self.logger.info(f"注释文件|Annotation: {self.config.annotation_gtf}")
        self.logger.info(f"输出目录|Output: {self.config.output_dir}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"链特异性|Strandness: {self.config.strandness}")
        if self.config.sr_dir:
            self.logger.info(f"二代数据|SR data: {self.config.sr_dir}")
        if self.config.lr_dir:
            self.logger.info(f"三代数据|LR data: {self.config.lr_dir} ({self.config.lr_platform})")
        self.logger.info(f"执行步骤|Steps: {', '.join(self.steps)}")
        self.logger.info("=" * 60)

    def _detect_samples(self):
        """检测输入样本|Detect input samples"""
        # 解析二代样本|Parse SR samples
        if self.config.sr_dir:
            self.config.samples_sr = self.sample_parser.parse_sr_samples(
                self.config.sr_dir,
                self.config.sr_pattern,
            )
            if not self.config.samples_sr:
                self.logger.warning("未检测到二代样本|No SR samples detected")

        # 解析三代样本|Parse LR samples
        if self.config.lr_dir:
            self.config.samples_lr = self.sample_parser.parse_lr_samples(
                self.config.lr_dir,
            )
            if not self.config.samples_lr:
                self.logger.warning("未检测到三代样本|No LR samples detected")

    def _execute_steps(self, start_time: float):
        """按顺序执行步骤|Execute steps in order"""
        steps_map = {
            "align_2nd": self._step_align_2nd,
            "align_3rd": self._step_align_3rd,
            "assemble_2nd": self._step_assemble_2nd,
            "assemble_3rd": self._step_assemble_3rd,
            "compare": self._step_compare,
            "correct": self._step_correct,
            "report": self._step_report,
        }

        for step_name in self.steps:
            if step_name not in steps_map:
                self.logger.warning(f"未知步骤，跳过|Unknown step, skipping: {step_name}")
                continue

            # 检查前置条件|Check prerequisites
            if step_name == "align_2nd" and not self.config.samples_sr:
                self.logger.info("无二代数据，跳过比对|No SR data, skipping align_2nd")
                continue
            if step_name == "align_3rd" and not self.config.samples_lr:
                self.logger.info("无三代数据，跳过比对|No LR data, skipping align_3rd")
                continue
            if step_name == "assemble_2nd" and not self.sr_bam_files:
                self.logger.info("无二代 BAM，跳过组装|No SR BAM, skipping assemble_2nd")
                continue
            if step_name == "assemble_3rd" and not self.lr_bam_files:
                self.logger.info("无三代 BAM，跳过组装|No LR BAM, skipping assemble_3rd")
                continue
            if step_name in ("compare", "correct", "report") and not self.sr_merged_gtf and not self.lr_gtf:
                self.logger.info("无组装结果，跳过后续步骤|No assembly results, skipping downstream steps")
                continue

            steps_map[step_name]()

        # 最终汇总|Final summary
        elapsed = time.time() - start_time
        self.logger.info("=" * 60)
        self.logger.info(f"Pipeline 完成|Pipeline completed in {elapsed:.1f}秒|seconds")
        self.logger.info(f"输出目录|Output: {self.config.output_dir}")
        self.logger.info("=" * 60)

    # ============================================================
    # 步骤实现|Step Implementations
    # ============================================================

    def _step_align_2nd(self):
        """步骤: 二代比对|Step: 2nd gen alignment"""
        self.logger.step("步骤: 二代转录组比对 (HISAT2)|Step: 2nd gen alignment (HISAT2)")

        indexer = HISAT2Indexer(self.config, self.logger, self.cmd_runner)
        index_prefix = indexer.build_index()

        aligner = HISAT2Aligner(self.config, self.logger, self.cmd_runner)
        self.sr_bam_files = aligner.align_all_samples(self.config.samples_sr, index_prefix)

    def _step_align_3rd(self):
        """步骤: 三代比对|Step: 3rd gen alignment"""
        self.logger.step("步骤: 三代转录组比对 (minimap2)|Step: 3rd gen alignment (minimap2)")

        aligner = Minimap2Aligner(self.config, self.logger, self.cmd_runner)
        self.lr_bam_files = aligner.align_all_samples(self.config.samples_lr)

    def _step_assemble_2nd(self):
        """步骤: 二代组装|Step: 2nd gen assembly"""
        self.logger.step("步骤: 二代转录本组装 (StringTie)|Step: 2nd gen assembly (StringTie)")

        assembler = StringTieAssembler(self.config, self.logger, self.cmd_runner)
        self.sr_merged_gtf = assembler.run_all(self.sr_bam_files)

    def _step_assemble_3rd(self):
        """步骤: 三代组装|Step: 3rd gen assembly"""
        self.logger.step("步骤: 三代转录本组装 (FLAIR)|Step: 3rd gen assembly (FLAIR)")

        assembler = FlairAssembler(self.config, self.logger, self.cmd_runner)

        # 提取二代 junction 用于辅助校正 ONT 数据
        sr_junctions = None
        if self.config.lr_platform == "ont" and self.sr_bam_files:
            sr_junctions = self._extract_sr_junctions()

        self.lr_gtf = assembler.run_all(
            samples=self.config.samples_lr,
            bam_files=self.lr_bam_files,
            sr_junctions=sr_junctions,
        )

    def _extract_sr_junctions(self) -> Optional[str]:
        """从二代 BAM 提取 junction 文件用于辅助 ONT 校正|Extract SR junctions for ONT correction

        Returns:
            Optional[str]: junction BED 文件路径|junction BED file path
        """
        if not self.sr_bam_files:
            return None

        out_dir = os.path.join(self.config.output_dir, "05_assemble_3rd")
        os.makedirs(out_dir, exist_ok=True)
        junction_file = os.path.join(out_dir, "sr_junctions.bed")

        if os.path.exists(junction_file):
            return junction_file

        # 使用第一个 BAM 的 junction（多样本时合并会更好，这里简化处理）
        bam = self.sr_bam_files[0]
        cmd = (
            f"regtools junctions extract -a 8 -m 50 -M 500000 "
            f"{bam} > {junction_file}"
        )
        cmd = build_conda_command(cmd, self.config.conda_env)
        success = self.cmd_runner.run(cmd, "提取二代 junction|Extracting SR junctions")
        if success:
            self.logger.info(f"二代 junction 文件|SR junction file: {junction_file}")
            return junction_file

        self.logger.warning("二代 junction 提取失败|SR junction extraction failed")
        return None

    def _step_compare(self):
        """步骤: GFFcompare 结构比较|Step: GFFcompare structure comparison"""
        self.logger.step("步骤: 转录本结构比较 (GFFcompare)|Step: Transcript comparison (GFFcompare)")

        query_gtfs = []
        if self.sr_merged_gtf:
            query_gtfs.append(self.sr_merged_gtf)
        if self.lr_gtf:
            query_gtfs.append(self.lr_gtf)

        if not query_gtfs:
            self.logger.error("无 query GTF 可比较|No query GTFs to compare")
            return

        runner = GFFCompareRunner(self.config, self.logger, self.cmd_runner)
        self.gffcmp_prefix = runner.run(query_gtfs)

    def _step_correct(self):
        """步骤: 自动校正决策|Step: Automatic correction decision"""
        self.logger.step("步骤: 自动校正决策|Step: Automatic correction decision")

        if not self.gffcmp_prefix:
            self.logger.error("无 GFFcompare 结果，跳过校正|No GFFcompare results, skipping correction")
            return

        annotated_gtf = f"{self.gffcmp_prefix}.annotated.gtf"
        if not os.path.exists(annotated_gtf):
            self.logger.error(f"annotated GTF 不存在|annotated GTF not found: {annotated_gtf}")
            return

        # 读取 annotated GTF|Read annotated GTF
        records = read_gffcompare_annotated(annotated_gtf)
        if not records:
            self.logger.error("annotated GTF 无有效记录|No valid records in annotated GTF")
            return

        self.logger.info(f"读取 {format_number(len(records))} 条转录本记录|Read {format_number(len(records))} transcript records")

        # 分级决策|Classification
        thresholds = {
            "min_cov": self.config.min_cov,
            "min_tpm": self.config.min_tpm,
            "min_tpm_sr_only": self.config.min_tpm_sr_only,
            "min_junction_ont": self.config.min_junction_ont,
        }

        classified = classify_transcripts(
            records=records,
            sr_assembly_gtf=self.sr_merged_gtf,
            lr_assembly_gtf=self.lr_gtf,
            thresholds=thresholds,
        )

        # 导出结果|Export results
        exporter = CorrectionExporter(self.config, self.logger)
        exporter.export_results(classified)

        # 写入校正后 GTF|Write corrected GTF
        exporter.write_corrected_gtf(classified, self.config.annotation_gtf)

    def _step_report(self):
        """步骤: 统计汇总|Step: Summary report"""
        self.logger.step("步骤: 统计汇总|Step: Summary report")

        reporter = SummaryReporter(self.config, self.logger, self.cmd_runner)
        reporter.generate_summary(
            gffcmp_prefix=self.gffcmp_prefix,
            correction_dir=os.path.join(self.config.output_dir, "07_correction"),
        )


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="转录组验证注释 Pipeline|Transcriptome Validation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -g genome.fa -a annotation.gtf --sr-dir ./clean_reads/ -o rnaseq_val_out
  %(prog)s -g genome.fa -a annotation.gtf --lr-dir ./lr_reads/ --lr-platform pacbio -o out
  %(prog)s -g genome.fa -a annotation.gtf --sr-dir ./sr/ --lr-dir ./lr/ --steps compare,correct,report -o out
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument("-g", "--genome", required=True,
                          help="基因组 FASTA 文件|Genome FASTA file")
    required.add_argument("-a", "--annotation", required=True,
                          help="参考注释 GTF 文件|Reference annotation GTF file")
    required.add_argument("-o", "--output", required=True,
                          help="输出目录|Output directory")

    # 数据输入|Data input
    data = parser.add_argument_group('数据输入|Data input')
    data.add_argument("--sr-dir", default=None,
                      help="二代 clean reads 目录（自动检测配对 fastq）|SR clean reads directory (auto-detect paired fastq)")
    data.add_argument("--sr-pattern", default=None,
                      help="二代 fastq 自定义命名模式|SR fastq custom pattern (e.g. *_1.clean.fq.gz)")
    data.add_argument("--lr-dir", default=None,
                      help="三代 clean reads 目录（自动检测单端 fastq）|LR clean reads directory (auto-detect single fastq)")
    data.add_argument("--lr-platform", default="pacbio",
                      choices=["pacbio", "ont"],
                      help="三代测序平台|LR sequencing platform")

    # 参数设置|Parameter settings
    params = parser.add_argument_group('参数设置|Parameters')
    params.add_argument("-t", "--threads", type=int, default=12,
                        help="线程数|Threads")
    params.add_argument("--strandness", default="RF",
                        choices=["RF", "FR", "unstranded"],
                        help="文库链特异性|Library strandness")
    params.add_argument("--max-intron", type=int, default=500000,
                        help="最大内含子长度|Max intron length")
    params.add_argument("--steps", default="all",
                        help="执行步骤（逗号分隔或 all）|Steps to run (comma-separated or all)")
    params.add_argument("--sample-timeout", type=int, default=21600,
                        help="单样本超时时间(秒)|Sample timeout in seconds")

    # StringTie 参数|StringTie parameters
    st = parser.add_argument_group('StringTie 参数|StringTie parameters')
    st.add_argument("--st-min-cov", type=float, default=5.0,
                    help="StringTie 最低覆盖度|StringTie min coverage")
    st.add_argument("--st-min-junction", type=int, default=3,
                    help="StringTie 最低 junction reads|StringTie min junction reads")
    st.add_argument("--st-min-isoform", type=float, default=0.1,
                    help="StringTie 最低亚型丰度比|StringTie min isoform fraction")

    # 校正阈值|Correction thresholds
    th = parser.add_argument_group('校正阈值|Correction thresholds')
    th.add_argument("--min-cov", type=float, default=5.0,
                    help="最低覆盖度阈值|Min coverage threshold")
    th.add_argument("--min-tpm", type=float, default=0.5,
                    help="最低 TPM 阈值|Min TPM threshold")
    th.add_argument("--min-tpm-sr-only", type=float, default=1.0,
                    help="仅二代支持时最低 TPM|Min TPM for SR-only support")
    th.add_argument("--min-junction-ont", type=int, default=5,
                    help="ONT 最低 junction reads|Min junction reads for ONT")

    # 日志选项|Logging options
    logging_group = parser.add_argument_group('日志选项|Logging options')
    logging_group.add_argument("-v", "--verbose", action="count", default=0,
                               help="详细模式 (-v, -vv)|Verbose mode")
    logging_group.add_argument("--quiet", action="store_true",
                               help="静默模式|Quiet mode")

    # 高级选项|Advanced options
    advanced = parser.add_argument_group('高级选项|Advanced options')
    advanced.add_argument("--force", action="store_true",
                          help="强制重新运行|Force re-run")
    advanced.add_argument("--dry-run", action="store_true",
                          help="试运行模式|Dry run mode")

    args = parser.parse_args()

    try:
        pipeline = RnaseqValPipeline(
            genome_fa=args.genome,
            annotation_gtf=args.annotation,
            output_dir=args.output,
            sr_dir=args.sr_dir,
            sr_pattern=args.sr_pattern,
            lr_dir=args.lr_dir,
            lr_platform=args.lr_platform,
            threads=args.threads,
            strandness=args.strandness,
            max_intron=args.max_intron,
            steps=args.steps,
            sample_timeout=args.sample_timeout,
            stringtie_min_cov=args.st_min_cov,
            stringtie_min_junction_reads=args.st_min_junction,
            stringtie_min_isoform_fraction=args.st_min_isoform,
            min_cov=args.min_cov,
            min_tpm=args.min_tpm,
            min_tpm_sr_only=args.min_tpm_sr_only,
            min_junction_ont=args.min_junction_ont,
            verbose=(args.verbose > 0),
            quiet=args.quiet,
            dry_run=args.dry_run,
            force=args.force,
        )
        pipeline.run()

    except ValueError as e:
        print(f"参数错误|Parameter error: {e}", file=sys.stderr)
        sys.exit(1)
    except KeyboardInterrupt:
        print("用户中断|Interrupted by user", file=sys.stderr)
        sys.exit(130)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
