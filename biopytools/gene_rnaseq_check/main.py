"""候选基因RNA-seq转录验证主程序模块|Candidate Gene RNA-seq Validation Main Module"""

import argparse
import os
import sys
import time
from typing import Dict, List, Optional

from .config import GeneRnaseqCheckConfig
from .utils import (
    GeneRnaseqCheckLogger, CommandRunner,
    generate_software_versions, format_number,
    build_conda_command, get_chrom_lengths,
)


class GeneRnaseqCheckPipeline:
    """候选基因RNA-seq转录验证Pipeline|Gene RNA-seq Check Pipeline"""

    def __init__(self, **kwargs):
        self.config = GeneRnaseqCheckConfig(**kwargs)
        self.config.validate()

        self.logger_manager = GeneRnaseqCheckLogger(
            self.config.log_dir,
            verbose=self.config.verbose,
            quiet=self.config.quiet,
        )
        self.logger = self.logger_manager.get_logger()
        self.cmd_runner = CommandRunner(self.logger)

        self.steps = self.config.get_steps_list()

        # 模块间共享数据|Inter-module shared data
        self.records: Dict = {}
        self.bam_files: List[str] = []
        self.chrom_lengths: Dict[str, int] = {}
        self.coverage_data: Dict = {}
        self.junction_data: Dict = {}
        self.class_codes: Dict = {}
        self.classifications: List = []
        self.boundary_suggestions: List = []

    def run(self):
        """运行Pipeline|Run the pipeline"""
        start_time = time.time()
        try:
            self._print_header(start_time)
            self._execute_steps()
            self._print_summary(start_time)
        except KeyboardInterrupt:
            self.logger.warning("用户中断|Interrupted by user")
            sys.exit(130)
        except Exception as e:
            self.logger.error(f"Pipeline异常终止|Pipeline terminated unexpectedly: {e}",
                             exc_info=True)
            sys.exit(1)

    def _print_header(self, start_time):
        """打印Pipeline头部信息|Print pipeline header"""
        self.logger.info("=" * 60)
        self.logger.info("候选基因RNA-seq转录验证|Gene RNA-seq Transcriptional Validation Pipeline")
        self.logger.info("=" * 60)
        self.logger.info(f"基因组|Genome: {self.config.genome_fa}")
        self.logger.info(f"注释文件|Annotation: {self.config.annotation_gff}")
        self.logger.info(f"目标基因数|Target genes: {len(self.config._gene_ids)}")
        self.logger.info(f"Reads目录|Reads dir: {self.config.reads_dir}")
        self.logger.info(f"检测到样本|Detected samples: {len(self.config._samples)}")
        self.logger.info(f"输出目录|Output: {self.config.output_dir}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"执行步骤|Steps: {', '.join(self.steps)}")
        self.logger.info("=" * 60)

        # 生成软件版本信息|Generate software versions
        tools = {
            'hisat2': self.config.hisat2_path,
            'hisat2-build': self.config.hisat2_build_path,
            'samtools': self.config.samtools_path,
            'bedtools': self.config.bedtools_path,
            'stringtie': self.config.stringtie_path,
            'gffcompare': self.config.gffcompare_path,
            'infer_experiment.py': self.config.infer_experiment_path,
        }
        params = {
            'threads': self.config.threads,
            'strandness_confidence': self.config.strandness_confidence,
            'flanking_window': self.config.flanking_window,
            'junction_tolerance': self.config.junction_tolerance,
        }
        generate_software_versions(
            self.config.output_dir, tools, params, start_time
        )

    def _execute_steps(self):
        """按顺序执行步骤|Execute steps in order"""
        steps_map = {
            'parse_gff': self._step_parse_gff,
            'align': self._step_align,
            'coverage': self._step_coverage,
            'junction': self._step_junction,
            'stringtie': self._step_stringtie,
            'classify': self._step_classify,
            'boundary': self._step_boundary,
        }
        for step_name in self.steps:
            if step_name not in steps_map:
                self.logger.warning(f"未知步骤，跳过|Unknown step, skipping: {step_name}")
                continue

            # 前置条件检查|Check prerequisites
            if step_name in ('coverage', 'junction', 'stringtie') and not self.bam_files:
                self.logger.info(f"无BAM文件，跳过 {step_name}|No BAM files, skipping {step_name}")
                continue
            if step_name in ('classify',) and not self.coverage_data:
                self.logger.info(f"无覆盖度数据，跳过 {step_name}|No coverage data, skipping {step_name}")
                continue
            if step_name == 'boundary' and not self.classifications:
                self.logger.info(f"无分类结果，跳过 {step_name}|No classification, skipping {step_name}")
                continue

            steps_map[step_name]()

    # ============================================================
    # 步骤实现|Step Implementations
    # ============================================================

    def _step_parse_gff(self):
        """步骤: 解析GFF3注释|Step: Parse GFF3 annotation"""
        self.logger_manager.step("步骤: 解析GFF3注释|Step: Parse GFF3 annotation")

        from .parse_gff import parse_gff3_effector_genes

        self.records = parse_gff3_effector_genes(
            self.config.annotation_gff,
            self.config._gene_ids,
            self.logger,
        )

        if not self.records:
            self.logger.error("未解析到任何目标基因|No target genes parsed from GFF3")
            sys.exit(1)

        self.logger.info(f"成功解析|Successfully parsed {len(self.records)} 个目标基因|target genes")

    def _step_align(self):
        """步骤: 建索引与比对|Step: Index and align"""
        self.logger_manager.step("步骤: 建索引与比对 (HISAT2)|Step: HISAT2 indexing and alignment")

        from .align import HISAT2Indexer, HISAT2Aligner, StrandnessDetector

        # 建索引|Build index
        indexer = HISAT2Indexer(self.config, self.logger, self.cmd_runner)
        index_prefix = indexer.build_index()
        if not index_prefix:
            self.logger.error("HISAT2索引构建失败|HISAT2 index building failed")
            sys.exit(1)

        # 初始比对（无strandness）|Initial alignment (no strandness)
        aligner = HISAT2Aligner(self.config, self.logger, self.cmd_runner)
        self.bam_files = aligner.align_all_samples(
            self.config._samples, index_prefix
        )

        if not self.bam_files:
            self.logger.error("所有样本比对失败|All samples failed alignment")
            sys.exit(1)

        self.logger.info(f"成功比对|Successfully aligned {len(self.bam_files)} 个样本|samples")

        # 获取染色体长度|Get chromosome lengths
        self.chrom_lengths = get_chrom_lengths(self.bam_files[0])

        # 链特异性检测|Strandness detection
        self.logger.info("检测建库链特异性|Detecting library strandness")
        detector = StrandnessDetector(self.config, self.logger, self.cmd_runner)
        strandness, confidence = detector.detect(self.bam_files[0], self.chrom_lengths)
        self.config.strandness = strandness

        # 如果需要重新比对|Re-align if strand-specific with high confidence
        if strandness != 'unstranded' and confidence >= self.config.strandness_confidence:
            self.logger.info(
                f"检测到链特异性建库（{strandness}，置信度 {confidence:.1f}%），重新比对|"
                f"Strand-specific ({strandness}, {confidence:.1f}%), re-aligning"
            )
            # 清理旧BAM|Clean old BAMs
            for bam in self.bam_files:
                for ext in ['', '.bai']:
                    p = bam + ext
                    if os.path.exists(p):
                        os.remove(p)
            # 用正确的strandness重新比对|Re-align with correct strandness
            self.bam_files = aligner.align_all_samples(
                self.config._samples, index_prefix, strandness=strandness
            )
        else:
            self.logger.info(f"使用unstranded模式|Using unstranded mode (confidence: {confidence:.1f}%)")

    def _step_coverage(self):
        """步骤: 覆盖度分析|Step: Coverage analysis"""
        self.logger_manager.step("步骤: 覆盖度分析 (bedtools)|Step: Coverage analysis (bedtools)")

        from .coverage import CoverageAnalyzer

        analyzer = CoverageAnalyzer(self.config, self.logger, self.cmd_runner)
        self.coverage_data = analyzer.compute_coverage(
            self.records, self.bam_files, self.chrom_lengths
        )

    def _step_junction(self):
        """步骤: Junction reads分析|Step: Junction reads analysis"""
        self.logger_manager.step("步骤: Junction reads分析|Step: Junction reads analysis")

        from .junction import JunctionAnalyzer

        analyzer = JunctionAnalyzer(self.config, self.logger)
        self.junction_data = analyzer.analyze_gene_junctions(
            self.bam_files, self.records, self.config.genome_fa
        )

    def _step_stringtie(self):
        """步骤: StringTie组装与gffcompare|Step: StringTie assembly and gffcompare"""
        self.logger_manager.step("步骤: StringTie组装与gffcompare|Step: StringTie assembly and gffcompare")

        from .stringtie_run import StringTieRunner, GFFCompareRunner, extract_effector_class_codes

        # StringTie组装|StringTie assembly
        runner = StringTieRunner(self.config, self.logger, self.cmd_runner)
        merged_gtf = runner.assemble_all(self.config._samples, self.bam_files)

        if not merged_gtf:
            self.logger.warning("StringTie组装失败，所有class_code设为'u'|StringTie failed, all class_code set to 'u'")
            self.class_codes = {gid: 'u' for gid in self.records}
            return

        # gffcompare|gffcompare
        gff_runner = GFFCompareRunner(self.config, self.logger, self.cmd_runner)
        gffcmp_prefix = gff_runner.run(merged_gtf)

        if not gffcmp_prefix:
            self.logger.warning("gffcompare失败，所有class_code设为'u'|gffcompare failed, all class_code set to 'u'")
            self.class_codes = {gid: 'u' for gid in self.records}
            return

        # 提取class_code|Extract class codes
        self.class_codes = extract_effector_class_codes(
            gffcmp_prefix, self.config._gene_ids, self.logger
        )

    def _step_classify(self):
        """步骤: 分类与报告生成|Step: Classification and report generation"""
        self.logger_manager.step("步骤: 分类与报告生成|Step: Classification and report generation")

        from .classify import classify_all_genes, write_validation_report, write_flagged_bed

        self.classifications = classify_all_genes(
            self.records, self.coverage_data, self.junction_data,
            self.class_codes, self.config, self.logger,
        )

        # 写入报告|Write reports
        results_dir = os.path.join(self.config.output_dir, '05_results')
        os.makedirs(results_dir, exist_ok=True)

        report_file = os.path.join(results_dir, 'gene_validation_report.tsv')
        write_validation_report(self.classifications, report_file)
        self.logger.info(f"主结果表格已生成|Main report generated: {report_file}")

        write_flagged_bed(self.classifications, results_dir)

        # 写入所有基因的BED（用于IGV批量查看）|Write BED for all genes
        all_bed = os.path.join(results_dir, 'all_target_genes.bed')
        with open(all_bed, 'w') as f:
            for r in self.classifications:
                f.write(f"{r.chrom}\t{r.start - 1}\t{r.end}\t{r.gene_id}\t0\t{r.strand}\n")

    def _step_boundary(self):
        """步骤: 边界修改建议|Step: Boundary suggestions"""
        self.logger_manager.step("步骤: 边界修改建议|Step: Boundary suggestions")

        from .boundary import BoundaryAnalyzer, write_boundary_suggestions

        # 找StringTie GTF
        stringtie_gtf = ''
        gtf_candidates = [
            os.path.join(self.config.output_dir, '04_stringtie', 'merged.gtf'),
        ]
        # 如果多样本，也检查第一个样本的GTF
        if self.config._samples:
            gtf_candidates.insert(0, os.path.join(
                self.config.output_dir, '04_stringtie',
                f"{self.config._samples[0]['name']}.gtf"
            ))
        for gtf_path in gtf_candidates:
            if os.path.exists(gtf_path):
                stringtie_gtf = gtf_path
                break

        if not stringtie_gtf:
            self.logger.warning("未找到StringTie GTF，跳过边界分析|No StringTie GTF found, skipping boundary analysis")
            return

        analyzer = BoundaryAnalyzer(self.config, self.logger)
        self.boundary_suggestions = analyzer.suggest_boundaries(
            self.classifications, self.records,
            self.coverage_data, self.junction_data,
            stringtie_gtf,
        )

        if self.boundary_suggestions:
            results_dir = os.path.join(self.config.output_dir, '05_results')
            boundary_file = os.path.join(results_dir, 'boundary_suggestions.tsv')
            write_boundary_suggestions(self.boundary_suggestions, boundary_file)
            self.logger.info(f"边界建议表已生成|Boundary suggestions generated: {boundary_file}")
            self.logger.info(f"共 {len(self.boundary_suggestions)} 个基因有边界建议|"
                             f"{len(self.boundary_suggestions)} genes have boundary suggestions")
        else:
            self.logger.info("无边界问题需要建议|No boundary issues need suggestions")

    def _print_summary(self, start_time):
        """打印汇总信息|Print summary"""
        elapsed = time.time() - start_time
        self.logger.info("=" * 60)
        self.logger.info(f"Pipeline 完成|Pipeline completed in {elapsed:.1f}秒|seconds")
        self.logger.info(f"输出目录|Output: {self.config.output_dir}")
        self.logger.info(f"主要结果|Main results:")
        self.logger.info(f"  - 05_results/gene_validation_report.tsv")
        if self.boundary_suggestions:
            self.logger.info(f"  - 05_results/boundary_suggestions.tsv")
        self.logger.info("=" * 60)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="候选基因RNA-seq转录验证|Gene RNA-seq Transcriptional Validation Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -g genome.fa -a annotation.gff -e gene_ids.txt -r ./reads/ -o output/
  %(prog)s -g genome.fa -a annotation.gff -e gene_ids.txt -r ./reads/ -o output/ --steps parse_gff,coverage,classify
        ''',
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument('-g', '--genome', required=True,
                          help='基因组FASTA文件|Genome FASTA file')
    required.add_argument('-a', '--annotation', required=True,
                          help='GFF3注释文件|GFF3 annotation file')
    required.add_argument('-e', '--gene-list', required=True,
                          help='目标基因ID列表文件|Target gene ID list file')
    required.add_argument('-r', '--reads-dir', required=True,
                          help='RNA-seq reads目录|RNA-seq reads directory')
    required.add_argument('-o', '--output-dir', default='./gene_rnaseq_check_output',
                          help='输出目录|Output directory')

    # 参数设置|Parameters
    params = parser.add_argument_group('参数设置|Parameters')
    params.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Threads')
    params.add_argument('--reads-pattern', default='',
                        help='FASTQ命名模式|FASTQ naming pattern (e.g. *_1.fq.gz)')
    params.add_argument('--steps', default='all',
                        help='执行步骤(逗号分隔或all)|Steps (comma-separated or all)')
    params.add_argument('--sample-timeout', type=int, default=21600,
                        help='单样本超时时间(秒)|Sample timeout in seconds')
    params.add_argument('--flanking-window', type=int, default=500,
                        help='上下游分析窗口(bp)|Flanking analysis window (bp)')
    params.add_argument('--junction-tolerance', type=int, default=5,
                        help='Junction容差(bp)|Junction tolerance (bp)')

    # 分类阈值|Classification thresholds
    th = parser.add_argument_group('分类阈值|Classification thresholds')
    th.add_argument('--strandness-confidence', type=float, default=70.0,
                    help='链特异性判定置信度(%%)|Strandness confidence threshold (%%)')

    # 日志选项|Logging options
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument('-v', '--verbose', action='store_true',
                           help='详细模式|Verbose mode')
    log_group.add_argument('--quiet', action='store_true',
                          help='静默模式|Quiet mode')

    # 高级选项|Advanced options
    adv = parser.add_argument_group('高级选项|Advanced options')
    adv.add_argument('--force', action='store_true',
                    help='强制重新运行|Force re-run')

    args = parser.parse_args()

    try:
        pipeline = GeneRnaseqCheckPipeline(
            genome_fa=args.genome,
            annotation_gff=args.annotation,
            gene_list=args.gene_list,
            reads_dir=args.reads_dir,
            output_dir=args.output_dir,
            reads_pattern=args.reads_pattern,
            threads=args.threads,
            steps=args.steps,
            sample_timeout=args.sample_timeout,
            flanking_window=args.flanking_window,
            junction_tolerance=args.junction_tolerance,
            strandness_confidence=args.strandness_confidence,
            verbose=args.verbose,
            quiet=args.quiet,
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
