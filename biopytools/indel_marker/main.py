"""INDEL分子标记主程序|INDEL Marker Main"""

import argparse
import sys
import os
import datetime
from pathlib import Path

from .config import IndelMarkerConfig
from .utils import IndelMarkerLogger, CommandRunner, check_dependencies
from .samplesheet import SamplesheetParser, SampleInfo
from .vcf_extractor import VCFExtractor
from .genotype_analyzer import classify_indel, IndelRecord, Candidate
from .coverage_validator import CoverageValidator, candidate_id
from .sequence_extractor import SequenceExtractor
from .primer_designer import PrimerDesigner
from .results_writer import ResultsWriter


class IndelMarkerRunner:
    """INDEL分子标记编排器|INDEL Marker Orchestrator"""

    def __init__(self, config: IndelMarkerConfig):
        self.config = config
        self.logger_manager = IndelMarkerLogger(config.logs_dir)
        self.logger = self.logger_manager.get_logger()
        self.cmd = CommandRunner(self.logger, Path.cwd())
        self.samples: dict = {}
        self.r_samples: list = []
        self.s_samples: list = []
        self.start_time = datetime.datetime.now()

    def _is_step_completed(self, output_file: str) -> bool:
        """断点续传：输出存在则跳过|Checkpoint: skip if output exists"""
        return bool(output_file) and os.path.exists(output_file)

    # ---- 各步骤|steps ----
    def _check_deps(self) -> bool:
        return check_dependencies(self.config, self.logger)

    def _parse_samplesheet(self):
        self.samples = SamplesheetParser.parse(self.config.samplesheet, logger=self.logger)
        self.r_samples, self.s_samples = SamplesheetParser.split_groups(self.samples)
        # 校验每组样品数|min samples per group
        errors = []
        if len(self.r_samples) < self.config.min_samples_per_group:
            errors.append(f"抗病组样品数|resistant count {len(self.r_samples)} < "
                          f"{self.config.min_samples_per_group}")
        if len(self.s_samples) < self.config.min_samples_per_group:
            errors.append(f"感病组样品数|susceptible count {len(self.s_samples)} < "
                          f"{self.config.min_samples_per_group}")
        # 校验bam存在|check bam existence
        for s in self.samples.values():
            if not os.path.exists(s.bam_path):
                errors.append(f"bam不存在|bam not found: {s.bam_path} ({s.name})")
        if errors:
            raise ValueError("\n".join(errors))
        self.logger.info(f"样品|Samples: 抗病|R={len(self.r_samples)} 感病|S={len(self.s_samples)}")

    def _check_samples_in_vcf(self, vcf_samples):
        missing = set(self.samples) - set(vcf_samples)
        if missing:
            raise ValueError(f"samplesheet样品不在VCF中|samples not in VCF: {missing}")

    def _extract_indels(self) -> list:
        out_tsv = str(self.config.vcf_extract_dir / 'indels.gt_matrix.tsv')
        ve = VCFExtractor(self.config, self.logger, self.cmd)
        # 样本名只取一次，复用给extract避免重复调用get_sample_names
        # fetch sample names once, reuse for extract (avoid double get_sample_names)
        vcf_samples = ve.get_sample_names()
        self._check_samples_in_vcf(vcf_samples)
        if self._is_step_completed(out_tsv):
            self.logger.info("跳过已完成步骤|Skipping completed step: vcf_extract")
            return ve._parse_matrix(out_tsv, vcf_samples)
        return ve.extract(out_tsv, samples=vcf_samples)

    def _classify(self, indels: list) -> list:
        candidates = []
        for ind in indels:
            c = classify_indel(ind, self.r_samples, self.s_samples, self.config)
            if c is not None:
                candidates.append(c)
        self.logger.info(f"群体判定候选|genotype candidates: {len(candidates)} / {len(indels)}")
        # 候选数限制（0=不限）|candidate cap (0=no limit)
        if self.config.max_candidates > 0 and len(candidates) > self.config.max_candidates:
            self.logger.info(
                f"候选数限制|limiting to max_candidates={self.config.max_candidates}")
            candidates = candidates[:self.config.max_candidates]
        return candidates

    def _compute_coverage(self, candidates):
        cv = CoverageValidator(self.config, self.logger, self.cmd, self.samples)
        return cv.compute(candidates)

    def _validate_coverage(self, candidates, coverage):
        cv = CoverageValidator(self.config, self.logger, self.cmd, self.samples)
        cv.validate(candidates, coverage)

    def _extract_sequences(self, candidates):
        se = SequenceExtractor(self.config, self.logger, self.cmd)
        return se.extract_flank(candidates)

    def _design_primers(self, candidates, sequences):
        pd = PrimerDesigner(self.config, self.logger)
        primers = {}
        for c in candidates:
            cid = candidate_id(c)
            primers[cid] = pd.design(cid, sequences.get(cid, ''))
        return primers

    def _gather_tool_versions(self) -> dict:
        """收集工具版本|Gather tool versions (bcftools/samtools/primer3-py)"""
        from .utils import build_conda_command
        versions = {}
        for name, path in [('bcftools', self.config.bcftools_path),
                           ('samtools', self.config.samtools_path)]:
            try:
                out = self.cmd.run_capture(
                    build_conda_command(path, ['--version']),
                    description=f"{name} --version")
                versions[name] = out.strip().splitlines()[0] if out else 'unavailable'
            except Exception:
                versions[name] = 'unavailable'
        try:
            import primer3
            versions['primer3'] = getattr(primer3, '__version__', 'unavailable')
        except Exception:
            versions['primer3'] = 'unavailable'
        return versions

    def _write_results(self, candidates, coverage, sequences, primers,
                       n_total_indels):
        rw = ResultsWriter(self.config, self.logger)
        rw.write_candidates_tsv(candidates, primers)
        rw.write_bed(candidates)
        rw.write_fasta(sequences)
        rw.write_coverage_tsv(coverage, sorted(self.samples.keys()))
        n_pass_qc = sum(1 for c in candidates if getattr(c, 'passes_coverage_qc', False))
        n_with_primer = sum(1 for p in primers.values() if p.get('primer_status') == 'ok')
        rw.write_summary(n_total_indels, len(candidates), n_pass_qc, n_with_primer)
        # 流程元数据|pipeline metadata (CLAUDE.md §6 step6 + §12.5)
        end_time = datetime.datetime.now()
        try:
            rw.write_pipeline_info(self.config, self._gather_tool_versions(),
                                   self.start_time, end_time)
        except Exception as e:
            self.logger.warning(f"流程元数据写出失败|pipeline info write failed: {e}")

    # ---- 主流程|pipeline ----
    def run(self):
        """运行完整流程|Run full pipeline"""
        self.logger.info("开始INDEL分子标记流程|Starting INDEL marker pipeline")
        try:
            if not self._check_deps():
                raise RuntimeError("依赖检查失败|Dependency check failed")

            self._parse_samplesheet()

            indels = self._extract_indels()
            if not indels:
                self.logger.warning("无合格INDEL，流程结束|No qualified INDELs")
                self._write_results([], {}, {}, [], 0)
                return

            candidates = self._classify(indels)
            if not candidates:
                self.logger.warning("无群体判定候选，流程结束|No candidates")
                self._write_results([], {}, {}, [], len(indels))
                return

            coverage = self._compute_coverage(candidates)
            self._validate_coverage(candidates, coverage)
            sequences = self._extract_sequences(candidates)
            primers = self._design_primers(candidates, sequences)
            self._write_results(candidates, coverage, sequences, primers, len(indels))

            self.logger.info("流程完成|Pipeline completed")

        except Exception as e:
            self.logger.error(f"流程终止|Pipeline aborted: {e}")
            raise


def main():
    """命令行入口|CLI entry"""
    parser = argparse.ArgumentParser(
        description='INDEL分子标记开发（抗病/感病共显性标记）|'
                    'INDEL marker development (R/S codominant markers)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="示例|Examples: biopytools indel-marker -v v.vcf.gz -s ss.tsv -g ref.fa -o out/")
    parser.add_argument('-v', '--vcf', required=True, help='多样本合并VCF|Multi-sample VCF')
    parser.add_argument('-s', '--samplesheet', required=True,
                        help='样品分组TSV(sample_name/group/bam_path)|Samplesheet TSV')
    parser.add_argument('-g', '--genome-fasta', required=True, help='参考基因组|Reference FASTA')
    parser.add_argument('-o', '--output-dir', default='./indel_marker_output',
                        help='输出目录|Output directory')
    parser.add_argument('-t', '--threads', type=int, default=12, help='线程数|Threads')
    parser.add_argument('--min-indel-size', type=int, default=10, help='最小INDEL长度|Min INDEL size')
    parser.add_argument('--max-indel-size', type=int, default=100, help='最大INDEL长度|Max INDEL size')
    parser.add_argument('--min-quality', type=float, default=20.0,
                        help='最低QUAL过滤(缺失QUAL保留)|Min QUAL filter (missing QUAL kept)')
    parser.add_argument('--max-candidates', type=int, default=0,
                        help='候选数上限(0=不限)|Candidate cap (0=no limit)')
    parser.add_argument('--min-group-consistency', type=float, default=0.9,
                        help='组内纯合一致比例阈值|Min within-group consistency (1.0=strict)')
    parser.add_argument('--min-samples-per-group', type=int, default=1,
                        help='每组最少样品数(默认1)|Min samples per group (default 1)')
    parser.add_argument('--min-depth', type=int, default=10, help='最低覆盖度|Min depth')
    parser.add_argument('--deletion-depth-ratio', type=float, default=0.3,
                        help='deletion骤降阈值|deletion drop threshold')
    parser.add_argument('--flank-length', type=int, default=300, help='侧翼长度|Flank length')

    args = parser.parse_args()

    config = IndelMarkerConfig(
        vcf_file=args.vcf, samplesheet=args.samplesheet,
        genome_fasta=args.genome_fasta, output_dir=args.output_dir,
        threads=args.threads, min_indel_size=args.min_indel_size,
        max_indel_size=args.max_indel_size,
        min_quality=args.min_quality,
        max_candidates=args.max_candidates,
        min_group_consistency=args.min_group_consistency,
        min_samples_per_group=args.min_samples_per_group,
        min_depth=args.min_depth, deletion_depth_ratio=args.deletion_depth_ratio,
        flank_length=args.flank_length)
    try:
        config.validate()
    except ValueError as e:
        print(f"配置错误|Config error: {e}", file=sys.stderr)
        sys.exit(1)

    runner = IndelMarkerRunner(config)
    try:
        runner.run()
    except Exception as e:
        print(f"运行失败|Run failed: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
