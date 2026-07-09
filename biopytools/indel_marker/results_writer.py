"""结果输出模块|Results Writer

写出候选主表 / BED / 侧翼FASTA / 覆盖度矩阵 / 摘要报告。
依赖 coverage_validator.candidate_id 与 genotype_analyzer.Candidate。
"""

from pathlib import Path
from typing import List, Dict

import yaml

from .genotype_analyzer import Candidate
from .coverage_validator import candidate_id


class ResultsWriter:
    """写候选主表/BED/FASTA/覆盖度矩阵/摘要|Write master table / BED / FASTA / coverage / summary"""

    # 候选主表完整列顺序|Full column order of candidate master table
    CANDIDATE_HEADER = [
        'candidate_id', 'chrom', 'pos', 'end', 'ref', 'alt',
        'indel_type', 'indel_size', 'direction', 'present_group',
        'R_hom_alt_rate', 'S_hom_alt_rate',
        'R_hom_ref_rate', 'S_hom_ref_rate',
        'R_depth_mean', 'S_depth_mean',
        'deletion_depth_ratio', 'passes_coverage_qc', 'passes_deletion_drop',
        'left_primer', 'right_primer', 'product_size', 'tm_left', 'tm_right',
        'primer_status',
    ]

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.results_dir = Path(config.output_dir) / '06_results'
        self.sequence_dir = Path(config.output_dir) / '04_sequence'

    def write_candidates_tsv(self, candidates: List[Candidate],
                             primers: Dict[str, dict]) -> str:
        """候选主表|Candidate master table

        缺失的 Candidate 属性（如未跑覆盖度验证）用 getattr 默认值兜底，
        避免单步缺失导致整张表无法输出|Missing attrs fall back via getattr
        so a missing upstream step never blocks the whole table.
        """
        self.results_dir.mkdir(parents=True, exist_ok=True)
        out = self.results_dir / 'indel_marker.candidates.tsv'
        with open(out, 'w', encoding='utf-8') as f:
            f.write('\t'.join(self.CANDIDATE_HEADER) + '\n')
            for c in candidates:
                cid = candidate_id(c)
                p = primers.get(cid, {})
                row = [
                    cid, c.indel.chrom, c.indel.pos, c.indel.end, c.indel.ref, c.indel.alt,
                    c.indel.indel_type, c.indel.indel_size, c.direction, c.present_group,
                    f"{c.r_stats.get('hom_alt_rate', 0):.3f}",
                    f"{c.s_stats.get('hom_alt_rate', 0):.3f}",
                    f"{c.r_stats.get('hom_ref_rate', 0):.3f}",
                    f"{c.s_stats.get('hom_ref_rate', 0):.3f}",
                    f"{getattr(c, 'r_depth_mean', 0):.1f}",
                    f"{getattr(c, 's_depth_mean', 0):.1f}",
                    str(getattr(c, 'deletion_depth_ratio', 'NA')),
                    str(getattr(c, 'passes_coverage_qc', '')),
                    str(getattr(c, 'passes_deletion_drop', '')),
                    p.get('left_primer', ''), p.get('right_primer', ''),
                    p.get('product_size', ''), p.get('tm_left', ''), p.get('tm_right', ''),
                    p.get('primer_status', 'fail'),
                ]
                f.write('\t'.join(str(x) for x in row) + '\n')
        self.logger.info(
            f"候选主表写出|master table: {out}（{len(candidates)} 候选|candidates）")
        return str(out)

    def write_bed(self, candidates: List[Candidate]) -> str:
        """候选区间BED（起始0-based）|Candidate BED (0-based start)"""
        self.results_dir.mkdir(parents=True, exist_ok=True)
        out = self.results_dir / 'indel_marker.candidates.bed'
        with open(out, 'w', encoding='utf-8') as f:
            for c in candidates:
                cid = candidate_id(c)
                # BED起始=0-based pos-1|BED start is 0-based
                f.write(f"{c.indel.chrom}\t{c.indel.pos - 1}\t{c.indel.end}\t"
                        f"{cid}\t{c.indel.indel_size}\t.\n")
        self.logger.info(f"BED写出|BED: {out}（{len(candidates)} 候选|candidates）")
        return str(out)

    def write_fasta(self, sequences: Dict[str, str]) -> str:
        """侧翼序列FASTA（header=candidate_id）|Flank FASTA (header=candidate_id)"""
        self.sequence_dir.mkdir(parents=True, exist_ok=True)
        out = self.sequence_dir / 'indels.flank.fa'
        with open(out, 'w', encoding='utf-8') as f:
            for cid, seq in sequences.items():
                f.write(f">{cid}\n{seq}\n")
        self.logger.info(
            f"FASTA写出|FASTA: {out}（{len(sequences)} 序列|sequences）")
        return str(out)

    def write_coverage_tsv(self, coverage: Dict[str, Dict[str, float]],
                           sample_names: List[str]) -> str:
        """覆盖度矩阵（行=candidate_id，列=样品）|Coverage matrix (row=cid, col=sample)"""
        out_dir = Path(self.config.output_dir) / '03_coverage'
        out_dir.mkdir(parents=True, exist_ok=True)
        out = out_dir / 'indels.coverage.tsv'
        with open(out, 'w', encoding='utf-8') as f:
            f.write('\t'.join(['candidate_id'] + sample_names) + '\n')
            for cid, depths in coverage.items():
                row = [cid] + [f"{depths.get(s, 0.0):.1f}" for s in sample_names]
                f.write('\t'.join(row) + '\n')
        self.logger.info(
            f"覆盖度矩阵写出|coverage matrix: {out}（{len(coverage)} 候选|candidates）")
        return str(out)

    def write_summary(self, n_total_indels: int, n_candidates: int,
                      n_pass_qc: int, n_with_primer: int) -> str:
        """摘要报告|Summary report"""
        self.results_dir.mkdir(parents=True, exist_ok=True)
        out = self.results_dir / 'indel_marker.summary.txt'
        with open(out, 'w', encoding='utf-8') as f:
            f.write("INDEL分子标记筛选摘要|INDEL Marker Summary\n")
            f.write("=" * 50 + "\n")
            f.write(f"合格INDEL总数|Qualified INDELs: {n_total_indels}\n")
            f.write(f"群体判定候选|Candidates after genotype: {n_candidates}\n")
            f.write(f"通过覆盖度质控|Passed coverage QC: {n_pass_qc}\n")
            f.write(f"成功设计引物|With primers: {n_with_primer}\n")
        self.logger.info(f"摘要写出|summary: {out}")
        return str(out)

    def write_pipeline_info(self, config, tools_versions: Dict[str, str],
                            start_time, end_time) -> str:
        """
        写出流程元数据|Write pipeline metadata

        生成两个文件|Produces two files:
        - 00_pipeline_info/software_versions.yml (CLAUDE.md §12.5.2 格式)
        - 00_pipeline_info/pipeline_params.yaml (关键运行参数|key run params)

        Args:
            config: IndelMarkerConfig
            tools_versions: {tool_name: version_str}，缺失版本用 'unavailable'|
                            {tool_name: version_str}; missing -> 'unavailable'
            start_time: datetime.datetime，流程开始时间|pipeline start time
            end_time: datetime.datetime，流程结束时间|pipeline end time

        Returns:
            software_versions.yml 路径|path to software_versions.yml
        """
        info_dir = Path(config.output_dir) / '00_pipeline_info'
        info_dir.mkdir(parents=True, exist_ok=True)

        runtime_seconds = int((end_time - start_time).total_seconds())

        # tools 子表：CLI工具记录 version/path/command；Python包(如primer3)仅version/path留空
        # tools sub-table: CLI tools record version/path/command; Python pkgs (primer3) version only
        cli_tools = {'bcftools', 'samtools'}
        tools_table = {}
        for name, ver in tools_versions.items():
            entry = {'version': ver}
            if name in cli_tools:
                entry['path'] = getattr(config, f'{name}_path', '')
                entry['command'] = f"{name} --version"
            else:
                entry['path'] = ''
                entry['command'] = f"import {name}"
            tools_table[name] = entry

        info = {
            'pipeline': {
                'name': 'biopytools indel_marker',
                'version': '1.0.0',
            },
            'tools': tools_table,
            'parameters': {
                'vcf_file': config.vcf_file,
                'samplesheet': config.samplesheet,
                'genome_fasta': config.genome_fasta,
                'min_indel_size': config.min_indel_size,
                'max_indel_size': config.max_indel_size,
                'min_quality': config.min_quality,
                'max_candidates': config.max_candidates,
                'min_group_consistency': config.min_group_consistency,
                'min_samples_per_group': config.min_samples_per_group,
                'max_missing_rate': config.max_missing_rate,
                'min_depth': config.min_depth,
                'min_mapq': config.min_mapq,
                'min_baseq': config.min_baseq,
                'deletion_depth_ratio': config.deletion_depth_ratio,
                'flank_length': config.flank_length,
                'threads': config.threads,
            },
            'execution': {
                'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': runtime_seconds,
            },
        }
        versions_out = info_dir / 'software_versions.yml'
        with open(versions_out, 'w', encoding='utf-8') as f:
            yaml.safe_dump(info, f, default_flow_style=False, allow_unicode=True)

        params_out = info_dir / 'pipeline_params.yaml'
        with open(params_out, 'w', encoding='utf-8') as f:
            yaml.safe_dump(info['parameters'], f, default_flow_style=False,
                           allow_unicode=True)

        self.logger.info(f"流程元数据写出|pipeline info: {versions_out}")
        return str(versions_out)
