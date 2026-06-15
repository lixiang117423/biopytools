"""
GFFcompare两两比较主程序模块|GFFcompare Pairwise Comparison Main Module
"""

import argparse
import os
import sys
from datetime import datetime
from pathlib import Path
from typing import List, Tuple

import yaml

from .config import GffCompareConfig
from .utils import GffCompareLogger, CommandRunner, build_conda_command, check_dependencies


class GffComparePairwise:
    """GFFcompare两两双向比较分析|GFFcompare Pairwise Bidirectional Comparison"""

    GFFCMP_OUTPUTS = ('.stats', '.tracking', '.annotated.gtf', '.tmap', '.refmap', '.loci')

    def __init__(self, **kwargs):
        self.config = GffCompareConfig(**kwargs)
        self.config.validate()

        self.logger_manager = GffCompareLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

    def generate_pairs(self) -> List[Tuple[str, str, str]]:
        """
        生成所有两两双向比较对|Generate all pairwise bidirectional comparison pairs

        Returns:
            [(query_file, ref_file, pair_name), ...]
        """
        pairs = []
        files = self.config.resolved_files
        names = self.config.sample_names

        for i in range(len(files)):
            for j in range(len(files)):
                if i == j:
                    continue
                pair_name = f"{names[i]}_vs_{names[j]}"
                pairs.append((files[i], files[j], pair_name))

        return pairs

    def _is_pair_completed(self, pair_dir: str) -> bool:
        """检查比较对是否已完成|Check if a comparison pair is completed"""
        annotated = os.path.join(pair_dir, "gffcmp.annotated.gtf")
        return os.path.exists(annotated) and os.path.getsize(annotated) > 0

    def _build_gffcompare_args(self, query_file: str, ref_file: str, prefix: str) -> List[str]:
        """构建gffcompare命令参数|Build gffcompare command arguments"""
        args = ["-r", ref_file, "-o", prefix]

        if self.config.exon_range is not None:
            args.extend(["-e", str(self.config.exon_range)])
        if self.config.tss_distance is not None:
            args.extend(["-d", str(self.config.tss_distance)])
        if self.config.discard_single_exon_query:
            args.append("-M")
        if self.config.discard_single_exon_ref:
            args.append("-N")
        if self.config.ref_overlap_only:
            args.append("-R")
        if self.config.query_overlap_only:
            args.append("-Q")
        if self.config.no_tmap_refmap:
            args.append("-T")
        if self.config.strict_match:
            args.append("--strict-match")
        if self.config.cds_match:
            args.append("--cds-match")
        if self.config.genome_seq:
            args.extend(["-s", self.config.genome_seq])
        if self.config.cprefix:
            args.extend(["-p", self.config.cprefix])
        if self.config.verbose_mode:
            args.append("-V")

        args.append(query_file)
        return args

    def run_pair(self, query_file: str, ref_file: str, pair_name: str) -> bool:
        """
        运行单次gffcompare比较|Run a single gffcompare comparison

        Args:
            query_file: 查询文件路径|Query file path
            ref_file: 参考文件路径|Reference file path
            pair_name: 比较对名称|Pair name (e.g., "sampleA_vs_sampleB")

        Returns:
            是否成功|Whether succeeded
        """
        pair_dir = os.path.join(self.config.output_dir, "01_gffcompare", pair_name)
        os.makedirs(pair_dir, exist_ok=True)

        # 断点续传检查|Checkpoint resume check
        if not self.config.force and self._is_pair_completed(pair_dir):
            self.logger.info(f"跳过已完成比较|Skipping completed comparison: {pair_name}")
            return True

        prefix = os.path.join(pair_dir, "gffcmp")
        args = self._build_gffcompare_args(query_file, ref_file, prefix)
        cmd = build_conda_command(self.config.gffcompare_path, args)

        description = f"{pair_name} ({Path(query_file).name} → {Path(ref_file).name})"
        success, stdout, stderr = self.cmd_runner.run(cmd, description)

        if not success:
            self.logger.error(f"比较失败|Comparison failed: {pair_name}")
            return False

        # 验证输出文件|Verify output files
        expected_files = ['gffcmp.stats', 'gffcmp.tracking', 'gffcmp.annotated.gtf']
        if not self.config.no_tmap_refmap:
            expected_files.extend(['gffcmp.tmap', 'gffcmp.refmap'])

        for fname in expected_files:
            fpath = os.path.join(pair_dir, fname)
            if not os.path.exists(fpath):
                self.logger.warning(f"预期输出文件未生成|Expected output file not generated: {fname}")
            else:
                self.logger.debug(f"输出文件已生成|Output file generated: {fname}")

        return True

    def generate_summary(self):
        """
        合并所有.stats文件到汇总TSV|Merge all .stats files into a summary TSV

        写入 02_summary/all_stats.tsv，每行前缀为 pair_name, query, reference
        """
        summary_dir = self.config.output_path / "02_summary"
        summary_file = summary_dir / "all_stats.tsv"

        pairs = self.generate_pairs()
        all_lines = []
        header_written = False

        for _, _, pair_name in pairs:
            stats_file = os.path.join(
                self.config.output_dir, "01_gffcompare", pair_name, "gffcmp.stats"
            )

            if not os.path.exists(stats_file):
                continue

            with open(stats_file, 'r') as f:
                for line in f:
                    line = line.rstrip('\n')
                    if not line.strip():
                        continue
                    if line.startswith('#'):
                        continue
                    # 拆分 pair_name 获取 query 和 reference
                    parts = pair_name.split('_vs_', 1)
                    query_name = parts[0] if len(parts) == 2 else ''
                    ref_name = parts[1] if len(parts) == 2 else ''

                    if not header_written:
                        all_lines.append(f"pair_name\tquery\treference\t{line}")
                        header_written = True
                    else:
                        all_lines.append(f"{pair_name}\t{query_name}\t{ref_name}\t{line}")

        if all_lines:
            with open(summary_file, 'w') as f:
                f.write('\n'.join(all_lines) + '\n')
            self.logger.info(f"汇总统计已生成|Summary statistics generated: {summary_file}")
        else:
            self.logger.warning("未找到.stats文件，跳过汇总|No .stats files found, skipping summary")

    def _generate_software_versions_yml(self, start_time: datetime):
        """生成software_versions.yml|Generate software_versions.yml"""
        info_dir = self.config.output_path / "00_pipeline_info"
        info_dir.mkdir(parents=True, exist_ok=True)
        output_file = info_dir / "software_versions.yml"

        # 获取gffcompare版本|Get gffcompare version
        try:
            cmd = build_conda_command(self.config.gffcompare_path, ['--version'])
            import subprocess as sp
            result = sp.run(cmd, capture_output=True, text=True, timeout=10)
            version = result.stdout.strip()
        except Exception:
            version = "unknown"

        end_time = datetime.now()
        runtime = int((end_time - start_time).total_seconds())

        # 收集gffcompare参数|Collect gffcompare parameters
        gc_params = {}
        param_map = {
            'exon_range': '-e', 'tss_distance': '-d',
            'discard_single_exon_query': '-M', 'discard_single_exon_ref': '-N',
            'ref_overlap_only': '-R', 'query_overlap_only': '-Q',
            'no_tmap_refmap': '-T', 'strict_match': '--strict-match',
            'cds_match': '--cds-match', 'verbose_mode': '-V',
        }
        for attr, flag in param_map.items():
            val = getattr(self.config, attr)
            if val:
                gc_params[flag] = val

        if self.config.genome_seq:
            gc_params['-s'] = self.config.genome_seq
        if self.config.cprefix:
            gc_params['-p'] = self.config.cprefix

        info = {
            'pipeline': {
                'name': 'biopytools gffcompare',
                'version': '1.0.0'
            },
            'tools': {
                'gffcompare': {
                    'version': version,
                    'path': self.config.gffcompare_path
                }
            },
            'parameters': {
                'input_files': self.config.sample_names,
                'total_pairs': len(self.generate_pairs()),
                'gffcompare_options': gc_params,
            },
            'execution': {
                'start_time': start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': runtime
            }
        }

        with open(output_file, 'w') as f:
            yaml.dump(info, f, default_flow_style=False, allow_unicode=True)

        self.logger.info(f"版本信息已生成|Software versions generated: {output_file}")

    def run_all(self):
        """运行完整的两两比较流程|Run complete pairwise comparison pipeline"""
        start_time = datetime.now()

        try:
            self.logger.info("=" * 60)
            self.logger.info("GFFcompare 两两双向比较流程|GFFcompare Pairwise Bidirectional Comparison")
            self.logger.info("=" * 60)
            self.logger.info(f"输入文件数|Input file count: {len(self.config.resolved_files)}")
            self.logger.info(f"样本名称|Sample names: {', '.join(self.config.sample_names)}")
            pairs = self.generate_pairs()
            self.logger.info(f"比较对数|Total comparison pairs: {len(pairs)}")

            # Step 1: 检查依赖|Check dependencies
            check_dependencies(self.config, self.logger)

            # Step 2: 两两比较|Pairwise comparisons
            success_count = 0
            fail_count = 0
            skip_count = 0

            for idx, (query, ref, pair_name) in enumerate(pairs, 1):
                self.logger.info("-" * 40)
                self.logger.info(f"比较进度|Progress: [{idx}/{len(pairs)}] {pair_name}")

                if not self.config.force and self._is_pair_completed(
                    os.path.join(self.config.output_dir, "01_gffcompare", pair_name)
                ):
                    skip_count += 1
                    self.logger.info(f"跳过已完成|Skipping completed: {pair_name}")
                    success_count += 1
                    continue

                if self.run_pair(query, ref, pair_name):
                    success_count += 1
                else:
                    fail_count += 1

            # Step 3: 汇总统计|Generate summary
            self.logger.info("-" * 40)
            self.logger.info("生成汇总统计|Generating summary statistics")
            self.generate_summary()

            # Step 4: 版本信息|Software versions
            self._generate_software_versions_yml(start_time)

            # 完成报告|Completion report
            self.logger.info("=" * 60)
            self.logger.info(f"比较流程完成|Comparison pipeline completed")
            self.logger.info(f"成功|Success: {success_count}, 失败|Failed: {fail_count}, 跳过|Skipped: {skip_count}")
            self.logger.info(f"结果目录|Results directory: {self.config.output_dir}")
            self.logger.info("=" * 60)

            if fail_count > 0:
                self.logger.warning(f"有{fail_count}个比较对失败|{fail_count} comparison pair(s) failed")
                sys.exit(1)

        except Exception as e:
            self.logger.error(f"分析流程异常终止|Pipeline terminated unexpectedly: {e}")
            sys.exit(1)


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='GFF/GTF文件两两双向比较分析工具|GFF/GTF Pairwise Bidirectional Comparison Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', nargs='+', required=True,
                        help='输入GFF/GTF文件或文件夹，支持混用(自动识别.gff/.gtf/.gff3)|'
                             'Input GFF/GTF file(s) or directory, supports mixing (auto-detect .gff/.gtf/.gff3)')

    # 输出参数|Output arguments
    parser.add_argument('-o', '--output-dir', default='./gffcompare_output',
                        help='输出目录|Output directory')

    # GFFcompare选项|GFFcompare options
    gc_group = parser.add_argument_group('GFFcompare选项|GFFcompare options')
    gc_group.add_argument('-e', '--exon-range', type=int, default=None,
                          help='端部外显子最大允许变异范围(默认100)|Max terminal exon range (default 100)')
    gc_group.add_argument('-d', '--tss-distance', type=int, default=None,
                          help='转录本起始位点分组距离(默认100)|TSS grouping distance (default 100)')
    gc_group.add_argument('-M', '--discard-single-exon-query', action='store_true',
                          help='丢弃单外显子query转录本|Discard single-exon query transcripts')
    gc_group.add_argument('-N', '--discard-single-exon-ref', action='store_true',
                          help='丢弃单外显子reference转录本|Discard single-exon reference transcripts')
    gc_group.add_argument('-R', '--ref-overlap-only', action='store_true',
                          help='仅考虑与query重叠的reference|Only consider reference overlapping query')
    gc_group.add_argument('-Q', '--query-overlap-only', action='store_true',
                          help='仅考虑与reference重叠的query|Only consider query overlapping reference')
    gc_group.add_argument('-T', '--no-tmap-refmap', action='store_true',
                          help='不生成.tmap和.refmap文件|Skip .tmap and .refmap files')
    gc_group.add_argument('--strict-match', action='store_true',
                          help='严格匹配模式(考虑端部外显子范围)|Strict match mode (consider -e range)')
    gc_group.add_argument('--cds-match', action='store_true',
                          help='启用CDS链匹配验证|Enable CDS chain matching validation')
    gc_group.add_argument('-s', '--genome-seq', default=None,
                          help='基因组序列路径(FASTA)|Genome sequence path (FASTA)')
    gc_group.add_argument('-p', '--cprefix', default=None,
                          help='合并GTF中转录本前缀(默认TCONS)|Transcript prefix in combined GTF (default TCONS)')
    gc_group.add_argument('-V', '--verbose-mode', action='store_true',
                          help='详细处理模式|Verbose processing mode')

    # 运行参数|Runtime arguments
    parser.add_argument('--force', '-f', action='store_true',
                        help='强制重新运行(覆盖断点续传)|Force re-run (override checkpoint resume)')
    parser.add_argument('--gffcompare-path', default=None,
                        help='gffcompare软件路径|gffcompare software path')

    args = parser.parse_args()

    # 构建config参数|Build config kwargs
    config_kwargs = {
        'input_files': args.input,
        'output_dir': args.output_dir,
        'exon_range': args.exon_range,
        'tss_distance': args.tss_distance,
        'discard_single_exon_query': args.discard_single_exon_query,
        'discard_single_exon_ref': args.discard_single_exon_ref,
        'ref_overlap_only': args.ref_overlap_only,
        'query_overlap_only': args.query_overlap_only,
        'no_tmap_refmap': args.no_tmap_refmap,
        'strict_match': args.strict_match,
        'cds_match': args.cds_match,
        'genome_seq': args.genome_seq,
        'cprefix': args.cprefix,
        'verbose_mode': args.verbose_mode,
        'force': args.force,
    }

    if args.gffcompare_path:
        config_kwargs['gffcompare_path'] = args.gffcompare_path

    analyzer = GffComparePairwise(**config_kwargs)
    analyzer.run_all()


if __name__ == "__main__":
    main()
