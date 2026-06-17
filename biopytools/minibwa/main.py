"""
Minibwa主程序模块|Minibwa Main Program Module

编排 index → align → stats → coverage 完整流程。
|Orchestrates the full index → align → stats → coverage pipeline.
"""

import argparse
import subprocess
import sys
from datetime import datetime
from pathlib import Path

from .config import MinibwaConfig, VALID_MODES, MODE_STANDARD
from .utils import (
    MinibwaLogger, CommandRunner, check_dependencies,
    find_fastq_pairs, build_conda_command,
)
from .indexer import GenomeIndexer
from .aligner import MinibwaAligner
from .postprocess import AlignmentStatsGenerator, CoverageAnalyzer


class MinibwaRunner:
    """Minibwa运行器|Minibwa Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = MinibwaConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = MinibwaLogger(self.config.log_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger)

        # 初始化各组件|Initialize components
        self.indexer = GenomeIndexer(self.config, self.logger, self.cmd_runner)
        self.aligner = None   # 延迟到索引构建后（需要index_prefix）|Delayed until index built
        self.stats_generator = AlignmentStatsGenerator(
            self.config, self.logger, self.cmd_runner
        )
        self.coverage_analyzer = CoverageAnalyzer(
            self.config, self.logger, self.cmd_runner
        )

        # 流程起始时间|Pipeline start time
        self.start_time = datetime.now()

    def _is_step_completed(self, output_file: str) -> bool:
        """断点续传检查|Checkpoint completion check"""
        return Path(output_file).exists()

    def run_analysis(self) -> bool:
        """运行完整比对分析流程|Run complete alignment pipeline"""
        self.logger.info("=" * 80)
        self.logger.info(
            f"开始Minibwa比对分析流程|Starting Minibwa pipeline "
            f"(mode={self.config.mode})"
        )
        self.logger.info("=" * 80)

        try:
            # Step 1: 检查依赖|Check dependencies
            self.logger.info("步骤1/5: 检查依赖软件|Step 1/5: Checking dependencies")
            check_dependencies(self.config, self.logger)

            # Step 2: 检查/构建索引|Check/build genome index
            self.logger.info("步骤2/5: 检查基因组索引|Step 2/5: Checking genome index")
            if not self.indexer.check_and_build_index():
                self.logger.error("基因组索引构建失败|Genome index build failed")
                return False
            self.aligner = MinibwaAligner(
                self.config, self.logger, self.cmd_runner,
                index_prefix=self.indexer.get_index_prefix()
            )

            # Step 3: 查找FASTQ文件对|Find FASTQ pairs
            self.logger.info("步骤3/5: 查找FASTQ文件|Step 3/5: Finding FASTQ files")
            fastq_pairs = find_fastq_pairs(
                self.config.input_dir, self.config.pattern, self.logger
            )
            if not fastq_pairs:
                self.logger.error("未找到任何FASTQ文件对|No FASTQ pairs found")
                return False

            # Step 4: 逐样品比对+后处理|Per-sample align + post-processing
            total = len(fastq_pairs)
            self.logger.info(
                f"步骤4/5: 批量比对与统计|Step 4/5: Batch alignment & stats "
                f"({total} samples)"
            )

            processed = []
            for i, (sample_name, r1, r2) in enumerate(fastq_pairs, 1):
                self.logger.info("-" * 80)
                self.logger.info(
                    f"处理样品|Processing sample {i}/{total}: {sample_name}"
                )

                # 比对|Alignment
                bam_file = self.aligner.align_sample(sample_name, r1, r2)
                if not bam_file:
                    self.logger.error(
                        f"样品比对失败|Sample alignment failed: {sample_name}"
                    )
                    continue

                # 统计|Stats
                if not self.stats_generator.generate_stats(sample_name, bam_file):
                    self.logger.warning(
                        f"统计生成失败，继续|Stats failed, continuing: {sample_name}"
                    )

                # 覆盖度|Coverage
                if not self.config.skip_coverage:
                    if not self.coverage_analyzer.analyze_sample_coverage(
                        sample_name, bam_file
                    ):
                        self.logger.warning(
                            f"覆盖度分析失败，继续|Coverage failed, continuing: {sample_name}"
                        )

                processed.append(sample_name)

            # Step 5: 汇总|Summary
            self.logger.info("步骤5/5: 生成汇总报告|Step 5/5: Generating summary")
            self.stats_generator.generate_summary_report(processed)

            # 元数据|Metadata
            self._write_pipeline_info(processed)

            self.logger.info("=" * 80)
            self.logger.info(
                f"Minibwa分析流程完成|Minibwa pipeline completed: "
                f"{len(processed)}/{total} samples succeeded"
            )
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info("=" * 80)
            return True

        except Exception as e:
            self.logger.error(f"流程异常|Pipeline error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _write_pipeline_info(self, processed_samples: list):
        """生成 software_versions.yml 和 pipeline_params.yaml|Generate metadata"""
        self.logger.info("写入流程元数据|Writing pipeline metadata")
        self._write_software_versions()
        self._write_pipeline_params(processed_samples)

    def _write_software_versions(self):
        """生成 software_versions.yml|Generate software_versions.yml"""
        end_time = datetime.now()
        runtime = int((end_time - self.start_time).total_seconds())

        versions = {}
        for tool_name, tool_path, version_args in [
            ('minibwa', self.config.minibwa_path, ['version']),
            ('samtools', self.config.samtools_path, ['--version']),
        ]:
            try:
                cmd = build_conda_command(tool_path, version_args)
                result = subprocess.run(
                    cmd, capture_output=True, text=True, timeout=30
                )
                first_line = (result.stdout or result.stderr).strip().split('\n')[0]
                versions[tool_name] = {
                    'version': first_line,
                    'path': tool_path,
                }
            except Exception as e:
                versions[tool_name] = {
                    'version': 'unknown',
                    'path': tool_path,
                    'error': str(e),
                }

        info = {
            'pipeline': {
                'name': 'biopytools.minibwa',
                'version': '1.0.0',
            },
            'tools': versions,
            'parameters': self.config.to_param_dict(),
            'execution': {
                'start_time': self.start_time.strftime('%Y-%m-%d %H:%M:%S'),
                'end_time': end_time.strftime('%Y-%m-%d %H:%M:%S'),
                'runtime_seconds': runtime,
            },
        }

        # 简单YAML输出（避免依赖PyYAML）|Minimal YAML without PyYAML dependency
        out_file = self.config.pipeline_info_dir / "software_versions.yml"
        with open(out_file, 'w') as f:
            self._dump_yaml(info, f, 0)

        self.logger.info(f"软件版本|Software versions: {out_file}")

    def _dump_yaml(self, data, f, indent: int):
        """简易YAML序列化|Minimal YAML serializer"""
        pad = '  ' * indent
        if isinstance(data, dict):
            for k, v in data.items():
                if isinstance(v, (dict, list)):
                    f.write(f"{pad}{k}:\n")
                    self._dump_yaml(v, f, indent + 1)
                else:
                    f.write(f"{pad}{k}: {self._yaml_scalar(v)}\n")
        elif isinstance(data, list):
            for item in data:
                if isinstance(item, (dict, list)):
                    f.write(f"{pad}-\n")
                    self._dump_yaml(item, f, indent + 1)
                else:
                    f.write(f"{pad}- {self._yaml_scalar(item)}\n")

    @staticmethod
    def _yaml_scalar(v) -> str:
        """简单YAML标量格式化|Simple YAML scalar formatter"""
        if v is None:
            return 'null'
        if isinstance(v, bool):
            return 'true' if v else 'false'
        if isinstance(v, (int, float)):
            return str(v)
        s = str(v)
        # 含特殊字符则加引号|Quote if contains special chars
        if any(c in s for c in ':#{}[]&*!|>"\'@\n'):
            return f'"{s}"'
        return s

    def _write_pipeline_params(self, processed_samples: list):
        """生成 pipeline_params.yaml|Generate pipeline_params.yaml"""
        out_file = self.config.pipeline_info_dir / "pipeline_params.yaml"
        params = self.config.to_param_dict()
        params['processed_samples'] = processed_samples
        params['skipped_coverage'] = self.config.skip_coverage

        with open(out_file, 'w') as f:
            self._dump_yaml(params, f, 0)

        self.logger.info(f"流程参数|Pipeline params: {out_file}")


def main():
    """命令行入口|Command-line entry"""
    parser = argparse.ArgumentParser(
        description='Minibwa短读长比对分析工具|Minibwa Short-read Alignment Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required')
    required.add_argument('-g', '--genome', required=True,
                          help='参考基因组FASTA|Reference genome FASTA')
    required.add_argument('-i', '--input', required=True,
                          help='FASTQ输入目录|FASTQ input directory')
    required.add_argument('-p', '--pattern', default='_1.fq.gz',
                          help='R1匹配模式|R1 pattern')

    # 输出参数|Output
    out_g = parser.add_argument_group('输出|Output')
    out_g.add_argument('-o', '--output-dir', default='./minibwa_output',
                       help='输出目录|Output directory')

    # 模式参数|Mode
    mode_g = parser.add_argument_group('比对模式|Alignment mode')
    mode_g.add_argument('--mode', choices=list(VALID_MODES), default=MODE_STANDARD,
                        help='比对模式|Alignment mode: '
                             'standard=短读PE/SE, hic=Hi-C, meth=BS-seq, long=长读')

    # 性能参数|Performance
    perf_g = parser.add_argument_group('性能|Performance')
    perf_g.add_argument('-t', '--threads', type=int, default=12,
                        help='线程数|Number of threads')

    # Minibwa map 参数|Minibwa map parameters
    map_g = parser.add_argument_group('Minibwa map 参数|Minibwa map parameters')
    map_g.add_argument('--preset', default='adap',
                       help='-x 预设: sr/lr/adap|preset')
    map_g.add_argument('-k', '--min-seed', type=int, default=19,
                       help='-k 最小种子长度|min seed length')
    map_g.add_argument('-c', '--max-occ', type=int, default=250,
                       help='-c 最大种子出现次数|max seed occurrences')
    map_g.add_argument('--max-gap', type=int, default=100,
                       help='minibwa -g 最大gap|max gap size')
    map_g.add_argument('-w', '--bandwidth', type=int, default=100,
                       help='-w 带宽|bandwidth')
    map_g.add_argument('-W', '--bandwidth-long', type=int, default=20000,
                       help='-W 长读带宽|long bandwidth')
    map_g.add_argument('-m', '--min-chain-score', type=int, default=25,
                       help='minibwa -m 最小链接分数|min chaining score')
    map_g.add_argument('--sec-ratio', type=float, default=0.5,
                       help='minibwa -p 次要/主要得分比|secondary-to-primary ratio')
    map_g.add_argument('-N', '--max-sec', type=int, default=50,
                       help='-N 保留次要比对数|retain N secondary alignments')
    map_g.add_argument('-s', '--min-dp-score', type=int, default=30,
                       help='-s 最小DP得分*{-A}|min DP score')

    # 打分参数|Scoring
    score_g = parser.add_argument_group('打分|Scoring')
    score_g.add_argument('-A', '--match-score', type=int, default=2,
                         help='-A 匹配得分|matching score')
    score_g.add_argument('-B', '--mismatch-penalty', type=int, default=8,
                         help='-B 错配罚分|mismatch penalty')
    score_g.add_argument('-O', '--gap-open', default='12,23',
                         help='-O gap开放罚分|gap open penalty')
    score_g.add_argument('-E', '--gap-ext', default='2,1',
                         help='-E gap延伸罚分|gap extension penalty')

    # IO参数|I/O
    io_g = parser.add_argument_group('IO选项|I/O options')
    io_g.add_argument('-R', '--read-group', default=None,
                      help='-R SAM read group line|read group line')
    io_g.add_argument('-u', '--no-unmap', action='store_true',
                      help='-u 不输出未比对read|do not output unmapped')
    io_g.add_argument('--outn', type=int, default=0,
                      help='--outn 输出次要比对上限|max secondary to output')
    io_g.add_argument('-y', '--copy-comment', action='store_true',
                      help='-y 复制FASTA/Q注释|copy comments')
    io_g.add_argument('-Y', '--soft-clip-supp', action='store_true',
                      help='-Y 软剪切补充比对|soft clip supplementary')
    io_g.add_argument('-K', '--batch-size', default='100m,1g',
                      help='-K 批处理大小|batch size')

    # 后处理|Post-processing
    post_g = parser.add_argument_group('后处理|Post-processing')
    post_g.add_argument('--markdup', action='store_true',
                        help='标记重复|Mark duplicates')
    post_g.add_argument('--remove-dup', action='store_true',
                        help='移除重复（需--markdup）|Remove duplicates (requires --markdup)')

    # 覆盖度|Coverage
    cov_g = parser.add_argument_group('覆盖度|Coverage')
    cov_g.add_argument('--skip-coverage', action='store_true',
                       help='跳过覆盖度分析|Skip coverage analysis')
    cov_g.add_argument('--min-base-quality', type=int, default=0,
                       help='最小碱基质量|Min base quality')
    cov_g.add_argument('--min-mapping-quality', type=int, default=0,
                       help='最小比对质量|Min mapping quality')
    cov_g.add_argument('--max-depth', type=int, default=0,
                       help='最大深度限制(0=无限)|Max depth (0=unlimited)')
    cov_g.add_argument('--window-size', type=int, default=1000000,
                       help='窗口大小|Window size bp')
    cov_g.add_argument('--step-size', type=int, default=100000,
                       help='步长|Step size bp')

    # 工具路径|Tool paths
    tool_g = parser.add_argument_group('工具路径|Tool paths')
    tool_g.add_argument('--minibwa-path', default='~/software/minibwa/minibwa',
                        help='minibwa二进制路径|minibwa binary path')
    tool_g.add_argument('--samtools-path', default='~/.local/bin/samtools',
                        help='samtools二进制路径|samtools binary path')

    # 运行控制|Run control
    rc_g = parser.add_argument_group('运行控制|Run control')
    rc_g.add_argument('--resume', action='store_true',
                      help='断点续传（跳过已完成样品）|Resume (skip completed samples)')

    args = parser.parse_args()

    try:
        runner = MinibwaRunner(
            genome=args.genome,
            input_dir=args.input,
            pattern=args.pattern,
            output_dir=args.output_dir,
            mode=args.mode,
            threads=args.threads,
            preset=args.preset,
            min_seed=args.min_seed,
            max_occ=args.max_occ,
            max_gap=args.max_gap,
            bandwidth=args.bandwidth,
            bandwidth_long=args.bandwidth_long,
            min_chain_score=args.min_chain_score,
            sec_ratio=args.sec_ratio,
            max_sec=args.max_sec,
            min_dp_score=args.min_dp_score,
            match_score=args.match_score,
            mismatch_penalty=args.mismatch_penalty,
            gap_open=args.gap_open,
            gap_ext=args.gap_ext,
            read_group=args.read_group,
            no_unmap=args.no_unmap,
            out_n=args.outn,
            copy_comment=args.copy_comment,
            soft_clip_supp=args.soft_clip_supp,
            batch_size=args.batch_size,
            markdup=args.markdup,
            remove_dup=args.remove_dup,
            skip_coverage=args.skip_coverage,
            min_base_quality=args.min_base_quality,
            min_mapping_quality=args.min_mapping_quality,
            max_depth=args.max_depth,
            window_size=args.window_size,
            step_size=args.step_size,
            minibwa_path=args.minibwa_path,
            samtools_path=args.samtools_path,
            resume=args.resume,
        )
        success = runner.run_analysis()
        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
