"""
TGS-GapCloser主程序模块|TGS-GapCloser Main Module
"""

import os
import sys
from pathlib import Path
from .config import TGSGapCloserConfig
from .utils import TGSGapCloserLogger, CommandRunner
from .quartet_filler import QuartetGapFiller


class TGSGapCloser:
    """TGS-GapCloser主类|TGS-GapCloser Main Class"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        # 初始化配置|Initialize configuration
        self.config = TGSGapCloserConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        output_dir = os.path.dirname(self.config.output_prefix) or '.'
        log_file = os.path.join(output_dir, 'tgsgapcloser.log')
        self.logger_manager = TGSGapCloserLogger(log_file)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, output_dir, dry_run=self.config.dry_run)

        # 输出文件路径|Output file paths
        self.output_gapclosed = f"{self.config.output_prefix}.gapcloser.fa"
        # TGS-GapCloser v1输出gapcloser.fa，v2输出scaff_seqs
        self.v2_output_file = f"{self.config.output_prefix}.scaff_seqs"

        self._log_initialization_info()

    def _log_initialization_info(self):
        """记录初始化信息|Log initialization information"""
        self.logger.info("="*80)
        self.logger.info("TGS-GapCloser Gap填充流程|TGS-GapCloser Gap Filling Pipeline")
        self.logger.info("="*80)
        self.logger.info(f"Scaffold文件|Scaffold file: {self.config.scaff_file}")
        self.logger.info(f"TGS reads文件|TGS reads file: {self.config.reads_file}")
        self.logger.info(f"输出前缀|Output prefix: {self.config.output_prefix}")
        self.logger.info(f"TGS类型|TGS type: {self.config.tgstype}")
        self.logger.info(f"纠错模式|Error correction mode: {self.config.mode}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"最小同一性|Min identity: {self.config.min_idy}")
        self.logger.info(f"最小匹配长度|Min match length: {self.config.min_match}")

        # 第2轮填充参数|Round 2 filling parameters
        if self.config.unitig_file:
            self.logger.info("="*80)
            self.logger.info("第2轮填充参数|Round 2 Filling Parameters:")
            self.logger.info(f"  Unitig文件|Unitig file: {self.config.unitig_file}")
            self.logger.info(f"  Flanking长度|Flanking length: {self.config.flanking_len} bp")
            self.logger.info(f"  最小比对长度|Min alignment length: {self.config.min_align_len} bp")
            self.logger.info(f"  最小比对同一性|Min alignment identity: {self.config.min_identity}%")
            self.logger.info(f"  最大填充长度|Max filling length: {self.config.max_filling_len} bp")

    def run(self):
        """运行Gap填充流程|Run gap filling pipeline"""
        try:
            # 提前创建 filler,供 gap 统计与前后对比报告复用|create filler early (reused for stats & report)
            filler = QuartetGapFiller(self.logger)
            # 检查断点续传|Check for resume capability
            round1_backup = f"{self.output_gapclosed}.round1"
            round2_tmp_file = f"{self.config.output_prefix}_round2_filled.fasta"

            # 如果最终输出已存在且有round1备份(且未强制重跑),说明流程已完成
            # Completed if final output + round1 backup exist (and not forced)
            if (os.path.exists(self.output_gapclosed) and os.path.exists(round1_backup)
                    and not self.config.force):
                self.logger.info("="*80)
                self.logger.info("检测到已完成的流程，跳过|Pipeline already completed, skipping")
                self.logger.info("="*80)
                self.logger.info(f"输出文件已存在|Output file exists: {self.output_gapclosed}")
                self._maybe_write_gap_report(filler)
                return True

            # 第1轮：运行TGS-GapCloser|Round 1: Run TGS-GapCloser
            self.logger.info("="*80)
            self.logger.info("开始TGS-GapCloser Gap填充流程|Starting TGS-GapCloser gap filling pipeline")
            self.logger.info("="*80)

            # 如果第1轮输出不存在(或强制重跑),运行第1轮
            # Run round 1 if output missing (or forced)
            if not os.path.exists(self.output_gapclosed) or self.config.force:
                self.logger.info("第1轮：TGS-GapCloser2 Gap填充|Round 1: TGS-GapCloser2 Gap Filling")
                self._run_tgsgapcloser()
                if self.config.dry_run:
                    self.logger.info("dry-run 模式,仅预览命令,跳过输出检查与第2轮|"
                                     "dry-run, preview only; skip output check & round 2")
                    return True
                self._check_output()
            else:
                self.logger.info("检测到第1轮已完成，跳过|Round 1 already completed, skipping")
                self.logger.info(f"第1轮输出文件|Round 1 output: {self.output_gapclosed}")

            # 检查是否还有剩余gap(复用开头创建的 filler)|check remaining gaps (reuse filler)
            has_gaps = filler.has_gaps(self.output_gapclosed)
            gap_count, _ = filler.count_gaps(self.output_gapclosed)

            self.logger.info(f"第1轮填充结果|Round 1 Results:")
            self.logger.info(f"  剩余gap数量|Remaining gaps: {gap_count}")

            # 第2轮：如果还有gap且提供了unitig文件，运行quarTeT gapfiller
            if has_gaps and self.config.unitig_file:
                # 第2轮已完成(且未强制重跑)则跳过|skip round 2 if done (and not forced)
                if os.path.exists(round1_backup) and not self.config.force:
                    self.logger.info("检测到第2轮已完成，跳过|Round 2 already completed, skipping")
                else:
                    self.logger.info("="*80)
                    self.logger.info("第2轮：quarTeT Gap Filler（使用unitig/contig）|Round 2: quarTeT Gap Filler (using unitigs/contigs)")
                    self.logger.info("="*80)

                    # 备份第1轮结果
                    os.rename(self.output_gapclosed, round1_backup)
                    self.logger.info(f"第1轮结果已备份|Round 1 result backed up to: {round1_backup}")

                    # 运行第2轮填充
                    round2_output = filler.fill_gaps(
                        draft_genome=round1_backup,
                        unitig_file=self.config.unitig_file,
                        output_prefix=f"{self.config.output_prefix}_round2",
                        flanking_len=self.config.flanking_len,
                        min_align_len=self.config.min_align_len,
                        min_identity=self.config.min_identity / 100.0,
                        max_filling_len=self.config.max_filling_len,
                        threads=self.config.threads,
                        minimap_option='-x asm20' if self.config.tgstype == 'hifi' else '-x asm5',
                        overwrite=self.config.force,
                        min_gap_length=self.config.min_gap_length,
                        minimap2_path=self.config.minimap2_path
                    )

                    # 将第2轮结果重命名为最终输出
                    os.rename(round2_output, self.output_gapclosed)
                    self.logger.info(f"第2轮填充结果|Round 2 Results:")

                    # 最终统计
                    final_gap_count, _ = filler.count_gaps(self.output_gapclosed)
                    total_filled = gap_count - final_gap_count
                    fill_rate = (total_filled / gap_count * 100) if gap_count > 0 else 0

                    self.logger.info(f"  第2轮填充gap数量|Filled in round 2: {total_filled}")
                    self.logger.info(f"  最终剩余gap|Final remaining gaps: {final_gap_count}")
                    self.logger.info(f"  第2轮填充率|Round 2 fill rate: {fill_rate:.1f}%")

            elif has_gaps and not self.config.unitig_file:
                self.logger.info("提示|Tip: 仍有剩余gap，可提供unitig文件（-ug）进行第2轮填充")
                self.logger.info("  Remaining gaps detected. You can provide unitig file (-ug) for round 2 filling")
            else:
                self.logger.info("所有gap已填充完成|All gaps filled successfully!")

            self.logger.info("="*80)
            self.logger.info("TGS-GapCloser Gap填充流程完成|TGS-GapCloser gap filling pipeline completed")
            self.logger.info(f"最终输出文件|Final output file: {self.output_gapclosed}")
            self.logger.info("="*80)

            # 清理done文件，避免影响批量处理|Clean up done files to avoid affecting batch processing
            output_dir = os.path.dirname(self.config.output_prefix) or '.'
            import glob
            done_files = glob.glob(os.path.join(output_dir, 'done_step*'))
            if done_files:
                self.logger.info(f"清理done文件|Cleaning up done files: {len(done_files)} 个|files")
                for done_file in done_files:
                    try:
                        os.remove(done_file)
                        self.logger.debug(f"  已删除|Deleted: {done_file}")
                    except Exception as e:
                        self.logger.warning(f"  删除失败|Failed to delete {done_file}: {e}")

            # 生成 gap 处理前后对比报告|generate before/after gap report
            self._maybe_write_gap_report(filler)
            return True

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {str(e)}")
            raise

    def _run_tgsgapcloser(self):
        """运行TGS-GapCloser命令|Run TGS-GapCloser command"""
        self.logger.info("步骤1: 运行TGS-GapCloser|Step 1: Running TGS-GapCloser")

        # 计算输出目录和相对前缀，用于隔离并发任务的done文件
        # Calculate output dir and relative prefix to isolate done files across concurrent jobs
        output_dir = os.path.dirname(self.config.output_prefix) or '.'
        output_basename = os.path.basename(self.config.output_prefix)

        # 清理输出目录中的done文件|Clean up done files in output directory
        import glob
        cleaned_count = 0
        output_done_files = glob.glob(os.path.join(output_dir, 'done_step*'))
        if output_done_files:
            self.logger.info(f"清理done文件|Cleaning up done files: {len(output_done_files)} 个|files")
            for done_file in output_done_files:
                try:
                    os.remove(done_file)
                    self.logger.debug(f"  已删除|Deleted: {done_file}")
                    cleaned_count += 1
                except Exception as e:
                    self.logger.warning(f"  删除失败|Failed to delete {done_file}: {e}")

        if cleaned_count > 0:
            self.logger.info(f"总共清理了|Total cleaned: {cleaned_count} 个done文件|done files")

        # 构建命令（使用绝对路径的scaff/reads，--output使用相对于输出目录的basename）
        # Build command (scaff/reads use absolute paths, --output uses basename relative to output dir)
        cmd = self.cmd_runner.build_tgsgapcloser_command(self.config)

        # 替换--output参数为仅文件名|Replace --output arg with basename only
        for i in range(len(cmd)):
            if cmd[i] == '--output' and i + 1 < len(cmd):
                cmd[i + 1] = output_basename
                break

        # 在输出目录中执行,隔离done文件;流式输出避免 capture 缓冲 OOM(§13.2.0)|
        # Execute in output dir (isolate done files); stream to avoid capture-buffer OOM
        ok = self.cmd_runner.run_with_progress(
            cmd, cwd=output_dir,
            description="TGS-GapCloser2 Gap填充|TGS-GapCloser2 gap filling")

        if not ok:
            raise RuntimeError(f"TGS-GapCloser执行失败|TGS-GapCloser execution failed")

        self.logger.info("TGS-GapCloser执行完成|TGS-GapCloser execution completed")

    def _check_output(self):
        """检查输出文件|Check output files"""
        self.logger.info("步骤2: 检查输出文件|Step 2: Checking output files")

        # TGS-GapCloser v1输出gapcloser.fa，v2输出scaff_seqs
        # 检查v1输出|Check v1 output
        if os.path.exists(self.output_gapclosed):
            self.logger.info(f"检测到TGS-GapCloser v1输出格式|Detected TGS-GapCloser v1 output format")
            actual_output = self.output_gapclosed
        # 检查v2输出|Check v2 output
        elif os.path.exists(self.v2_output_file):
            self.logger.info(f"检测到TGS-GapCloser v2输出格式|Detected TGS-GapCloser v2 output format")
            self.logger.info(f"重命名|Renaming: {self.v2_output_file} -> {self.output_gapclosed}")
            # 重命名为统一格式|Rename to unified format
            os.rename(self.v2_output_file, self.output_gapclosed)
            actual_output = self.output_gapclosed
        else:
            raise FileNotFoundError(
                f"输出文件未找到|Output file not found: "
                f"{self.output_gapclosed} 或|or {self.v2_output_file}"
            )

        # 检查文件大小|Check file size
        file_size = os.path.getsize(actual_output)
        if file_size == 0:
            raise ValueError(f"输出文件为空|Output file is empty: {actual_output}")

        self.logger.info(f"输出文件检查通过|Output file check passed: {actual_output} ({file_size} bytes)")

    def _maybe_write_gap_report(self, filler):
        """非 dry_run 且最终输出存在时,生成 gap 前后对比报告|
        Write before/after gap report when not dry-run and output exists."""
        if self.config.dry_run:
            return
        if not os.path.exists(self.output_gapclosed):
            return
        try:
            self._write_gap_report(self.config.scaff_file, self.output_gapclosed,
                                   self.config.min_gap_length, filler)
            self._write_gap_table(self.config.scaff_file, self.output_gapclosed,
                                  self.config.min_gap_length, filler)
        except Exception as e:
            self.logger.warning(f"生成Gap报告失败|Failed to write gap report: {e}")

    def _write_gap_report(self, before_file: str, after_file: str,
                          min_gap_length: int, filler) -> str:
        """生成 gap 处理前后对比报告(数量/长度分布/top/填充效果)|
        Write before/after gap comparison (counts / length dist / top / fill effect)."""
        report_file = f"{self.config.output_prefix}_gap_report.txt"
        before = filler.analyze_gaps(before_file, min_gap_length)
        after = filler.analyze_gaps(after_file, min_gap_length)

        lines = []
        lines.append("Gap填充报告|Gap Filling Report")
        lines.append("=" * 60)
        lines.append(f"处理前文件|Before (input): {before_file}")
        lines.append(f"处理后文件|After (output): {after_file}")
        lines.append(f"最小gap长度阈值|Min gap length threshold: {min_gap_length} bp")
        lines.append("")
        for title, stats in (("处理前|BEFORE", before), ("处理后|AFTER", after)):
            lines.append(f"--- {title} ---")
            lines.append(f"  总gap数量|Total gaps: {stats['total']}")
            lines.append(f"  总gap长度|Total gap length: {stats['total_length']} bp")
            lines.append(f"  含gap序列数|Sequences with gaps: {stats['seq_count']}")
            lines.append(f"  平均gap长度|Avg gap length: {stats['avg']} bp")
            lines.append(f"  最大gap长度|Max gap length: {stats['max']} bp")
            lines.append(f"  长度分布(bp)|Length distribution:")
            for _, _, label in stats['bins']:
                lines.append(f"    {label}: {stats['by_bin'].get(label, 0)}")
            lines.append(f"  最大10个gap|Top 10 largest gaps:")
            if stats['top']:
                for length, sid, start, end in stats['top']:
                    lines.append(f"    {sid}:{start}-{end} ({length} bp)")
            else:
                lines.append("    (无|none)")
            lines.append("")
        lines.append("--- 填充效果|FILLING EFFECT ---")
        filled_n = before['total'] - after['total']
        filled_bp = before['total_length'] - after['total_length']
        rate = (filled_n / before['total'] * 100) if before['total'] else 0.0
        lines.append(f"  填充gap数量|Gaps filled: {filled_n}")
        lines.append(f"  填充碱基数|Bases filled: {filled_bp} bp")
        lines.append(f"  填充率(按数量)|Fill rate (by count): {rate:.1f}%")
        lines.append(f"  剩余gap|Remaining gaps: {after['total']}")

        with open(report_file, 'w', encoding='utf-8') as f:
            f.write('\n'.join(lines) + '\n')
        self.logger.info(f"Gap报告已写|Gap report written: {report_file}")
        return report_file

    def _write_gap_table(self, before_file: str, after_file: str,
                         min_gap_length: int, filler) -> str:
        """生成 per-gap 明细 TSV(每个 gap: 处理前长度 + 填充状态 + 处理后残留长度)|
        Per-gap detail TSV (length + fill status + residual length)."""
        table_file = f"{self.config.output_prefix}_gap_table.tsv"
        rows = filler.track_gaps(before_file, after_file, min_gap_length)
        with open(table_file, 'w', encoding='utf-8') as f:
            f.write("seq\tgap_idx\tbefore_start\tbefore_end\tbefore_length\tstatus\tafter_length\n")
            for r in rows:
                f.write(f"{r['seq']}\t{r['gap_idx']}\t{r['before_start']}\t{r['before_end']}\t"
                        f"{r['before_length']}\t{r['status']}\t{r['after_length']}\n")
        filled = sum(1 for r in rows if r['status'] == 'Filled')
        remaining = sum(1 for r in rows if r['status'] == 'Remaining')
        unknown = sum(1 for r in rows if r['status'] == 'Unknown')
        self.logger.info(f"Gap明细表已写|Gap table written: {table_file} "
                         f"(填充|filled {filled}, 残留|remaining {remaining}, 未知|unknown {unknown})")
        return table_file


def main():
    """命令行入口函数|Command line entry function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='TGS-GapCloser Gap填充工具|TGS-GapCloser Gap Filling Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-s', '--scaff_file', required=True,
                       help='输入scaffold文件|Input scaffold file')
    parser.add_argument('-t', '--tgstype', required=True, choices=['ont', 'pb', 'hifi'],
                       help='TGS类型|TGS type (ont/pb/hifi)')
    parser.add_argument('-ir', '--reads_file', required=True,
                       help='输入TGS reads文件|Input TGS reads file')
    parser.add_argument('-o', '--output_prefix', required=True,
                       help='输出前缀|Output prefix')

    # 可选参数|Optional parameters
    parser.add_argument('-m', '--mode', default='none', choices=['none', 'racon', 'pilon'],
                       help='纠错模式|Error correction mode (default: none)')
    parser.add_argument('--tgsgapcloser_path',
                       default='~/software/TGS-GapCloser2/TGS-GapCloser2-master/tgsgapcloser2',
                       help='TGS-GapCloser路径|TGS-GapCloser path')
    parser.add_argument('-idy', '--min_idy', type=float,
                       help='最小同一性|Min identity (auto-set if not specified)')
    parser.add_argument('-l', '--min_match', type=int,
                       help='最小匹配长度|Min match length (auto-set if not specified)')
    parser.add_argument('-threads', '--threads', type=int, default=12,
                       help='线程数|Threads (default: 12)')
    parser.add_argument('-chunk', type=int, default=3,
                       help='分块数量|Chunk count (default: 3)')
    parser.add_argument('-g_check', action='store_true',
                       help='启用Gap大小差异检查|Enable gap size difference check')
    parser.add_argument('-min_nread', type=int, default=1,
                       help='最小reads数量|Min read count (default: 1)')
    parser.add_argument('-max_nread', type=int, default=-1,
                       help='最大reads数量|Max read count (default: -1)')
    parser.add_argument('-max_candidate', type=int, default=200,
                       help='最大候选数|Max candidates (default: 200)')

    # Racon参数|Racon parameters
    parser.add_argument('-racon', '--racon_path',
                       help='Racon路径|Racon path')
    parser.add_argument('-racon_round', type=int, default=3,
                       help='Racon轮数|Racon rounds (default: 3)')

    # Pilon参数|Pilon parameters
    parser.add_argument('-pilon', '--pilon_path',
                       help='Pilon路径|Pilon path')
    parser.add_argument('-ngs', '--ngs_file',
                       help='NGS reads文件|NGS reads file')
    parser.add_argument('-java', '--java_path',
                       help='Java路径|Java path')
    parser.add_argument('-samtools', '--samtools_path',
                       help='Samtools路径|Samtools path')
    parser.add_argument('-pilon_mem', default='300G',
                       help='Pilon内存|Pilon memory (default: 300G)')
    parser.add_argument('-pilon_round', type=int, default=3,
                       help='Pilon轮数|Pilon rounds (default: 3)')

    # 其他参数|Other parameters
    parser.add_argument('-minmap_arg',
                       help='自定义minimap2参数|Custom minimap2 arguments')

    # quarTeT GapFiller参数（第2轮填充）|quarTeT GapFiller parameters (round 2)
    parser.add_argument('-unitig_file',
                       help='hifiasm unitig/contig文件（第2轮填充）|hifiasm unitig/contig file (round 2)')
    parser.add_argument('-flanking_len', type=int, default=5000,
                       help='Flanking序列长度（bp）|Flanking sequence length (default: 5000)')
    parser.add_argument('-min_align_len', type=int, default=1000,
                       help='最小比对长度（bp）|Min alignment length (default: 1000)')
    parser.add_argument('-min_identity', type=int, default=40,
                       help='最小比对同一性（%%）|Min alignment identity %% (default: 40)')
    parser.add_argument('-max_filling_len', type=int, default=1000000,
                       help='最大填充长度（bp）|Max filling length (default: 1000000)')
    parser.add_argument('-min_gap_length', type=int, default=100,
                       help='第2轮识别/填充的最小gap长度(bp)|Min gap length (bp) for round-2 (default: 100)')

    # 流程控制|Pipeline control
    parser.add_argument('-f', '--force', action='store_true',
                       help='忽略断点续传强制重跑|Force rerun, ignore checkpoint')
    parser.add_argument('--dry-run', action='store_true',
                       help='只打印命令不执行|Dry run, print commands only')

    args = parser.parse_args()

    # 创建并运行TGSGapCloser|Create and run TGSGapCloser
    try:
        gapcloser = TGSGapCloser(**vars(args))
        success = gapcloser.run()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"错误|Error: {e}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
