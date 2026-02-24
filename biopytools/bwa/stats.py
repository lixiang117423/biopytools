"""
比对统计模块|Alignment Statistics Module
"""

from pathlib import Path

class AlignmentStatsGenerator:
    """比对统计生成器|Alignment Statistics Generator"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def generate_stats(self, sample_name: str, bam_file: str):
        """生成比对统计|Generate alignment statistics"""
        self.logger.info(f"生成比对统计|Generating alignment statistics: {sample_name}")
        
        # Step 1: flagstat统计|Generate flagstat
        if not self._run_flagstat(sample_name, bam_file):
            return False
        
        # Step 2: 详细统计|Generate detailed stats
        if not self._run_stats(sample_name, bam_file):
            return False
        
        self.logger.info(f"统计生成完成|Statistics generation completed: {sample_name}")
        return True
    
    def _run_flagstat(self, sample_name: str, bam_file: str) -> bool:
        """运行samtools flagstat|Run samtools flagstat"""
        output_file = self.config.stats_dir / f"{sample_name}.flagstat.txt"

        cmd = f"{self.config.samtools_path} flagstat {Path(bam_file).absolute()} > {output_file.absolute()}"

        return self.cmd_runner.run(cmd, f"生成flagstat统计|Generating flagstat: {sample_name}")

    def _run_stats(self, sample_name: str, bam_file: str) -> bool:
        """运行samtools stats|Run samtools stats"""
        output_file = self.config.stats_dir / f"{sample_name}.stats.txt"

        cmd = f"{self.config.samtools_path} stats {Path(bam_file).absolute()} > {output_file.absolute()}"

        return self.cmd_runner.run(cmd, f"生成详细统计|Generating detailed stats: {sample_name}")
    
    def generate_summary_report(self, samples):
        """生成汇总报告|Generate summary report"""
        self.logger.info("生成汇总报告|Generating summary report")
        
        report_file = self.config.output_path / "alignment_summary.txt"
        
        with open(report_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("BWA全基因组比对分析汇总报告|BWA Whole Genome Alignment Summary\n")
            f.write("=" * 80 + "\n\n")

            f.write(f"分析日期|Analysis date: {Path(report_file).stat().st_mtime}\n")
            f.write(f"参考基因组|Reference genome: {self.config.genome}\n")
            f.write(f"输入目录|Input directory: {self.config.input_dir}\n")
            f.write(f"匹配模式|Pattern: {self.config.pattern}\n")
            f.write(f"样品数量|Number of samples: {len(samples)}\n\n")

            f.write("-" * 80 + "\n")
            f.write("分析参数|Analysis Parameters\n")
            f.write("-" * 80 + "\n")
            f.write(f"线程数|Threads: {self.config.threads}\n")
            f.write(f"标记重复|Mark duplicates: {'Yes' if self.config.markdup else 'No'}\n")
            f.write(f"窗口大小|Window size: {self.config.window_size:,} bp\n")
            f.write(f"步长|Step size: {self.config.step_size:,} bp\n\n")

            f.write("-" * 80 + "\n")
            f.write("输出文件|Output Files\n")
            f.write("-" * 80 + "\n")
            f.write(f"  BAM文件目录|BAM files: {self.config.bam_dir}\n")
            f.write(f"  覆盖度文件|Coverage files: {self.config.coverage_dir}\n")
            f.write(f"  滑窗文件|Window files: {self.config.window_dir}\n")
            f.write(f"  统计文件|Stats files: {self.config.stats_dir}\n")
            f.write(f"  日志文件|Log files: {self.config.log_dir}\n\n")

            f.write("-" * 80 + "\n")
            f.write("已处理样品列表|Processed Samples\n")
            f.write("-" * 80 + "\n")
            for i, sample in enumerate(samples, 1):
                f.write(f"  {i}. {sample}\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("分析完成|Analysis Completed\n")
            f.write("=" * 80 + "\n")
        
        self.logger.info(f"汇总报告已生成|Summary report generated: {report_file}")
