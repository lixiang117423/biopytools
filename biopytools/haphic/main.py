"""
HapHiC主程序模块 | HapHiC Main Module
"""

import argparse
import subprocess
import sys
import time
import os
from pathlib import Path
from typing import Optional

from .config import HapHiCConfig
from .utils import HapHiCLogger, QualityController, ResultValidator, BWAAligner
from .pipeline import HapHiCPipeline


class HapHiCProcessor:
    """HapHiC处理器主类 | Main HapHiC Processor Class"""

    def __init__(self, **kwargs):
        """初始化HapHiC处理器 | Initialize HapHiC processor"""
        # 初始化配置 | Initialize configuration
        self.config = HapHiCConfig(**kwargs)
        self.config.validate()

        # 初始化日志 | Initialize logging
        self.logger_manager = HapHiCLogger(
            log_file=self.config.log_file,
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化组件 | Initialize components
        self.quality_controller = QualityController(self.logger, self.config)
        self.result_validator = ResultValidator(self.logger, self.config)

        # 初始化BWA比对器 | Initialize BWA aligner
        if self.config.hic_file_type == "fastq":
            self.bwa_aligner = BWAAligner(self.config, self.logger)
        else:
            self.bwa_aligner = None

        self.logger.info("🧬 HapHiC处理器已初始化 | HapHiC processor initialized")
        self._log_configuration()

    def _log_configuration(self):
        """记录配置信息 | Log configuration information"""
        self.logger.info("⚙️ 配置参数 | Configuration Parameters:")

        summary = self.config.get_summary()
        for key, value in summary.items():
            self.logger.info(f"  {key}: {value}")

    def run_pipeline(self) -> bool:
        """运行完整流程 | Run complete pipeline"""
        try:
            self.logger.info("🚀 开始HapHiC流程 | Starting HapHiC pipeline")
            start_time = time.time()

            # 预检查 | Pre-flight checks
            if not self._pre_flight_checks():
                return False

            # BWA比对 (如果需要) | BWA alignment (if needed)
            if self.config.hic_file_type == "fastq":
                # 检查是否已存在比对结果
                alignment_dir = os.path.join(self.config.output_dir, "00.mapping")
                potential_bam = os.path.join(alignment_dir, "HiC.bam")
                potential_filtered_bam = os.path.join(alignment_dir, "HiC.filtered.bam")

                # 添加调试信息
                self.logger.info(f"🔍 检查已存在的比对结果 | Checking for existing alignment results")
                self.logger.info(f"📂 输出目录: {self.config.output_dir}")
                self.logger.info(f"📂 比对目录: {alignment_dir}")
                self.logger.info(f"📄 潜在BAM文件: {potential_bam}")
                self.logger.info(f"📄 潜在过滤BAM文件: {potential_filtered_bam}")

                # 检查目录是否存在
                if os.path.exists(alignment_dir):
                    self.logger.info(f"✅ 比对目录存在 | Alignment directory exists")
                    # 列出目录内容
                    try:
                        dir_contents = os.listdir(alignment_dir)
                        self.logger.info(f"📋 比对目录内容: {dir_contents}")
                    except Exception as e:
                        self.logger.warning(f"⚠️ 无法列出比对目录内容: {e}")
                else:
                    self.logger.info(f"❌ 比对目录不存在 | Alignment directory does not exist")

                # 更宽松的检查：即使没有索引文件，如果BAM文件存在且大小合理，也跳过比对
                if os.path.exists(potential_filtered_bam):
                    try:
                        size = os.path.getsize(potential_filtered_bam)
                        if size > 1000000:  # 至少1MB
                            self.logger.info(f"✅ 发现已存在的过滤BAM文件，跳过BWA比对步骤 | Found existing filtered BAM file, skipping BWA alignment: {potential_filtered_bam} ({size:,} bytes)")
                            # HapHiC需要按read name排序的BAM文件，使用原始的过滤文件
                            # HapHiC requires read name-sorted BAM, use the original filtered file
                            self.logger.info(f"✅ 使用按read name排序的过滤BAM文件以满足HapHiC要求 | Using read name-sorted filtered BAM file to meet HapHiC requirements")
                            self.config.bam_file = potential_filtered_bam
                            # HapHiC不需要BAM索引，跳过索引创建以避免错误
                            # HapHiC does not need BAM index, skip index creation to avoid errors
                        else:
                            self.logger.info(f"❌ 过滤BAM文件太小，重新执行比对 | Filtered BAM file too small, re-running alignment: {size:,} bytes")
                    except Exception as e:
                        self.logger.warning(f"⚠️ 无法检查过滤BAM文件大小: {e}")
                elif os.path.exists(potential_bam):
                    try:
                        size = os.path.getsize(potential_bam)
                        if size > 1000000:  # 至少1MB
                            self.logger.info(f"✅ 发现已存在的BAM文件，跳过BWA比对步骤 | Found existing BAM file, skipping BWA alignment: {potential_bam} ({size:,} bytes)")
                            # HapHiC需要按read name排序的BAM文件，不是coordinate排序的
                            # coordinate排序的文件通常是HiC.sorted.bam，我们需要使用原始的HiC.bam
                            # HapHiC requires read name-sorted BAM, not coordinate-sorted
                            # coordinate-sorted files are usually HiC.sorted.bam, we need to use the original HiC.bam
                            self.logger.info(f"✅ 使用按read name排序的BAM文件以满足HapHiC要求 | Using read name-sorted BAM file to meet HapHiC requirements")
                            self.config.bam_file = potential_bam
                            # HapHiC不需要BAM索引，跳过索引创建以避免错误
                            # HapHiC does not need BAM index, skip index creation to avoid errors
                        else:
                            self.logger.info(f"❌ BAM文件太小，重新执行比对 | BAM file too small, re-running alignment: {size:,} bytes")
                    except Exception as e:
                        self.logger.warning(f"⚠️ 无法检查BAM文件大小: {e}")
                else:
                    self.logger.info("🔗 执行BWA比对 | Executing BWA alignment")
                    bam_file = self.bwa_aligner.run_alignment()
                    # 更新配置中的bam_file用于后续步骤 | Update bam_file in config for subsequent steps
                    self.config.bam_file = bam_file
            else:
                # 对于BAM输入，检查索引文件是否存在
                if os.path.exists(self.config.hic_file + ".bai"):
                    self.logger.info(f"✅ BAM索引已存在，跳过索引创建 | BAM index already exists, skipping index creation")
                self.config.bam_file = self.config.hic_file

            # 创建流程管理器 | Create pipeline manager
            pipeline = HapHiCPipeline(self.config, self.logger)

            # 运行流程 | Run pipeline
            success = pipeline.run_complete_pipeline()

            # 检查是否有纠错后的contig，需要重新比对 | Check if corrected contigs exist and need realignment
            if success and self.config.correct_nrounds > 0:
                corrected_asm = os.path.join(self.config.output_dir, "01.cluster", "corrected_asm.fa")
                if os.path.exists(corrected_asm):
                    self.logger.info("🔄 检测到纠错后的contig，需要重新比对Hi-C数据 | Detected corrected contigs, need to realign Hi-C data")

                    # 重新比对Hi-C数据到纠错后的组装 | Realign Hi-C data to corrected assembly
                    self.logger.info("🔗 开始重新比对Hi-C数据到纠错后的组装 | Starting Hi-C realignment to corrected assembly")

                    # 更新组装文件引用 | Update assembly file reference
                    original_asm = self.config.asm_file
                    self.config.asm_file = corrected_asm

                    # 重新运行比对 | Re-run alignment
                    self.bwa_aligner = BWAAligner(self.config, self.logger)
                    corrected_bam_file = self.bwa_aligner.run_alignment(is_realignment=True)

                    # 更新BAM文件引用 | Update BAM file reference
                    self.config.bam_file = corrected_bam_file

                    # 重新运行后续步骤（从reassign开始）| Rerun subsequent steps (starting from reassign)
                    self.logger.info("🔄 重新运行HapHiC后续步骤 | Rerunning HapHiC subsequent steps")
                    success = pipeline.run_step("reassign") and pipeline.run_step("sort") and pipeline.run_step("build")

                    # 恢复原始组装文件引用 | Restore original assembly file reference
                    self.config.asm_file = original_asm

            # 结果验证 | Result validation
            if success:
                self._post_processing()

            # 统计信息 | Statistics
            elapsed_time = time.time() - start_time
            self._log_statistics(elapsed_time)

            return success

        except Exception as e:
            self.logger.error(f"💥 HapHiC流程执行失败 | HapHiC pipeline execution failed: {e}")
            return False

    
    def _pre_flight_checks(self) -> bool:
        """预检查 | Pre-flight checks"""
        self.logger.info("🔍 系统预检查 | System pre-flight checks")

        # 验证输入文件 | Validate input files
        if not self.quality_controller.validate_inputs():
            return False

        # 检查系统资源 | Check system resources
        if not self.quality_controller.check_system_resources():
            self.logger.warning("⚠️ 系统资源检查有警告，但将继续执行 | System resource check has warnings, but will continue")

        # 检查HapHiC工具 | Check HapHiC tool
        if not os.path.exists(self.config.haphic_bin):
            self.logger.error(f"❌ HapHiC工具不存在 | HapHiC tool not found: {self.config.haphic_bin}")
            return False

        # 测试HapHiC工具 | Test HapHiC tool
        try:
            result = subprocess.run(
                [self.config.haphic_bin, "--help"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode != 0:
                self.logger.error("❌ HapHiC工具无法正常运行 | HapHiC tool cannot run properly")
                return False
        except Exception as e:
            self.logger.error(f"❌ HapHiC工具测试失败 | HapHiC tool test failed: {e}")
            return False

        self.logger.info("✅ 预检查通过 | Pre-flight checks passed")
        return True

    def _post_processing(self):
        """后处理 | Post processing"""
        self.logger.info("📊 后处理和验证 | Post-processing and validation")

        # 验证输出文件 | Validate output files
        self.result_validator.validate_output_files()

        # 获取统计信息 | Get statistics
        stats = self.result_validator.get_statistics()
        if stats:
            self.logger.info("📈 输出统计 | Output Statistics:")
            for key, value in stats.items():
                self.logger.info(f"  {key}: {value}")

        # 生成报告 | Generate report
        self._generate_report()

    def _log_statistics(self, elapsed_time: float):
        """记录统计信息 | Log statistics"""
        self.logger.info("📊 处理统计 | Processing Statistics:")
        self.logger.info(f"  ⏱️ 总耗时 | Total time: {elapsed_time:.2f}秒")

        output_files = self.config.get_output_files()
        for file_type, file_path in output_files.items():
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.logger.info(f"  📄 {file_type}: {os.path.basename(file_path)} ({size:,} bytes)")

    def _try_create_bam_index(self, bam_file: str):
        """已弃用：HapHiC不需要BAM索引 | Deprecated: HapHiC does not need BAM index"""
        self.logger.info("ℹ️ HapHiC不需要BAM索引，跳过索引创建 | HapHiC does not need BAM index, skipping index creation")
        pass

    def _generate_report(self):
        """生成报告 | Generate report"""
        report_file = os.path.join(self.config.output_dir, f"{self.config.prefix}_haphic_report.txt")

        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("HapHiC基因组Scaffolding分析报告 | HapHiC Genome Scaffolding Report\n")
                f.write("=" * 60 + "\n\n")
                f.write(f"分析时间 | Analysis Time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n\n")

                f.write("配置参数 | Configuration Parameters:\n")
                f.write("-" * 30 + "\n")
                summary = self.config.get_summary()
                for key, value in summary.items():
                    f.write(f"{key}: {value}\n")

                f.write("\n输出文件 | Output Files:\n")
                f.write("-" * 30 + "\n")
                output_files = self.config.get_output_files()
                for file_type, file_path in output_files.items():
                    if os.path.exists(file_path):
                        size = os.path.getsize(file_path)
                        f.write(f"{file_type}: {file_path} ({size:,} bytes)\n")

                f.write("\n" + "=" * 60 + "\n")
                f.write("报告生成完成 | Report generation completed\n")

            self.logger.info(f"📄 报告已生成: {report_file}")

        except Exception as e:
            self.logger.error(f"❌ 报告生成失败 | Report generation failed: {e}")


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="HapHiC基因组Scaffolding工具 - Pipeline模式 | HapHiC Genome Scaffolding Tool - Pipeline Mode\n\n"
                   "该工具使用HapHiC pipeline模式一步完成所有scaffolding步骤，包括：\n"
                   "1. BWA比对（如果输入为FASTQ）\n"
                   "2. HapHiC聚类 (cluster)\n"
                   "3. 重新分配 (reassign)\n"
                   "4. 排序和定向 (sort)\n"
                   "5. 构建scaffolds (build)\n"
                   "6. 生成Juicebox文件（可选）\n\n"
                   "This tool uses HapHiC pipeline mode to complete all scaffolding steps in one go, including:\n"
                   "1. BWA alignment (if input is FASTQ)\n"
                   "2. HapHiC clustering\n"
                   "3. Reassignment\n"
                   "4. Ordering and orientation\n"
                   "5. Building scaffolds\n"
                   "6. Generating Juicebox files (optional)",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数 | Required parameters
    parser.add_argument('asm_file',
                       help='基因组组装文件路径 | Genome assembly file path (FASTA)')
    parser.add_argument('bam_file',
                       help='Hi-C BAM文件路径 | Hi-C BAM file path')
    parser.add_argument('nchrs', type=int,
                       help='染色体数量 | Number of chromosomes')

    # 工具路径 | Tool path
    parser.add_argument('--haphic-bin',
                       default="/share/org/YZWL/yzwl_lixg/miniforge3/envs/haphic/bin/haphic",
                       help='HapHiC可执行文件路径 | HapHiC executable path')

    # 输出配置 | Output configuration
    parser.add_argument('-o', '--output-dir',
                       help='输出目录路径 | Output directory path')
    parser.add_argument('--prefix',
                       help='输出文件前缀 | Output file prefix')
    parser.add_argument('--force-rerun', action='store_true',
                       help='强制重新运行所有步骤 | Force rerun all steps (disable resume mode)')

    # Hi-C数据处理参数 | Hi-C data processing parameters
    parser.add_argument('--mapq-threshold', type=int, default=1,
                       help='MAPQ阈值 | MAPQ threshold (default: 1)')
    parser.add_argument('--edit-distance', type=int, default=3,
                       help='编辑距离阈值 | Edit distance threshold (default: 3)')

    # 聚类参数 | Clustering parameters
    parser.add_argument('--min-inflation', type=float, default=1.0,
                       help='最小膨胀参数 | Min inflation (default: 1.0)')
    parser.add_argument('--max-inflation', type=float, default=3.0,
                       help='最大膨胀参数 | Max inflation (default: 3.0)')
    parser.add_argument('--inflation-step', type=float, default=0.2,
                       help='膨胀参数步长 | Inflation step (default: 0.2)')
    parser.add_argument('--Nx', type=int, default=80,
                       help='Nx参数 | Nx parameter (default: 80)')
    parser.add_argument('--min-group-len', type=int, default=0,
                       help='最小分组长度 | Min group length (default: 0)')

    # 排序和定向参数 | Ordering and orientation parameters
    parser.add_argument('--processes', type=int, default=8,
                       help='并行进程数 | Number of parallel processes (default: 8)')
    parser.add_argument('--no-fast-sorting', action='store_true',
                       help='禁用快速排序 | Disable fast sorting')
    parser.add_argument('--no-allhic-optimization', action='store_true',
                       help='禁用ALLHiC优化 | Disable ALLHiC optimization')

    # 组装校正参数 | Assembly correction parameters
    parser.add_argument('--correct-nrounds', type=int, default=0,
                       help='组装校正轮数 | Assembly correction rounds (default: 0)')
    parser.add_argument('--correct-min-coverage', type=float, default=10.0,
                       help='校正最小覆盖度 | Correction min coverage (default: 10.0)')

    # 单倍型分相参数 | Haplotype phasing parameters
    parser.add_argument('--remove-allelic-links', type=int,
                       help='移除等位基因连接数 | Remove allelic links count')
    parser.add_argument('--phasing-weight', type=float, default=1.0,
                       help='分相权重 | Phasing weight (default: 1.0)')
    parser.add_argument('--gfa-files',
                       help='GFA文件路径 | GFA files path (comma-separated)')

    # 可视化参数 | Visualization parameters
    parser.add_argument('--generate-plots', action='store_true',
                       help='生成可视化图表 | Generate visualization plots')
    parser.add_argument('--bin-size', type=int, default=500,
                       help='接触图装箱大小 | Contact map bin size (default: 500)')
    parser.add_argument('--min-len', type=float, default=1.0,
                       help='最小scaffold长度 | Min scaffold length (default: 1.0)')
    parser.add_argument('--separate-plots', action='store_true',
                       help='生成单独图表 | Generate separate plots')

    # 性能参数 | Performance parameters
    parser.add_argument('--threads', type=int, default=8,
                       help='线程数 | Number of threads (default: 8)')
    parser.add_argument('--memory-limit',
                       help='内存限制 | Memory limit (e.g., 64G)')

    # 高级选项 | Advanced options
    parser.add_argument('--quick-view', action='store_true',
                       help='快速查看模式 | Quick view mode')
    parser.add_argument('--RE', default="GATC",
                       help='限制性内切酶位点 | Restriction enzyme sites (default: GATC)')
    parser.add_argument('--skip-clustering', action='store_true',
                       help='跳过聚类步骤 | Skip clustering step')

    # 输出格式选项 | Output format options
    parser.add_argument('--no-agp', action='store_true',
                       help='不输出AGP文件 | Don\'t output AGP file')
    parser.add_argument('--no-fasta', action='store_true',
                       help='不输出FASTA文件 | Don\'t output FASTA file')
    parser.add_argument('--no-juicebox', action='store_true',
                       help='不生成Juicebox脚本 | Don\'t generate Juicebox script')

    # 日志配置 | Logging configuration
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='详细输出模式 | Verbose output mode')
    parser.add_argument('--log-file',
                       help='日志文件路径 | Log file path')

    args = parser.parse_args()

    try:
        # 创建处理器 | Create processor
        processor = HapHiCProcessor(
            asm_file=args.asm_file,
            bam_file=args.bam_file,
            nchrs=args.nchrs,
            haphic_bin=args.haphic_bin,
            output_dir=args.output_dir,
            prefix=args.prefix,
            force_rerun=args.force_rerun,
            mapq_threshold=args.mapq_threshold,
            edit_distance=args.edit_distance,
            min_re_sites=args.min_RE_sites,
            min_inflation=args.min_inflation,
            max_inflation=args.max_inflation,
            inflation_step=args.inflation_step,
            nx=args.Nx,
            min_group_len=args.min_group_len,
            processes=args.processes,
            fast_sorting=not args.no_fast_sorting,
            allhic_optimization=not args.no_allhic_optimization,
            correct_nrounds=args.correct_nrounds,
            correct_min_coverage=args.correct_min_coverage,
            remove_allelic_links=args.remove_allelic_links,
            phasing_weight=args.phasing_weight,
            gfa_files=args.gfa_files,
            generate_plots=args.generate_plots,
            bin_size=args.bin_size,
            min_len=args.min_len,
            separate_plots=args.separate_plots,
            threads=args.threads,
            memory_limit=args.memory_limit,
            quick_view=args.quick_view,
            re_sites=args.RE,
            skip_clustering=args.skip_clustering,
            output_agp=not args.no_agp,
            output_fasta=not args.no_fasta,
            output_juicebox=not args.no_juicebox,
            verbose=args.verbose,
            log_file=args.log_file,
            dry_run=False  # main.py不支持dry_run模式，默认为False
        )

        # 运行分析 | Run analysis
        success = processor.run_pipeline()

        if not success:
            sys.exit(1)

    except Exception as e:
        print(f"❌ 程序执行失败 | Program execution failed: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()