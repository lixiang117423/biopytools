"""
基因组组装主程序模块|Genome Assembly Main Module
"""

import argparse
import sys
import os
from pathlib import Path
from datetime import datetime

# 添加当前目录到sys.path，解决相对导入问题
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

try:
    from config import AssemblyConfig
    from logger import AssemblyLogger
    from utils import check_dependencies
    from assembler import HifiasmAssembler
    from report import ReportGenerator
    from ngs_polisher import NGSPolisher
    from purge_dups_wrapper import PurgeDupsWrapper
except ImportError:
    # 如果相对导入失败，尝试绝对导入
    sys.path.append(str(current_dir.parent))
    from biopytools.genome_assembler.config import AssemblyConfig
    from biopytools.genome_assembler.logger import AssemblyLogger
    from biopytools.genome_assembler.utils import check_dependencies
    from biopytools.genome_assembler.assembler import HifiasmAssembler
    from biopytools.genome_assembler.report import ReportGenerator
    from biopytools.genome_assembler.ngs_polisher import NGSPolisher
    from biopytools.genome_assembler.purge_dups_wrapper import PurgeDupsWrapper

class GenomeAssembler:
    """基因组组装主类|Main Genome Assembler Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = AssemblyConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        logger_path = Path(self.config.log_dir)
        self.logger_manager = AssemblyLogger(logger_path)
        self.logger = self.logger_manager.get_logger()

        # 初始化组件|Initialize components
        self.assembler = HifiasmAssembler(self.config, self.logger)
        self.report_generator = ReportGenerator(self.config, self.logger)

        # 如果有NGS数据，初始化NGS polisher|Initialize NGS polisher if NGS data provided
        if self.config.has_ngs:
            self.ngs_polisher = NGSPolisher(self.config, self.logger)
        else:
            self.ngs_polisher = None

        # 如果启用了Purge_Dups，初始化去冗余器|Initialize Purge_Dups if enabled
        if self.config.has_purge_dups:
            self.purge_dups_wrapper = PurgeDupsWrapper(self.config, self.logger)
        else:
            self.purge_dups_wrapper = None

    def run_assembly(self):
        """运行完整的组装流程|Run complete assembly pipeline"""
        start_time = datetime.now()

        try:
            self.logger.info("=" * 80)
            self.logger.info("开始基因组组装流程|Starting Genome Assembly Pipeline")
            self.logger.info("=" * 80)

            # 1. 检查依赖|Check dependencies
            check_dependencies(self.logger)

            # 2. 显示配置信息|Display configuration
            self._print_config()

            # 检查断点续传状态|Check resume status
            if self.config.resume:
                completed_steps = self.config.get_completed_steps()
                self._log_resume_status(completed_steps)
            else:
                completed_steps = {'assembly': False, 'fasta_converted': False, 'reads_mapped': False}

            # 3. 执行初次组装|Execute initial assembly
            if not completed_steps['assembly']:
                if not self.assembler.run_assembly():
                    sys.exit(1)
            else:
                self.logger.info("检测到组装已完成，跳过|Assembly already completed, skipping")

            # 4. 转换格式|Convert formats
            if not completed_steps['fasta_converted']:
                fasta_results = self.assembler.convert_gfa_to_fasta()
            else:
                self.logger.info("检测到FASTA转换已完成，跳过|FASTA conversion already completed, skipping")
                # 重新获取统计信息|Re-fetch statistics
                fasta_results = self.assembler.convert_gfa_to_fasta()

            # 5. 生成contig-reads映射|Generate contig-reads mapping
            if not completed_steps['reads_mapped']:
                self.assembler.generate_contig_reads_mapping()
            else:
                self.logger.info("检测到contig-reads映射已生成，跳过|Contig-reads mapping already generated, skipping")

            # 6. 如果有NGS数据，执行NGS polish|If NGS data provided, run NGS polish
            if self.config.has_ngs:
                self.logger.info("=" * 80)
                self.logger.info("检测到NGS数据，开始NGS polish流程|NGS data detected, starting NGS polish")
                self.logger.info("=" * 80)

                if not self.ngs_polisher.run_ngs_polish(completed_steps):
                    self.logger.error("NGS polish流程失败|NGS polish pipeline failed")
                    sys.exit(1)

                # 更新fasta_results为polish后的结果|Update fasta_results to polished results
                polished_genome = os.path.join(self.config.ngs_polish_dir,
                                             f"{self.config.prefix}.polished.fa")
                if os.path.exists(polished_genome):
                    # 可以选择使用polish后的结果|Optionally use polished results
                    self.logger.info(f"Polished基因组|Polished genome: {polished_genome}")

            # 7. 如果启用了Purge_Dups，执行去冗余|If Purge_Dups enabled, run deduplication
            final_genome = None
            if self.config.has_purge_dups:
                self.logger.info("=" * 80)
                self.logger.info("检测到Purge_Dups启用，开始去冗余流程|Purge_Dups enabled, starting deduplication")
                self.logger.info("=" * 80)

                # 确定去冗余的输入文件|Determine input file for deduplication
                # 优先使用NGS polish后的结果，否则使用组装结果
                # Prefer NGS polished result, otherwise use assembly result
                if self.config.has_ngs:
                    polished_genome = os.path.join(self.config.ngs_polish_dir,
                                                 f"{self.config.prefix}.polished.fa")
                    if os.path.exists(polished_genome):
                        purge_input = polished_genome
                        self.logger.info(f"使用Polished基因组进行去冗余|Using Polished genome for deduplication: {purge_input}")
                    else:
                        # 使用primary.fa|Use primary.fa
                        if self.config.has_hic:
                            purge_input = os.path.join(self.config.fasta_dir, f"{self.config.prefix}.primary.fa")
                        else:
                            purge_input = os.path.join(self.config.fasta_dir, f"{self.config.prefix}.primary.fa")
                        self.logger.info(f"使用Assembly基因组进行去冗余|Using Assembly genome for deduplication: {purge_input}")
                else:
                    # 没有NGS polish，直接使用组装结果|No NGS polish, use assembly result directly
                    if self.config.has_hic:
                        purge_input = os.path.join(self.config.fasta_dir, f"{self.config.prefix}.primary.fa")
                    else:
                        purge_input = os.path.join(self.config.fasta_dir, f"{self.config.prefix}.primary.fa")
                    self.logger.info(f"使用Assembly基因组进行去冗余|Using Assembly genome for deduplication: {purge_input}")

                # 运行去冗余|Run deduplication
                final_genome = self.purge_dups_wrapper.run_purge_dups(purge_input)

                if not final_genome:
                    self.logger.error("Purge_Dups去冗余流程失败|Purge_Dups deduplication pipeline failed")
                    sys.exit(1)

            end_time = datetime.now()

            # 7. 生成报告|Generate report
            self.report_generator.generate_summary(fasta_results, start_time, end_time)

            # 8. 输出结果|Output results
            duration = int((end_time - start_time).total_seconds())
            self._print_summary(duration, fasta_results)

        except Exception as e:
            self.logger.error(f"组装流程在执行过程中意外终止|Assembly pipeline terminated unexpectedly: {e}")
            sys.exit(1)

    def _log_resume_status(self, completed_steps: dict):
        """
        输出断点续传状态|Log resume status

        Args:
            completed_steps: 已完成步骤|Completed steps
        """
        self.logger.info("=" * 80)
        self.logger.info("断点续传模式|Resume mode: enabled")
        self.logger.info("=" * 80)
        self.logger.info("检查已完成步骤|Checking completed steps:")

        # 初次组装步骤|Initial assembly steps
        self.logger.info("初次组装步骤|Initial assembly steps:")
        step_names = {
            'assembly': 'Hifiasm组装|Hifiasm assembly',
            'fasta_converted': 'GFA转FASTA|GFA to FASTA conversion',
            'reads_mapped': 'Contig-reads映射生成|Contig-reads mapping generated'
        }

        for step, name in step_names.items():
            status = "已完成|completed" if completed_steps[step] else "未完成|pending"
            self.logger.info(f"  [{status}] {name}")

        # NGS polish步骤|NGS polish steps
        if self.config.has_ngs:
            self.logger.info("NGS polish步骤|NGS polish steps:")
            ngs_step_names = {
                'bwa_alignment': 'BWA比对|BWA alignment',
                'coverage_filter': '覆盖度过滤|Coverage filter',
                'filtered_reads': '筛选reads|Filter reads',
                'reassembly': '重新组装|Reassembly'
            }

            for step, name in ngs_step_names.items():
                if step in completed_steps:
                    status = "已完成|completed" if completed_steps[step] else "未完成|pending"
                    self.logger.info(f"  [{status}] {name}")

        # Purge_Dups去冗余步骤|Purge_Dups deduplication steps
        if self.config.has_purge_dups:
            self.logger.info("Purge_Dups去冗余步骤|Purge_Dups deduplication steps:")
            purge_step_names = {
                'purge_dups_coverage': '覆盖度计算|Coverage calculation',
                'purge_dups_cutoffs': '阈值计算|Cutoffs calculation',
                'purge_dups_split_align': '分割和自比对|Split and alignment',
                'purge_dups_purge': '去冗余识别|Deduplication identification',
                'purge_dups_seqs': '序列提取|Sequence extraction'
            }

            for step, name in purge_step_names.items():
                if step in completed_steps:
                    status = "已完成|completed" if completed_steps[step] else "未完成|pending"
                    self.logger.info(f"  [{status}] {name}")

        self.logger.info("=" * 80)

    def _print_config(self):
        """打印配置信息|Print configuration"""
        self.logger.info("组装配置|Assembly Configuration:")
        self.logger.info(f"  样本名称|Sample Name: {self.config.prefix}")
        self.logger.info(f"  基因组大小|Genome Size: {self.config.genome_size}")
        self.logger.info(f"  倍性|Ploidy: {self.config.n_hap}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  HiFi数据|HiFi Data: {self.config.hifi_data}")

        if self.config.has_hic:
            self.logger.info(f"  Hi-C R1|Hi-C R1: {self.config.hic_r1}")
            self.logger.info(f"  Hi-C R2|Hi-C R2: {self.config.hic_r2}")
            self.logger.info(f"  组装模式|Assembly Mode: HiFi + Hi-C")
        else:
            self.logger.info(f"  组装模式|Assembly Mode: HiFi only")

        if self.config.has_ngs:
            self.logger.info(f"  NGS数据目录|NGS Data Directory: {self.config.ngs_data}")
            self.logger.info(f"  高质量覆盖度阈值|High Coverage Threshold: {self.config.high_cov}%")
            self.logger.info(f"  NGS polish|NGS Polish: 启用|Enabled")

        if self.config.has_purge_dups:
            self.logger.info(f"  Purge_Dups去冗余|Purge_Dups Deduplication: 启用|Enabled")
            self.logger.info(f"  去冗余线程数|Deduplication Threads: {self.config.purge_dups_threads}")
            self.logger.info(f"  去冗余reads类型|Deduplication Reads Type: {self.config.purge_dups_read_type}")

        self.logger.info(f"  工作目录|Work Directory: {self.config.work_dir}")

    def _print_summary(self, duration: int, fasta_results: dict):
        """打印摘要|Print summary"""
        self.logger.info("=" * 80)
        self.logger.info("组装流程完成|Assembly Pipeline Completed")
        self.logger.info("=" * 80)
        self.logger.info(f"总耗时|Total Duration: {duration} 秒")
        self.logger.info(f"生成FASTA文件|Generated FASTA files: {len(fasta_results)}")

        # 打印文件结构|Print file structure
        self.report_generator.print_file_tree()

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='HiFi基因组组装工具(hifiasm)|HiFi Genome Assembly Tool (hifiasm)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s --hifi hifi.fq -p sample1
  %(prog)s --hifi hifi.fq --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz -p sample1
  %(prog)s --hifi hifi.fq --ngs ngs_data_dir -p sample1 --high-cov 90
  %(prog)s --hifi hifi.fq --hic-r1 hic_R1.fq.gz --hic-r2 hic_R2.fq.gz -p sample1 --no-purge-dups
        """
    )

    # 必需参数|Required arguments
    parser.add_argument('--hifi', '-i', required=True,
                       help='HiFi数据文件路径|Path to HiFi data file')

    # 可选参数：Hi-C数据|Optional parameters: Hi-C data
    parser.add_argument('--hic-r1', '-1',
                       help='Hi-C R1文件路径|Path to Hi-C R1 file')
    parser.add_argument('--hic-r2', '-2',
                       help='Hi-C R2文件路径|Path to Hi-C R2 file')

    # 样本参数|Sample parameters
    parser.add_argument('--prefix', '-p', default="genome_sample",
                       help='样本前缀|Sample prefix')

    # 组装参数|Assembly parameters
    parser.add_argument('--threads', '-t', type=int, default=88,
                       help='线程数|Number of threads')

    parser.add_argument('--genome-size', '-g', default="1.45g",
                       help='预估基因组大小|Estimated genome size (e.g., 1.45g, 250m)')

    parser.add_argument('--n-hap', type=int, default=2,
                       help='倍性|Ploidy (haploid count)')

    parser.add_argument('--purge-level', '-l', type=int, default=None,
                       help='Purge level (0=no purging, 1=light, 2/3=aggressive) [default: 3 for unzip]')

    parser.add_argument('--hom-cov', type=int, default=None,
                       help='Homozygous read coverage (--hom-cov) [default: auto]')

    # 路径参数|Path parameters
    parser.add_argument('--output', '-o', default="./assembly_output",
                       help='输出目录|Output directory')

    # NGS polish参数|NGS polish parameters
    parser.add_argument('--ngs',
                       help='NGS二代数据目录|NGS second-generation data directory (optional)')
    parser.add_argument('--ngs-pattern', default="_1.clean.fq.gz",
                       help='NGS文件匹配模式|NGS file matching pattern (default: _1.clean.fq.gz)')
    parser.add_argument('--high-cov', type=float, default=95.0,
                       help='高质量contig覆盖度阈值|High quality contig coverage threshold (default: 95.0)')
    parser.add_argument('--medium-cov-min', type=float, default=30.0,
                       help='中等质量contig最小覆盖度|Medium quality contig minimum coverage (default: 30.0)')

    # Purge_Dups去冗余参数|Purge_Dups deduplication parameters
    parser.add_argument('--no-purge-dups', action='store_true',
                       help='禁用Purge_Dups去冗余|Disable Purge_Dups deduplication (enabled by default)')
    parser.add_argument('--purge-dups-path', default='~/miniforge3/envs/purge_dups_v.1.2.6',
                       help='Purge_Dups软件路径|Purge_Dups software path (default: ~/miniforge3/envs/purge_dups_v.1.2.6)')
    parser.add_argument('--purge-dups-threads', type=int, default=None,
                       help='去冗余线程数|Deduplication threads (default: same as assembly threads)')
    parser.add_argument('--purge-dups-read-type', choices=['pacbio', 'hifi', 'illumina'],
                       default='hifi', help='去冗余reads类型|Deduplication reads type (default: hifi)')

    # 执行控制|Execution control
    parser.add_argument('--no-resume',
                       action='store_true',
                       help='禁用断点续传（强制重新运行所有步骤）|Disable resume mode (force rerun all steps)')
    parser.add_argument('--resume',
                       action='store_true',
                       help='启用断点续传（默认已启用）|Enable resume mode (enabled by default)')

    args = parser.parse_args()

    # 确定resume参数值|Determine resume parameter value
    # 默认启用断点续传，使用--no-resume可禁用|Resume enabled by default, use --no-resume to disable
    resume_value = not args.no_resume

    # 确定purge_dups参数值|Determine purge_dups parameter value
    # 默认启用去冗余，使用--no-purge-dups可禁用|Purge_dups enabled by default, use --no-purge-dups to disable
    purge_dups_value = not args.no_purge_dups

    # 创建组装器并运行|Create assembler and run
    assembler = GenomeAssembler(
        hifi_data=args.hifi,
        hic_r1=args.hic_r1,
        hic_r2=args.hic_r2,
        prefix=args.prefix,
        threads=args.threads,
        genome_size=args.genome_size,
        n_hap=args.n_hap,
        purge_level=args.purge_level,
        hom_cov=args.hom_cov,
        base_dir=args.output,
        resume=resume_value,
        ngs_data=args.ngs,
        ngs_pattern=args.ngs_pattern,
        high_cov=args.high_cov,
        medium_cov_min=args.medium_cov_min,
        enable_purge_dups=purge_dups_value,
        purge_dups_path=args.purge_dups_path,
        purge_dups_threads=args.purge_dups_threads,
        purge_dups_read_type=args.purge_dups_read_type
    )

    assembler.run_assembly()

if __name__ == "__main__":
    main()
