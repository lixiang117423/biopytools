"""
HiCanu计算器模块|HiCanu Calculator Module
"""

import os
import glob
from typing import Optional, Dict
from .config import HiCanuConfig
from .utils import CommandRunner, parse_fasta_stats, format_size, generate_contig_reads_map


class HiCanuCalculator:
    """HiCanu计算器类|HiCanu Calculator Class"""

    def __init__(self, config: HiCanuConfig, logger, cmd_runner: CommandRunner):
        """
        初始化HiCanu计算器|Initialize HiCanu calculator

        Args:
            config: HiCanu配置对象|HiCanu config object
            logger: 日志对象|Logger object
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_assembly(self) -> bool:
        """
        运行HiCanu组装|Run HiCanu assembly

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始HiCanu组装|Starting HiCanu assembly")
        self.logger.info("=" * 60)

        # 输出配置信息|Output configuration info
        self._log_configuration()

        # 检查断点续传状态|Check resume status
        if self.config.resume:
            completed_steps = self.config.get_completed_steps()
            self._log_resume_status(completed_steps)
        else:
            completed_steps = {'assembly': False, 'fasta_copied': False, 'reads_mapped': False}

        # 步骤1: 执行Canu组装|Step 1: Run Canu assembly
        if not completed_steps['assembly']:
            cmd = self.config.get_canu_command()
            success = self.cmd_runner.run_command(
                cmd,
                dry_run=self.config.dry_run
            )

            if not success:
                self.logger.error("HiCanu组装失败|HiCanu assembly failed")
                return False
        else:
            self.logger.info("检测到组装已完成，跳过|Assembly already completed, skipping")

        # 检查输出文件|Check output files
        if not self.config.dry_run:
            if not self._check_outputs():
                self.logger.error("输出文件检查失败|Output files check failed")
                return False

            # 解析组装结果|Parse assembly results
            if not completed_steps.get('assembly'):
                self._parse_assembly_results()

            # 步骤2: 复制FASTA文件|Step 2: Copy FASTA files
            if not completed_steps['fasta_copied']:
                self._copy_fasta_files()
            else:
                self.logger.info("检测到FASTA文件已复制，跳过|FASTA files already copied, skipping")

            # 步骤3: 生成contig-reads映射文件|Step 3: Generate contig-reads mapping
            if not completed_steps['reads_mapped']:
                self._generate_contig_reads_mapping()
            else:
                self.logger.info("检测到contig-reads映射已生成，跳过|Contig-reads mapping already generated, skipping")

        self.logger.info("=" * 60)
        self.logger.info("HiCanu组装完成|HiCanu assembly completed successfully")
        self.logger.info("=" * 60)

        return True

    def _log_resume_status(self, completed_steps: dict):
        """
        输出断点续传状态|Log resume status

        Args:
            completed_steps: 已完成步骤|Completed steps
        """
        self.logger.info("=" * 60)
        self.logger.info("断点续传模式|Resume mode: enabled")
        self.logger.info("=" * 60)
        self.logger.info("检查已完成步骤|Checking completed steps:")

        step_names = {
            'assembly': 'Canu组装|Canu assembly',
            'fasta_copied': 'FASTA文件复制|FASTA files copied',
            'reads_mapped': 'Contig-reads映射生成|Contig-reads mapping generated'
        }

        for step, name in step_names.items():
            status = "已完成|completed" if completed_steps[step] else "未完成|pending"
            self.logger.info(f"  [{status}] {name}")

        self.logger.info("=" * 60)

    def _log_configuration(self):
        """输出配置信息|Log configuration information"""
        self.logger.info("组装配置|Assembly configuration:")
        self.logger.info(f"  Reads文件|Reads file: {self.config.reads_file}")
        self.logger.info(f"  基因组大小|Genome size: {self.config.genome_size}")
        self.logger.info(f"  工作目录|Work directory: {self.config.work_dir}")
        self.logger.info(f"  前缀|Prefix: {self.config.prefix}")
        self.logger.info(f"  组装阶段|Assembly stage: {self.config.stage}")
        self.logger.info(f"  最小reads长度|Min read length: {self.config.min_read_length}")
        self.logger.info(f"  最小重叠长度|Min overlap length: {self.config.min_overlap_length}")

        if self.config.corrected_error_rate:
            self.logger.info(f"  纠正错误率|Corrected error rate: {self.config.corrected_error_rate}")

        if self.config.raw_error_rate:
            self.logger.info(f"  原始错误率|Raw error rate: {self.config.raw_error_rate}")

        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  内存限制|Memory limit: {self.config.memory}")

        if self.config.dry_run:
            self.logger.info("  模拟运行模式|Dry run mode: True")

    def _check_outputs(self) -> bool:
        """
        检查输出文件|Check output files

        Returns:
            bool: 所有文件是否都存在|Whether all files exist
        """
        self.logger.info("检查输出文件|Checking output files")

        output_files = self.config.get_output_files()
        all_exist = True

        # 检查contigs文件|Check contigs file
        contigs_file = output_files['contigs']
        if self.cmd_runner.check_output_exists(contigs_file, error_on_missing=True):
            self.logger.info(f"  Contigs文件|Contigs file: {contigs_file}")
        else:
            self.logger.warning(f"  Contigs文件不存在|Contigs file not found: {contigs_file}")
            # contigs文件可能不存在（取决于stage），继续检查其他文件

        # 检查unitigs文件|Check unitigs file (unitigs对于HiFi组装是可选的)|Check unitigs file (optional for HiFi assembly)
        unitigs_file = output_files['unitigs']
        if self.cmd_runner.check_output_exists(unitigs_file, error_on_missing=False, silent=True):
            self.logger.info(f"  Unitigs文件|Unitigs file: {unitigs_file}")
        else:
            self.logger.debug(f"  Unitigs文件不存在（对于HiFi组装这是正常的）|Unitigs file not found (normal for HiFi assembly): {unitigs_file}")

        # 检查报告文件|Check report file
        report_file = output_files['report']
        if os.path.exists(report_file):
            self.logger.info(f"  报告文件|Report file: {report_file}")
        else:
            self.logger.warning(f"  报告文件不存在|Report file not found: {report_file}")

        # 至少要有一个组装结果文件|At least one assembly result file must exist
        if not os.path.exists(contigs_file) and not os.path.exists(unitigs_file):
            self.logger.error("未找到任何组装结果文件|No assembly result files found")
            return False

        return True

    def _parse_assembly_results(self):
        """解析组装结果|Parse assembly results"""
        self.logger.info("=" * 60)
        self.logger.info("解析组装结果|Parsing assembly results")
        self.logger.info("=" * 60)

        output_files = self.config.get_output_files()

        # 解析contigs|Parse contigs
        if os.path.exists(output_files['contigs']):
            self.logger.info("解析Contigs|Parsing contigs:")
            contig_stats = parse_fasta_stats(output_files['contigs'], self.logger)

            if contig_stats:
                self.logger.info(f"  组装大小|Assembly size: {format_size(contig_stats['total_bp'])}")
                self.logger.info(f"  Contigs数量|Number of contigs: {contig_stats['total_sequences']:,}")
                self.logger.info(f"  N50: {format_size(contig_stats['n50'])}")
                self.logger.info(f"  L50: {contig_stats['l50']:,}")

        # 解析unitigs|Parse unitigs
        if os.path.exists(output_files['unitigs']):
            self.logger.info("解析Unitigs|Parsing unitigs:")
            unitig_stats = parse_fasta_stats(output_files['unitigs'], self.logger)

            if unitig_stats:
                self.logger.info(f"  组装大小|Assembly size: {format_size(unitig_stats['total_bp'])}")
                self.logger.info(f"  Unitigs数量|Number of unitigs: {unitig_stats['total_sequences']:,}")
                self.logger.info(f"  N50: {format_size(unitig_stats['n50'])}")

        # 解析Canu报告|Parse Canu report
        if os.path.exists(output_files['report']):
            self.logger.info("解析Canu报告|Parsing Canu report:")
            self._parse_canu_report(output_files['report'])

    def _parse_canu_report(self, report_file: str):
        """
        解析Canu报告文件|Parse Canu report file

        Args:
            report_file: 报告文件路径|Report file path
        """
        try:
            with open(report_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    # 提取关键信息|Extract key information
                    if 'genomeSize' in line or 'total length' in line or \
                       'number of contigs' in line or 'N50' in line:
                        self.logger.info(f"  {line}")

        except Exception as e:
            self.logger.warning(f"解析报告文件失败|Failed to parse report file: {str(e)}")

    def get_assembly_summary(self) -> Dict:
        """
        获取组装结果摘要|Get assembly results summary

        Returns:
            dict: 组装摘要|Assembly summary
        """
        summary = {
            'assembler': 'Canu',
            'version': '2.3',
            'config': {
                'genome_size': self.config.genome_size,
                'min_read_length': self.config.min_read_length,
                'min_overlap_length': self.config.min_overlap_length,
                'stage': self.config.stage,
            },
            'output_files': self.config.get_output_files(),
        }

        # 添加统计信息|Add statistics
        output_files = self.config.get_output_files()

        if os.path.exists(output_files['contigs']):
            contig_stats = parse_fasta_stats(output_files['contigs'], self.logger)
            summary['contigs'] = contig_stats

        if os.path.exists(output_files['unitigs']):
            unitig_stats = parse_fasta_stats(output_files['unitigs'], self.logger)
            summary['unitigs'] = unitig_stats

        return summary

    def _copy_fasta_files(self):
        """复制fasta文件到02.fasta目录|Copy fasta files to 02.fasta directory"""
        import shutil

        self.logger.info("=" * 60)
        self.logger.info("复制FASTA文件|Copying FASTA files")
        self.logger.info("=" * 60)

        # 查找01.raw_output中的所有fasta文件|Find all fasta files in 01.raw_output
        fasta_files = glob.glob(os.path.join(self.config.raw_dir, "*.fasta"))

        if not fasta_files:
            self.logger.warning(f"未找到FASTA文件|No FASTA files found in {self.config.raw_dir}")
            return

        copied_count = 0
        for fasta_file in fasta_files:
            filename = os.path.basename(fasta_file)
            dest_path = os.path.join(self.config.fasta_dir, filename)
            shutil.copy2(fasta_file, dest_path)
            copied_count += 1
            self.logger.info(f"  已复制|Copied: {filename}")

        self.logger.info(f"共复制 {copied_count} 个文件到|Copied {copied_count} files to: {self.config.fasta_dir}")

    def _generate_contig_reads_mapping(self):
        """生成contig-reads映射文件|Generate contig-reads mapping file"""
        generate_contig_reads_map(
            raw_dir=self.config.raw_dir,
            prefix=self.config.prefix,
            fasta_dir=self.config.fasta_dir,
            logger=self.logger
        )
