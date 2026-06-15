"""
叶绿体基因组组装主程序模块|Plastome Assembly Main Module
"""

import os
import sys
from pathlib import Path
from typing import Optional
from .config import PlastomeConfig
from .utils import PlastomeLogger, CommandRunner, detect_reads_files, detect_samples_and_reads


class PlastomeAssembler:
    """叶绿体基因组组装主类|Main Plastome Assembly Class"""

    def __init__(self, **kwargs):
        """
        初始化叶绿体组装器|Initialize plastome assembler

        Args:
            **kwargs: 配置参数|Configuration parameters
        """
        # 初始化配置|Initialize configuration
        self.config = PlastomeConfig(**kwargs)

        # 自动检测reads文件|Auto-detect reads files
        r1, r2, unpaired = detect_reads_files(
            self.config.input_dir,
            None,  # logger - 还未初始化|not yet initialized
            self.config.read1_suffix,
            self.config.read2_suffix
        )
        self.config.r1_file = r1
        self.config.r2_file = r2
        self.config.unpaired_files = unpaired

        # 验证配置|Validate configuration
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = PlastomeLogger(
            self.config.output_dir,
            False,  # verbose
            os.path.join(self.config.output_dir, 'plastome_assembly.log')
        )
        self.logger = self.logger_manager.get_logger()

        # 获取conda环境路径|Get conda environment path
        self.conda_env_path = os.path.dirname(os.path.dirname(self.config.getorganelle_path))
        self.logger.debug(f"Conda环境路径|Conda environment path: {self.conda_env_path}")

        # 获取当前目录作为工作目录（避免中文路径问题）|Get current directory as working directory (to avoid Chinese path issues)
        self.working_dir = os.getcwd()
        self.logger.debug(f"工作目录|Working directory: {self.working_dir}")

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir, self.conda_env_path, self.working_dir)

        # 输出配置信息|Output configuration information
        self._log_configuration()

    def update_reads_files(self, r1_file: str, r2_file: str):
        """
        更新reads文件并重新记录配置信息|Update reads files and re-log configuration

        Args:
            r1_file: R1文件路径|R1 file path
            r2_file: R2文件路径|R2 file path
        """
        self.config.r1_file = r1_file
        self.config.r2_file = r2_file
        # 重新记录配置以显示正确的文件信息|Re-log configuration to show correct file info
        self._log_configuration()

    @classmethod
    def run_batch(cls, input_dir: str, output_dir: str, **kwargs) -> dict:
        """
        批量处理多个样品|Batch process multiple samples

        Args:
            input_dir: 输入目录|Input directory
            output_dir: 输出目录|Output directory
            **kwargs: 其他配置参数|Other configuration parameters

        Returns:
            dict: 每个样品的处理结果|Processing results for each sample
                  {
                      'sample_name': True/False  # 是否成功|Whether successful
                  }
        """
        # 创建临时logger用于检测|Create temporary logger for detection
        import logging
        detection_logger = logging.getLogger('batch_detection')
        detection_logger.setLevel(logging.INFO)
        detection_logger.handlers.clear()  # 清除已有的handlers|Clear existing handlers
        handler = logging.StreamHandler(sys.stdout)
        handler.setLevel(logging.INFO)

        # 设置时间格式|Set time format
        formatter = logging.Formatter(
            '%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S'
        )
        handler.setFormatter(formatter)
        detection_logger.addHandler(handler)

        # 检测所有样品|Detect all samples
        # 从 kwargs 中获取后缀配置，如果没有则使用默认值|Get suffix config from kwargs, use defaults if not specified
        read1_suffix = kwargs.get('read1_suffix', '_1.clean.fq.gz')
        read2_suffix = kwargs.get('read2_suffix', '_2.clean.fq.gz')
        samples = detect_samples_and_reads(input_dir, detection_logger, read1_suffix, read2_suffix)

        if not samples:
            detection_logger.error("未检测到任何样品|No samples detected")
            return {}

        results = {}

        # 处理每个样品|Process each sample
        for sample_name, files in sorted(samples.items()):
            detection_logger.info("")
            detection_logger.info("=" * 60)
            detection_logger.info(f"处理样品|Processing sample: {sample_name} ({len(results) + 1}/{len(samples)})")
            detection_logger.info("=" * 60)

            # 为每个样品创建独立的输出目录|Create separate output directory for each sample
            # GetOrganelle会自动在输出目录下创建以prefix命名的子目录
            # 所以这里直接使用output_dir作为基础目录，让GetOrganelle创建 sample_name 子目录
            # GetOrganelle will auto-create subdirectory named with prefix under output directory
            # So use output_dir directly as base directory, let GetOrganelle create sample_name subdirectory
            sample_output_dir = output_dir

            try:
                # 创建组装器并运行|Create assembler and run
                assembler = cls(
                    input_dir=input_dir,
                    output_dir=sample_output_dir,
                    output_prefix=sample_name,
                    **kwargs
                )

                # 手动设置reads文件并重新记录配置|Manually set reads files and re-log configuration
                assembler.update_reads_files(files['r1'], files['r2'])

                # 运行组装|Run assembly
                success = assembler.run()
                results[sample_name] = success

                if success:
                    detection_logger.info(f"[OK] 样品 {sample_name} 组装成功|Sample {sample_name} assembly successful")
                else:
                    detection_logger.error(f"[FAIL] 样品 {sample_name} 组装失败|Sample {sample_name} assembly failed")

            except Exception as e:
                detection_logger.error(f"[FAIL] 样品 {sample_name} 处理异常|Sample {sample_name} processing error: {str(e)}")
                results[sample_name] = False

        # 输出总结|Output summary
        detection_logger.info("")
        detection_logger.info("=" * 60)
        detection_logger.info("批量处理总结|Batch Processing Summary")
        detection_logger.info("=" * 60)
        detection_logger.info(f"总样品数|Total samples: {len(samples)}")
        successful = sum(1 for v in results.values() if v)
        failed = len(results) - successful
        detection_logger.info(f"成功|Successful: {successful}")
        detection_logger.info(f"失败|Failed: {failed}")

        if successful == len(samples):
            detection_logger.info("[OK] 所有样品组装成功|All samples assembled successfully")
        elif successful > 0:
            detection_logger.warning(f"[WARNING] 部分样品组装失败|Some samples failed: {failed}")

        # 显示最终结果文件路径|Show final result file paths
        detection_logger.info("")
        detection_logger.info("=" * 60)
        detection_logger.info("最终叶绿体基因组文件|Final Plastome Genome Files")
        detection_logger.info("=" * 60)

        for sample_name in sorted(results.keys()):
            if results[sample_name]:
                # 查找最终文件|Find final file
                final_file = os.path.join(output_dir, sample_name, f"{sample_name}.plastome.fasta")
                if os.path.exists(final_file):
                    detection_logger.info(f"[OK] {sample_name}: {final_file}")
                else:
                    detection_logger.warning(f"[WARNING] {sample_name}: 未找到最终文件|Final file not found")
            else:
                detection_logger.error(f"[FAIL] {sample_name}: 组装失败|Assembly failed")

        detection_logger.info("")
        detection_logger.info(f"所有结果文件保存在|All result files saved in: {output_dir}")
        detection_logger.info("")

        return results

    def _log_configuration(self):
        """记录配置信息|Log configuration information"""
        self.logger.info("=" * 60)
        self.logger.info("叶绿体基因组组装|Plastome Assembly")
        self.logger.info("=" * 60)
        self.logger.info(f"输入目录|Input directory: {self.config.input_dir}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"输出前缀|Output prefix: {self.config.output_prefix}")
        self.logger.info(f"Organelle类型|Organelle type: {self.config.organelle_type}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        if self.config.r1_file:
            self.logger.info(f"R1文件|R1 file: {os.path.basename(self.config.r1_file)}")
        if self.config.r2_file:
            self.logger.info(f"R2文件|R2 file: {os.path.basename(self.config.r2_file)}")
        if self.config.unpaired_files:
            files = self.config.unpaired_files.split(',')
            self.logger.info(f"单端文件数|Single-end files: {len(files)}")
        self.logger.info("=" * 60)

    def _build_getorganelle_command(self) -> Optional[list]:
        """
        构建GetOrganelle命令|Build GetOrganelle command

        Returns:
            list: 命令列表，如果无法构建则返回None|Command list, or None if cannot build
        """
        command = [
            self.config.getorganelle_path,
            '-o', os.path.join(self.config.output_dir, self.config.output_prefix)
        ]

        # 添加配对reads|Add paired reads
        if self.config.r1_file and self.config.r2_file:
            command.extend(['-1', self.config.r1_file])
            command.extend(['-2', self.config.r2_file])
        # 添加单端reads|Add single-end reads
        elif self.config.unpaired_files:
            command.extend(['-u', self.config.unpaired_files])
        else:
            self.logger.error("没有可用的reads文件|No available reads files")
            return None

        # 添加参数|Add parameters
        command.extend(['-F', self.config.organelle_type])
        command.extend(['-R', str(self.config.max_rounds)])
        command.extend(['-k', self.config.kmer_list])
        command.extend(['-t', str(self.config.threads)])

        return command

    def run_assembly(self) -> bool:
        """
        运行叶绿体基因组组装|Run plastome assembly

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("开始叶绿体基因组组装|Starting plastome assembly")
        self.logger.info("=" * 60)
        self.logger.info("")

        try:
            # 构建命令|Build command
            command = self._build_getorganelle_command()
            if command is None:
                return False

            # 运行GetOrganelle|Run GetOrganelle
            success = self.cmd_runner.run_command(
                command,
                description="叶绿体基因组组装|Plastome assembly with GetOrganelle"
            )

            if not success:
                self.logger.error("叶绿体基因组组装失败|Plastome assembly failed")
                return False

            # 检查输出文件|Check output files
            self._check_output_files()

            # 后处理：整理最终结果文件|Post-processing: organize final result files
            final_fasta = self._organize_final_results()

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("叶绿体基因组组装完成|Plastome assembly completed successfully")
            self.logger.info("=" * 60)
            self.logger.info("")
            self.logger.info("最终结果文件|Final result file:")
            if final_fasta:
                self.logger.info(f"  [OK] {final_fasta}")
            self.logger.info("")
            self.logger.info(f"完整输出目录|Complete output directory: {os.path.join(self.config.output_dir, self.config.output_prefix)}")
            self.logger.info("")

            return True

        except Exception as e:
            self.logger.error(f"叶绿体基因组组装异常|Plastome assembly exception: {str(e)}")
            return False

    def _check_output_files(self):
        """检查并报告输出文件|Check and report output files"""
        output_dir = Path(self.config.output_dir) / self.config.output_prefix

        if not output_dir.exists():
            self.logger.warning(f"输出目录不存在|Output directory not found: {output_dir}")
            return

        # 使用通配符查找结果文件|Use wildcards to find result files
        fasta_files = sorted(output_dir.glob("*.path_sequence.fasta"))
        gfa_files = sorted(output_dir.glob("*.selected_graph.gfa"))
        log_files = sorted(output_dir.glob("*.log.txt"))

        self.logger.info("")
        self.logger.info("检查输出文件|Checking output files:")

        # 报告FASTA文件|Report FASTA files
        if fasta_files:
            self.logger.info(f"  [OK] 叶绿体基因组序列文件|Plastome sequence files: {len(fasta_files)}")
            for f in fasta_files:
                size = f.stat().st_size
                self.logger.info(f"    - {f.name} ({size:,} bytes)")
        else:
            self.logger.warning(f"  [FAIL] 未找到FASTA文件|No FASTA files found")

        # 报告GFA文件|Report GFA files
        if gfa_files:
            self.logger.info(f"  [OK] 组装图文件|Assembly graph files: {len(gfa_files)}")
            for f in gfa_files:
                size = f.stat().st_size
                self.logger.info(f"    - {f.name} ({size:,} bytes)")
        else:
            self.logger.warning(f"  [FAIL] 未找到GFA文件|No GFA files found")

        # 报告日志文件|Report log files
        if log_files:
            self.logger.info(f"  [OK] 日志文件|Log files: {len(log_files)}")
            for f in log_files:
                size = f.stat().st_size
                self.logger.info(f"    - {f.name} ({size:,} bytes)")
        else:
            self.logger.warning(f"  [FAIL] 未找到日志文件|No log files found")

    def _organize_final_results(self) -> Optional[str]:
        """
        整理最终结果文件，重命名为易于识别的名称，并格式化序列
        Organize final result files, rename to easy-to-recognize names, and format sequences

        Returns:
            str: 最终文件路径|Final file path
        """
        output_dir = Path(self.config.output_dir) / self.config.output_prefix

        if not output_dir.exists():
            self.logger.warning(f"输出目录不存在|Output directory not found: {output_dir}")
            return None

        # 查找所有 path_sequence.fasta 文件|Find all path_sequence.fasta files
        fasta_files = sorted(output_dir.glob("*.path_sequence.fasta"))

        if not fasta_files:
            self.logger.warning("未找到叶绿体基因组序列文件|No plastome sequence files found")
            return None

        # 通常 graph1.1 是主要结果，优先选择|Usually graph1.1 is the main result, prioritize it
        final_file = None
        for f in fasta_files:
            if 'graph1.1' in f.name:
                final_file = f
                break

        # 如果没有 graph1.1，使用第一个文件|If no graph1.1, use first file
        if not final_file:
            final_file = fasta_files[0]

        # 在样品目录内创建易于识别的文件名|Create easy-to-recognize filename within sample directory
        final_output = output_dir / f"{self.config.output_prefix}.plastome.fasta"
        final_gfa = None

        # 查找 GFA 文件|Find GFA file
        gfa_files = sorted(output_dir.glob("*.selected_graph.gfa"))
        if gfa_files:
            final_gfa = output_dir / f"{self.config.output_prefix}.graph.gfa"

        try:
            # 处理叶绿体基因组序列文件|Process plastome sequence file
            self._process_and_write_fasta(final_file, final_output)

            # 复制组装图文件|Copy assembly graph file
            if final_gfa and gfa_files and not final_gfa.exists():
                import shutil
                shutil.copy2(gfa_files[0], final_gfa)
                self.logger.info(f"已重命名组装图|Renamed assembly graph:")
                self.logger.info(f"  {gfa_files[0].name} -> {final_gfa.name}")

            return str(final_output)
        except Exception as e:
            self.logger.error(f"重命名文件失败|Failed to rename file: {str(e)}")
            return str(final_file)  # 返回原始文件路径|Return original file path

    def _process_and_write_fasta(self, input_file: Path, output_file: Path):
        """
        处理FASTA文件：替换header为样品名，格式化序列为每行60个核苷酸
        Process FASTA file: replace header with sample name, format sequence to 60 nt per line

        Args:
            input_file: 输入FASTA文件|Input FASTA file
            output_file: 输出FASTA文件|Output FASTA file
        """
        # 读取序列|Read sequence
        sequences = []
        current_header = None
        current_seq = []

        with open(input_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # 保存前一个序列|Save previous sequence
                    if current_header:
                        sequences.append((current_header, ''.join(current_seq)))
                    current_header = line[1:]  # 去掉 > 符号|Remove > symbol
                    current_seq = []
                else:
                    current_seq.append(line)

            # 保存最后一个序列|Save last sequence
            if current_header:
                sequences.append((current_header, ''.join(current_seq)))

        if not sequences:
            self.logger.warning(f"FASTA文件中没有序列|No sequences found in FASTA file: {input_file}")
            return

        # 写入处理后的序列|Write processed sequences
        with open(output_file, 'w') as f:
            for idx, (original_header, sequence) in enumerate(sequences, 1):
                # 使用样品名作为header|Use sample name as header
                # 如果有多个序列，添加数字后缀|Add numeric suffix if multiple sequences
                if len(sequences) > 1:
                    header = f"{self.config.output_prefix}_{idx}"
                else:
                    header = self.config.output_prefix

                f.write(f">{header}\n")

                # 格式化序列为每行60个字符|Format sequence to 60 characters per line
                for i in range(0, len(sequence), 60):
                    f.write(sequence[i:i+60] + '\n')

        # 统计信息|Statistics
        total_length = sum(len(seq) for _, seq in sequences)
        self.logger.info(f"已重命名叶绿体基因组序列|Renamed plastome sequence:")
        self.logger.info(f"  {input_file.name} -> {output_file.name}")
        self.logger.info(f"  序列数量|Sequence count: {len(sequences)}")
        self.logger.info(f"  总长度|Total length: {total_length:,} bp")

    def run(self) -> bool:
        """
        运行叶绿体基因组组装分析|Run plastome assembly analysis

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            return self.run_assembly()
        except Exception as e:
            self.logger.error(f"分析失败|Analysis failed: {str(e)}")
            return False


def main():
    """主函数入口|Main function entry point"""
    import argparse

    parser = argparse.ArgumentParser(
        description='叶绿体基因组组装自动化脚本|Plastome Assembly Automation Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i reads_folder -o plastome_output
        '''
    )

    # 必需参数|Required parameters
    required = parser.add_argument_group('必需参数|Required parameters')
    required.add_argument('-i', '--input', required=True,
                         help='输入目录(包含reads文件)|Input directory containing reads files')

    # 输出参数|Output parameters
    output_group = parser.add_argument_group('输出参数|Output parameters')
    output_group.add_argument('-o', '--output-dir', default='./plastome_output',
                             help='输出目录|Output directory')
    output_group.add_argument('-p', '--prefix',
                             help='输出前缀|Output prefix')

    # GetOrganelle参数|GetOrganelle parameters
    params_group = parser.add_argument_group('GetOrganelle参数|GetOrganelle parameters')
    params_group.add_argument('--organelle-type',
                              default='embplant_pt',
                              choices=['embplant_pt', 'embplant_mt', 'embplant_nr', 'other_pt',
                                      'animal_mt', 'fungus_mt', 'fungus_nr'],
                              help='Organelle类型|Organelle type')
    params_group.add_argument('-R', '--max-rounds', type=int, default=15,
                             help='最大扩展轮数|Maximum extension rounds')
    params_group.add_argument('-k', '--kmer-list', default='21,45,65,85,105',
                             help='Kmer列表(逗号分隔)|Kmer list comma-separated')
    params_group.add_argument('-t', '--threads', type=int, default=12,
                             help='线程数|Threads')

    # Reads文件识别参数|Reads file detection parameters
    reads_group = parser.add_argument_group('Reads文件识别|Reads file detection')
    reads_group.add_argument('--read1-suffix', default='_1.clean.fq.gz',
                             help='R1文件后缀模式|R1 file suffix pattern (default: %(default)s)')
    reads_group.add_argument('--read2-suffix', default='_2.clean.fq.gz',
                             help='R2文件后缀模式|R2 file suffix pattern (default: %(default)s)')

    # 软件配置|Software configuration
    software_group = parser.add_argument_group('软件配置|Software configuration')
    software_group.add_argument('--getorganelle-path',
                               default='~/miniforge3/envs/getorganelle_v.1.7.71/bin/get_organelle_from_reads.py',
                               help='GetOrganelle脚本路径|GetOrganelle script path')

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument('-v', '--verbose', action='store_true',
                          help='详细输出模式|Verbose output mode')
    log_group.add_argument('--log-file',
                          help='日志文件路径|Log file path')

    args = parser.parse_args()

    # 构建配置字典|Build configuration dictionary
    config = {
        'input_dir': args.input,
        'output_dir': args.output_dir,
        'output_prefix': args.prefix,
        'getorganelle_path': args.getorganelle_path,
        'organelle_type': args.organelle_type,
        'max_rounds': args.max_rounds,
        'kmer_list': args.kmer_list,
        'threads': args.threads,
        'read1_suffix': args.read1_suffix,
        'read2_suffix': args.read2_suffix,
    }

    # 创建组装器并运行|Create assembler and run
    try:
        assembler = PlastomeAssembler(**config)
        success = assembler.run()
        sys.exit(0 if success else 1)
    except Exception as e:
        print(f"错误|Error: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == '__main__':
    main()
