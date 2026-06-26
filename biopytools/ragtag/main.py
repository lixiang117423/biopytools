"""
RagTag主程序模块|RagTag Main Module
"""

import os
import sys
from pathlib import Path
from .config import RagTagConfig
from .utils import RagTagLogger, CommandRunner, parse_fasta_header, rename_sequence_id, is_scaffolded_sequence


class RagTagScaffolder:
    """RagTag Scaffolder主类|RagTag Scaffolder Main Class"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        # 初始化配置|Initialize configuration
        self.config = RagTagConfig(**kwargs)
        self.config.validate()

        # 确保输出目录存在|Ensure output directory exists
        # (已在config的__post_init__中创建，这里再次确认|Already created in config's __post_init__, confirm here)
        os.makedirs(self.config.output_dir, exist_ok=True)

        # 初始化日志|Initialize logging
        log_file = os.path.join(self.config.output_dir, 'ragtag.log')
        self.logger_manager = RagTagLogger(log_file)
        self.logger = self.logger_manager.get_logger()

        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir)

        # 输出文件路径|Output file paths
        self.scaffolded_output = os.path.join(
            self.config.output_dir,
            f"{self.config.sample_name}_RagTag_scaffolded.fa"
        )
        self.unscaffolded_output = os.path.join(
            self.config.output_dir,
            f"{self.config.sample_name}_RagTag_unscaffolded.fa"
        )
        self.combined_output = os.path.join(
            self.config.output_dir,
            f"{self.config.sample_name}_RagTag_combined.fa"
        )

    def run(self):
        """运行scaffolding流程|Run scaffolding pipeline"""
        self.logger.info("开始RagTag scaffolding流程|Starting RagTag scaffolding pipeline")
        self.logger.info(f"参考基因组|Reference genome: {self.config.reference}")
        self.logger.info(f"查询基因组|Query genome: {self.config.query}")
        self.logger.info(f"样品名称|Sample name: {self.config.sample_name}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")

        try:
            # 步骤1: 运行RagTag scaffold|Step 1: Run RagTag scaffold
            self._run_ragtag_scaffold()

            # 步骤2: 处理输出文件，重命名序列ID并分类|Step 2: Process output files, rename sequence IDs and categorize
            self._process_output()

            self.logger.info("RagTag scaffolding流程完成|RagTag scaffolding pipeline completed")
            self.logger.info(f"Scaffolded输出|Scaffolded output: {self.scaffolded_output}")
            self.logger.info(f"Unscaffolded输出|Unscaffolded output: {self.unscaffolded_output}")
            self.logger.info(f"合并输出|Combined output: {self.combined_output}")

            return True

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {str(e)}")
            return False

    def _run_ragtag_scaffold(self):
        """运行RagTag scaffold命令|Run RagTag scaffold command"""
        self.logger.info("步骤1: 运行RagTag scaffold|Step 1: Running RagTag scaffold")

        result = self.cmd_runner.run_ragtag_scaffold(self.config)

        if result.returncode != 0:
            raise RuntimeError(f"RagTag scaffold执行失败|RagTag scaffold execution failed")

        # 检查输出文件|Check output files
        expected_output = os.path.join(self.config.output_dir, 'ragtag.scaffold.fasta')
        if not os.path.exists(expected_output):
            raise FileNotFoundError(f"RagTag输出文件未找到|RagTag output file not found: {expected_output}")

        self.logger.info("RagTag scaffold执行完成|RagTag scaffold execution completed")

    def _process_output(self):
        """处理输出文件，重命名序列ID并分类|Process output files, rename sequence IDs and categorize"""
        self.logger.info("步骤2: 处理输出文件，重命名序列ID并分类|Step 2: Processing output files, renaming sequence IDs and categorizing")

        if self.config.prefix:
            self.logger.info(f"使用序列ID前缀|Using sequence ID prefix: {self.config.prefix}")

        ragtag_output = os.path.join(self.config.output_dir, 'ragtag.scaffold.fasta')

        if not os.path.exists(ragtag_output):
            raise FileNotFoundError(f"RagTag输出文件不存在|RagTag output file not found: {ragtag_output}")

        # 读取并处理FASTA文件|Read and process FASTA file
        scaffolded_seqs = []
        unscaffolded_seqs = []
        unscaffolded_counter = 1

        with open(ragtag_output, 'r') as f:
            original_seq_id = None  # 保存原始序列ID|Save original sequence ID
            current_seq = []

            for line in f:
                line = line.strip()

                if line.startswith('>'):
                    # 保存前一个序列|Save previous sequence
                    if original_seq_id is not None:
                        # 使用原始ID判断是否scaffolded|Use original ID to check if scaffolded
                        if self._is_scaffolded(original_seq_id):
                            # Scaffolded序列：重命名并应用前缀|Scaffolded sequence: rename and apply prefix
                            final_seq_id = rename_sequence_id(original_seq_id, self.config.sample_name)
                            final_seq_id = self._apply_prefix(final_seq_id)
                            scaffolded_seqs.append((final_seq_id, ''.join(current_seq)))
                        else:
                            # Unscaffolded序列：重命名为scaffold_N并应用前缀|Unscaffolded sequence: rename to scaffold_N and apply prefix
                            final_seq_id = f"scaffold_{unscaffolded_counter}"
                            final_seq_id = self._apply_prefix(final_seq_id)
                            unscaffolded_seqs.append((final_seq_id, ''.join(current_seq)))
                            unscaffolded_counter += 1

                    # 解析新的序列头|Parse new sequence header
                    header = line[1:]  # 移除>符号|Remove > symbol
                    original_seq_id, description = parse_fasta_header(header)
                    current_seq = []

                else:
                    current_seq.append(line)

            # 保存最后一个序列|Save last sequence
            if original_seq_id is not None:
                # 使用原始ID判断是否scaffolded|Use original ID to check if scaffolded
                if self._is_scaffolded(original_seq_id):
                    # Scaffolded序列：重命名并应用前缀|Scaffolded sequence: rename and apply prefix
                    final_seq_id = rename_sequence_id(original_seq_id, self.config.sample_name)
                    final_seq_id = self._apply_prefix(final_seq_id)
                    scaffolded_seqs.append((final_seq_id, ''.join(current_seq)))
                else:
                    # Unscaffolded序列：重命名为scaffold_N并应用前缀|Unscaffolded sequence: rename to scaffold_N and apply prefix
                    final_seq_id = f"scaffold_{unscaffolded_counter}"
                    final_seq_id = self._apply_prefix(final_seq_id)
                    unscaffolded_seqs.append((final_seq_id, ''.join(current_seq)))
                    unscaffolded_counter += 1

        # 写入分类后的FASTA文件|Write categorized FASTA files
        self._write_fasta(self.scaffolded_output, scaffolded_seqs)
        self.logger.info(f"写入scaffolded序列|Wrote scaffolded sequences: {len(scaffolded_seqs)}条序列|sequences")

        self._write_fasta(self.unscaffolded_output, unscaffolded_seqs)
        self.logger.info(f"写入unscaffolded序列|Wrote unscaffolded sequences: {len(unscaffolded_seqs)}条序列|sequences")

        # 合并scaffolded和unscaffolded文件|Merge scaffolded and unscaffolded files
        self._merge_fasta_files(self.scaffolded_output, self.unscaffolded_output, self.combined_output)
        total_seqs = len(scaffolded_seqs) + len(unscaffolded_seqs)
        self.logger.info(f"写入合并序列文件|Wrote combined sequences file: {total_seqs}条序列|sequences")

        # 重命名RagTag原始输出文件|Rename RagTag original output files
        self._rename_ragtag_outputs()

    def _apply_prefix(self, seq_id):
        """
        应用前缀到序列ID|Apply prefix to sequence ID

        Args:
            seq_id: 原始序列ID|Original sequence ID

        Returns:
            str: 添加前缀后的序列ID|Sequence ID with prefix
        """
        if self.config.prefix:
            return f"{self.config.prefix}_{seq_id}"
        return seq_id

    def _is_scaffolded(self, seq_id):
        """
        判断序列是否已被scaffold|Check if sequence is scaffolded

        Args:
            seq_id: 序列ID|Sequence ID

        Returns:
            bool: 如果是scaffolded序列返回True，否则返回False
        """
        return is_scaffolded_sequence(seq_id)

    def _write_fasta(self, output_file, sequences):
        """
        写入FASTA文件|Write FASTA file

        Args:
            output_file: 输出文件路径|Output file path
            sequences: 序列列表|List of sequences, each is a tuple of (seq_id, sequence)
        """
        with open(output_file, 'w') as f:
            for seq_id, sequence in sequences:
                f.write(f">{seq_id}\n")
                # 每行写80个碱基|Write 80 bases per line
                for i in range(0, len(sequence), 80):
                    f.write(sequence[i:i+80] + '\n')

    def _merge_fasta_files(self, file1, file2, output_file):
        """
        合并两个FASTA文件|Merge two FASTA files

        将scaffolded和unscaffolded序列合并为一个完整的FASTA文件
        Merge scaffolded and unscaffolded sequences into one complete FASTA file

        Args:
            file1: 第一个FASTA文件路径|First FASTA file path
            file2: 第二个FASTA文件路径|Second FASTA file path
            output_file: 合并输出文件路径|Merged output file path
        """
        import shutil
        with open(output_file, 'w') as outfile:
            for input_file in (file1, file2):
                if os.path.exists(input_file):
                    with open(input_file, 'r') as infile:
                        shutil.copyfileobj(infile, outfile)
        self.logger.debug(f"合并FASTA文件完成|Merged FASTA files: {file1} + {file2} -> {output_file}")

    def _rename_ragtag_outputs(self):
        """重命名RagTag原始输出文件|Rename RagTag original output files"""
        self.logger.info("步骤3: 重命名RagTag原始输出文件|Step 3: Renaming RagTag original output files")

        # RagTag输出文件列表|RagTag output files list
        ragtag_files = [
            'ragtag.scaffold.fasta',
            'ragtag.scaffold.agp',
            'ragtag.scaffold.asm.paf',
            'ragtag.scaffold.asm.paf.log',
            'ragtag.scaffold.confidence.txt',
            'ragtag.scaffold.err',
            'ragtag.scaffold.stats',
            'ragtag.log'
        ]

        renamed_count = 0
        for filename in ragtag_files:
            old_path = os.path.join(self.config.output_dir, filename)
            new_filename = f"{self.config.sample_name}_{filename}"
            new_path = os.path.join(self.config.output_dir, new_filename)

            if os.path.exists(old_path):
                try:
                    os.rename(old_path, new_path)
                    self.logger.debug(f"重命名|Renamed: {filename} -> {new_filename}")
                    renamed_count += 1
                except Exception as e:
                    self.logger.warning(f"重命名失败|Failed to rename {filename}: {str(e)}")

        self.logger.info(f"完成重命名RagTag输出文件|Completed renaming RagTag output files: {renamed_count}个文件|files")


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='RagTag基因组scaffolding工具|RagTag Genome Scaffolding Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-r', '--reference', required=True,
                       help='参考基因组FASTA文件|Reference genome FASTA file')
    parser.add_argument('-q', '--query', required=True,
                       help='查询基因组FASTA文件|Query genome FASTA file')
    parser.add_argument('-s', '--sample-name', required=True,
                       help='样品名称，用于输出文件命名|Sample name for output file naming')

    # 可选参数|Optional parameters
    parser.add_argument('-p', '--prefix', default=None,
                       help='序列ID前缀，会添加到所有输出序列的ID前面|Sequence ID prefix to add to all output sequence IDs')
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')
    parser.add_argument('-o', '--output-dir', default='./ragtag_output',
                       help='输出目录|Output directory')

    # RagTag选项|RagTag options
    parser.add_argument('--aligner', default='minimap2',
                       choices=['minimap2', 'unimap', 'nucmer'],
                       help='比对器|Aligner to use')
    parser.add_argument('-C', '--concatenate-unplaced', action='store_true',
                       help='将未定位的contigs合并为chr0|Concatenate unplaced contigs into chr0')
    parser.add_argument('-R', '--infer-gaps', action='store_true',
                       help='推断gap大小|Infer gap sizes')

    args = parser.parse_args()

    # 创建并运行scaffolder|Create and run scaffolder
    scaffolder = RagTagScaffolder(
        reference=args.reference,
        query=args.query,
        sample_name=args.sample_name,
        prefix=args.prefix,
        threads=args.threads,
        output_dir=args.output_dir,
        aligner=args.aligner,
        concatenate_unplaced=args.concatenate_unplaced,
        infer_gaps=args.infer_gaps
    )

    success = scaffolder.run()

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
