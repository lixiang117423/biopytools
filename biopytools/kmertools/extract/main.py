"""Kmer提取器|Kmer Extractor

从FASTA文件中提取kmer并生成canonical kmer
Extract kmers from FASTA file and generate canonical kmers

支持两种提取方法：
Supports two extraction methods:
1. unikmer：使用unikmer工具（推荐，速度快）
2. pyfastx：使用pyfastx库（纯Python实现）
"""

import os
import subprocess
import pyfastx
from typing import Dict, Optional


class KmerExtractor:
    """Kmer提取器|Kmer Extractor

    从FASTA文件提取kmer，支持canonical kmer和位置信息
    Extract kmers from FASTA file, supports canonical kmer and position information
    """

    def __init__(self, config, logger):
        """初始化Kmer提取器|Initialize Kmer extractor

        Args:
            config: KmerPAVConfig配置对象|KmerPAVConfig object
            logger: 日志器|Logger instance
        """
        self.config = config
        self.logger = logger
        self.complement_dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'N': 'N'}

    def get_reverse_complement(self, sequence: str) -> str:
        """获取反向互补序列|Get reverse complement sequence

        Args:
            sequence: DNA序列|DNA sequence

        Returns:
            str: 反向互补序列|Reverse complement sequence
        """
        sequence = sequence.upper()
        reverse_sequence = sequence[::-1]
        complement = ''.join(self.complement_dict.get(base, base) for base in reverse_sequence)
        return complement

    def get_canonical_kmer(self, kmer: str) -> str:
        """获取canonical kmer|Get canonical kmer

        取kmer和其反向互补序列中字典序较小的作为canonical kmer
        Take the lexicographically smaller one between kmer and its reverse complement as canonical kmer

        Args:
            kmer: kmer序列|Kmer sequence

        Returns:
            str: canonical kmer|Canonical kmer
        """
        rc_kmer = self.get_reverse_complement(kmer)
        return min(kmer, rc_kmer)

    def get_kmers_with_positions(self, sequence: str, k: int) -> Dict[str, str]:
        """获取kmer及其位置|Get kmers with positions

        Args:
            sequence: DNA序列|DNA sequence
            k: kmer大小|Kmer size

        Returns:
            Dict[str, str]: kmer到位置的映射|Kmer to position mapping
        """
        sequence = sequence.upper()
        kmers = {}

        for i in range(len(sequence) - k + 1):
            kmer = sequence[i:i+k]
            if kmer not in kmers:
                kmers[kmer] = str(i + 1)  # 1-based position

        return kmers

    def extract_kmers_from_fasta(self) -> bool:
        """从FASTA文件提取kmer|Extract kmers from FASTA file

        根据配置选择提取方法：
        - unikmer：使用unikmer工具（默认，速度快，支持BED输出）
        - pyfastx：使用pyfastx库（纯Python，无需外部工具）

        Returns:
            bool: 成功返回True|Return True if successful
        """
        method = self.config.extract_method

        if method == "unikmer":
            self.logger.info(f"使用unikmer方法提取kmer|Using unikmer method for extraction")
            return self.extract_kmers_with_unikmer()
        elif method == "pyfastx":
            self.logger.info(f"使用pyfastx方法提取kmer|Using pyfastx method for extraction")
            return self.extract_kmers_with_pyfastx()
        else:
            self.logger.error(f"不支持的提取方法|Unsupported extraction method: {method}")
            return False

    def extract_kmers_with_unikmer(self) -> bool:
        """使用unikmer工具从FASTA文件提取kmer|Extract kmers from FASTA using unikmer

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info(f"从FASTA文件提取kmer (unikmer)|Extracting kmers from FASTA using unikmer: {self.config.fasta_file}")

        fasta_file = self.config.fasta_file
        output_dir = self.config.output_dir
        kmer_size = self.config.kmer_size
        threads = self.config.threads
        unikmer_path = self.config.unikmer_path

        # 准备输出文件路径|Prepare output file paths
        base_name = os.path.basename(fasta_file).rsplit('.', 1)[0]
        unik_output_prefix = os.path.join(output_dir, base_name)
        unik_file = f"{unik_output_prefix}.unik"

        # 构建unikmer count命令|Build unikmer count command
        cmd_parts = [
            unikmer_path,
            "count",
            "-k", str(kmer_size),
            "-j", str(threads),
            "--canonical",
            "-o", base_name,  # unikmer会自动添加.unik后缀
            fasta_file
        ]

        cmd = " ".join(cmd_parts)
        self.logger.info(f"执行unikmer命令|Running unikmer command: {cmd}")

        try:
            # 在输出目录中执行命令|Run command in output directory
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=output_dir,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"unikmer执行失败|unikmer execution failed: {result.stderr}")
                return False

            self.logger.info(f"unikmer执行成功|unikmer execution completed: {unik_file}")

            # 如果需要，生成BED文件|Generate BED file if requested
            if self.config.extract_output_bed:
                bed_output = f"{unik_output_prefix}.bed"
                success = self.generate_bed_with_unikmer(unik_file, fasta_file, bed_output)
                if not success:
                    self.logger.warning(f"BED文件生成失败|BED file generation failed")
                    return False
                self.logger.info(f"BED文件已生成|BED file generated: {bed_output}")

            # 如果需要，生成kmer列表文件（用于向后兼容）|Generate kmer list file (for backward compatibility)
            kmer_list_output = self.config.kmer_output
            if kmer_list_output and self.config.extract_output_bed:
                # 从BED文件生成FASTA格式的kmer列表
                bed_output = f"{unik_output_prefix}.bed"
                success = self.generate_fasta_from_bed(bed_output, kmer_list_output)
                if success:
                    self.logger.info(f"Kmer FASTA文件已生成|Kmer FASTA file generated: {kmer_list_output}")

            self.logger.info(f"kmer提取完成|Kmer extraction completed")
            return True

        except Exception as e:
            self.logger.error(f"kmer提取失败|Kmer extraction failed: {e}")
            return False

    def generate_bed_with_unikmer(self, unik_file: str, fasta_file: str, bed_output: str) -> bool:
        """使用unikmer locate生成BED文件|Generate BED file using unikmer locate

        Args:
            unik_file: unikmer输出文件|unikmer output file
            fasta_file: 输入FASTA文件|Input FASTA file
            bed_output: 输出BED文件|Output BED file

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info(f"生成BED文件|Generating BED file: {bed_output}")

        unikmer_path = self.config.unikmer_path
        output_dir = self.config.output_dir

        # 构建unikmer locate命令|Build unikmer locate command
        # 使用相对路径，因为工作目录已设为output_dir
        unik_basename = os.path.basename(unik_file)
        temp_bed = os.path.join(output_dir, "temp_locate")
        temp_bed_basename = os.path.basename(temp_bed)

        cmd_parts = [
            unikmer_path,
            "locate",
            unik_basename,
            "-g", fasta_file,
            "-o", temp_bed_basename
        ]

        cmd = " ".join(cmd_parts)

        try:
            result = subprocess.run(
                cmd,
                shell=True,
                cwd=output_dir,
                capture_output=True,
                text=True
            )

            if result.returncode != 0:
                self.logger.error(f"unikmer locate执行失败|unikmer locate failed: {result.stderr}")
                return False

            # 查找实际生成的BED文件|Find the actual generated BED file
            possible_files = [
                temp_bed,
                f"{temp_bed}.bed",
                f"{temp_bed}.txt",
                os.path.join(output_dir, os.path.basename(fasta_file).rsplit('.', 1)[0])
            ]

            actual_bed_file = None
            for possible_file in possible_files:
                if os.path.exists(possible_file):
                    actual_bed_file = possible_file
                    break

            if not actual_bed_file:
                self.logger.error(f"找不到unikmer locate输出文件|Cannot find unikmer locate output file")
                return False

            # 如果需要，重命名文件|Rename file if needed
            if actual_bed_file != bed_output:
                import shutil
                shutil.move(actual_bed_file, bed_output)

            # 对BED文件排序|Sort BED file
            self._sort_bed_file(bed_output)

            self.logger.info(f"BED文件生成完成|BED file generation completed")
            return True

        except Exception as e:
            self.logger.error(f"BED文件生成失败|BED file generation failed: {e}")
            return False

    def _sort_bed_file(self, bed_file: str) -> bool:
        """对BED文件按位置排序|Sort BED file by position"""
        try:
            with open(bed_file, 'r') as f:
                lines = f.readlines()

            bed_entries = []
            for line in lines:
                if line.strip():
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        try:
                            seq_name = parts[0]
                            start_pos = int(parts[1])
                            bed_entries.append((seq_name, start_pos, line.strip()))
                        except ValueError:
                            continue

            # 按序列名称和起始位置排序|Sort by sequence name and start position
            bed_entries.sort(key=lambda x: (x[0], x[1]))

            # 写回文件|Write back to file
            with open(bed_file, 'w') as f:
                for _, _, line in bed_entries:
                    f.write(f"{line}\n")

            return True

        except Exception as e:
            self.logger.warning(f"BED文件排序失败|BED file sorting failed: {e}")
            return False

    def generate_fasta_from_bed(self, bed_file: str, fasta_output: str) -> bool:
        """从BED文件生成FASTA格式的kmer文件|Generate FASTA format kmer file from BED file

        使用与kmer_query模块相同的逻辑，保留所有位置的kmer（不去重）
        Use the same logic as kmer_query module, keep all kmers at all positions (no deduplication)

        Args:
            bed_file: 输入BED文件|Input BED file
            fasta_output: 输出FASTA文件|Output FASTA file

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info(f"从BED文件生成FASTA|Generating FASTA from BED: {bed_file}")

        try:
            # 读取BED文件并生成FASTA条目|Read BED file and generate FASTA entries
            fasta_entries = []

            with open(bed_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 4:
                            seq_name = parts[0]
                            start_pos = int(parts[1]) + 1  # BED格式是0-based，转为1-based
                            end_pos = int(parts[2])
                            kmer_seq = parts[3]

                            # 生成FASTA ID：>seq_name_start_end
                            # Generate FASTA ID: >seq_name_start_end
                            fasta_id = f">{seq_name}_{start_pos}_{end_pos}"
                            fasta_entries.append((seq_name, start_pos, fasta_id, kmer_seq))

            # 按序列名称和位置排序|Sort by sequence name and position
            fasta_entries.sort(key=lambda x: (x[0], x[1]))

            # 写入FASTA文件|Write FASTA file
            with open(fasta_output, 'w') as f:
                for _, _, fasta_id, kmer_seq in fasta_entries:
                    f.write(f"{fasta_id}\n{kmer_seq}\n")

            self.logger.info(f"FASTA文件生成完成|FASTA file generation completed: {len(fasta_entries)} k-mers")
            return True

        except Exception as e:
            self.logger.error(f"FASTA文件生成失败|FASTA file generation failed: {e}")
            return False

    def extract_kmers_with_pyfastx(self) -> bool:
        """使用pyfastx从FASTA文件提取kmer|Extract kmers from FASTA using pyfastx

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info(f"从FASTA文件提取kmer (pyfastx)|Extracting kmers from FASTA using pyfastx: {self.config.fasta_file}")

        kmer_size = self.config.kmer_size
        kmer_output = self.config.kmer_output
        kmer_pos_output = self.config.kmer_pos_output

        # 打开输出文件|Open output files
        try:
            outfile = open(kmer_output, 'w')
            outfile.write('kmer\n')
            outfile2 = open(kmer_pos_output, 'w')
            outfile2.write('kmer\tsample\tPos\n')
        except IOError as e:
            self.logger.error(f"无法创建输出文件|Cannot create output files: {e}")
            return False

        # 统计|Statistics
        total_kmers = 0
        unique_kmers = set()

        try:
            # 遍历FASTA序列|Iterate over FASTA sequences
            for name, seq in pyfastx.Fasta(self.config.fasta_file, build_index=False):
                self.logger.debug(f"处理序列|Processing sequence: {name}")

                kmer_positions = self.get_kmers_with_positions(seq, kmer_size)

                for kmer, positions in kmer_positions.items():
                    canonical_kmer = self.get_canonical_kmer(kmer)
                    outfile2.write(f"{canonical_kmer}\t{name}\t{positions}\n")
                    outfile.write(f"{canonical_kmer}\n")

                    unique_kmers.add(canonical_kmer)
                    total_kmers += 1

            outfile.close()
            outfile2.close()

            self.logger.info(f"kmer提取完成|Kmer extraction completed")
            self.logger.info(f"总kmer数|Total kmers: {total_kmers}")
            self.logger.info(f"唯一kmer数|Unique kmers: {len(unique_kmers)}")
            self.logger.info(f"kmer列表|Kmer list: {kmer_output}")
            self.logger.info(f"kmer位置|Kmer positions: {kmer_pos_output}")

            return True

        except Exception as e:
            self.logger.error(f"kmer提取失败|Kmer extraction failed: {e}")
            outfile.close()
            outfile2.close()
            return False
