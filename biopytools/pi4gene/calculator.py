"""
Pi4Gene核心计算模块|Pi4Gene Core Calculation Module
功能: 序列提取、MAFFT比对、pi计算|
Features: Sequence extraction, MAFFT alignment, pi calculation
"""

import os
from itertools import combinations
from pathlib import Path
from typing import Dict, List, Tuple
from Bio import SeqIO

from .config import Pi4GeneConfig
from .utils import CommandRunner, build_conda_command, format_number


class Pi4GeneCalculator:
    """Pi4Gene计算器|Pi4Gene Calculator"""

    def __init__(self, config: Pi4GeneConfig, logger):
        self.config = config
        self.logger = logger
        self.runner = CommandRunner(logger)

    def extract_sequences_by_group(
        self, group_to_ids: Dict[str, List[str]]
    ) -> Dict[str, str]:
        """
        按分组提取序列并写入FASTA文件
        Extract sequences by group and write to FASTA files

        Args:
            group_to_ids: 分组到序列ID的映射|Group to sequence ID mapping

        Returns:
            {分组名: 输出FASTA路径}|{group_name: output_fasta_path}
        """
        self.logger.info("开始按分组提取序列|Starting sequence extraction by group")

        # 加载所有序列|Load all sequences
        seq_dict = SeqIO.to_dict(SeqIO.parse(self.config.input_file, "fasta"))
        self.logger.info(
            f"序列文件中共{len(seq_dict)}条序列|"
            f"Total {len(seq_dict)} sequences in input file"
        )

        group_fasta_map: Dict[str, str] = {}
        mafft_dir = self.config.output_path / '01_mafft'
        not_found_ids = []

        for group_name, seq_ids in sorted(group_to_ids.items()):
            output_fasta = mafft_dir / f"{group_name}.fasta"

            # 断点续传|Checkpoint resume
            if output_fasta.exists() and output_fasta.stat().st_size > 0:
                self.logger.info(
                    f"跳过已完成步骤|Skipping completed step: "
                    f"序列提取|sequence extraction for {group_name}"
                )
                group_fasta_map[group_name] = str(output_fasta)
                continue

            records = []
            for seq_id in seq_ids:
                if seq_id in seq_dict:
                    records.append(seq_dict[seq_id])
                else:
                    not_found_ids.append(seq_id)

            if not records:
                self.logger.warning(
                    f"分组{group_name}中没有找到任何有效序列，跳过|"
                    f"No valid sequences found for group {group_name}, skipping"
                )
                continue

            SeqIO.write(records, str(output_fasta), "fasta")
            group_fasta_map[group_name] = str(output_fasta)
            self.logger.info(
                f"  {group_name}: 提取{len(records)}条序列|"
                f"{group_name}: extracted {len(records)} sequences"
            )

        if not_found_ids:
            self.logger.warning(
                f"共{len(not_found_ids)}个序列ID未在序列文件中找到|"
                f"Total {len(not_found_ids)} sequence IDs not found in input file"
            )

        self.logger.info(
            f"序列提取完成，共{len(group_fasta_map)}个分组|"
            f"Sequence extraction completed for {len(group_fasta_map)} groups"
        )
        return group_fasta_map

    def run_mafft(self, group_fasta_map: Dict[str, str]) -> Dict[str, str]:
        """
        对每个分组的序列运行MAFFT比对
        Run MAFFT alignment for each group

        Args:
            group_fasta_map: {分组名: fasta文件路径}|{group_name: fasta_path}

        Returns:
            {分组名: 比对结果路径}|{group_name: aligned_fasta_path}
        """
        self.logger.info("开始MAFFT比对|Starting MAFFT alignment")

        aligned_map: Dict[str, str] = {}
        mafft_dir = self.config.output_path / '01_mafft'

        for group_name, fasta_path in sorted(group_fasta_map.items()):
            aligned_path = mafft_dir / f"{group_name}.aligned.fasta"

            # 断点续传|Checkpoint resume
            if aligned_path.exists() and aligned_path.stat().st_size > 0:
                self.logger.info(
                    f"跳过已完成步骤|Skipping completed step: "
                    f"MAFFT比对|alignment for {group_name}"
                )
                aligned_map[group_name] = str(aligned_path)
                continue

            # 检查序列数量|Check sequence count
            records = list(SeqIO.parse(fasta_path, "fasta"))
            n_seq = len(records)

            if n_seq < 2:
                self.logger.warning(
                    f"分组{group_name}只有{n}条序列，无法计算pi，跳过比对|"
                    f"Group {group_name} has only {n} sequence(s), "
                    f"cannot calculate pi, skipping alignment"
                )
                continue

            # 构建MAFFT命令|Build MAFFT command
            args = [
                '--auto',
                '--thread', str(self.config.threads),
                '--quiet',
                fasta_path
            ]

            cmd = build_conda_command(self.config.mafft_path, args)
            success = self.runner.run(
                cmd,
                f"MAFFT比对|MAFFT alignment for {group_name} ({n_seq} seqs)",
                stdout_file=str(aligned_path)
            )

            if success and aligned_path.exists():
                aligned_map[group_name] = str(aligned_path)
            else:
                self.logger.error(
                    f"分组{group_name}的MAFFT比对失败|"
                    f"MAFFT alignment failed for group {group_name}"
                )

        self.logger.info(
            f"MAFFT比对完成，共{len(aligned_map)}个分组|"
            f"MAFFT alignment completed for {len(aligned_map)} groups"
        )
        return aligned_map

    def calculate_pi(self, aligned_map: Dict[str, str]) -> List[Tuple[str, float, int]]:
        """
        从MAFFT比对结果直接计算pi (Nei & Li 1979)
        Calculate pi directly from MAFFT alignment (Nei & Li 1979)

        pi = 平均每对序列的核苷酸差异数 / 有效比对长度
        gap处理: 全gap位点跳过，部分gap的序列对在该位点不计入
        """
        self.logger.info("开始pi计算|Starting pi calculation")

        results: List[Tuple[str, float, int]] = []

        for group_name, aligned_path in sorted(aligned_map.items()):
            records = list(SeqIO.parse(aligned_path, "fasta"))
            n_seq = len(records)

            if n_seq < 2:
                self.logger.warning(
                    f"分组{group_name}只有{n_seq}条序列，跳过|"
                    f"Group {group_name} has only {n_seq} sequence(s), skipping"
                )
                continue

            seqs = [str(rec.seq).upper() for rec in records]
            aln_len = len(seqs[0])

            total_diffs = 0.0
            total_comparisons = 0

            for i, j in combinations(range(n_seq), 2):
                diffs = 0
                compared = 0
                for pos in range(aln_len):
                    b_i = seqs[i][pos]
                    b_j = seqs[j][pos]
                    if b_i in ('-', '.') or b_j in ('-', '.'):
                        continue
                    compared += 1
                    if b_i != b_j:
                        diffs += 1
                if compared > 0:
                    total_diffs += diffs / compared
                    total_comparisons += 1

            if total_comparisons > 0:
                pi_value = total_diffs / total_comparisons
                results.append((group_name, pi_value, n_seq))
                self.logger.info(
                    f"  {group_name}: pi={pi_value:.6f}, "
                    f"序列{n_seq}条, 比对长度{aln_len}bp"
                )

        self.logger.info(
            f"pi计算完成，共{len(results)}个分组|"
            f"Pi calculation completed for {len(results)} groups"
        )
        return results

    def write_results(self, results: List[Tuple[str, float, int]]) -> str:
        """
        写入pi计算结果到TSV文件
        Write pi results to TSV file

        Args:
            results: [(group_name, pi_value, n_seq), ...]

        Returns:
            输出文件路径|Output file path
        """
        output_file = self.config.output_path / 'pi_results.tsv'

        with open(output_file, 'w') as f:
            f.write("Group\tPi\tN_seq\n")
            for group_name, pi_value, n_seq in sorted(results):
                f.write(f"{group_name}\t{pi_value:.6f}\t{n_seq}\n")

        self.logger.info(f"结果已保存|Results saved: {output_file}")
        return str(output_file)
