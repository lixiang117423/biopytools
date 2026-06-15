"""
Ka/Ks Calculator增强版输入验证模块|Ka/Ks Calculator Enhanced Input Validation Module
功能: 详细的错误提示和用户友好的验证信息|Detailed error messages and user-friendly validation
"""

import re
import pandas as pd
from typing import Dict, List, Tuple, Set
from Bio import SeqIO
from .config import KaKsConfig
from .utils import KaKsLogger


class SequenceValidator:
    """序列验证工具类|Sequence validation utility"""

    def __init__(self, logger: KaKsLogger):
        """
        初始化验证器|Initialize validator

        Args:
            logger: 日志器实例|Logger instance
        """
        self.logger = logger
        self.config = KaKsConfig()
        self.validation_stats = {
            'character_issues': 0,
            'n_ratio_issues': 0,
            'length_issues': 0,
            'codon_issues': 0,
            'total_processed': 0,
            'total_valid': 0,
            'failed_sequences': []
        }

    def validate_fasta_file(self, fasta_file: str) -> Tuple[bool, Dict[str, str]]:
        """
        验证FASTA文件格式和内容|Validate FASTA file format and content

        Args:
            fasta_file: FASTA文件路径|FASTA file path

        Returns:
            (is_valid, sequences_dict): 验证结果和序列字典|Validation result and sequences dictionary
        """
        try:
            self.logger.info(f"验证FASTA文件|Validating FASTA file: {fasta_file}")

            sequences = {}
            self._reset_validation_stats()

            for record in SeqIO.parse(fasta_file, "fasta"):
                self.validation_stats['total_processed'] += 1
                seq_str = str(record.seq).upper()

                is_valid, processed_seq = self._detailed_sequence_validation(record.id, seq_str)

                if is_valid:
                    sequences[record.id] = processed_seq
                    self.validation_stats['total_valid'] += 1
                else:
                    self.validation_stats.setdefault('total_invalid', 0)
                    self.validation_stats['total_invalid'] += 1

            self._display_validation_summary(fasta_file)

            if len(sequences) == 0:
                self.logger.error("文件中没有有效序列，请根据上述建议修复后重试|No valid sequences in file, please fix and retry")
                return False, {}

            return True, sequences

        except Exception as e:
            self.logger.error(f"FASTA文件验证失败|FASTA validation failed: {fasta_file} - {e}")
            return False, {}

    def _detailed_sequence_validation(self, seq_id: str, seq_str: str) -> Tuple[bool, str]:
        """
        详细的单序列验证|Detailed single sequence validation

        Args:
            seq_id: 序列ID|Sequence ID
            seq_str: 序列字符串（已转大写）|Sequence string (already uppercased)

        Returns:
            (is_valid, processed_sequence): 是否有效和处理后的序列|Whether valid and processed sequence
        """
        original_seq = seq_str
        issues_found = []

        # 1. 字符验证|Character validation
        char_valid, char_issues = self._validate_characters(seq_str)
        if not char_valid:
            self.validation_stats['character_issues'] += 1
            issues_found.extend(char_issues)
            self._log_character_validation_failure(seq_id, char_issues)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'invalid_characters', 'details': char_issues
            })
            return False, ""

        # 2. N字符比例验证|N character ratio validation
        n_valid, n_ratio = self._validate_n_ratio(seq_str)
        if not n_valid:
            self.validation_stats['n_ratio_issues'] += 1
            self._log_n_ratio_validation_failure(seq_id, n_ratio)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'high_n_ratio', 'details': f'{n_ratio:.2%}'
            })
            return False, ""

        # 3. 长度验证|Length validation
        length_valid, length_issue = self._validate_length(seq_str)
        if not length_valid:
            self.validation_stats['length_issues'] += 1
            self._log_length_validation_failure(seq_id, len(seq_str), length_issue)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'invalid_length', 'details': f'{len(seq_str)}bp'
            })
            return False, ""

        # 4. 密码子完整性验证|Codon completeness validation
        codon_valid, remainder = self._validate_codon_completeness(seq_str)
        if not codon_valid:
            self.validation_stats['codon_issues'] += 1
            self._log_codon_validation_failure(seq_id, len(seq_str), remainder)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'incomplete_codon', 'details': f'余数|remainder {remainder}'
            })
            return False, ""

        # 5. 最终处理|Final processing
        final_seq = seq_str.replace('N', '')

        return True, final_seq

    def _validate_characters(self, seq_str: str) -> Tuple[bool, List[str]]:
        """验证字符并返回详细问题|Validate characters and return detailed issues"""
        if re.match(r'^[ATCGN]+$', seq_str):
            return True, []

        invalid_chars = set(char for char in seq_str if char not in 'ATCGN')
        issues = []

        lowercase_chars = set(char for char in invalid_chars if char.lower() in 'atcg')
        iupac_chars = set(char for char in invalid_chars if char in 'RYSWKMDBHV')
        other_chars = invalid_chars - lowercase_chars - iupac_chars

        if lowercase_chars:
            issues.append(f"包含小写核苷酸|Contains lowercase nucleotides: {sorted(lowercase_chars)}")
        if iupac_chars:
            issues.append(f"包含IUPAC模糊代码|Contains IUPAC ambiguity codes: {sorted(iupac_chars)}")
        if other_chars:
            issues.append(f"包含其他无效字符|Contains other invalid characters: {sorted(other_chars)}")

        return False, issues

    def _validate_n_ratio(self, seq_str: str) -> Tuple[bool, float]:
        """验证N字符比例|Validate N character ratio"""
        n_ratio = seq_str.count('N') / len(seq_str) if len(seq_str) > 0 else 1
        return n_ratio <= self.config.max_missing_ratio, n_ratio

    def _validate_length(self, seq_str: str) -> Tuple[bool, str]:
        """验证序列长度|Validate sequence length"""
        length = len(seq_str)
        if length < self.config.min_sequence_length:
            return False, "too_short"
        elif length > self.config.max_sequence_length:
            return False, "too_long"
        return True, "valid"

    def _validate_codon_completeness(self, seq_str: str) -> Tuple[bool, int]:
        """验证密码子完整性|Validate codon completeness"""
        remainder = len(seq_str) % 3
        return remainder == 0, remainder

    def _log_character_validation_failure(self, seq_id: str, issues: List[str]):
        """记录字符验证失败的详细信息|Log character validation failure details"""
        self.logger.error(f"       字符验证失败|Character validation failed!")
        self.logger.error(f"         序列ID|Sequence ID: {seq_id}")
        for issue in issues:
            self.logger.error(f"         问题|Issue: {issue}")

        self.logger.info("       修复建议|Fix suggestions:")
        if any("小写" in issue or "lowercase" in issue for issue in issues):
            self.logger.info("         - 转换小写字母为大写|Convert lowercase to uppercase: seq = seq.upper()")
        if any("IUPAC" in issue for issue in issues):
            self.logger.info("         - 替换IUPAC模糊代码为标准字符或N|Replace IUPAC codes with standard chars or N")
            self.logger.info("           例如|Example: R->A, Y->C, S->G, W->A, K->G, M->A")
        if any("其他" in issue or "other" in issue for issue in issues):
            self.logger.info("         - 移除或替换非核苷酸字符|Remove or replace non-nucleotide characters")
            self.logger.info("         - 确保序列只包含 A, T, C, G, N 字符|Ensure only A, T, C, G, N characters")

    def _log_n_ratio_validation_failure(self, seq_id: str, n_ratio: float):
        """记录N字符比例验证失败的详细信息|Log N ratio validation failure details"""
        self.logger.error(f"       N字符比例过高|N character ratio too high!")
        self.logger.error(f"         序列ID|Sequence ID: {seq_id}")
        self.logger.error(f"         N字符比例|N ratio: {n_ratio:.2%} (阈值|threshold: <={self.config.max_missing_ratio:.1%})")

        self.logger.info("       修复建议|Fix suggestions:")
        self.logger.info("         - 检查序列质量，可能包含低质量区域|Check sequence quality, may contain low-quality regions")
        self.logger.info("         - 使用序列清理工具移除低质量片段|Use sequence cleaning tools to remove low-quality segments")
        self.logger.info("         - 如果是组装序列，检查gap填充质量|If assembly sequence, check gap filling quality")
        self.logger.info("         - 考虑从序列两端修剪N字符|Consider trimming N characters from both ends")

    def _log_length_validation_failure(self, seq_id: str, length: int, issue_type: str):
        """记录长度验证失败的详细信息|Log length validation failure details"""
        self.logger.error(f"       序列长度不合格|Sequence length invalid!")
        self.logger.error(f"         序列ID|Sequence ID: {seq_id}")
        self.logger.error(f"         当前长度|Current length: {length} bp")

        if issue_type == "too_short":
            self.logger.error(f"         要求最小长度|Required minimum length: {self.config.min_sequence_length} bp")
            self.logger.info("       修复建议|Fix suggestions:")
            self.logger.info("         - 检查序列是否完整（可能缺失部分片段）|Check if sequence is complete (may be missing fragments)")
            self.logger.info("         - 确认这是完整的CDS序列|Confirm this is a complete CDS sequence")
            self.logger.info("         - 如果是部分序列，考虑获取完整序列|If partial sequence, consider obtaining complete sequence")
        else:
            self.logger.error(f"         允许最大长度|Allowed maximum length: {self.config.max_sequence_length} bp")
            self.logger.info("       修复建议|Fix suggestions:")
            self.logger.info("         - 确认这是CDS序列而非完整基因或基因组片段|Confirm this is CDS, not full gene or genome fragment")
            self.logger.info("         - 移除内含子、UTR等非编码区域|Remove introns, UTRs and other non-coding regions")
            self.logger.info("         - 如果确实是超长CDS，可调整max_sequence_length参数|If truly long CDS, adjust max_sequence_length parameter")

    def _log_codon_validation_failure(self, seq_id: str, length: int, remainder: int):
        """记录密码子验证失败的详细信息|Log codon validation failure details"""
        self.logger.error(f"       密码子不完整|Incomplete codon!")
        self.logger.error(f"         序列ID|Sequence ID: {seq_id}")
        self.logger.error(f"         序列长度|Sequence length: {length} bp")
        self.logger.error(f"         余数|Remainder: {remainder} bp (长度必须是3的倍数|length must be multiple of 3)")

        self.logger.info("       修复建议|Fix suggestions:")
        if remainder == 1:
            self.logger.info("         - 序列末尾多了1个核苷酸，建议移除最后1个字符|Sequence has 1 extra nucleotide at end, remove last character")
            self.logger.info("         - 检查序列是否包含终止密码子后的额外字符|Check if sequence includes characters after stop codon")
        elif remainder == 2:
            self.logger.info("         - 序列末尾多了2个核苷酸，建议移除最后2个字符|Sequence has 2 extra nucleotides at end, remove last 2 characters")
            self.logger.info("         - 检查序列起始位置是否正确|Check if sequence start position is correct")

        self.logger.info("         - 自动修复|Auto-fix: seq = seq[:len(seq)//3*3]")
        self.logger.info("         - 手动检查序列边界，确保正确的起始和终止位置|Manually check sequence boundaries for correct start/end positions")

    def _display_validation_summary(self, fasta_file: str):
        """显示验证汇总信息和建议|Display validation summary and suggestions"""
        stats = self.validation_stats

        self.logger.separator("验证结果汇总|Validation Summary")
        self.logger.info(f"文件|File: {fasta_file}")
        self.logger.info(f"总序列数|Total sequences: {stats['total_processed']}")
        self.logger.info(f"有效序列|Valid sequences: {stats['total_valid']}")
        self.logger.info(f"无效序列|Invalid sequences: {stats['total_processed'] - stats['total_valid']}")

        if stats['total_valid'] > 0:
            success_rate = (stats['total_valid'] / stats['total_processed']) * 100
            self.logger.success(f"成功率|Success rate: {success_rate:.1f}%")

        if stats['total_processed'] > stats['total_valid']:
            self.logger.warning("验证失败统计|Validation failure statistics:")

            if stats['character_issues'] > 0:
                self.logger.warning(f"   字符问题|Character issues: {stats['character_issues']} 个序列|sequences")
            if stats['n_ratio_issues'] > 0:
                self.logger.warning(f"   N字符过多|High N ratio: {stats['n_ratio_issues']} 个序列|sequences")
            if stats['length_issues'] > 0:
                self.logger.warning(f"   长度问题|Length issues: {stats['length_issues']} 个序列|sequences")
            if stats['codon_issues'] > 0:
                self.logger.warning(f"   密码子问题|Codon issues: {stats['codon_issues']} 个序列|sequences")

            self.logger.info("失败序列详情|Failed sequence details:")
            for failed_seq in stats['failed_sequences'][:5]:
                self.logger.info(f"   - {failed_seq['id']}: {failed_seq['reason']} ({failed_seq['details']})")

            if len(stats['failed_sequences']) > 5:
                remaining = len(stats['failed_sequences']) - 5
                self.logger.info(f"   ... 还有 {remaining} 个失败序列|... {remaining} more failed sequences")

            self._provide_overall_suggestions(stats)

    def _provide_overall_suggestions(self, stats: Dict):
        """提供总体修复建议|Provide overall fix suggestions"""
        self.logger.info("总体修复建议|Overall fix suggestions:")

        if stats['character_issues'] > 0:
            self.logger.info("字符问题修复|Character issue fixes:")
            self.logger.info("   1. 批量转换为大写|Batch convert to uppercase: awk '{if(/^>/){print}else{print toupper($0)}}' input.fasta > output.fasta")
            self.logger.info("   2. 替换IUPAC字符|Replace IUPAC chars: sed 's/[RYSWKMDBHV]/N/g' input.fasta > output.fasta")

        if stats['n_ratio_issues'] > 0:
            self.logger.info("N字符问题修复|N character issue fixes:")
            self.logger.info("   1. 移除首尾N|Remove leading/trailing N: use seqkit tool")
            self.logger.info("   2. 过滤高N比例序列|Filter high N ratio sequences: use bioinformatics QC tools")

        if stats['length_issues'] > 0:
            self.logger.info("长度问题修复|Length issue fixes:")
            self.logger.info("   1. 确认序列来源和完整性|Confirm sequence source and completeness")
            self.logger.info("   2. 调整配置文件中的长度阈值|Adjust length thresholds in config")

        if stats['codon_issues'] > 0:
            self.logger.info("密码子问题修复|Codon issue fixes:")
            self.logger.info("   1. 自动修剪|Auto-trim: keep length as multiple of 3")
            self.logger.info("   2. 检查序列注释和起止位置|Check sequence annotation and start/end positions")

        self.logger.info("推荐工具|Recommended tools:")
        self.logger.info("   - SeqKit: 序列处理瑞士军刀|Sequence processing swiss army knife")
        self.logger.info("   - FASTX-Toolkit: 序列质控工具|Sequence QC tool")

    def _reset_validation_stats(self):
        """重置验证统计|Reset validation statistics"""
        self.validation_stats = {
            'character_issues': 0,
            'n_ratio_issues': 0,
            'length_issues': 0,
            'codon_issues': 0,
            'total_processed': 0,
            'total_valid': 0,
            'failed_sequences': []
        }

    def parse_pair_file(self, pair_file: str) -> List[Tuple[str, str, str]]:
        """
        解析配对文件（不验证ID是否存在）|Parse pair file without validating IDs

        Args:
            pair_file: 配对文件路径|Pair file path

        Returns:
            配对列表|Pairs list: [(seq1_id, seq2_id, pair_name), ...]
        """
        self.logger.info(f"解析配对文件|Parsing pair file: {pair_file}")

        if pair_file.endswith('.csv'):
            df = pd.read_csv(pair_file, comment='#', header=None)
        else:
            df = pd.read_csv(pair_file, sep='\t', comment='#', header=None)

        if len(df.columns) < 2:
            raise ValueError("配对文件至少需要2列|Pair file must have at least 2 columns")

        id1_col, id2_col = df.columns[0], df.columns[1]
        name_col = df.columns[2] if len(df.columns) > 2 else None

        pairs = []
        for idx, row in df.iterrows():
            seq1_id = str(row[id1_col]).strip()
            seq2_id = str(row[id2_col]).strip()
            pair_name = str(row[name_col]).strip() if name_col and pd.notna(row[name_col]) else f"pair_{idx+1:04d}"
            pairs.append((seq1_id, seq2_id, pair_name))

        self.logger.info(f"解析完成|Parsing complete: {len(pairs)} 对|pairs")
        return pairs

    def validate_pair_file(self, pair_file: str, seq1_ids: Set[str], seq2_ids: Set[str]) -> Tuple[bool, List[Tuple[str, str, str]]]:
        """
        验证序列配对文件|Validate sequence pair file

        Args:
            pair_file: 配对文件路径|Pair file path
            seq1_ids: 第一个FASTA文件的序列ID集合|Sequence IDs from first FASTA
            seq2_ids: 第二个FASTA文件的序列ID集合|Sequence IDs from second FASTA

        Returns:
            (is_valid, valid_pairs): 验证结果和有效配对列表|Validation result and valid pairs list
        """
        try:
            self.logger.info(f"验证配对文件|Validating pair file: {pair_file}")
            valid_pairs = []
            invalid_count = 0
            total_count = 0

            if pair_file.endswith('.csv'):
                df = pd.read_csv(pair_file, comment='#', header=None)
            else:
                df = pd.read_csv(pair_file, sep='\t', comment='#', header=None)

            total_count = len(df)

            if len(df.columns) < 2:
                raise ValueError("配对文件至少需要2列|Pair file must have at least 2 columns")

            id1_col, id2_col = df.columns[0], df.columns[1]
            name_col = df.columns[2] if len(df.columns) > 2 else None

            for idx, row in df.iterrows():
                seq1_id = str(row[id1_col]).strip()
                seq2_id = str(row[id2_col]).strip()
                pair_name = str(row[name_col]).strip() if name_col and pd.notna(row[name_col]) else f"pair_{idx+1:04d}"

                if seq1_id not in seq1_ids:
                    self.logger.warning(f"序列ID '{seq1_id}' 在第一个FASTA文件中未找到|Sequence ID '{seq1_id}' not found in first FASTA")
                    invalid_count += 1
                    continue

                if seq2_id not in seq2_ids:
                    self.logger.warning(f"序列ID '{seq2_id}' 在第二个FASTA文件中未找到|Sequence ID '{seq2_id}' not found in second FASTA")
                    invalid_count += 1
                    continue

                pair_tuple = (seq1_id, seq2_id, pair_name)
                if pair_tuple in valid_pairs:
                    self.logger.warning(f"发现重复配对|Duplicate pair found: {seq1_id} - {seq2_id}")
                    invalid_count += 1
                    continue

                valid_pairs.append(pair_tuple)

            if invalid_count > 0:
                self.logger.warning(f"发现|Found {invalid_count}/{total_count} 个无效配对|invalid pairs")

            if not valid_pairs:
                self.logger.error(f"文件中没有有效配对|No valid pairs in file: {pair_file}")
                return False, []

            self.logger.success(f"验证完成|Validation complete: {len(valid_pairs)}/{total_count} 个有效配对|valid pairs")
            return True, valid_pairs

        except Exception as e:
            self.logger.error(f"配对文件验证失败|Pair file validation failed: {pair_file} - {e}")
            return False, []
