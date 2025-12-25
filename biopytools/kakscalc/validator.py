# """
# 🔬 Ka/Ks Calculator输入验证模块
# 功能: 验证FASTA文件和配对文件的格式与内容 | Validate FASTA and pair file formats and contents
# """

# import re
# import pandas as pd
# from typing import Dict, List, Tuple, Set
# from Bio import SeqIO
# from .config import KaKsConfig
# from .logger import Logger

# class SequenceValidator:
#     """🔬 序列验证工具类 | Sequence validation utility"""
    
#     def __init__(self, logger: Logger):
#         """
#         🏗️ 初始化验证器 | Initialize validator
        
#         Args:
#             logger: 日志器实例 | Logger instance
#         """
#         self.logger = logger
#         self.config = KaKsConfig()
        
#     def validate_fasta_file(self, fasta_file: str) -> Tuple[bool, Dict[str, str]]:
#         """
#         🧪 验证FASTA文件格式和内容 | Validate FASTA file format and content
        
#         Args:
#             fasta_file: FASTA文件路径 | FASTA file path
            
#         Returns:
#             (is_valid, sequences_dict): 验证结果和序列字典 | Validation result and sequences dictionary
#         """
#         try:
#             self.logger.info(f"验证FASTA文件 | Validating FASTA file: {fasta_file}", "🧪")
#             sequences = {}
#             invalid_count = 0
#             total_count = 0
            
#             for record in SeqIO.parse(fasta_file, "fasta"):
#                 total_count += 1
#                 seq_str = str(record.seq).upper()
                
#                 # 🔤 检查核苷酸字符 | Check nucleotide characters
#                 if not re.match(r'^[ATCGN]+$', seq_str):
#                     self.logger.warning(f"序列 {record.id} 包含非标准字符 | Sequence {record.id} contains non-standard characters")
#                     invalid_count += 1
#                     continue
                
#                 # 🚫 检查N字符比例 | Check N character ratio
#                 n_ratio = seq_str.count('N') / len(seq_str) if len(seq_str) > 0 else 1
#                 if n_ratio > self.config.MAX_MISSING_RATIO:
#                     self.logger.warning(f"序列 {record.id} N字符比例过高 ({n_ratio:.2%}) | Sequence {record.id} has high N ratio ({n_ratio:.2%})")
#                     invalid_count += 1
#                     continue
                
#                 # 📏 检查序列长度 | Check sequence length
#                 if len(seq_str) < self.config.MIN_SEQUENCE_LENGTH:
#                     self.logger.warning(f"序列 {record.id} 长度过短 ({len(seq_str)} bp) | Sequence {record.id} too short ({len(seq_str)} bp)")
#                     invalid_count += 1
#                     continue
                
#                 if len(seq_str) > self.config.MAX_SEQUENCE_LENGTH:
#                     self.logger.warning(f"序列 {record.id} 长度过长 ({len(seq_str)} bp) | Sequence {record.id} too long ({len(seq_str)} bp)")
#                     invalid_count += 1
#                     continue
                
#                 # 🧬 检查密码子完整性 | Check codon completeness
#                 if len(seq_str) % 3 != 0:
#                     self.logger.warning(f"序列 {record.id} 长度不是3的倍数 | Sequence {record.id} length not multiple of 3")
#                     invalid_count += 1
#                     continue
                
#                 # ✅ 有效序列 | Valid sequence
#                 sequences[record.id] = seq_str.replace('N', '')  # 移除N字符 | Remove N characters
            
#             if invalid_count > 0:
#                 self.logger.warning(f"发现 {invalid_count}/{total_count} 个无效序列 | Found {invalid_count}/{total_count} invalid sequences")
            
#             if not sequences:
#                 self.logger.error(f"文件中没有有效序列 | No valid sequences in file: {fasta_file}")
#                 return False, {}
            
#             self.logger.success(f"验证完成：{len(sequences)}/{total_count} 个有效序列 | Validated: {len(sequences)}/{total_count} valid sequences")
#             return True, sequences
            
#         except Exception as e:
#             self.logger.error(f"FASTA文件验证失败 | FASTA validation failed: {fasta_file} - {e}")
#             return False, {}
    
#     def validate_pair_file(self, pair_file: str, seq1_ids: Set[str], seq2_ids: Set[str]) -> Tuple[bool, List[Tuple[str, str, str]]]:
#         """
#         🔗 验证序列配对文件 | Validate sequence pair file
        
#         Args:
#             pair_file: 配对文件路径 | Pair file path
#             seq1_ids: 第一个FASTA文件的序列ID集合 | Sequence IDs from first FASTA
#             seq2_ids: 第二个FASTA文件的序列ID集合 | Sequence IDs from second FASTA
            
#         Returns:
#             (is_valid, valid_pairs): 验证结果和有效配对列表 | Validation result and valid pairs list
#         """
#         try:
#             self.logger.info(f"验证配对文件 | Validating pair file: {pair_file}", "🔗")
#             valid_pairs = []
#             invalid_count = 0
#             total_count = 0
            
#             # 📖 读取配对文件 | Read pair file
#             if pair_file.endswith('.csv'):
#                 df = pd.read_csv(pair_file, comment='#')
#             else:
#                 df = pd.read_csv(pair_file, sep='\t', comment='#')
            
#             total_count = len(df)
            
#             # 🏷️ 处理不同列格式 | Handle different column formats
#             if len(df.columns) < 2:
#                 raise ValueError("配对文件至少需要2列 | Pair file must have at least 2 columns")
            
#             id1_col, id2_col = df.columns[0], df.columns[1]
#             name_col = df.columns[2] if len(df.columns) > 2 else None
            
#             for idx, row in df.iterrows():
#                 seq1_id = str(row[id1_col]).strip()
#                 seq2_id = str(row[id2_col]).strip()
#                 pair_name = str(row[name_col]).strip() if name_col and pd.notna(row[name_col]) else f"pair_{idx+1:04d}"
                
#                 # ✅ 检查ID是否存在 | Check if IDs exist
#                 if seq1_id not in seq1_ids:
#                     self.logger.warning(f"序列ID '{seq1_id}' 在第一个FASTA文件中未找到 | Sequence ID '{seq1_id}' not found in first FASTA")
#                     invalid_count += 1
#                     continue
                
#                 if seq2_id not in seq2_ids:
#                     self.logger.warning(f"序列ID '{seq2_id}' 在第二个FASTA文件中未找到 | Sequence ID '{seq2_id}' not found in second FASTA")
#                     invalid_count += 1
#                     continue
                
#                 # 🔄 检查重复配对 | Check duplicate pairs
#                 pair_tuple = (seq1_id, seq2_id, pair_name)
#                 if pair_tuple in valid_pairs:
#                     self.logger.warning(f"发现重复配对 | Duplicate pair found: {seq1_id} - {seq2_id}")
#                     invalid_count += 1
#                     continue
                
#                 valid_pairs.append(pair_tuple)
            
#             if invalid_count > 0:
#                 self.logger.warning(f"发现 {invalid_count}/{total_count} 个无效配对 | Found {invalid_count}/{total_count} invalid pairs")
            
#             if not valid_pairs:
#                 self.logger.error(f"文件中没有有效配对 | No valid pairs in file: {pair_file}")
#                 return False, []
            
#             self.logger.success(f"验证完成：{len(valid_pairs)}/{total_count} 个有效配对 | Validated: {len(valid_pairs)}/{total_count} valid pairs")
#             return True, valid_pairs
            
#         except Exception as e:
#             self.logger.error(f"配对文件验证失败 | Pair file validation failed: {pair_file} - {e}")
#             return False, []
    
#     def estimate_sequence_similarity(self, seq1: str, seq2: str) -> float:
#         """
#         🔍 估算序列相似度 | Estimate sequence similarity
        
#         Args:
#             seq1: 第一个序列 | First sequence
#             seq2: 第二个序列 | Second sequence
            
#         Returns:
#             相似度分数 (0-1) | Similarity score (0-1)
#         """
#         try:
#             if len(seq1) == 0 or len(seq2) == 0:
#                 return 0.0
            
#             min_len = min(len(seq1), len(seq2))
#             max_len = max(len(seq1), len(seq2))
            
#             # 长度差异太大 | Length difference too large
#             if min_len / max_len < 0.5:
#                 return 0.0
            
#             # 计算匹配字符数 | Calculate matching characters
#             matches = sum(1 for i in range(min_len) if seq1[i] == seq2[i])
#             similarity = matches / max_len
            
#             return similarity
            
#         except Exception:
#             return 0.0

"""
🔬 Ka/Ks Calculator增强版输入验证模块
功能: 详细的错误提示和用户友好的验证信息
"""

import re
import pandas as pd
from typing import Dict, List, Tuple, Set
from Bio import SeqIO
from .config import KaKsConfig
from .logger import Logger

class SequenceValidator:
    """🔬 增强版序列验证工具类 - 详细用户提醒版"""
    
    def __init__(self, logger: Logger):
        """
        🏗️ 初始化验证器
        
        Args:
            logger: 日志器实例
        """
        self.logger = logger
        self.config = KaKsConfig()
        # 记录各种验证失败的统计
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
        🧪 增强版FASTA文件验证 - 带详细用户提醒
        
        Args:
            fasta_file: FASTA文件路径
            
        Returns:
            (is_valid, sequences_dict): 验证结果和序列字典
        """
        try:
            self.logger.info(f"验证FASTA文件 | Validating FASTA file: {fasta_file}", "🧪")
            self.logger.info("开始逐个验证序列，详细信息如下：", "📋")
            
            sequences = {}
            self._reset_validation_stats()
            
            for record in SeqIO.parse(fasta_file, "fasta"):
                self.validation_stats['total_processed'] += 1
                seq_str = str(record.seq).upper()
                
                # 显示当前处理的序列信息
                self.logger.info(f"\n🔍 验证序列 #{self.validation_stats['total_processed']}: {record.id}")
                self.logger.info(f"   原始长度: {len(str(record.seq))} bp")
                
                # 执行详细验证
                is_valid, processed_seq = self._detailed_sequence_validation(record.id, seq_str)
                
                if is_valid:
                    sequences[record.id] = processed_seq
                    self.validation_stats['total_valid'] += 1
                    self.logger.success(f"   ✅ 验证通过，序列已接受")
                else:
                    self.logger.warning(f"   ❌ 验证失败，序列被跳过")
            
            # 显示验证汇总和用户建议
            self._display_validation_summary(fasta_file)
            
            if len(sequences) == 0:
                self.logger.error("❌ 文件中没有有效序列！请根据上述建议修复后重试。")
                return False, {}
            
            return True, sequences
            
        except Exception as e:
            self.logger.error(f"FASTA文件验证失败 | FASTA validation failed: {fasta_file} - {e}")
            return False, {}
    
    def _detailed_sequence_validation(self, seq_id: str, seq_str: str) -> Tuple[bool, str]:
        """
        🔍 详细的单序列验证，包含用户友好的错误信息
        
        Args:
            seq_id: 序列ID
            seq_str: 序列字符串（已转大写）
            
        Returns:
            (is_valid, processed_sequence): 是否有效和处理后的序列
        """
        original_seq = seq_str
        issues_found = []
        
        # 1. 🔤 字符验证
        self.logger.info("   🔤 检查1: 核苷酸字符验证...")
        char_valid, char_issues = self._validate_characters(seq_str)
        if not char_valid:
            self.validation_stats['character_issues'] += 1
            issues_found.extend(char_issues)
            self._log_character_validation_failure(seq_id, char_issues)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'invalid_characters', 'details': char_issues
            })
            return False, ""
        else:
            self.logger.info("      ✓ 字符验证通过：只包含标准核苷酸字符 (A,T,C,G,N)")
        
        # 2. 🚫 N字符比例验证
        self.logger.info("   🚫 检查2: N字符比例验证...")
        n_valid, n_ratio = self._validate_n_ratio(seq_str)
        if not n_valid:
            self.validation_stats['n_ratio_issues'] += 1
            self._log_n_ratio_validation_failure(seq_id, n_ratio)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'high_n_ratio', 'details': f'{n_ratio:.2%}'
            })
            return False, ""
        else:
            self.logger.info(f"      ✓ N字符比例合格：{n_ratio:.2%} (阈值: ≤{self.config.MAX_MISSING_RATIO:.1%})")
        
        # 3. 📏 长度验证
        self.logger.info("   📏 检查3: 序列长度验证...")
        length_valid, length_issue = self._validate_length(seq_str)
        if not length_valid:
            self.validation_stats['length_issues'] += 1
            self._log_length_validation_failure(seq_id, len(seq_str), length_issue)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'invalid_length', 'details': f'{len(seq_str)}bp'
            })
            return False, ""
        else:
            self.logger.info(f"      ✓ 长度合格：{len(seq_str)} bp (范围: {self.config.MIN_SEQUENCE_LENGTH}-{self.config.MAX_SEQUENCE_LENGTH} bp)")
        
        # 4. 🧬 密码子完整性验证
        self.logger.info("   🧬 检查4: 密码子完整性验证...")
        codon_valid, remainder = self._validate_codon_completeness(seq_str)
        if not codon_valid:
            self.validation_stats['codon_issues'] += 1
            self._log_codon_validation_failure(seq_id, len(seq_str), remainder)
            self.validation_stats['failed_sequences'].append({
                'id': seq_id, 'reason': 'incomplete_codon', 'details': f'余数{remainder}'
            })
            return False, ""
        else:
            codon_count = len(seq_str) // 3
            self.logger.info(f"      ✓ 密码子完整：{len(seq_str)} bp = {codon_count} 个完整密码子")
        
        # 5. ✨ 最终处理
        self.logger.info("   ✨ 最终处理: 移除N字符...")
        final_seq = seq_str.replace('N', '')
        removed_n_count = len(seq_str) - len(final_seq)
        if removed_n_count > 0:
            self.logger.info(f"      ✓ 移除了 {removed_n_count} 个N字符，最终长度: {len(final_seq)} bp")
        else:
            self.logger.info(f"      ✓ 无需移除N字符，最终长度: {len(final_seq)} bp")
        
        return True, final_seq
    
    def _validate_characters(self, seq_str: str) -> Tuple[bool, List[str]]:
        """验证字符并返回详细问题"""
        if re.match(r'^[ATCGN]+$', seq_str):
            return True, []
        
        # 找出所有无效字符
        invalid_chars = set(char for char in seq_str if char not in 'ATCGN')
        issues = []
        
        # 分类无效字符
        lowercase_chars = set(char for char in invalid_chars if char.lower() in 'atcg')
        iupac_chars = set(char for char in invalid_chars if char in 'RYSWKMDBHV')
        other_chars = invalid_chars - lowercase_chars - iupac_chars
        
        if lowercase_chars:
            issues.append(f"包含小写核苷酸: {sorted(lowercase_chars)}")
        if iupac_chars:
            issues.append(f"包含IUPAC模糊代码: {sorted(iupac_chars)}")
        if other_chars:
            issues.append(f"包含其他无效字符: {sorted(other_chars)}")
        
        return False, issues
    
    def _validate_n_ratio(self, seq_str: str) -> Tuple[bool, float]:
        """验证N字符比例"""
        n_ratio = seq_str.count('N') / len(seq_str) if len(seq_str) > 0 else 1
        return n_ratio <= self.config.MAX_MISSING_RATIO, n_ratio
    
    def _validate_length(self, seq_str: str) -> Tuple[bool, str]:
        """验证序列长度"""
        length = len(seq_str)
        if length < self.config.MIN_SEQUENCE_LENGTH:
            return False, "too_short"
        elif length > self.config.MAX_SEQUENCE_LENGTH:
            return False, "too_long"
        return True, "valid"
    
    def _validate_codon_completeness(self, seq_str: str) -> Tuple[bool, int]:
        """验证密码子完整性"""
        remainder = len(seq_str) % 3
        return remainder == 0, remainder
    
    def _log_character_validation_failure(self, seq_id: str, issues: List[str]):
        """记录字符验证失败的详细信息"""
        self.logger.error(f"      ❌ 字符验证失败！")
        self.logger.error(f"         序列ID: {seq_id}")
        for issue in issues:
            self.logger.error(f"         问题: {issue}")
        
        self.logger.info("      💡 修复建议：")
        if any("小写" in issue for issue in issues):
            self.logger.info("         • 转换小写字母为大写：seq = seq.upper()")
        if any("IUPAC" in issue for issue in issues):
            self.logger.info("         • 替换IUPAC模糊代码为标准字符或N")
            self.logger.info("           例如: R→A, Y→C, S→G, W→A, K→G, M→A")
        if any("其他" in issue for issue in issues):
            self.logger.info("         • 移除或替换非核苷酸字符")
            self.logger.info("         • 确保序列只包含 A, T, C, G, N 字符")
    
    def _log_n_ratio_validation_failure(self, seq_id: str, n_ratio: float):
        """记录N字符比例验证失败的详细信息"""
        self.logger.error(f"      ❌ N字符比例过高！")
        self.logger.error(f"         序列ID: {seq_id}")
        self.logger.error(f"         N字符比例: {n_ratio:.2%} (阈值: ≤{self.config.MAX_MISSING_RATIO:.1%})")
        
        self.logger.info("      💡 修复建议：")
        self.logger.info("         • 检查序列质量，可能包含低质量区域")
        self.logger.info("         • 使用序列清理工具移除低质量片段")
        self.logger.info("         • 如果是组装序列，检查gap填充质量")
        self.logger.info("         • 考虑从序列两端修剪N字符")
    
    def _log_length_validation_failure(self, seq_id: str, length: int, issue_type: str):
        """记录长度验证失败的详细信息"""
        self.logger.error(f"      ❌ 序列长度不合格！")
        self.logger.error(f"         序列ID: {seq_id}")
        self.logger.error(f"         当前长度: {length} bp")
        
        if issue_type == "too_short":
            self.logger.error(f"         要求最小长度: {self.config.MIN_SEQUENCE_LENGTH} bp")
            self.logger.info("      💡 修复建议：")
            self.logger.info("         • 检查序列是否完整（可能缺失部分片段）")
            self.logger.info("         • 确认这是完整的CDS序列")
            self.logger.info("         • 如果是部分序列，考虑获取完整序列")
        else:  # too_long
            self.logger.error(f"         允许最大长度: {self.config.MAX_SEQUENCE_LENGTH} bp")
            self.logger.info("      💡 修复建议：")
            self.logger.info("         • 确认这是CDS序列而非完整基因或基因组片段")
            self.logger.info("         • 移除内含子、UTR等非编码区域")
            self.logger.info("         • 如果确实是超长CDS，可调整MAX_SEQUENCE_LENGTH参数")
    
    def _log_codon_validation_failure(self, seq_id: str, length: int, remainder: int):
        """记录密码子验证失败的详细信息"""
        self.logger.error(f"      ❌ 密码子不完整！")
        self.logger.error(f"         序列ID: {seq_id}")
        self.logger.error(f"         序列长度: {length} bp")
        self.logger.error(f"         余数: {remainder} bp (长度必须是3的倍数)")
        
        self.logger.info("      💡 修复建议：")
        if remainder == 1:
            self.logger.info("         • 序列末尾多了1个核苷酸，建议移除最后1个字符")
            self.logger.info("         • 检查序列是否包含终止密码子后的额外字符")
        elif remainder == 2:
            self.logger.info("         • 序列末尾多了2个核苷酸，建议移除最后2个字符")
            self.logger.info("         • 检查序列起始位置是否正确")
        
        self.logger.info("         • 自动修复: seq = seq[:len(seq)//3*3]")
        self.logger.info("         • 手动检查序列边界，确保正确的起始和终止位置")
    
    def _display_validation_summary(self, fasta_file: str):
        """显示验证汇总信息和建议"""
        stats = self.validation_stats
        
        self.logger.separator("📊 验证结果汇总")
        self.logger.info(f"文件: {fasta_file}")
        self.logger.info(f"总序列数: {stats['total_processed']}")
        self.logger.info(f"有效序列: {stats['total_valid']}")
        self.logger.info(f"无效序列: {stats['total_processed'] - stats['total_valid']}")
        
        if stats['total_valid'] > 0:
            success_rate = (stats['total_valid'] / stats['total_processed']) * 100
            self.logger.success(f"成功率: {success_rate:.1f}%")
        
        # 显示失败统计
        if stats['total_processed'] > stats['total_valid']:
            self.logger.warning("\n❌ 验证失败统计:")
            
            if stats['character_issues'] > 0:
                self.logger.warning(f"   字符问题: {stats['character_issues']} 个序列")
            if stats['n_ratio_issues'] > 0:
                self.logger.warning(f"   N字符过多: {stats['n_ratio_issues']} 个序列")
            if stats['length_issues'] > 0:
                self.logger.warning(f"   长度问题: {stats['length_issues']} 个序列")
            if stats['codon_issues'] > 0:
                self.logger.warning(f"   密码子问题: {stats['codon_issues']} 个序列")
            
            # 显示具体失败的序列
            self.logger.info("\n📋 失败序列详情:")
            for failed_seq in stats['failed_sequences'][:5]:  # 只显示前5个
                self.logger.info(f"   • {failed_seq['id']}: {failed_seq['reason']} ({failed_seq['details']})")
            
            if len(stats['failed_sequences']) > 5:
                remaining = len(stats['failed_sequences']) - 5
                self.logger.info(f"   ... 还有 {remaining} 个失败序列")
            
            # 提供总体修复建议
            self._provide_overall_suggestions(stats)
    
    def _provide_overall_suggestions(self, stats: Dict):
        """提供总体修复建议"""
        self.logger.info("\n💡 总体修复建议:")
        
        if stats['character_issues'] > 0:
            self.logger.info("🔤 字符问题修复:")
            self.logger.info("   1. 批量转换为大写: awk '{if(/^>/){print}else{print toupper($0)}}' input.fasta > output.fasta")
            self.logger.info("   2. 替换IUPAC字符: sed 's/[RYSWKMDBHV]/N/g' input.fasta > output.fasta")
        
        if stats['n_ratio_issues'] > 0:
            self.logger.info("🚫 N字符问题修复:")
            self.logger.info("   1. 移除首尾N: 可使用seqkit工具")
            self.logger.info("   2. 过滤高N比例序列: 使用生物信息学质控工具")
        
        if stats['length_issues'] > 0:
            self.logger.info("📏 长度问题修复:")
            self.logger.info("   1. 确认序列来源和完整性")
            self.logger.info("   2. 调整配置文件中的长度阈值")
        
        if stats['codon_issues'] > 0:
            self.logger.info("🧬 密码子问题修复:")
            self.logger.info("   1. 自动修剪: 保留3的倍数长度")
            self.logger.info("   2. 检查序列注释和起止位置")
        
        self.logger.info("\n🛠️ 推荐工具:")
        self.logger.info("   • SeqKit: 序列处理瑞士军刀")
        self.logger.info("   • FASTX-Toolkit: 序列质控工具")
        self.logger.info("   • 我提供的诊断和修复工具")
        
        self.logger.info("\n📞 获取帮助:")
        self.logger.info("   • 运行详细诊断: python fasta_diagnostic_tool.py your_file.fasta")
        self.logger.info("   • 启用详细日志: 添加 --verbose 参数")
    
    def _reset_validation_stats(self):
        """重置验证统计"""
        self.validation_stats = {
            'character_issues': 0,
            'n_ratio_issues': 0,
            'length_issues': 0, 
            'codon_issues': 0,
            'total_processed': 0,
            'total_valid': 0,
            'failed_sequences': []
        }
    
    def validate_pair_file(self, pair_file: str, seq1_ids: Set[str], seq2_ids: Set[str]) -> Tuple[bool, List[Tuple[str, str, str]]]:
        """
        🔗 验证序列配对文件 - 保持原有逻辑
        """
        try:
            self.logger.info(f"验证配对文件 | Validating pair file: {pair_file}", "🔗")
            valid_pairs = []
            invalid_count = 0
            total_count = 0
            
            # 📖 读取配对文件
            if pair_file.endswith('.csv'):
                df = pd.read_csv(pair_file, comment='#', header=None)
            else:
                df = pd.read_csv(pair_file, sep='\t', comment='#', header=None)
            
            total_count = len(df)
            
            # 🏷️ 处理不同列格式
            if len(df.columns) < 2:
                raise ValueError("配对文件至少需要2列 | Pair file must have at least 2 columns")
            
            id1_col, id2_col = df.columns[0], df.columns[1]
            name_col = df.columns[2] if len(df.columns) > 2 else None
            
            for idx, row in df.iterrows():
                seq1_id = str(row[id1_col]).strip()
                seq2_id = str(row[id2_col]).strip()
                pair_name = str(row[name_col]).strip() if name_col and pd.notna(row[name_col]) else f"pair_{idx+1:04d}"
                
                # ✅ 检查ID是否存在
                if seq1_id not in seq1_ids:
                    self.logger.warning(f"序列ID '{seq1_id}' 在第一个FASTA文件中未找到")
                    invalid_count += 1
                    continue
                
                if seq2_id not in seq2_ids:
                    self.logger.warning(f"序列ID '{seq2_id}' 在第二个FASTA文件中未找到")
                    invalid_count += 1
                    continue
                
                # 🔄 检查重复配对
                pair_tuple = (seq1_id, seq2_id, pair_name)
                if pair_tuple in valid_pairs:
                    self.logger.warning(f"发现重复配对: {seq1_id} - {seq2_id}")
                    invalid_count += 1
                    continue
                
                valid_pairs.append(pair_tuple)
            
            if invalid_count > 0:
                self.logger.warning(f"发现 {invalid_count}/{total_count} 个无效配对")
            
            if not valid_pairs:
                self.logger.error(f"文件中没有有效配对: {pair_file}")
                return False, []
            
            self.logger.success(f"验证完成：{len(valid_pairs)}/{total_count} 个有效配对")
            return True, valid_pairs
            
        except Exception as e:
            self.logger.error(f"配对文件验证失败: {pair_file} - {e}")
            return False, []