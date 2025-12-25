# """
# 🧬 序列处理模块 | Sequence Processing Module
# """

# import re
# from pathlib import Path
# from Bio import SeqIO
# from collections import Counter

# class SequenceProcessor:
#     """序列处理器 | Sequence Processor"""
    
#     def __init__(self, config, logger):
#         self.config = config
#         self.logger = logger
    
#     def detect_sequence_type(self, input_file: Path) -> str:
#         """自动检测序列类型 | Auto-detect sequence type"""
#         self.logger.info("🔬 检测序列类型 | Detecting sequence type")
        
#         nucleotides = set('ATGCNatgcn')
#         amino_acids = set('ACDEFGHIKLMNPQRSTVWYacdefghiklmnpqrstvwy')
        
#         total_chars = 0
#         nucleotide_chars = 0
        
#         try:
#             with open(input_file, 'r') as f:
#                 for record in SeqIO.parse(f, 'fasta'):
#                     seq_str = str(record.seq)
#                     total_chars += len(seq_str)
#                     nucleotide_chars += sum(1 for c in seq_str if c in nucleotides)
                    
#                     # 只检查前几条序列即可
#                     if total_chars > 1000:
#                         break
            
#             if total_chars == 0:
#                 raise ValueError("❌ 输入文件为空或格式错误 | Input file is empty or invalid format")
            
#             nucleotide_ratio = nucleotide_chars / total_chars
            
#             if nucleotide_ratio > 0.9:
#                 seq_type = 'nucleotide'
#                 self.logger.info(f"✅ 检测到核酸序列 | Detected nucleotide sequences (ratio: {nucleotide_ratio:.2%})")
#             else:
#                 seq_type = 'protein'
#                 self.logger.info(f"✅ 检测到蛋白序列 | Detected protein sequences (ratio: {nucleotide_ratio:.2%})")
            
#             return seq_type
            
#         except Exception as e:
#             self.logger.error(f"❌ 序列类型检测失败 | Sequence type detection failed: {e}")
#             raise
    
#     def clean_sequences(self, input_file: Path, output_file: Path) -> dict:
#         """清理序列并处理ID | Clean sequences and process IDs"""
#         self.logger.info("🧹 清理序列和处理ID | Cleaning sequences and processing IDs")
        
#         # 读取原始序列
#         records = list(SeqIO.parse(input_file, 'fasta'))
#         self.logger.info(f"📊 读取到 {len(records)} 条序列 | Read {len(records)} sequences")
        
#         # 清理序列中的特殊字符
#         cleaned_records = []
#         for record in records:
#             # 替换序列中的 . 和 * 等特殊字符为空
#             original_seq = str(record.seq)
#             cleaned_seq = re.sub(r'[.*]', '', original_seq)
            
#             if len(cleaned_seq) != len(original_seq):
#                 removed_chars = len(original_seq) - len(cleaned_seq)
#                 self.logger.debug(f"🔧 序列 {record.id}: 移除了 {removed_chars} 个特殊字符 | Removed {removed_chars} special characters")
            
#             record.seq = cleaned_seq
#             cleaned_records.append(record)
        
#         # 处理ID重复问题
#         id_mapping = self._process_duplicate_ids(cleaned_records)
        
#         # 写入清理后的序列
#         SeqIO.write(cleaned_records, output_file, 'fasta')
#         self.logger.info(f"✅ 清理后的序列已保存 | Cleaned sequences saved: {output_file}")
        
#         # 保存ID映射关系
#         if id_mapping:
#             self._save_id_mapping(id_mapping)
        
#         return id_mapping
    
#     def _process_duplicate_ids(self, records: list) -> dict:
#         """处理重复ID | Process duplicate IDs"""
#         self.logger.info("🔄 处理重复ID | Processing duplicate IDs")
        
#         # 首先将冒号替换为下划线
#         for record in records:
#             if ':' in record.id:
#                 original_id = record.id
#                 record.id = record.id.replace(':', '_')
#                 record.description = record.description.replace(original_id, record.id, 1)
#                 self.logger.debug(f"🔄 ID中的冒号已替换 | Replaced colon in ID: {original_id} -> {record.id}")
        
#         # 检查是否还有重复ID
#         id_counts = Counter([record.id for record in records])
#         duplicates = {id_: count for id_, count in id_counts.items() if count > 1}
        
#         id_mapping = {}
        
#         if duplicates:
#             self.logger.warning(f"⚠️ 发现 {len(duplicates)} 个重复ID | Found {len(duplicates)} duplicate IDs")
            
#             # 为重复的ID添加后缀
#             id_counter = {}
#             for i, record in enumerate(records):
#                 if record.id in duplicates:
#                     if record.id not in id_counter:
#                         id_counter[record.id] = 0
#                     else:
#                         id_counter[record.id] += 1
#                         original_id = record.id
#                         new_id = f"{record.id}_{id_counter[record.id]}"
#                         id_mapping[new_id] = original_id
#                         record.id = new_id
#                         record.description = record.description.replace(original_id, new_id, 1)
#                         self.logger.debug(f"🔄 重复ID已重命名 | Renamed duplicate ID: {original_id} -> {new_id}")
#         else:
#             self.logger.info("✅ 没有重复ID | No duplicate IDs found")
        
#         return id_mapping
    
#     def _save_id_mapping(self, id_mapping: dict):
#         """保存ID映射关系 | Save ID mapping"""
#         if not id_mapping:
#             return
        
#         mapping_file = self.config.id_mapping_file
        
#         with open(mapping_file, 'w') as f:
#             f.write("# ID映射关系 | ID Mapping\n")
#             f.write("# New_ID\tOriginal_ID\n")
#             for new_id, original_id in id_mapping.items():
#                 f.write(f"{new_id}\t{original_id}\n")
        
#         self.logger.info(f"📝 ID映射文件已保存 | ID mapping file saved: {mapping_file}")
#         self.logger.info(f"📊 共重命名 {len(id_mapping)} 个ID | Renamed {len(id_mapping)} IDs")

"""
序列处理模块 | Sequence Processor Module
"""

from Bio import SeqIO
from pathlib import Path
import re
from collections import Counter

class SequenceProcessor:
    """序列处理器类 | Sequence Processor Class"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def detect_sequence_type(self, input_file):
        """
        检测序列类型 | Detect sequence type
        
        Args:
            input_file: 输入序列文件路径
            
        Returns:
            str: 'protein' 或 'nucleotide'
        """
        self.logger.info(f"🔍 检测序列类型 | Detecting sequence type: {input_file}")
        
        # 读取前100条序列或所有序列（取较小值）
        sequences = []
        with open(input_file) as f:
            for i, record in enumerate(SeqIO.parse(f, "fasta")):
                sequences.append(str(record.seq).upper())
                if i >= 99:  # 读取前100条
                    break
        
        if not sequences:
            raise ValueError("输入文件中没有找到序列 | No sequences found in input file")
        
        # 统计所有序列的字符
        all_chars = ''.join(sequences)
        char_counter = Counter(all_chars)
        
        # 核酸特征字符
        nucleotide_chars = set('ATCGUN')
        # 蛋白质特征字符（不在核酸中）
        protein_specific_chars = set('EFILPQZ')
        
        # 计算核酸字符比例
        nucleotide_count = sum(char_counter[c] for c in nucleotide_chars if c in char_counter)
        total_count = sum(char_counter.values())
        nucleotide_ratio = nucleotide_count / total_count if total_count > 0 else 0
        
        # 检查是否有蛋白质特异字符
        has_protein_chars = any(c in protein_specific_chars for c in char_counter)
        
        # 判断逻辑
        if has_protein_chars:
            seq_type = 'protein'
        elif nucleotide_ratio > 0.85:  # 如果85%以上是核酸字符
            seq_type = 'nucleotide'
        else:
            seq_type = 'protein'
        
        self.logger.info(f"✅ 检测到序列类型 | Detected sequence type: {seq_type}")
        self.logger.debug(f"   核酸字符比例 | Nucleotide ratio: {nucleotide_ratio:.2%}")
        self.logger.debug(f"   蛋白质特异字符 | Protein-specific chars: {has_protein_chars}")
        
        return seq_type
    
    def clean_sequences(self, input_file, output_file):
        """
        清理序列并记录所有ID映射关系
        
        Args:
            input_file: 输入序列文件
            output_file: 输出清理后的序列文件
            
        Returns:
            dict: ID映射字典 {new_id: original_id}
        """
        self.logger.info(f"🧹 开始清理序列 | Starting sequence cleaning")
        self.logger.info(f"   输入文件 | Input: {input_file}")
        self.logger.info(f"   输出文件 | Output: {output_file}")
        
        id_mapping = {}  # 存储所有ID映射关系
        seen_ids = set()  # 用于检测重复ID
        cleaned_count = 0
        modified_count = 0
        
        # 创建输出目录
        Path(output_file).parent.mkdir(parents=True, exist_ok=True)
        
        with open(input_file) as f_in, open(output_file, 'w') as f_out:
            for record in SeqIO.parse(f_in, "fasta"):
                original_id = record.id
                new_id = original_id
                modifications = []  # 记录所有修改
                
                # 1. 替换特殊字符
                special_chars = {
                    ':': '_',
                    '|': '_',
                    ' ': '_',
                    ';': '_',
                    ',': '_',
                    '(': '_',
                    ')': '_',
                    '[': '_',
                    ']': '_',
                    '{': '_',
                    '}': '_',
                }
                
                for char, replacement in special_chars.items():
                    if char in new_id:
                        new_id = new_id.replace(char, replacement)
                        modifications.append(f"替换'{char}'为'{replacement}'")
                
                # 2. 移除连续的下划线
                if '__' in new_id:
                    new_id = re.sub(r'_+', '_', new_id)
                    modifications.append("移除连续下划线")
                
                # 3. 移除首尾的下划线
                new_id_stripped = new_id.strip('_')
                if new_id_stripped != new_id:
                    new_id = new_id_stripped
                    modifications.append("移除首尾下划线")
                
                # 4. 处理重复ID
                if new_id in seen_ids:
                    counter = 1
                    base_id = new_id
                    while f"{base_id}_{counter}" in seen_ids:
                        counter += 1
                    new_id = f"{base_id}_{counter}"
                    modifications.append(f"重复ID重命名（添加后缀_{counter}）")
                
                # 5. 记录修改
                if modifications:
                    self.logger.debug(f"🔄 ID修改 | ID modified: {original_id} -> {new_id}")
                    self.logger.debug(f"   修改内容 | Modifications: {'; '.join(modifications)}")
                    modified_count += 1
                
                # 6. 记录映射关系（所有ID都记录，无论是否修改）
                id_mapping[new_id] = original_id
                seen_ids.add(new_id)
                
                # 7. 写入清理后的序列
                record.id = new_id
                record.description = ""  # 清空description避免重复
                SeqIO.write(record, f_out, "fasta")
                cleaned_count += 1
        
        # 8. 写入完整的ID映射文件
        self._write_complete_mapping(id_mapping)
        
        self.logger.info(f"✅ 序列清理完成 | Sequence cleaning completed")
        self.logger.info(f"   总序列数 | Total sequences: {cleaned_count}")
        self.logger.info(f"   修改ID数 | Modified IDs: {modified_count}")
        self.logger.info(f"   未修改ID数 | Unchanged IDs: {cleaned_count - modified_count}")
        
        return id_mapping
    
    def _write_complete_mapping(self, id_mapping):
        """
        写入完整的ID映射文件
        
        Args:
            id_mapping: ID映射字典 {new_id: original_id}
        """
        mapping_file = self.config.id_mapping_file
        
        # 创建输出目录
        Path(mapping_file).parent.mkdir(parents=True, exist_ok=True)
        
        # 统计修改和未修改的数量
        modified_ids = [(new_id, orig_id) for new_id, orig_id in id_mapping.items() if new_id != orig_id]
        unchanged_ids = [(new_id, orig_id) for new_id, orig_id in id_mapping.items() if new_id == orig_id]
        
        with open(mapping_file, 'w') as f:
            # 写入文件头
            f.write("# ID映射关系 | ID Mapping\n")
            f.write(f"# 总计: {len(id_mapping)} 个ID | Total: {len(id_mapping)} IDs\n")
            f.write(f"# 修改: {len(modified_ids)} 个 | Modified: {len(modified_ids)}\n")
            f.write(f"# 未修改: {len(unchanged_ids)} 个 | Unchanged: {len(unchanged_ids)}\n")
            f.write("#\n")
            f.write("# New_ID\tOriginal_ID\tStatus\n")
            
            # 先写入修改的ID
            if modified_ids:
                f.write("\n# === 修改的ID | Modified IDs ===\n")
                for new_id, orig_id in sorted(modified_ids):
                    f.write(f"{new_id}\t{orig_id}\tModified\n")
            
            # 再写入未修改的ID
            if unchanged_ids:
                f.write("\n# === 未修改的ID | Unchanged IDs ===\n")
                for new_id, orig_id in sorted(unchanged_ids):
                    f.write(f"{new_id}\t{orig_id}\tUnchanged\n")
        
        self.logger.info(f"✅ ID映射文件已生成 | ID mapping file created: {mapping_file}")
        self.logger.info(f"📊 映射统计 | Mapping statistics:")
        self.logger.info(f"   总计 | Total: {len(id_mapping)} IDs")
        self.logger.info(f"   修改 | Modified: {len(modified_ids)} IDs")
        self.logger.info(f"   未修改 | Unchanged: {len(unchanged_ids)} IDs")
    
    def validate_sequences(self, input_file, seq_type=None):
        """
        验证序列文件
        
        Args:
            input_file: 输入序列文件
            seq_type: 序列类型 ('protein' 或 'nucleotide')
            
        Returns:
            bool: 验证是否通过
        """
        self.logger.info(f"🔍 验证序列文件 | Validating sequence file: {input_file}")
        
        try:
            seq_count = 0
            for record in SeqIO.parse(input_file, "fasta"):
                seq_count += 1
                
                # 检查序列ID
                if not record.id:
                    self.logger.error(f"❌ 发现空序列ID | Found empty sequence ID")
                    return False
                
                # 检查序列长度
                if len(record.seq) == 0:
                    self.logger.error(f"❌ 发现空序列 | Found empty sequence: {record.id}")
                    return False
            
            if seq_count == 0:
                self.logger.error(f"❌ 文件中没有序列 | No sequences found in file")
                return False
            
            self.logger.info(f"✅ 序列验证通过 | Sequence validation passed")
            self.logger.info(f"   序列数量 | Sequence count: {seq_count}")
            return True
            
        except Exception as e:
            self.logger.error(f"❌ 序列验证失败 | Sequence validation failed: {e}")
            return False