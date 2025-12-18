"""
序列子集提取核心模块 | Sequence Subsequence Extraction Core Module
"""

import os
from pathlib import Path
from Bio import SeqIO
from .utils import validate_file_exists, validate_id_list


class SequenceExtractor:
    """序列提取器 | Sequence Extractor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def extract_sequences_by_id_list(self):
        """根据ID列表提取序列 | Extract sequences by ID list"""
        input_fasta = self.config.input_fasta
        id_list_file = self.config.id_list_file
        output_fasta = self.config.output_fasta
        keep_order = self.config.keep_order

        # 验证输入文件 | Validate input files
        validate_file_exists(input_fasta, "输入FASTA文件 | Input FASTA file")
        validate_file_exists(id_list_file, "ID列表文件 | ID list file")

        # 验证ID列表 | Validate ID list
        id_list = validate_id_list(id_list_file, self.logger)
        if not id_list:
            self.logger.error("❌ 无法获取有效的ID列表 | Cannot get valid ID list")
            return False

        # 确保输出目录存在 | Ensure output directory exists
        output_dir = Path(output_fasta).parent
        os.makedirs(output_dir, exist_ok=True)

        self.logger.info(f"📄 输入FASTA文件 | Input FASTA file: {input_fasta}")
        self.logger.info(f"📋 ID列表文件 | ID list file: {id_list_file}")
        self.logger.info(f"📤 输出FASTA文件 | Output FASTA file: {output_fasta}")
        self.logger.info(f"📋 待提取ID数量 | Number of IDs to extract: {len(id_list)}")
        self.logger.info(f"🔄 保持原始顺序 | Keep original order: {keep_order}")

        try:
            # 将FASTA读入字典 | Load FASTA into dictionary
            self.logger.info("📖 正在读取FASTA文件 | Reading FASTA file...")
            seq_dict = SeqIO.to_dict(SeqIO.parse(input_fasta, 'fasta'))
            self.logger.info(f"✅ 成功读取 {len(seq_dict)} 条序列 | Successfully loaded {len(seq_dict)} sequences")

            # 检查哪些ID在FASTA文件中存在 | Check which IDs exist in FASTA file
            found_ids = []
            missing_ids = []

            for seq_id in id_list:
                if seq_id in seq_dict:
                    found_ids.append(seq_id)
                else:
                    missing_ids.append(seq_id)

            self.logger.info(f"✅ 找到 {len(found_ids)} 个匹配的ID | Found {len(found_ids)} matching IDs")

            if missing_ids:
                self.logger.warning(f"⚠️ {len(missing_ids)} 个ID在FASTA文件中未找到 | {len(missing_ids)} IDs not found in FASTA file:")
                for missing_id in missing_ids[:10]:  # 只显示前10个 | Only show first 10
                    self.logger.warning(f"    - {missing_id}")
                if len(missing_ids) > 10:
                    self.logger.warning(f"    ... 还有 {len(missing_ids) - 10} 个未显示 | ... and {len(missing_ids) - 10} more not shown")

            if not found_ids:
                self.logger.error("❌ 没有找到任何匹配的ID，请检查ID列表和FASTA文件 | No matching IDs found, please check ID list and FASTA file")
                return False

            # 按照ID列表顺序输出 | Output in ID list order
            self.logger.info("💾 正在写入输出文件 | Writing output file...")
            extracted_count = 0

            with open(output_fasta, 'w') as out:
                for seq_id in found_ids:
                    try:
                        SeqIO.write(seq_dict[seq_id], out, 'fasta')
                        extracted_count += 1
                    except Exception as e:
                        self.logger.error(f"❌ 写入序列 {seq_id} 失败 | Failed to write sequence {seq_id}: {e}")
                        continue

            self.logger.info(f"✅ 成功提取 {extracted_count} 条序列到 | Successfully extracted {extracted_count} sequences to: {output_fasta}")

            # 统计信息 | Statistics
            total_sequences = len(seq_dict)
            extraction_rate = (extracted_count / len(found_ids)) * 100 if found_ids else 0
            coverage_rate = (len(found_ids) / total_sequences) * 100 if total_sequences > 0 else 0

            self.logger.info("📊 统计信息 | Statistics:")
            self.logger.info(f"    📄 输入序列总数 | Total input sequences: {total_sequences}")
            self.logger.info(f"    📋 请求提取ID数 | Requested IDs: {len(id_list)}")
            self.logger.info(f"    ✅ 成功匹配ID数 | Successfully matched IDs: {len(found_ids)}")
            self.logger.info(f"    📤 成功提取序列数 | Successfully extracted sequences: {extracted_count}")
            self.logger.info(f"    📈 提取成功率 | Extraction success rate: {extraction_rate:.2f}%")

            return True

        except Exception as e:
            self.logger.error(f"❌ 序列提取过程中发生错误 | Error occurred during sequence extraction: {e}")
            return False

    def extract_sequences_by_pattern(self):
        """根据模式匹配提取序列 | Extract sequences by pattern matching"""
        input_fasta = self.config.input_fasta
        output_fasta = self.config.output_fasta
        pattern = self.config.pattern
        pattern_type = getattr(self.config, 'pattern_type', 'contains')  # 'contains', 'startswith', 'endswith', 'regex'
        case_sensitive = getattr(self.config, 'case_sensitive', True)

        # 验证输入文件 | Validate input files
        validate_file_exists(input_fasta, "输入FASTA文件 | Input FASTA file")

        # 确保输出目录存在 | Ensure output directory exists
        output_dir = Path(output_fasta).parent
        os.makedirs(output_dir, exist_ok=True)

        self.logger.info(f"📄 输入FASTA文件 | Input FASTA file: {input_fasta}")
        self.logger.info(f"📤 输出FASTA文件 | Output FASTA file: {output_fasta}")
        self.logger.info(f"🔍 匹配模式 | Pattern: {pattern}")
        self.logger.info(f"📝 模式类型 | Pattern type: {pattern_type}")
        self.logger.info(f"🔤 区分大小写 | Case sensitive: {case_sensitive}")

        try:
            self.logger.info("📖 正在读取FASTA文件 | Reading FASTA file...")
            matching_sequences = []

            # 根据模式类型进行匹配 | Match based on pattern type
            for record in SeqIO.parse(input_fasta, 'fasta'):
                seq_id = record.id

                # 应用大小写敏感性 | Apply case sensitivity
                search_id = seq_id if case_sensitive else seq_id.lower()
                search_pattern = pattern if case_sensitive else pattern.lower()

                match = False
                if pattern_type == 'contains':
                    match = search_pattern in search_id
                elif pattern_type == 'startswith':
                    match = search_id.startswith(search_pattern)
                elif pattern_type == 'endswith':
                    match = search_id.endswith(search_pattern)
                elif pattern_type == 'regex':
                    import re
                    flags = 0 if case_sensitive else re.IGNORECASE
                    match = re.search(search_pattern, search_id, flags) is not None

                if match:
                    matching_sequences.append(record)

            self.logger.info(f"✅ 找到 {len(matching_sequences)} 条匹配序列 | Found {len(matching_sequences)} matching sequences")

            if not matching_sequences:
                self.logger.warning("⚠️ 没有找到匹配的序列 | No matching sequences found")
                return False

            # 写入匹配的序列 | Write matching sequences
            self.logger.info("💾 正在写入输出文件 | Writing output file...")
            with open(output_fasta, 'w') as out:
                SeqIO.write(matching_sequences, out, 'fasta')

            self.logger.info(f"✅ 成功提取 {len(matching_sequences)} 条序列到 | Successfully extracted {len(matching_sequences)} sequences to: {output_fasta}")
            return True

        except Exception as e:
            self.logger.error(f"❌ 模式匹配提取过程中发生错误 | Error occurred during pattern matching extraction: {e}")
            return False

    def extract_sequences_by_length(self):
        """根据长度范围提取序列 | Extract sequences by length range"""
        input_fasta = self.config.input_fasta
        output_fasta = self.config.output_fasta
        min_length = getattr(self.config, 'min_length', 0)
        max_length = getattr(self.config, 'max_length', float('inf'))

        # 验证输入文件 | Validate input files
        validate_file_exists(input_fasta, "输入FASTA文件 | Input FASTA file")

        # 确保输出目录存在 | Ensure output directory exists
        output_dir = Path(output_fasta).parent
        os.makedirs(output_dir, exist_ok=True)

        self.logger.info(f"📄 输入FASTA文件 | Input FASTA file: {input_fasta}")
        self.logger.info(f"📤 输出FASTA文件 | Output FASTA file: {output_fasta}")
        self.logger.info(f"📏 最小长度 | Minimum length: {min_length}")
        self.logger.info(f"📏 最大长度 | Maximum length: {max_length}")

        try:
            self.logger.info("📖 正在读取FASTA文件 | Reading FASTA file...")
            matching_sequences = []

            for record in SeqIO.parse(input_fasta, 'fasta'):
                seq_length = len(record.seq)
                if min_length <= seq_length <= max_length:
                    matching_sequences.append(record)

            self.logger.info(f"✅ 找到 {len(matching_sequences)} 条符合长度范围的序列 | Found {len(matching_sequences)} sequences within length range")

            if not matching_sequences:
                self.logger.warning("⚠️ 没有找到符合长度范围的序列 | No sequences found within the specified length range")
                return False

            # 写入匹配的序列 | Write matching sequences
            self.logger.info("💾 正在写入输出文件 | Writing output file...")
            with open(output_fasta, 'w') as out:
                SeqIO.write(matching_sequences, out, 'fasta')

            self.logger.info(f"✅ 成功提取 {len(matching_sequences)} 条序列到 | Successfully extracted {len(matching_sequences)} sequences to: {output_fasta}")
            return True

        except Exception as e:
            self.logger.error(f"❌ 长度筛选提取过程中发生错误 | Error occurred during length filtering extraction: {e}")
            return False