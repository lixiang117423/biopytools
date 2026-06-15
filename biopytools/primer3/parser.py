"""
Primer3输入输出解析模块|Primer3 Input/Output Parser Module
"""

import re
from pathlib import Path
from typing import List, Dict, Tuple
from Bio import SeqIO
from .utils import format_sequence_id


class FastaParser:
    """FASTA文件解析器|FASTA File Parser"""

    def __init__(self, fasta_file: str):
        """
        初始化FASTA解析器|Initialize FASTA parser

        Args:
            fasta_file: FASTA文件路径|FASTA file path
        """
        self.fasta_file = fasta_file

    def parse(self) -> List[Tuple[str, str]]:
        """
        解析FASTA文件|Parse FASTA file

        Returns:
            序列列表，每个元素为(seq_id, sequence)|List of sequences, each element is (seq_id, sequence)
        """
        sequences = []
        try:
            for record in SeqIO.parse(self.fasta_file, "fasta"):
                seq_id = format_sequence_id(record.id)
                sequence = str(record.seq).upper()
                sequences.append((seq_id, sequence))
        except Exception as e:
            raise RuntimeError(f"FASTA文件解析失败|Failed to parse FASTA file: {e}")

        return sequences


class Primer3InputGenerator:
    """Primer3输入格式生成器|Primer3 Input Format Generator"""

    def __init__(self, settings: dict):
        """
        初始化生成器|Initialize generator

        Args:
            settings: Primer3参数字典|Primer3 parameters dictionary
        """
        self.settings = settings

    def generate_single_sequence_input(self, seq_id: str, sequence: str, seq_specific_settings: dict = None) -> str:
        """
        为单条序列生成Primer3输入格式|Generate Primer3 input format for single sequence

        Args:
            seq_id: 序列ID|Sequence ID
            sequence: 序列|Sequence
            seq_specific_settings: 序列特定的设置|Sequence-specific settings

        Returns:
            Primer3输入格式字符串|Primer3 input format string
        """
        lines = []

        # 序列ID|Sequence ID
        lines.append(f"SEQUENCE_ID={seq_id}")

        # 序列模板|Sequence template
        lines.append(f"SEQUENCE_TEMPLATE={sequence}")

        # 添加全局设置参数|Add global settings
        for key, value in self.settings.items():
            # 如果该参数在序列特定设置中，跳过（让序列特定设置优先）
            # Skip if this parameter is in sequence-specific settings (let seq-specific override)
            if seq_specific_settings and key in seq_specific_settings:
                continue
            lines.append(f"{key}={value}")

        # 添加序列特定的设置（在全局设置之后，会覆盖全局设置）|Add sequence-specific settings
        # Sequence-specific settings are added AFTER global settings, so they override global settings
        if seq_specific_settings:
            for key, value in seq_specific_settings.items():
                lines.append(f"{key}={value}")

        # 序列分隔符|Sequence separator
        lines.append("=")

        return "\n".join(lines) + "\n"

    def generate_batch_input(self, sequences: List[Tuple[str, str]], config=None) -> str:
        """
        为批量序列生成Primer3输入格式|Generate Primer3 input format for batch sequences

        Args:
            sequences: 序列列表|List of sequences
            config: 配置对象|Configuration object (for sequence-specific settings)

        Returns:
            Primer3输入格式字符串|Primer3 input format string
        """
        input_lines = []
        for seq_id, sequence in sequences:
            # 获取序列特定的设置|Get sequence-specific settings
            seq_specific_settings = None
            if config and hasattr(config, 'get_sequence_specific_settings'):
                seq_specific_settings = config.get_sequence_specific_settings(len(sequence))

            seq_input = self.generate_single_sequence_input(seq_id, sequence, seq_specific_settings)
            input_lines.append(seq_input)

        return "".join(input_lines)


class Primer3OutputParser:
    """Primer3输出解析器|Primer3 Output Parser"""

    def __init__(self):
        """初始化解析器|Initialize parser"""
        self.current_record = {}
        self.current_pair_num = -1  # 初始化为-1，确保引物对0能被正确创建|Initialize to -1 to ensure pair 0 is created

    def parse(self, output: str) -> Dict[str, List[Dict]]:
        """
        解析Primer3输出|Parse Primer3 output

        Args:
            output: Primer3输出字符串|Primer3 output string

        Returns:
            字典，键为序列ID，值为引物对列表|Dictionary with seq_id as key and primer pair list as value
        """
        results = {}
        lines = output.strip().split('\n')

        for line in lines:
            line = line.strip()
            if not line:
                continue

            # 处理序列ID|Handle sequence ID
            if line.startswith('SEQUENCE_ID='):
                if self.current_record:
                    seq_id = self.current_record.get('SEQUENCE_ID', '')
                    if seq_id not in results:
                        results[seq_id] = []
                    if self.current_record.get('primers'):
                        results[seq_id].extend(self.current_record['primers'])

                self.current_record = {'SEQUENCE_ID': line.split('=', 1)[1]}
                self.current_pair_num = -1  # 重置为-1，确保新序列的引物对0能被创建|Reset to -1 to ensure pair 0 is created for new sequence
                self.current_record['primers'] = []

            # 处理引物对惩罚值（表示新的一对引物）|Handle primer pair penalty (indicates new pair)
            elif line.startswith('PRIMER_PAIR_') and '_PENALTY=' in line:
                match = re.match(r'PRIMER_PAIR_(\d+)_PENALTY=(.*)', line)
                if match:
                    pair_num = int(match.group(1))
                    penalty = float(match.group(2))

                    if pair_num != self.current_pair_num:
                        self.current_pair_num = pair_num
                        primer_pair = {
                            'pair_num': pair_num,
                            'penalty': penalty
                        }
                        self.current_record['primers'].append(primer_pair)

            # 处理其他引物参数|Handle other primer parameters
            elif '=' in line:
                key, value = line.split('=', 1)

                # 解析引物序列相关参数|Parse primer sequence related parameters
                pair_match = re.match(r'PRIMER_(LEFT|RIGHT)_(\d+)_(.+)', key)
                if pair_match and self.current_record['primers']:
                    primer_type = pair_match.group(1).lower()  # left or right
                    pair_num = int(pair_match.group(2))
                    param_name = pair_match.group(3)

                    # 找到对应的引物对|Find corresponding primer pair
                    primer_pair = None
                    for p in self.current_record['primers']:
                        if p['pair_num'] == pair_num:
                            primer_pair = p
                            break

                    if primer_pair:
                        # 组织参数名|Organize parameter name
                        if param_name == 'SEQUENCE':
                            param_key = f'{primer_type}_primer'
                        elif param_name == 'TM':
                            param_key = f'{primer_type}_tm'
                        elif param_name == 'GC_PERCENT':
                            param_key = f'{primer_type}_gc_percent'
                        elif param_name == 'SELF_ANY_TH':
                            param_key = f'{primer_type}_self_any_th'
                        elif param_name == 'SELF_END_TH':
                            param_key = f'{primer_type}_self_end_th'
                        elif param_name == 'HAIRPIN_TH':
                            param_key = f'{primer_type}_hairpin_th'
                        elif param_name == 'END_STABILITY':
                            param_key = f'{primer_type}_end_stability'
                        else:
                            param_key = f'{primer_type}_{param_name.lower()}'

                        primer_pair[param_key] = value

                # 解析引物对参数|Parse primer pair parameters
                elif key.startswith('PRIMER_PAIR_') and self.current_record['primers']:
                    pair_match = re.match(r'PRIMER_PAIR_(\d+)_(.+)', key)
                    if pair_match:
                        pair_num = int(pair_match.group(1))
                        param_name = pair_match.group(2)

                        # 找到对应的引物对|Find corresponding primer pair
                        primer_pair = None
                        for p in self.current_record['primers']:
                            if p['pair_num'] == pair_num:
                                primer_pair = p
                                break

                        if primer_pair:
                            if param_name == 'PRODUCT_SIZE':
                                primer_pair['product_size'] = int(value)
                            elif param_name == 'PRODUCT_TM':
                                primer_pair['product_tm'] = float(value)
                            elif param_name == 'COMPL_ANY_TH':
                                primer_pair['compl_any_th'] = float(value)
                            elif param_name == 'COMPL_END_TH':
                                primer_pair['compl_end_th'] = float(value)
                            else:
                                primer_pair[param_name.lower()] = value

                # 保存序列级别参数|Save sequence level parameters
                elif key.startswith('SEQUENCE_'):
                    self.current_record[key] = value

        # 添加最后一条记录|Add last record
        if self.current_record:
            seq_id = self.current_record.get('SEQUENCE_ID', '')
            if seq_id not in results:
                results[seq_id] = []
            if self.current_record.get('primers'):
                results[seq_id].extend(self.current_record['primers'])

        return results


class ResultsFormatter:
    """结果格式化器|Results Formatter"""

    # 中英文列名映射|Chinese-English column name mapping
    COLUMN_NAMES = {
        'zh': {
            'ID': '序列ID',
            'Input_Seq': '输入序列',
            'Pair_Num': '引物对编号',
            'F_Primer': '正向引物',
            'R_Primer': '反向引物',
            'F_Primer_Length': '正向引物长度',
            'R_Primer_Length': '反向引物长度',
            'F_Primer_TM': '正向引物Tm值',
            'R_Primer_TM': '反向引物Tm值',
            'Annealing_Temp': '退火温度',
            'F_GC_Percent': '正向GC含量(%)',
            'R_GC_Percent': '反向GC含量(%)',
            'Product_Size': '产物大小(bp)',
            'Product_TM': '产物Tm值',
            'Penalty': '惩罚值',
            'F_Self_Any_TH': '正向自身二聚体Any',
            'R_Self_Any_TH': '反向自身二聚体Any',
            'Cross_Dimer': '引物间二聚体Any',
            'F_Self_End_TH': '正向自身二聚体End',
            'R_Self_End_TH': '反向自身二聚体End',
            'Cross_End_Dimer': '引物间二聚体End',
        },
        'en': {
            'ID': 'Sequence_ID',
            'Input_Seq': 'Input_Sequence',
            'Pair_Num': 'Pair_Number',
            'F_Primer': 'F_Primer',
            'R_Primer': 'R_Primer',
            'F_Primer_Length': 'F_Primer_Length',
            'R_Primer_Length': 'R_Primer_Length',
            'F_Primer_TM': 'F_Primer_TM',
            'R_Primer_TM': 'R_Primer_TM',
            'Annealing_Temp': 'Annealing_Temp',
            'F_GC_Percent': 'F_GC_Percent',
            'R_GC_Percent': 'R_GC_Percent',
            'Product_Size': 'Product_Size_bp',
            'Product_TM': 'Product_TM',
            'Penalty': 'Penalty',
            'F_Self_Any_TH': 'F_Self_Dimer_Any',
            'R_Self_Any_TH': 'R_Self_Dimer_Any',
            'Cross_Dimer': 'Cross_Dimer_Any',
            'F_Self_End_TH': 'F_Self_Dimer_End',
            'R_Self_End_TH': 'R_Self_Dimer_End',
            'Cross_End_Dimer': 'Cross_Dimer_End',
        }
    }

    @staticmethod
    def to_dataframe(results: Dict[str, List[Dict]], sequences: Dict[str, str], header_lang: str = 'zh') -> 'pd.DataFrame':
        """
        将结果转换为pandas DataFrame|Convert results to pandas DataFrame

        Args:
            results: Primer3解析结果|Primer3 parsed results
            sequences: 原始序列字典|Original sequences dictionary
            header_lang: 表头语言(zh/en)|Header language (zh/en)

        Returns:
            pandas DataFrame|pandas DataFrame
        """
        import pandas as pd

        rows = []
        for seq_id, primers in results.items():
            input_seq = sequences.get(seq_id, '')
            for primer in primers:
                row = {
                    'ID': seq_id,
                    'Input_Seq': input_seq,
                    'Pair_Num': primer.get('pair_num', '') + 1,  # 从1开始编号|Start numbering from 1
                    'F_Primer': primer.get('left_primer', ''),
                    'R_Primer': primer.get('right_primer', ''),
                    'F_Primer_Length': len(primer.get('left_primer', '')),
                    'R_Primer_Length': len(primer.get('right_primer', '')),
                    'F_Primer_TM': primer.get('left_tm', ''),
                    'R_Primer_TM': primer.get('right_tm', ''),
                    'Annealing_Temp': min(
                        float(primer.get('left_tm', 0)),
                        float(primer.get('right_tm', 0))
                    ) if primer.get('left_tm') and primer.get('right_tm') else '',
                    'F_GC_Percent': primer.get('left_gc_percent', ''),
                    'R_GC_Percent': primer.get('right_gc_percent', ''),
                    'Product_Size': primer.get('product_size', ''),
                    'Product_TM': primer.get('product_tm', ''),
                    'Penalty': primer.get('penalty', ''),
                    'F_Self_Any_TH': primer.get('left_self_any_th', ''),
                    'R_Self_Any_TH': primer.get('right_self_any_th', ''),
                    'Cross_Dimer': primer.get('compl_any_th', ''),
                    'F_Self_End_TH': primer.get('left_self_end_th', ''),
                    'R_Self_End_TH': primer.get('right_self_end_th', ''),
                    'Cross_End_Dimer': primer.get('compl_end_th', ''),
                }
                rows.append(row)

        df = pd.DataFrame(rows)

        # 重命名列为指定语言|Rename columns to specified language
        column_mapping = ResultsFormatter.COLUMN_NAMES.get(header_lang, ResultsFormatter.COLUMN_NAMES['zh'])
        df = df.rename(columns=column_mapping)

        return df

    @staticmethod
    def save_table(dataframe: 'pd.DataFrame', output_file: str, format: str = 'csv'):
        """
        保存结果表格|Save results table

        Args:
            dataframe: 结果DataFrame|Result DataFrame
            output_file: 输出文件路径|Output file path
            format: 输出格式（csv, tsv, xlsx）|Output format (csv, tsv, xlsx)
        """
        if format == 'csv':
            dataframe.to_csv(output_file, index=False)
        elif format == 'tsv':
            dataframe.to_csv(output_file, sep='\t', index=False)
        elif format == 'xlsx':
            dataframe.to_excel(output_file, index=False)
        else:
            raise ValueError(f"不支持的输出格式|Unsupported output format: {format}")
