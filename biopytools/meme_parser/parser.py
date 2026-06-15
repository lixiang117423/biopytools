"""
MEME结果解析模块|MEME Results Parser Module
"""

import pandas as pd
import xml.etree.ElementTree as ET
from pathlib import Path
from typing import List, Dict
import re


class MemeParser:
    """MEME结果解析器|MEME Results Parser"""

    def __init__(self, logger, config):
        self.logger = logger
        self.config = config

    def parse(self, format: str = None) -> pd.DataFrame:
        """
        解析MEME输出文件|Parse MEME output file

        Args:
            format: 文件格式 ('xml' 或 'txt')，None表示自动检测

        Returns:
            pd.DataFrame: 解析结果|Parsed results
        """
        # 自动检测格式|Auto-detect format
        if format is None:
            format = self._detect_format()

        self.logger.info(f"解析MEME输出|Parsing MEME output ({format.upper()} format)")

        if format == 'xml':
            df = self.parse_xml()
        elif format == 'txt':
            df = self.parse_txt()
        else:
            raise ValueError(f"不支持的格式|Unsupported format: {format}")

        self.logger.info(f"解析完成，共{len(df)}条motif记录|Parsing completed, {len(df)} motif records")
        return df

    def _detect_format(self) -> str:
        """
        自动检测文件格式|Auto-detect file format

        Returns:
            str: 'xml' 或 'txt'
        """
        file_path = self.config.input_file

        if file_path.endswith('.xml'):
            return 'xml'
        elif file_path.endswith('.txt'):
            return 'txt'
        else:
            # 默认尝试XML
            return 'xml'

    def parse_xml(self) -> pd.DataFrame:
        """
        解析MEME XML格式输出|Parse MEME XML format output

        Returns:
            pd.DataFrame: 解析结果|Parsed results
        """
        xml_file = self.config.get_meme_xml_path()

        self.logger.info(f"解析XML文件|Parsing XML file: {xml_file}")

        # 解析XML|Parse XML
        tree = ET.parse(xml_file)
        root = tree.getroot()

        # 提取序列信息|Extract sequence information
        seq_info = self._extract_sequence_info_xml(root)

        # 提取motif信息|Extract motif information
        motif_info = self._extract_motif_info_xml(root)

        # 提取序列上的motif位置信息|Extract motif positions on sequences
        motif_positions = self._extract_motif_positions_xml(root)

        # 合并信息|Merge information
        df = self._merge_xml_data(seq_info, motif_info, motif_positions)

        return df

    def _extract_sequence_info_xml(self, root: ET.Element) -> pd.DataFrame:
        """提取XML中的序列信息|Extract sequence information from XML"""
        sequences = []

        # MEME XML结构: sequence元素有id和name属性
        for seq_elem in root.findall('.//sequence'):
            seq_id = seq_elem.get('id')  # sequence_0, sequence_1, etc.
            seq_name = seq_elem.get('name')  # 实际序列ID，如Ccu05G1702-mRNA-1
            length = int(seq_elem.get('length'))
            sequences.append({
                'sequence_id': seq_id,  # MEME内部ID
                'sequence_name': seq_name,  # 实际序列名称
                'length': length
            })

        return pd.DataFrame(sequences)

    def _extract_motif_info_xml(self, root: ET.Element) -> pd.DataFrame:
        """提取XML中的motif信息|Extract motif information from XML"""
        motifs = []

        # MEME XML结构: root -> motifs -> motif
        for motif_elem in root.findall('.//motif'):
            motif_id = motif_elem.get('id')
            motif_name = motif_elem.get('name', '')

            # 提取motif统计信息
            evalue = float(motif_elem.get('e-value', 0))
            sites = int(motif_elem.get('sites', 0))
            width = int(motif_elem.get('width', 0))

            motifs.append({
                'motif_id': motif_id,
                'motif_name': motif_name,
                'evalue': evalue,
                'sites': sites,
                'width': width
            })

        return pd.DataFrame(motifs)

    def _extract_motif_positions_xml(self, root: ET.Element) -> List[Dict]:
        """提取XML中序列上的motif位置信息|Extract motif positions on sequences from XML"""
        positions = []

        # MEME XML结构: scanned_sites包含每个序列的所有motif位点
        for scanned_sites in root.findall('.//scanned_sites'):
            sequence_id = scanned_sites.get('sequence_id')  # sequence_0, sequence_1, etc.

            # 遍历该序列的所有motif位点
            for site_elem in scanned_sites.findall('scanned_site'):
                motif_id = site_elem.get('motif_id')  # motif_1, motif_2, etc.
                position = int(site_elem.get('position', 0))
                strand = site_elem.get('strand', '+')
                pvalue = float(site_elem.get('pvalue', 1.0))

                positions.append({
                    'sequence_id': sequence_id,
                    'motif_id': motif_id,
                    'position': position,
                    'strand': strand,
                    'pvalue': pvalue
                })

        return positions

    def _merge_xml_data(self, seq_info: pd.DataFrame,
                        motif_info: pd.DataFrame,
                        motif_positions: List[Dict]) -> pd.DataFrame:
        """合并XML解析的数据|Merge XML parsed data"""
        # 转换位置信息为DataFrame
        pos_df = pd.DataFrame(motif_positions)

        # 合并序列信息（通过sequence_id映射到实际序列名称）|Merge sequence info
        merged = pos_df.merge(
            seq_info[['sequence_id', 'sequence_name', 'length']],
            on='sequence_id',
            how='left'
        )

        # 合并motif信息（宽度、e-value等）|Merge motif information
        merged = merged.merge(
            motif_info[['motif_id', 'width', 'evalue', 'sites']],
            on='motif_id',
            how='left'
        )

        # 计算起始和结束位置|Calculate start and end positions
        merged['start'] = merged['position'] + 1
        merged['end'] = merged['position'] + merged['width']

        # 使用实际序列名称作为sequence_id|Use actual sequence name as sequence_id
        merged['sequence_id'] = merged['sequence_name']

        # 选择并排序列|Select and reorder columns
        final_columns = ['sequence_id', 'length', 'motif_id', 'start', 'end',
                        'width', 'strand', 'pvalue', 'evalue', 'sites']

        # 只保留存在的列|Only keep existing columns
        final_columns = [col for col in final_columns if col in merged.columns]

        result = merged[final_columns].copy()

        # 重命名sites为num_sites|Rename sites to num_sites
        if 'sites' in result.columns:
            result = result.rename(columns={'sites': 'num_sites'})

        return result

    def parse_txt(self) -> pd.DataFrame:
        """
        解析MEME TXT格式输出|Parse MEME TXT format output

        Returns:
            pd.DataFrame: 解析结果|Parsed results
        """
        txt_file = self.config.get_meme_txt_path()

        self.logger.info(f"解析TXT文件|Parsing TXT file: {txt_file}")

        with open(txt_file, 'r') as f:
            content = f.read()

        # 提取motif信息|Extract motif information
        motifs = self._parse_txt_motifs(content)

        # 提取序列信息|Extract sequence information
        sequences = self._parse_txt_sequences(content)

        # 合并数据|Merge data
        df = self._merge_txt_data(motifs, sequences)

        return df

    def _parse_txt_motifs(self, content: str) -> List[Dict]:
        """
        从TXT内容中提取motif信息|Extract motif information from TXT content

        MEME TXT格式示例:
        ********************************************************************************
        MOTIF 1 MEME-1
        ********************************************************************************
            ...
            Motif 1 sites sorted by position p-value
        Sequence name            Strand  Start    P-value  Best possible
        ...
        """
        motifs = []

        # 分割motif块|Split motif blocks
        motif_blocks = re.split(
            r'\n\s*\*+\s*\n\s*MOTIF\s+\d+\s+.*?\n\s*\*+\s*\n',
            content,
            flags=re.IGNORECASE
        )

        for i, block in enumerate(motif_blocks[1:], 1):  # 跳过第一个空块
            motif_id = f"Motif-{i}"

            # 查找sites部分
            sites_pattern = r'Motif \d+ sites sorted by position.*?\n(.*?)(?=\n\n|\n\s*$)'
            sites_match = re.search(sites_pattern, block, re.DOTALL | re.IGNORECASE)

            if sites_match:
                sites_section = sites_match.group(1)
                lines = sites_section.strip().split('\n')

                # 跳过表头|Skip header
                for line in lines[1:]:
                    if not line.strip():
                        continue

                    parts = line.split()
                    if len(parts) >= 4:
                        sequence_id = parts[0]
                        strand = parts[1]
                        start = int(parts[2])
                        pvalue = self._parse_scientific_notation(parts[3])

                        motifs.append({
                            'sequence_id': sequence_id,
                            'motif_id': motif_id,
                            'strand': strand,
                            'start': start,
                            'pvalue': pvalue
                        })

        return motifs

    def _parse_txt_sequences(self, content: str) -> Dict[str, int]:
        """从TXT内容中提取序列长度信息|Extract sequence length information from TXT content"""
        sequences = {}

        # 查找序列部分|Find sequences section
        seq_pattern = r'EXECUTION COMMAND.*?\n\n(.*?)(?=\n\n|\n\s*--+|\n\s*\*)'
        seq_match = re.search(seq_pattern, content, re.DOTALL)

        if seq_match:
            seq_section = seq_match.group(1)
            lines = seq_section.split('\n')

            for line in lines:
                if 'sequence' in line.lower():
                    # 提取序列名称和长度
                    parts = line.split()
                    if len(parts) >= 2:
                        seq_name = parts[0]
                        # 尝试提取长度
                        length_match = re.search(r'(\d+)\s+sequences?', line)
                        if length_match:
                            # 这是总序列数，不是单个序列长度
                            pass

        # 在MEME TXT中，序列长度在motif sites部分隐含
        # 我们需要从motif位置推断或者查找特定部分
        # 这里简化处理，返回空字典，由调用方处理

        return sequences

    def _parse_scientific_notation(self, s: str) -> float:
        """解析科学计数法字符串|Parse scientific notation string"""
        try:
            return float(s)
        except ValueError:
            # 处理形如 "1e-10" 或 "1.0e-10"
            s = s.replace('e', 'E')
            return float(s)

    def _merge_txt_data(self, motifs: List[Dict], sequences: Dict) -> pd.DataFrame:
        """合并TXT解析的数据|Merge TXT parsed data"""
        if not motifs:
            return pd.DataFrame()

        df = pd.DataFrame(motifs)

        # 计算end位置（需要motif宽度信息，这在TXT中可能需要额外解析）
        # 简化处理：从start和典型motif宽度推断
        # 实际宽度需要从motif的consensus或sites部分获取

        # 对于TXT格式，width可能需要额外解析
        # 这里暂时设为0或从其他信息推断
        if 'width' not in df.columns:
            df['width'] = 0  # 需要进一步解析

        if 'end' not in df.columns:
            df['end'] = df['start'] + df['width']

        return df
