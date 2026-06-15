"""
RxLR基序扫描器|RxLR Motif Scanner
"""

import re
from dataclasses import dataclass
from typing import List, Dict, Optional, Tuple
from .utils import SequenceExtractor


@dataclass
class MotifMatch:
    """基序匹配结果|Motif Match Result"""
    motif_type: str  # 基序类型|Motif type (RxLR/QxLR/GxLR/EER)
    position: int    # 在窗口中的位置|Position in window (0-based)
    context: str     # 上下文序列|Context sequence
    absolute_position: int  # 在原序列中的绝对位置|Absolute position in original sequence


@dataclass
class ScanResult:
    """扫描结果|Scan Result"""
    sequence_id: str
    sequence_length: int
    is_valid_length: bool
    window_sequence: Optional[str]
    motif_matches: List[MotifMatch]
    has_rxlr: bool
    has_eer: bool
    is_candidate: bool


class RxLRMotifScanner:
    """RxLR基序扫描器|RxLR Motif Scanner"""

    # 定义基序模式|Define motif patterns
    RXLR_PATTERNS = {
        'RxLR': re.compile(r'R.LR'),
        'QxLR': re.compile(r'Q.LR'),
        'GxLR': re.compile(r'G.LR'),
    }
    EER_PATTERN = re.compile(r'EER')

    def __init__(self, window_start: int = 20, window_end: int = 120, logger=None):
        """
        初始化扫描器|Initialize scanner

        Args:
            window_start: 窗口起始位置|Window start position (0-based)
            window_end: 窗口结束位置|Window end position
            logger: 日志记录器|Logger instance
        """
        self.window_start = window_start
        self.window_end = window_end
        self.logger = logger

    def scan_sequence(self, seq_id: str, sequence: str, min_length: int = 120) -> ScanResult:
        """
        扫描单个序列|Scan a single sequence

        Args:
            seq_id: 序列ID|Sequence ID
            sequence: 序列|Sequence
            min_length: 最小长度|Minimum length

        Returns:
            扫描结果|Scan result
        """
        seq_len = len(sequence)
        is_valid_length = seq_len >= min_length

        # 提取窗口序列|Extract window sequence
        window_seq = SequenceExtractor.extract_window(
            sequence,
            self.window_start,
            self.window_end
        )

        # 如果窗口提取失败（序列太短），仍然记录结果
        if window_seq is None:
            window_seq = ""

        # 搜索基序|Search motifs
        motif_matches = []

        if window_seq:
            # 搜索RxLR相关基序|Search RxLR-related motifs
            for motif_name, pattern in self.RXLR_PATTERNS.items():
                matches = self._search_motif_in_window(
                    window_seq,
                    pattern,
                    motif_name
                )
                motif_matches.extend(matches)

            # 搜索EER基序|Search EER motif
            eer_matches = self._search_motif_in_window(
                window_seq,
                self.EER_PATTERN,
                'EER'
            )
            motif_matches.extend(eer_matches)

        # 判定是否为候选|Determine if candidate
        has_rxlr = any(m.motif_type in ['RxLR', 'QxLR', 'GxLR'] for m in motif_matches)
        has_eer = any(m.motif_type == 'EER' for m in motif_matches)
        is_candidate = has_rxlr or has_eer

        return ScanResult(
            sequence_id=seq_id,
            sequence_length=seq_len,
            is_valid_length=is_valid_length,
            window_sequence=window_seq if window_seq else None,
            motif_matches=motif_matches,
            has_rxlr=has_rxlr,
            has_eer=has_eer,
            is_candidate=is_candidate
        )

    def _search_motif_in_window(self, window_seq: str, pattern: re.Pattern, motif_name: str) -> List[MotifMatch]:
        """
        在窗口中搜索基序|Search motif in window

        Args:
            window_seq: 窗口序列|Window sequence
            pattern: 正则表达式模式|Regex pattern
            motif_name: 基序名称|Motif name

        Returns:
            匹配结果列表|List of matches
        """
        matches = []

        for match in pattern.finditer(window_seq):
            start_pos = match.start()
            matched_text = match.group()

            # 提取上下文（前后各5个氨基酸）|Extract context (5 aa before and after)
            context_start = max(0, start_pos - 5)
            context_end = min(len(window_seq), start_pos + len(matched_text) + 5)
            context = window_seq[context_start:context_end]

            # 计算绝对位置|Calculate absolute position
            absolute_position = self.window_start + start_pos

            matches.append(MotifMatch(
                motif_type=motif_name,
                position=start_pos,
                context=context,
                absolute_position=absolute_position
            ))

        return matches

    def get_motif_summary(self, result: ScanResult) -> Dict:
        """
        获取基序匹配摘要|Get motif match summary

        Args:
            result: 扫描结果|Scan result

        Returns:
            摘要字典|Summary dictionary
        """
        motif_counts = {}
        for match in result.motif_matches:
            motif_counts[match.motif_type] = motif_counts.get(match.motif_type, 0) + 1

        return {
            'sequence_id': result.sequence_id,
            'length': result.sequence_length,
            'valid_length': result.is_valid_length,
            'has_rxlr': result.has_rxlr,
            'has_eer': result.has_eer,
            'is_candidate': result.is_candidate,
            'motif_counts': motif_counts,
            'total_matches': len(result.motif_matches)
        }

    def format_result_for_output(self, result: ScanResult) -> Dict:
        """
        格式化结果用于输出|Format result for output

        Args:
            result: 扫描结果|Scan result

        Returns:
            格式化的结果字典|Formatted result dictionary
        """
        # 构建基序匹配信息|Build motif match information
        motif_info = []
        for match in result.motif_matches:
            motif_info.append(
                f"{match.motif_type}@{match.absolute_position + 1}"  # 转为1-based|Convert to 1-based
            )

        return {
            'Sequence_ID': result.sequence_id,
            'Sequence_Length': result.sequence_length,
            'Valid_Length_(≥120)': 'Yes' if result.is_valid_length else 'No',
            'Has_RxLR_motif': 'Yes' if result.has_rxlr else 'No',
            'Has_EER_motif': 'Yes' if result.has_eer else 'No',
            'RxLR_Candidate': 'Yes' if result.is_candidate else 'No',
            'Motif_Types': ';'.join(set(m.motif_type for m in result.motif_matches)) if result.motif_matches else 'None',
            'Motif_Positions': ';'.join(motif_info) if motif_info else 'None',
            'Total_Motif_Count': len(result.motif_matches)
        }
