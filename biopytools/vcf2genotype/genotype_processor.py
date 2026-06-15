"""
基因型数据处理模块|Genotype Data Processing Module
"""

from typing import Dict, Any, Optional


def calculate_variant_length(ref: str, alt: str) -> int:
    """
    计算变异长度|Calculate variant length

    Args:
        ref: 参考序列|Reference sequence
        alt: 替换序列（多个ALT用逗号分隔）|Alternative sequence(s), comma-separated

    Returns:
        int: 变异长度|Variant length

    Note:
        - SNV: 长度为0（ref和alt长度相同）|SNV: length is 0 (ref and alt have same length)
        - Insertion: alt长度 - ref长度|Insertion: alt length - ref length
        - Deletion: ref长度 - alt长度|Deletion: ref length - alt length
        - 多个ALT时取第一个|For multiple ALTs, use the first one
    """
    if alt == '.' or not alt:
        return 0

    # 处理多个ALT的情况，取第一个|Handle multiple ALTs, use the first one
    first_alt = alt.split(',')[0] if ',' in alt else alt

    return abs(len(first_alt) - len(ref))


def is_biallelic(ref: str, alt: str) -> bool:
    """检查是否为双等位变异|Check if variant is biallelic"""
    return ',' not in alt


class GenotypeProcessor:
    """基因型处理器（纯过滤，无内存缓存）|Genotype Processor (filter only, no memory buffering)"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        # 预计算过滤条件，避免循环内重复判断|Pre-compute filter conditions
        self._need_length_filter = (config.min_length is not None
                                    or config.max_length is not None)
        self._min_length = config.min_length
        self._max_length = config.max_length

    def should_keep(self, row: Dict[str, Any]) -> bool:
        """判断该记录是否通过过滤条件|Check if record passes filter conditions"""
        # 检查是否只要双等位位点|Check if only biallelic sites are needed
        if self.config.biallelic_only and not is_biallelic(row['REF'], row['ALT']):
            return False

        # 检查长度过滤|Check length filter
        if self._need_length_filter:
            variant_length = calculate_variant_length(row['REF'], row['ALT'])
            if self._min_length is not None and variant_length < self._min_length:
                return False
            if self._max_length is not None and variant_length > self._max_length:
                return False

        return True
