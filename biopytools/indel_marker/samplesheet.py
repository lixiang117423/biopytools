"""samplesheet解析模块|Samplesheet Parser Module

契约|Contract:
    - 无需表头（若首行形如 `...<分隔>group<分隔>...` 会自动识别并跳过）|
      No header required (a first row whose 2nd column is literally "group" is
      auto-detected and skipped).
    - 分隔符默认 tab，行内无 tab 时回退为任意空白（空格/多空白）|
      Delimiter defaults to tab; falls back to arbitrary whitespace (spaces) when
      no tab is present in the line.
    - 分组不写死 resistant/susceptible：任意两种不同标签均可，自动映射为
      resistant / susceptible 两个角色（内部仍用规范名，下游无需改动）|
      Groups are not restricted to resistant/susceptible: any two distinct labels
      are accepted and auto-mapped to the resistant / susceptible roles (the
      canonical names are stored internally, so downstream code is unchanged).
"""

from dataclasses import dataclass
from typing import Dict, Tuple, List, Optional
import logging

from ..common.paths import expand_path


# 角色 -> 常见别名（小写匹配）|role -> common aliases (lowercase match)
_RESISTANT_ALIASES = {"resistant", "res", "r", "抗病", "抗"}
_SUSCEPTIBLE_ALIASES = {"susceptible", "sus", "s", "感病", "感"}
# 表头关键列名（小写）|header keyword columns (lowercase)
_HEADER_FIRST_COL = {"sample", "sample_name", "sample_id", "name"}


@dataclass
class SampleInfo:
    """样品信息|Sample Info"""
    name: str
    group: str       # 规范化为 resistant / susceptible|canonicalized to resistant/susceptible
    bam_path: str


def _split_cols(line: str) -> List[str]:
    """切列：默认tab，行内无tab时回退任意空白|Split cols: tab default, whitespace fallback"""
    # 默认tab：路径内含空格也安全；无tab时才按空白切（路径含空格需用tab）|
    # tab keeps spaces inside paths safe; whitespace split only when no tab present
    if '\t' in line:
        return line.split('\t')
    return line.split()


def _resolve_role(label: str) -> str:
    """将分组标签解析为角色|Resolve a group label to a role

    Returns: 'resistant' | 'susceptible' | 'unknown'
    """
    low = label.strip().lower()
    if low in _RESISTANT_ALIASES:
        return 'resistant'
    if low in _SUSCEPTIBLE_ALIASES:
        return 'susceptible'
    return 'unknown'


def _is_header_row(cols: List[str]) -> bool:
    """启发式识别表头行|Heuristic header detection（仅首行判断|first-row only）"""
    if len(cols) < 2:
        return False
    return cols[0].strip().lower() in _HEADER_FIRST_COL or cols[1].strip().lower() == 'group'


class SamplesheetParser:
    """samplesheet解析器|Samplesheet Parser"""

    @staticmethod
    def parse(path: str, logger: Optional[logging.Logger] = None) -> Dict[str, SampleInfo]:
        """
        解析samplesheet|Parse samplesheet

        格式|Format: ``sample_name <分隔> group <分隔> bam_path``
        - 无需表头（自动跳过形如 ``... group ...`` 的首行）|no header required
        - 分隔符默认tab，无tab时回退空白|delimiter: tab default, whitespace fallback
        - 分组：恰好两种不同标签，自动映射 resistant/susceptible|exactly 2 distinct groups

        Args:
            path: samplesheet路径|samplesheet path
            logger: 可选日志器，用于告警非标准分组标签的角色映射|
                    optional logger to warn on non-standard group role mapping

        Returns:
            sample_name -> SampleInfo（group已规范化）|sample_name -> SampleInfo (canonicalized)
        """
        path = expand_path(path)
        rows: List[List[str]] = []
        errors: List[str] = []
        header_checked = False

        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.rstrip('\n')
                if not line.strip() or line.lstrip().startswith('#'):
                    continue
                cols = _split_cols(line)
                if len(cols) < 3:
                    errors.append(f"列数不足|<3 columns: {line}")
                    continue
                # 仅对首条有效行做表头识别|header heuristic on the first valid row only
                if not header_checked:
                    header_checked = True
                    if _is_header_row(cols):
                        if logger is not None:
                            logger.info("检测到表头行并跳过|header row detected and skipped")
                        continue
                rows.append(cols)

        if errors:
            raise ValueError("\n".join(errors))
        if not rows:
            raise ValueError("samplesheet无有效数据行|no valid data rows in samplesheet")

        # 收集出现的分组标签（保持首次出现顺序）|distinct group labels in first-seen order
        seen_labels: List[str] = []
        for cols in rows:
            g = cols[1].strip()
            if g not in seen_labels:
                seen_labels.append(g)
        if len(seen_labels) != 2:
            raise ValueError(
                f"需要恰好2种不同分组，实际|need exactly 2 distinct groups, got "
                f"{len(seen_labels)}: {seen_labels}")

        # 两种标签 -> resistant/susceptible 角色|map 2 labels to resistant/susceptible roles
        label_to_role = SamplesheetParser._assign_roles(seen_labels, logger)

        # 构建samples|build samples
        samples: Dict[str, SampleInfo] = {}
        for cols in rows:
            name = cols[0].strip()
            group_label = cols[1].strip()
            bam = expand_path(cols[2].strip())
            if name in samples:
                errors.append(f"样品名重复|duplicate sample name: {name}")
                continue
            samples[name] = SampleInfo(
                name=name, group=label_to_role[group_label], bam_path=bam)

        if errors:
            raise ValueError("\n".join(errors))
        return samples

    @staticmethod
    def _assign_roles(labels: List[str],
                      logger: Optional[logging.Logger]) -> Dict[str, str]:
        """把两个标签映射为 resistant/susceptible|Map 2 labels to resistant/susceptible roles"""
        roles = [_resolve_role(l) for l in labels]
        # 两个标签都能识别为别名|both recognized as standard aliases
        if roles[0] != 'unknown' and roles[1] != 'unknown':
            if roles[0] == roles[1]:
                raise ValueError(
                    f"两个分组解析为同一角色|both groups map to same role "
                    f"'{roles[0]}': {labels}")
            return {l: r for l, r in zip(labels, roles)}
        # 否则按首次出现顺序赋角色并告警|otherwise assign by first-seen order with a warning
        if logger is not None:
            logger.warning(
                f"分组标签非标准别名，按出现顺序赋角色（第1组=抗病resistant，"
                f"第2组=感病susceptible）|non-standard group labels; roles assigned by "
                f"order (1st=resistant, 2nd=susceptible): "
                f"{labels[0]}->resistant, {labels[1]}->susceptible")
        return {labels[0]: 'resistant', labels[1]: 'susceptible'}

    @staticmethod
    def split_groups(samples: Dict[str, SampleInfo]) -> Tuple[List[str], List[str]]:
        """拆分抗/感两组样本名|Split into resistant/susceptible name lists"""
        r = sorted(s.name for s in samples.values() if s.group == "resistant")
        s = sorted(s.name for s in samples.values() if s.group == "susceptible")
        return r, s
