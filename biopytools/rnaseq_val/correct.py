"""
转录本校正决策模块|Transcript Correction Decision Module

基于 GFFcompare class code + 定量指标，自动分级决策并输出校正后 GTF
"""

import os
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple, TYPE_CHECKING

from .utils import CommandRunner, format_number

if TYPE_CHECKING:
    from .config import RnaseqValConfig


# ============================================================
# GTF 读写工具|GTF I/O Utilities
# ============================================================

def read_gffcompare_annotated(annotated_gtf: str) -> List[dict]:
    """读取 GFFcompare 输出的 annotated GTF|Read GFFcompare annotated GTF

    GFFcompare 会在 GTF 的 attribute 列添加 class_code 和其他字段

    Args:
        annotated_gtf: annotated GTF 文件路径|annotated GTF file path

    Returns:
        List[dict]: 每条记录包含 seqname, source, feature, start, end, score,
                    strand, frame, attributes (原始字符串), class_code, transcript_id 等
                    Each record contains seqname, source, feature, start, end, score,
                    strand, frame, attributes (raw string), class_code, transcript_id, etc.
    """
    records = []

    if not os.path.exists(annotated_gtf):
        return records

    with open(annotated_gtf, "r") as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue

            parts = line.split("\t")
            if len(parts) < 9:
                continue

            attrs_str = parts[8]

            # 提取 class_code|Extract class_code
            class_code = "-"
            if "class_code " in attrs_str:
                for attr in attrs_str.split(";"):
                    attr = attr.strip()
                    if attr.startswith("class_code "):
                        class_code = attr.split('"')[1] if '"' in attr else attr.split()[1]
                        break

            # 提取 transcript_id|Extract transcript_id
            transcript_id = "-"
            for attr in attrs_str.split(";"):
                attr = attr.strip()
                if attr.startswith("transcript_id "):
                    transcript_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break

            # 提取 gene_id|Extract gene_id
            gene_id = "-"
            for attr in attrs_str.split(";"):
                attr = attr.strip()
                if attr.startswith("gene_id "):
                    gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                    break

            records.append({
                "seqname": parts[0],
                "source": parts[1],
                "feature": parts[2],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "score": parts[5],
                "strand": parts[6],
                "frame": parts[7],
                "attributes": attrs_str,
                "class_code": class_code,
                "transcript_id": transcript_id,
                "gene_id": gene_id,
                "line": line,  # 保留原始行|Keep original line
            })

    return records


def filter_gtf_by_transcripts(original_gtf: str, keep_transcript_ids: set, output_gtf: str) -> bool:
    """根据转录本 ID 过滤原始 GTF，输出校正后 GTF|Filter original GTF by transcript IDs

    Args:
        original_gtf: 原始注释 GTF 文件路径|Original annotation GTF path
        keep_transcript_ids: 需要保留的转录本 ID 集合|Set of transcript IDs to keep
        output_gtf: 输出 GTF 文件路径|Output GTF file path

    Returns:
        bool: 成功返回 True|True if succeeded
    """
    if not os.path.exists(original_gtf):
        return False

    os.makedirs(os.path.dirname(output_gtf), exist_ok=True)

    written_lines = []
    with open(original_gtf, "r") as f:
        for line in f:
            stripped = line.strip()

            # 保留注释行和空行|Keep comment and empty lines
            if not stripped or stripped.startswith("#"):
                written_lines.append(line)
                continue

            parts = stripped.split("\t")
            if len(parts) < 9:
                written_lines.append(line)
                continue

            attrs_str = parts[8]

            # 对于 transcript / exon / CDS / UTR 等 feature，检查关联的 transcript_id
            if parts[2] in ("transcript", "exon", "CDS", "UTR", "three_prime_utr", "five_prime_utr", "start_codon", "stop_codon"):
                # 提取 transcript_id
                tid = None
                for attr in attrs_str.split(";"):
                    attr = attr.strip()
                    if attr.startswith("transcript_id "):
                        tid = attr.split('"')[1] if '"' in attr else attr.split()[1]
                        break

                if tid and tid in keep_transcript_ids:
                    written_lines.append(line)
                elif parts[2] == "gene":
                    # gene feature: 检查是否有任何子 transcript 被保留
                    # 简单处理：也保留 gene 行（gene_id 匹配）
                    gene_id = None
                    for attr in attrs_str.split(";"):
                        attr = attr.strip()
                        if attr.startswith("gene_id "):
                            gene_id = attr.split('"')[1] if '"' in attr else attr.split()[1]
                            break
                    written_lines.append(line)
                # 其他 feature 类型不匹配则跳过
            elif parts[2] == "gene":
                written_lines.append(line)
            else:
                written_lines.append(line)

    with open(output_gtf, "w") as f:
        f.writelines(written_lines)

    return True


# ============================================================
# 校正决策器|Correction Decision Maker
# ============================================================

# 决策类别常量|Decision category constants
KEEP_HIGH_CONFIDENCE = "KEEP_HIGH_CONFIDENCE"    # 二代+三代均支持，高置信
KEEP_LR_SUPPORTED = "KEEP_LR_SUPPORTED"          # 三代全长支持
KEEP_SR_SUPPORTED = "KEEP_SR_SUPPORTED"          # 仅二代支持，TPM 充足
KEEP_WITH_FLAG = "KEEP_WITH_FLAG"                # class j/k，有支持
MANUAL_REVIEW = "MANUAL_REVIEW"                  # 需人工复核


def classify_transcripts(
    records: List[dict],
    sr_assembly_gtf: Optional[str] = None,
    lr_assembly_gtf: Optional[str] = None,
    thresholds: Optional[Dict] = None,
) -> List[dict]:
    """对所有转录本进行分级决策|Classify all transcripts into decision categories

    决策逻辑：
    - KEEP_HIGH_CONFIDENCE: class = 或 c，二代+三代均有组装支持
    - KEEP_LR_SUPPORTED: class = 或 c，仅三代有支持
    - KEEP_SR_SUPPORTED: class = 或 c，仅二代有支持且 TPM >= min_tpm_sr_only
    - KEEP_WITH_FLAG: class j 或 k，有任一数据支持且覆盖度达标
    - MANUAL_REVIEW: 其余

    Args:
        records: annotated GTF 解析后的记录列表|Parsed annotated GTF records
        sr_assembly_gtf: 二代组装 GTF（用于检查 SR 支持）|SR assembly GTF
        lr_assembly_gtf: 三代组装 GTF（用于检查 LR 支持）|LR assembly GTF
        thresholds: 阈值字典|Threshold dict

    Returns:
        List[dict]: 每条记录附加 decision 字段|Each record with decision field appended
    """
    if thresholds is None:
        thresholds = {
            "min_cov": 5.0,
            "min_tpm": 0.5,
            "min_tpm_sr_only": 1.0,
            "min_junction_ont": 5,
        }

    min_cov = thresholds.get("min_cov", 5.0)
    min_tpm = thresholds.get("min_tpm", 0.5)
    min_tpm_sr_only = thresholds.get("min_tpm_sr_only", 1.0)

    # 收集二代和三代组装中的转录本 ID
    sr_transcripts = _get_transcript_ids_from_gtf(sr_assembly_gtf) if sr_assembly_gtf else set()
    lr_transcripts = _get_transcript_ids_from_gtf(lr_assembly_gtf) if lr_assembly_gtf else set()

    for record in records:
        code = record["class_code"]
        tid = record["transcript_id"]

        has_sr = tid in sr_transcripts if sr_transcripts else False
        has_lr = tid in lr_transcripts if lr_transcripts else False

        # 从 attributes 提取 cov 和 TPM（如果存在）
        cov = _extract_attr_float(record["attributes"], "cov", 0.0)
        tpm = _extract_attr_float(record["attributes"], "TPM", 0.0)

        # 分级决策|Classification decision
        if code in ("=", "c"):
            if has_sr and has_lr:
                record["decision"] = KEEP_HIGH_CONFIDENCE
            elif has_lr:
                record["decision"] = KEEP_LR_SUPPORTED
            elif has_sr and tpm >= min_tpm_sr_only:
                record["decision"] = KEEP_SR_SUPPORTED
            elif has_sr:
                # 二代支持但 TPM 不够，降为人工复核
                record["decision"] = MANUAL_REVIEW
            else:
                record["decision"] = MANUAL_REVIEW

        elif code in ("j", "k"):
            if (has_sr or has_lr) and cov >= min_cov:
                record["decision"] = KEEP_WITH_FLAG
            elif has_lr:
                # 三代全长支持，即使覆盖度不足也可保留
                record["decision"] = KEEP_WITH_FLAG
            elif has_sr and tpm >= min_tpm:
                record["decision"] = KEEP_WITH_FLAG
            else:
                record["decision"] = MANUAL_REVIEW

        else:
            # class u, e, i, x 等
            record["decision"] = MANUAL_REVIEW

        # 附加支持信息|Append support info
        record["sr_support"] = has_sr
        record["lr_support"] = has_lr
        record["cov"] = cov
        record["tpm"] = tpm

    return records


def _get_transcript_ids_from_gtf(gtf_path: str) -> set:
    """从 GTF 文件提取所有 transcript_id|Extract all transcript IDs from GTF

    Args:
        gtf_path: GTF 文件路径|GTF file path

    Returns:
        set: transcript ID 集合|transcript ID set
    """
    tids = set()
    if not os.path.exists(gtf_path):
        return tids

    with open(gtf_path, "r") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 9 and parts[2] == "transcript":
                for attr in parts[8].split(";"):
                    attr = attr.strip()
                    if attr.startswith("transcript_id "):
                        tid = attr.split('"')[1] if '"' in attr else attr.split()[1]
                        tids.add(tid)
                        break

    return tids


def _extract_attr_float(attrs_str: str, attr_name: str, default: float = 0.0) -> float:
    """从 GTF attribute 字符串提取浮点数值|Extract float value from GTF attribute string

    Args:
        attrs_str: GTF attribute 列|GTF attribute column
        attr_name: 属性名|Attribute name
        default: 默认值|Default value

    Returns:
        float: 提取的浮点值|Extracted float value
    """
    for attr in attrs_str.split(";"):
        attr = attr.strip()
        if attr.startswith(f"{attr_name} "):
            try:
                return float(attr.split('"')[1] if '"' in attr else attr.split()[1])
            except (ValueError, IndexError):
                return default
    return default


# ============================================================
# 结果导出器|Result Exporter
# ============================================================

class CorrectionExporter:
    """校正结果导出器|Correction Result Exporter"""

    def __init__(self, config: 'RnaseqValConfig', logger: logging.Logger):
        """初始化|Initialize

        Args:
            config: 配置对象|Configuration object
            logger: logger 实例|Logger instance
        """
        self.config = config
        self.logger = logger

    def export_results(self, records: List[dict]) -> Dict[str, str]:
        """导出各类别结果|Export results by category

        Args:
            records: 含 decision 字段的记录列表|Records with decision field

        Returns:
            Dict[str, str]: 各类别输出文件路径字典|Dict of output file paths by category
        """
        out_dir = os.path.join(self.config.output_dir, "07_correction")
        os.makedirs(out_dir, exist_ok=True)

        # 按类别分组|Group by decision
        categories = {
            KEEP_HIGH_CONFIDENCE: [],
            KEEP_LR_SUPPORTED: [],
            KEEP_SR_SUPPORTED: [],
            KEEP_WITH_FLAG: [],
            MANUAL_REVIEW: [],
        }

        for record in records:
            decision = record.get("decision", MANUAL_REVIEW)
            if decision in categories:
                categories[decision].append(record)

        # 导出各类别 TSV|Export each category as TSV
        output_files = {}
        tsv_columns = ["transcript_id", "gene_id", "seqname", "class_code", "decision",
                        "sr_support", "lr_support", "cov", "tpm", "start", "end", "strand"]

        for cat_name, cat_records in categories.items():
            if not cat_records:
                continue

            tsv_file = os.path.join(out_dir, f"{cat_name.lower()}.tsv")
            with open(tsv_file, "w") as f:
                f.write("\t".join(tsv_columns) + "\n")
                for rec in cat_records:
                    row = [str(rec.get(col, "")) for col in tsv_columns]
                    f.write("\t".join(row) + "\n")

            output_files[cat_name] = tsv_file
            self.logger.info(f"导出类别|Exported category {cat_name}: {len(cat_records)} 条|records -> {tsv_file}")

        # 统计汇总|Summary statistics
        total = len(records)
        self.logger.info("=" * 50)
        self.logger.info(f"校正决策统计|Correction decision summary (total: {total})")
        self.logger.info("=" * 50)
        for cat_name in [KEEP_HIGH_CONFIDENCE, KEEP_LR_SUPPORTED, KEEP_SR_SUPPORTED, KEEP_WITH_FLAG, MANUAL_REVIEW]:
            count = len(categories[cat_name])
            pct = (count / total * 100) if total > 0 else 0
            self.logger.info(f"  {cat_name}: {count} ({pct:.1f}%)")
        self.logger.info("=" * 50)

        return output_files

    def write_corrected_gtf(self, records: List[dict], original_gtf: str) -> Optional[str]:
        """根据保留的转录本 ID 过滤原始注释 GTF|Filter original GTF by kept transcript IDs

        Args:
            records: 含 decision 字段的记录列表|Records with decision field
            original_gtf: 原始注释 GTF 文件路径|Original annotation GTF path

        Returns:
            Optional[str]: 校正后 GTF 文件路径，失败返回 None|Corrected GTF path, None on failure
        """
        out_dir = os.path.join(self.config.output_dir, "07_correction")
        os.makedirs(out_dir, exist_ok=True)
        corrected_gtf = os.path.join(out_dir, "corrected_annotation.gtf")

        # 收集需要保留的转录本 ID（所有非 MANUAL_REVIEW 的参考转录本）
        keep_ids = set()
        for rec in records:
            if rec.get("decision", MANUAL_REVIEW) != MANUAL_REVIEW:
                # 只保留参考转录本（class_code 不是 "-" 或 "u" 的）
                tid = rec.get("transcript_id", "")
                if tid and tid != "-":
                    # 从 annotated GTF 的 transcript_id 是 query 侧的
                    # 需要从 tracking 或 ref_id 获取参考侧 ID
                    # 这里保留所有非 MANUAL_REVIEW 的 query transcript_id
                    keep_ids.add(tid)

        if not keep_ids:
            self.logger.warning("无转录本被保留|No transcripts to keep")
            return None

        self.logger.info(f"保留 {len(keep_ids)} 个转录本|Keeping {len(keep_ids)} transcripts")

        # 过滤原始 GTF|Filter original GTF
        success = filter_gtf_by_transcripts(original_gtf, keep_ids, corrected_gtf)

        if success:
            self.logger.info(f"校正后 GTF 已写入|Corrected GTF written: {corrected_gtf}")
            return corrected_gtf
        else:
            self.logger.error("GTF 过滤失败|GTF filtering failed")
            return None
