"""VCF→DeepBSA CSV 核心转换逻辑|VCF to DeepBSA CSV core conversion logic

CSV 输出格式(无 header,与原 deepbsa vcf2csv 逐字节一致):
CHROM,POS,REF,ALT,AD_ref_s1,AD_alt_s1,AD_ref_s2,AD_alt_s2,...
|CSV output (no header, byte-identical to the original deepbsa vcf2csv):
CHROM,POS,REF,ALT,AD_ref_s1,AD_alt_s1,AD_ref_s2,AD_alt_s2,...
"""
import csv
import logging
from pathlib import Path
from typing import List, Optional, Tuple

from .config import Vcf2DeepBsaConfig
from .utils import VcfReader, format_number

# 跳过原因 → (stats键, 中文|英文标签)|Skip reason → (stats key, label)
SKIP_REASONS = {
    "no_ad": "无AD字段|no AD field",
    "incomplete_ad": "AD字段不完整|incomplete AD",
    "bad_value": "AD/POS值非法|invalid AD/POS value",
    "column_mismatch": "样本列数不一致|inconsistent sample columns",
}


def make_output_name(vcf_path: str) -> str:
    """由输入名生成输出 CSV 名|Derive output CSV name from input name

    input.vcf / input.vcf.gz → input_deepbsa.csv
    """
    name = Path(vcf_path).name
    if name.endswith(".vcf.gz"):
        stem = name[:-len(".vcf.gz")]
    elif name.endswith(".vcf"):
        stem = name[:-len(".vcf")]
    else:
        stem = Path(name).stem
    return f"{stem}_deepbsa.csv"


class Vcf2DeepBsaConverter:
    """VCF→DeepBSA CSV 转换器|VCF to DeepBSA CSV converter"""

    def __init__(self, config: Vcf2DeepBsaConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.stats = {
            "total": 0,
            "written": 0,
            "sample_count": 0,
            "skipped_no_ad": 0,
            "skipped_incomplete_ad": 0,
            "skipped_bad_value": 0,
            "skipped_column_mismatch": 0,
        }

    def run(self) -> str:
        """执行转换,返回 CSV 路径|Run conversion, return CSV path"""
        csv_path = self.config.output_path / make_output_name(self.config.input_vcf)
        self.logger.info("开始 VCF→DeepBSA CSV 转换|Starting VCF to DeepBSA CSV conversion")
        self.logger.info(f"输入文件|Input: {self.config.input_vcf}")
        self.logger.info(f"输出文件|Output: {csv_path}")

        expected_cols: Optional[int] = None
        # lineterminator=\n 与原 pandas to_csv 在 Linux 上的输出逐字节一致
        # |lineterminator=\n keeps output byte-identical to pandas to_csv on Linux
        with VcfReader(self.config.input_vcf) as reader, \
                open(csv_path, "w", newline="") as out:
            writer = csv.writer(out, lineterminator="\n")
            if reader.samples:
                self.logger.info(
                    f"样本列表|Samples: {', '.join(reader.samples)}")
            for record in reader:
                self.stats["total"] += 1
                row, reason = self._record_to_row(record)
                if row is None:
                    self.stats[f"skipped_{reason}"] += 1
                    continue
                if expected_cols is None:
                    expected_cols = len(row)
                    self.stats["sample_count"] = (len(row) - 4) // 2
                elif len(row) != expected_cols:
                    self.stats["skipped_column_mismatch"] += 1
                    continue
                writer.writerow(row)
                self.stats["written"] += 1

        self._log_stats()
        return str(csv_path)

    def _record_to_row(self, record) -> Tuple[Optional[List], Optional[str]]:
        """单条记录转 CSV 行,失败返回 (None, 跳过原因)|Convert one record; (None, reason) on failure"""
        try:
            row = [record.CHROM, int(record.POS), record.REF, record.ALT]
        except ValueError:
            return None, "bad_value"
        for sample_gt in record.GT:
            ad = sample_gt.get("AD")
            if ad is None:
                return None, "no_ad"
            values = ad.split(",")
            if len(values) < 2:
                return None, "incomplete_ad"
            try:
                row.append(int(values[0]))
                row.append(int(values[1]))
            except ValueError:
                return None, "bad_value"
        return row, None

    def _log_stats(self):
        """输出转换统计|Log conversion statistics"""
        s = self.stats
        self.logger.info("=" * 60)
        self.logger.info(f"总记录数|Total records: {format_number(s['total'])}")
        self.logger.info(f"写入记录数|Written records: {format_number(s['written'])}")
        self.logger.info(f"样本数|Sample count: {s['sample_count']}")
        skipped = s["total"] - s["written"]
        self.logger.info(f"跳过总数|Skipped total: {format_number(skipped)}")
        for reason, label in SKIP_REASONS.items():
            count = s[f"skipped_{reason}"]
            if count:
                self.logger.warning(
                    f"跳过-{label}: {format_number(count)}")
        if s["total"] > 0 and s["written"] == 0:
            self.logger.warning(
                "无有效记录输出,请检查 VCF 的 FORMAT 是否含 AD 字段|"
                "No valid records written; check that VCF FORMAT contains AD")
