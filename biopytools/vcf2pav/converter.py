"""VCF→PAV矩阵核心转换逻辑|VCF→PAV matrix core conversion logic"""

import csv
import logging
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .config import Vcf2PavConfig
from .utils import decode_gt, make_sv_id


class Vcf2PavConverter:
    """VCF→PAV矩阵转换器|VCF→PAV matrix converter"""

    def __init__(self, config: Vcf2PavConfig, logger: logging.Logger):
        self.config = config
        self.logger = logger

    def run(self) -> str:
        """执行转换,返回输出文件路径|Run conversion, return output path.

        Returns:
            pav_matrix.tsv 的完整路径|Full path to pav_matrix.tsv
        """
        self.logger.info("开始 VCF→PAV 转换|Starting VCF→PAV conversion")
        self.logger.info(f"输入文件|Input: {self.config.input_vcf}")

        # 第一遍扫描: 收集样本列表|First pass: collect sample list
        samples = self._parse_header()
        self.logger.info(f"检测到样本数|Detected samples: {len(samples)}")

        # 第二遍扫描: 解析 SV 记录|Second pass: parse SV records
        sv_rows = self._parse_records(samples)
        self.logger.info(f"解析SV记录数|Parsed SV records: {len(sv_rows)}")

        # 写入矩阵 TSV|Write matrix TSV
        output_path = str(self.config.output_path / "pav_matrix.tsv")
        self._write_matrix(output_path, samples, sv_rows)
        self.logger.info(f"PAV矩阵已写入|PAV matrix written: {output_path}")

        # 写入统计摘要|Write summary
        summary_path = str(self.config.output_path / "pav_summary.tsv")
        self._write_summary(summary_path, samples, sv_rows)
        self.logger.info(f"统计摘要已写入|Summary written: {summary_path}")

        return output_path

    def _parse_header(self) -> List[str]:
        """解析 VCF header 获取样本名列表|Parse VCF header for sample names.

        Returns:
            样本名列表(按 VCF 列顺序)|List of sample names (VCF column order)
        """
        with open(self.config.input_vcf) as fh:
            for line in fh:
                line = line.rstrip("\n")
                if line.startswith("#CHROM"):
                    fields = line.split("\t")
                    # FORMAT 是第 9 列(索引8),样本从第10列(索引9)开始
                    samples = fields[9:]
                    self.logger.info(f"样本列表|Samples: {', '.join(samples)}")
                    return samples
        raise ValueError(f"VCF 文件中未找到 #CHROM header|No #CHROM header in VCF: {self.config.input_vcf}")

    def _parse_records(self, samples: List[str]) -> List[Dict]:
        """解析 VCF 数据行,构建 SV→样本 PAV 矩阵|Parse VCF records, build SV→sample PAV.

        Args:
            samples: 样本名列表|List of sample names

        Returns:
            [{SV_ID, SVTYPE, SVLEN, CHROM, START, END, S1: v, S2: v, ...}, ...]
        """
        rows = []
        unparseable_count = 0

        with open(self.config.input_vcf) as fh:
            for lineno, line in enumerate(fh, 1):
                line = line.rstrip("\n")
                if not line.strip() or line.startswith("#"):
                    continue

                fields = line.split("\t")
                if len(fields) < 10:
                    continue

                chrom = fields[0]
                pos = int(fields[1])
                info = fields[7]
                fmt = fields[8]
                sample_fields = fields[9:]

                # 解析 INFO|Parse INFO
                info_dict = {}
                for entry in info.split(";"):
                    if "=" in entry:
                        k, v = entry.split("=", 1)
                        info_dict[k] = v
                    else:
                        info_dict[entry] = True

                svtype = info_dict.get("SVTYPE", "UNK")
                end = int(info_dict.get("END", pos))
                svlen_str = info_dict.get("SVLEN", ".")
                chr2 = info_dict.get("CHR2", "")

                # 解析 GT 在 FORMAT 中的位置|Find GT position in FORMAT
                fmt_keys = fmt.split(":")
                try:
                    gt_idx = fmt_keys.index("GT")
                except ValueError:
                    # GT 不在 FORMAT 中,跳过|GT not in FORMAT, skip
                    continue

                # 构建行数据|Build row data
                sv_id = make_sv_id(chrom, pos, end, svtype, chr2)
                row = {
                    "SV_ID": sv_id,
                    "SVTYPE": svtype,
                    "SVLEN": svlen_str,
                    "CHROM": chrom,
                    "START": str(pos),
                    "END": str(end),
                }

                # 解析每个样本的 GT|Parse GT for each sample
                for i, sample_name in enumerate(samples):
                    if i >= len(sample_fields):
                        row[sample_name] = "NaN"
                        continue
                    sample_vals = sample_fields[i].split(":")
                    gt_val = sample_vals[gt_idx] if gt_idx < len(sample_vals) else "."
                    decoded = decode_gt(gt_val)
                    if decoded is None:
                        unparseable_count += 1
                        row[sample_name] = "NaN"
                    else:
                        row[sample_name] = str(decoded)

                rows.append(row)

        if unparseable_count > 0:
            self.logger.warning(
                f"无法解析的GT值数量|Unparseable GT count: {unparseable_count}")
        return rows

    def _write_matrix(self, output_path: str, samples: List[str],
                      rows: List[Dict]) -> None:
        """写入 PAV 矩阵 TSV|Write PAV matrix TSV.

        列顺序: SV_ID, SVTYPE, SVLEN, CHROM, START, END, Sample1, Sample2, ...
        """
        columns = ["SV_ID", "SVTYPE", "SVLEN", "CHROM", "START", "END"] + samples
        with open(output_path, "w", newline="") as fh:
            writer = csv.DictWriter(fh, fieldnames=columns, delimiter="\t",
                                    extrasaction="ignore")
            writer.writeheader()
            writer.writerows(rows)

    def _write_summary(self, summary_path: str, samples: List[str],
                       rows: List[Dict]) -> None:
        """写入 PAV 统计摘要|Write PAV summary.

        每样本各 SVTYPE 计数、存在/缺失总数|Per-sample SVTYPE counts, present/absent totals.
        """
        from collections import Counter

        header = ["sample", "total_sv", "present", "absent",
                  "DEL", "INS", "INV", "DUP", "BND", "TRA", "UNK"]
        with open(summary_path, "w", newline="") as fh:
            writer = csv.writer(fh, delimiter="\t")
            writer.writerow(header)
            for sample in samples:
                svtype_counts = Counter()
                present = 0
                absent = 0
                for row in rows:
                    svtype_counts[row["SVTYPE"]] += 1
                    val = row.get(sample, "NaN")
                    if val == "1":
                        present += 1
                    elif val == "0":
                        absent += 1
                writer.writerow([
                    sample,
                    len(rows),
                    present,
                    absent,
                    svtype_counts.get("DEL", 0),
                    svtype_counts.get("INS", 0),
                    svtype_counts.get("INV", 0),
                    svtype_counts.get("DUP", 0),
                    svtype_counts.get("BND", 0),
                    svtype_counts.get("TRA", 0),
                    svtype_counts.get("UNK", 0),
                ])
