"""
输出格式化模块|Output Formatting Module
"""

import csv
import os
from typing import List, Dict, Any, Optional
from pathlib import Path


class OutputFormatter:
    """输出格式化器（支持流式写入）|Output Formatter (supports streaming write)"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def _get_output_filename(self, suffix: str = "") -> str:
        """获取输出文件名|Get output filename"""
        if suffix:
            filename = f"{self.config.output_prefix}_{suffix}"
        else:
            filename = f"{self.config.output_prefix}_all"

        # 添加文件扩展名|Add file extension
        if self.config.output_type == "csv":
            filename += ".csv"
        elif self.config.output_type == "excel":
            filename += ".xlsx"
        else:  # txt
            filename += ".txt"

        return os.path.join(self.config.output_dir, filename)

    # ---- 流式写入接口|Streaming write interface ----

    def open_stream(self, suffix: str = "", fieldnames: Optional[List[str]] = None,
                    defer_header: bool = False):
        """
        打开流式写出文件|Open streaming output file

        Args:
            suffix: 文件名后缀（如染色体名）|Filename suffix (e.g., chromosome name)
            fieldnames: 列名列表|Column names
            defer_header: 延迟写header，需手动调用flush_header()|Defer header, call flush_header() manually
        """
        filename = self._get_output_filename(suffix)
        self._stream_file = open(filename, 'w', encoding='utf-8', newline='')
        self._stream_fieldnames = fieldnames
        self._header_written = False

        if fieldnames and not defer_header:
            if self.config.output_type == "csv":
                self._csv_writer = csv.DictWriter(self._stream_file, fieldnames=fieldnames)
                self._csv_writer.writeheader()
            else:
                self._stream_file.write('\t'.join(fieldnames) + '\n')
            self._header_written = True

        return self

    def update_fieldnames(self, fieldnames: List[str]):
        """更新列名并写入header（延迟模式用）|Update fieldnames and write header (for deferred mode)"""
        self._stream_fieldnames = fieldnames
        if self.config.output_type == "csv":
            self._csv_writer = csv.DictWriter(self._stream_file, fieldnames=fieldnames)
            self._csv_writer.writeheader()
        else:
            self._stream_file.write('\t'.join(fieldnames) + '\n')
        self._header_written = True

    def write_row(self, row: Dict[str, Any]):
        """写入单行数据|Write a single data row"""
        if self.config.output_type == "csv":
            self._csv_writer.writerow(row)
        else:
            # txt: 按 fieldnames 顺序输出|txt: output in fieldnames order
            values = [str(row.get(f, '')) for f in self._stream_fieldnames]
            self._stream_file.write('\t'.join(values) + '\n')

    def close_stream(self, total_rows: int):
        """关闭流式写出文件|Close streaming output file"""
        self._stream_file.close()
        filename = self._stream_file.name
        self.logger.info(f"已保存文件|Saved file: {filename} ({total_rows} variants)")

    # ---- 批量写入接口（保留兼容）|Batch write interface (kept for compatibility) ----

    def _write_txt_file(self, data: List[Dict[str, Any]], filename: str):
        """写入TXT文件|Write TXT file"""
        if not data:
            self.logger.warning(f"没有数据可写入|No data to write: {filename}")
            return

        fieldnames = list(data[0].keys())

        with open(filename, 'w', encoding='utf-8') as f:
            f.write('\t'.join(fieldnames) + '\n')
            for row in data:
                values = [str(row.get(field, '')) for field in fieldnames]
                f.write('\t'.join(values) + '\n')

        self.logger.info(f"已保存TXT文件|Saved TXT file: {filename} ({len(data)} variants)")

    def _write_csv_file(self, data: List[Dict[str, Any]], filename: str):
        """写入CSV文件|Write CSV file"""
        if not data:
            self.logger.warning(f"没有数据可写入|No data to write: {filename}")
            return

        fieldnames = list(data[0].keys())

        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)

        self.logger.info(f"已保存CSV文件|Saved CSV file: {filename} ({len(data)} variants)")

    def _write_excel_file(self, data: List[Dict[str, Any]], filename: str):
        """写入Excel文件|Write Excel file"""
        try:
            import pandas as pd
            if not data:
                self.logger.warning(f"没有数据可写入|No data to write: {filename}")
                return

            df = pd.DataFrame(data)

            with pd.ExcelWriter(filename, engine='openpyxl') as writer:
                df.to_excel(writer, index=False)
                worksheet = writer.book.active

                for row in worksheet.iter_rows():
                    for cell in row:
                        cell.number_format = '@'

            self.logger.info(f"已保存Excel文件|Saved Excel file: {filename} ({len(data)} variants)")

        except ImportError:
            self.logger.error("pandas未安装，无法导出Excel格式|pandas not installed, cannot export Excel format")
            csv_filename = filename.replace('.xlsx', '.csv')
            self._write_csv_file(data, csv_filename)

    def write_output(self, data: List[Dict[str, Any]], suffix: str = ""):
        """写入输出文件|Write output file"""
        filename = self._get_output_filename(suffix)

        if self.config.output_type == "csv":
            self._write_csv_file(data, filename)
        elif self.config.output_type == "excel":
            self._write_excel_file(data, filename)
        else:  # txt
            self._write_txt_file(data, filename)

    def write_summary(self, stats: Dict[str, Any]):
        """写入汇总文件|Write summary file"""
        summary_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_summary.txt")

        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("VCF基因型提取汇总|VCF Genotype Extraction Summary\n")
            f.write("=" * 50 + "\n\n")

            f.write(f"输入文件|Input file: {self.config.vcf_file}\n")
            f.write(f"输出格式|Output format: {self.config.output_type}\n")
            f.write(f"样本选择|Sample selection: {self.config.samples}\n")
            f.write(f"只保留双等位位点|Biallelic only: {self.config.biallelic_only}\n")
            f.write(f"按染色体拆分|Split by chromosome: {self.config.split_by_chromosome}\n\n")

            f.write(f"总变异数|Total variants: {stats['total_variants']}\n")
            f.write(f"染色体数|Number of chromosomes: {len(stats['chromosomes'])}\n\n")

            f.write("各染色体变异数|Variants per chromosome:\n")
            for chrom in sorted(stats['chromosome_counts'].keys()):
                count = stats['chromosome_counts'][chrom]
                f.write(f"  {chrom}: {count}\n")

        self.logger.info(f"已保存汇总文件|Saved summary file: {summary_file}")
