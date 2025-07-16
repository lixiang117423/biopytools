"""
输出格式化模块 | Output Formatting Module
"""

import csv
import os
from typing import List, Dict, Any
from pathlib import Path

class OutputFormatter:
    """输出格式化器 | Output Formatter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def _get_output_filename(self, suffix: str = "") -> str:
        """获取输出文件名 | Get output filename"""
        if suffix:
            filename = f"{self.config.output_prefix}_{suffix}"
        else:
            filename = f"{self.config.output_prefix}_all"
        
        # 添加文件扩展名 | Add file extension
        if self.config.output_type == "csv":
            filename += ".csv"
        elif self.config.output_type == "excel":
            filename += ".xlsx"
        else:  # txt
            filename += ".txt"
        
        return os.path.join(self.config.output_dir, filename)
    
    def _write_txt_file(self, data: List[Dict[str, Any]], filename: str):
        """写入TXT文件 | Write TXT file"""
        if not data:
            self.logger.warning(f"没有数据可写入 | No data to write: {filename}")
            return
        
        fieldnames = list(data[0].keys())
        
        with open(filename, 'w', encoding='utf-8') as f:
            # 写入头部 | Write header
            f.write('\t'.join(fieldnames) + '\n')
            
            # 写入数据 | Write data
            for row in data:
                values = [str(row.get(field, '')) for field in fieldnames]
                f.write('\t'.join(values) + '\n')
        
        self.logger.info(f"已保存TXT文件 | Saved TXT file: {filename} ({len(data)} variants)")
    
    def _write_csv_file(self, data: List[Dict[str, Any]], filename: str):
        """写入CSV文件 | Write CSV file"""
        if not data:
            self.logger.warning(f"没有数据可写入 | No data to write: {filename}")
            return
        
        fieldnames = list(data[0].keys())
        
        with open(filename, 'w', newline='', encoding='utf-8') as f:
            writer = csv.DictWriter(f, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(data)
        
        self.logger.info(f"已保存CSV文件 | Saved CSV file: {filename} ({len(data)} variants)")
    
    def _write_excel_file(self, data: List[Dict[str, Any]], filename: str):
        """写入Excel文件 | Write Excel file"""
        try:
            import pandas as pd
            
            if not data:
                self.logger.warning(f"没有数据可写入 | No data to write: {filename}")
                return
            
            df = pd.DataFrame(data)
            df.to_excel(filename, index=False)
            
            self.logger.info(f"已保存Excel文件 | Saved Excel file: {filename} ({len(data)} variants)")
            
        except ImportError:
            self.logger.error("pandas未安装，无法导出Excel格式 | pandas not installed, cannot export Excel format")
            # 回退到CSV格式 | Fallback to CSV format
            csv_filename = filename.replace('.xlsx', '.csv')
            self._write_csv_file(data, csv_filename)
    
    def write_output(self, data: List[Dict[str, Any]], suffix: str = ""):
        """写入输出文件 | Write output file"""
        filename = self._get_output_filename(suffix)
        
        if self.config.output_type == "csv":
            self._write_csv_file(data, filename)
        elif self.config.output_type == "excel":
            self._write_excel_file(data, filename)
        else:  # txt
            self._write_txt_file(data, filename)
    
    def write_summary(self, stats: Dict[str, Any]):
        """写入汇总文件 | Write summary file"""
        summary_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}_summary.txt")
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("VCF基因型提取汇总 | VCF Genotype Extraction Summary\n")
            f.write("=" * 50 + "\n\n")
            
            f.write(f"输入文件 | Input file: {self.config.vcf_file}\n")
            f.write(f"输出格式 | Output format: {self.config.output_type}\n")
            f.write(f"样本选择 | Sample selection: {self.config.samples}\n")
            f.write(f"只保留双等位位点 | Biallelic only: {self.config.biallelic_only}\n")
            f.write(f"按染色体拆分 | Split by chromosome: {self.config.split_by_chromosome}\n\n")
            
            f.write(f"总变异数 | Total variants: {stats['total_variants']}\n")
            f.write(f"染色体数 | Number of chromosomes: {len(stats['chromosomes'])}\n\n")
            
            f.write("各染色体变异数 | Variants per chromosome:\n")
            for chrom in sorted(stats['chromosome_counts'].keys()):
                count = stats['chromosome_counts'][chrom]
                f.write(f"  {chrom}: {count}\n")
        
        self.logger.info(f"已保存汇总文件 | Saved summary file: {summary_file}")
