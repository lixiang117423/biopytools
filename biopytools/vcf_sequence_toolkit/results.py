"""
序列提取结果处理模块 | Sequence Extraction Results Processing Module
"""

import os
import csv
from typing import Dict
from .utils import SequenceFormatter

class SequenceExporter:
    """序列导出器 | Sequence Exporter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.formatter = SequenceFormatter(logger)
    
    def export_sequences(self, sample_sequences: Dict[str, str], reference_seq: str = ""):
        """导出序列到文件 | Export sequences to file"""
        if not sample_sequences:
            self.logger.warning("没有序列数据需要导出 | No sequence data to export")
            return False
        
        try:
            if self.config.export_format == "tab":
                return self._export_tab_format(sample_sequences, reference_seq)
            elif self.config.export_format == "fasta":
                return self._export_fasta_format(sample_sequences, reference_seq)
            elif self.config.export_format == "csv":
                return self._export_csv_format(sample_sequences, reference_seq)
            else:
                self.logger.error(f"不支持的导出格式 | Unsupported export format: {self.config.export_format}")
                return False
        
        except Exception as e:
            self.logger.error(f"导出序列失败 | Failed to export sequences: {e}")
            return False
    
    def _export_tab_format(self, sample_sequences: Dict[str, str], reference_seq: str) -> bool:
        """导出为制表符分隔格式 | Export as tab-delimited format"""
        output_file = os.path.join(
            self.config.output_dir, 
            f"{self.config.chrom}_{self.config.start}_{self.config.end}_sequences.txt"
        )
        
        with open(output_file, 'w') as f:
            # 写入头部信息 | Write header information
            f.write(f"# Sequences for {self.config.chrom}:{self.config.start}-{self.config.end}\n")
            f.write(f"# Region length: {self.config.end - self.config.start + 1} bp\n")
            f.write(f"# Export format: tab-delimited\n")
            f.write(f"# Use first allele: {self.config.use_first_allele}\n")
            
            if reference_seq and self.config.include_reference:
                f.write(f"# Reference sequence: {reference_seq}\n")
            
            f.write("#\n")
            f.write("Sample\tSequence\tLength\n")
            
            # 写入参考序列 | Write reference sequence
            if reference_seq and self.config.include_reference:
                f.write(f"Reference\t{reference_seq}\t{len(reference_seq)}\n")
            
            # 写入样品序列 | Write sample sequences
            for sample_name in sorted(sample_sequences.keys()):
                sequence = sample_sequences[sample_name]
                f.write(f"{sample_name}\t{sequence}\t{len(sequence)}\n")
        
        self.logger.info(f"序列已导出为制表符格式 | Sequences exported in tab format: {output_file}")
        return True
    
    def _export_fasta_format(self, sample_sequences: Dict[str, str], reference_seq: str) -> bool:
        """导出为FASTA格式 | Export as FASTA format"""
        output_file = os.path.join(
            self.config.output_dir, 
            f"{self.config.chrom}_{self.config.start}_{self.config.end}_sequences.fasta"
        )
        
        with open(output_file, 'w') as f:
            # 写入参考序列 | Write reference sequence
            if reference_seq and self.config.include_reference:
                header = f">Reference_{self.config.chrom}:{self.config.start}-{self.config.end}"
                f.write(f"{header}\n")
                f.write(f"{self.formatter.wrap_sequence(reference_seq)}\n")
            
            # 写入样品序列 | Write sample sequences
            for sample_name in sorted(sample_sequences.keys()):
                sequence = sample_sequences[sample_name]
                header = f">{self.formatter.format_sequence_name(sample_name, self.config.chrom, self.config.start, self.config.end)}"
                f.write(f"{header}\n")
                f.write(f"{self.formatter.wrap_sequence(sequence)}\n")
        
        self.logger.info(f"序列已导出为FASTA格式 | Sequences exported in FASTA format: {output_file}")
        return True
    
    def _export_csv_format(self, sample_sequences: Dict[str, str], reference_seq: str) -> bool:
        """导出为CSV格式 | Export as CSV format"""
        output_file = os.path.join(
            self.config.output_dir, 
            f"{self.config.chrom}_{self.config.start}_{self.config.end}_sequences.csv"
        )
        
        with open(output_file, 'w', newline='') as f:
            writer = csv.writer(f)
            
            # 写入头部 | Write header
            writer.writerow(['Sample', 'Sequence', 'Length', 'Region'])
            
            region_info = f"{self.config.chrom}:{self.config.start}-{self.config.end}"
            
            # 写入参考序列 | Write reference sequence
            if reference_seq and self.config.include_reference:
                writer.writerow(['Reference', reference_seq, len(reference_seq), region_info])
            
            # 写入样品序列 | Write sample sequences
            for sample_name in sorted(sample_sequences.keys()):
                sequence = sample_sequences[sample_name]
                writer.writerow([sample_name, sequence, len(sequence), region_info])
        
        self.logger.info(f"序列已导出为CSV格式 | Sequences exported in CSV format: {output_file}")
        return True

class StatisticsGenerator:
    """统计信息生成器 | Statistics Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.formatter = SequenceFormatter(logger)
    
    def generate_statistics_report(self, sample_sequences: Dict[str, str], 
                                 reference_seq: str = "", variants_count: int = 0):
        """生成统计报告 | Generate statistics report"""
        if not sample_sequences:
            self.logger.warning("没有序列数据用于统计 | No sequence data for statistics")
            return False
        
        try:
            # 计算序列统计 | Calculate sequence statistics
            stats = self.formatter.calculate_sequence_stats(sample_sequences)
            
            # 生成报告文件 | Generate report file
            report_file = os.path.join(
                self.config.output_dir, 
                f"{self.config.chrom}_{self.config.start}_{self.config.end}_statistics.txt"
            )
            
            with open(report_file, 'w') as f:
                f.write("序列提取统计报告 | Sequence Extraction Statistics Report\n")
                f.write("=" * 60 + "\n\n")
                
                # 基本信息 | Basic information
                f.write("基本信息 | Basic Information:\n")
                f.write(f"  区间 | Region: {self.config.chrom}:{self.config.start}-{self.config.end}\n")
                f.write(f"  长度 | Length: {self.config.end - self.config.start + 1} bp\n")
                f.write(f"  样品数量 | Sample count: {len(sample_sequences)}\n")
                f.write(f"  变异数量 | Variant count: {variants_count}\n")
                f.write(f"  使用第一等位基因 | Use first allele: {self.config.use_first_allele}\n\n")
                
                # 参考序列信息 | Reference sequence information
                if reference_seq:
                    ref_stats = self.formatter.calculate_sequence_stats({"Reference": reference_seq})["Reference"]
                    f.write("参考序列统计 | Reference Sequence Statistics:\n")
                    f.write(f"  A: {ref_stats['A']} ({ref_stats['A']/ref_stats['total']*100:.1f}%)\n")
                    f.write(f"  T: {ref_stats['T']} ({ref_stats['T']/ref_stats['total']*100:.1f}%)\n")
                    f.write(f"  C: {ref_stats['C']} ({ref_stats['C']/ref_stats['total']*100:.1f}%)\n")
                    f.write(f"  G: {ref_stats['G']} ({ref_stats['G']/ref_stats['total']*100:.1f}%)\n")
                    f.write(f"  GC含量 | GC content: {ref_stats['GC_percent']}%\n\n")
                
                # 样品序列统计汇总 | Sample sequence statistics summary
                f.write("样品序列统计汇总 | Sample Sequence Statistics Summary:\n")
                
                total_a = sum(stat['A'] for stat in stats.values())
                total_t = sum(stat['T'] for stat in stats.values())
                total_c = sum(stat['C'] for stat in stats.values())
                total_g = sum(stat['G'] for stat in stats.values())
                total_n = sum(stat['N'] for stat in stats.values())
                total_missing = sum(stat['-'] for stat in stats.values())
                total_bases = sum(stat['total'] for stat in stats.values())
                
                f.write(f"  总碱基数 | Total bases: {total_bases}\n")
                f.write(f"  A: {total_a} ({total_a/total_bases*100:.1f}%)\n")
                f.write(f"  T: {total_t} ({total_t/total_bases*100:.1f}%)\n")
                f.write(f"  C: {total_c} ({total_c/total_bases*100:.1f}%)\n")
                f.write(f"  G: {total_g} ({total_g/total_bases*100:.1f}%)\n")
                f.write(f"  N: {total_n} ({total_n/total_bases*100:.1f}%)\n")
                f.write(f"  缺失 | Missing: {total_missing} ({total_missing/total_bases*100:.1f}%)\n\n")
                
                # 详细样品统计 | Detailed sample statistics
                f.write("详细样品统计 | Detailed Sample Statistics:\n")
                f.write("Sample\tA\tT\tC\tG\tN\tMissing\tGC%\n")
                
                for sample_name in sorted(stats.keys()):
                    stat = stats[sample_name]
                    f.write(f"{sample_name}\t{stat['A']}\t{stat['T']}\t{stat['C']}\t"
                           f"{stat['G']}\t{stat['N']}\t{stat['-']}\t{stat['GC_percent']}\n")
            
            self.logger.info(f"统计报告已生成 | Statistics report generated: {report_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成统计报告失败 | Failed to generate statistics report: {e}")
            return False

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, sample_count: int, variants_count: int, 
                              sequence_length: int, output_files: list):
        """生成总结报告 | Generate summary report"""
        report_file = os.path.join(self.config.output_dir, "extraction_summary.txt")
        
        try:
            with open(report_file, 'w', encoding='utf-8') as f:
                f.write("序列提取总结报告 | Sequence Extraction Summary Report\n")
                f.write("=" * 50 + "\n\n")
                
                # 输入信息 | Input information
                f.write("输入信息 | Input Information:\n")
                f.write(f"  VCF文件 | VCF file: {self.config.vcf_file}\n")
                f.write(f"  基因组文件 | Genome file: {self.config.genome_file}\n")
                f.write(f"  提取区间 | Extraction region: {self.config.chrom}:{self.config.start}-{self.config.end}\n")
                f.write(f"  区间长度 | Region length: {sequence_length} bp\n\n")
                
                # 处理结果 | Processing results
                f.write("处理结果 | Processing Results:\n")
                f.write(f"  样品数量 | Sample count: {sample_count}\n")
                f.write(f"  变异数量 | Variant count: {variants_count}\n")
                f.write(f"  使用等位基因 | Allele used: {'第一个 | First' if self.config.use_first_allele else '第二个 | Second'}\n")
                f.write(f"  包含参考序列 | Include reference: {'是 | Yes' if self.config.include_reference else '否 | No'}\n")
                f.write(f"  导出格式 | Export format: {self.config.export_format}\n\n")
                
                # 过滤参数 | Filtering parameters
                f.write("过滤参数 | Filtering Parameters:\n")
                if self.config.min_qual:
                    f.write(f"  最小质量值 | Minimum quality: {self.config.min_qual}\n")
                if self.config.sample_list:
                    f.write(f"  指定样品数量 | Specified samples: {len(self.config.sample_list)}\n")
                if self.config.exclude_samples:
                    f.write(f"  排除样品数量 | Excluded samples: {len(self.config.exclude_samples)}\n")
                f.write("\n")
                
                # 输出文件 | Output files
                f.write("输出文件 | Output Files:\n")
                for file_path in output_files:
                    f.write(f"  - {file_path}\n")
            
            self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
            
        except Exception as e:
            self.logger.error(f"生成总结报告失败 | Failed to generate summary report: {e}")
