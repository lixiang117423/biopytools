"""
HiFiasm数据处理模块 | HiFiasm Data Processing Module
"""

import os
import logging
import json
import pandas as pd
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import gzip
import subprocess

class StatisticsCalculator:
    """统计信息计算器 | Statistics Calculator"""
    
    def __init__(self, config, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.stats_dir = Path(config.output_dir) / 'statistics'
        self.stats_dir.mkdir(parents=True, exist_ok=True)
    
    def calculate_all_statistics(self) -> bool:
        """计算所有统计信息 | Calculate all statistics"""
        try:
            self.logger.info("计算组装统计信息 | Calculating assembly statistics")
            
            # 获取组装文件
            from .assembly import GFAConverter
            converter = GFAConverter(self.config, self.logger, None)
            assembly_files = converter.get_converted_files()
            
            if not assembly_files:
                self.logger.error("未找到组装文件 | No assembly files found")
                return False
            
            # 计算每个组装的统计信息
            all_stats = {}
            for assembly_type, file_path in assembly_files.items():
                stats = self._calculate_assembly_statistics(file_path, assembly_type)
                if stats:
                    all_stats[assembly_type] = stats
            
            # 生成统计报告
            self._generate_statistics_report(all_stats)
            
            # 生成比较表格
            self._generate_comparison_table(all_stats)
            
            # 如果启用绘图，生成统计图表
            if self.config.generate_plots:
                self._generate_statistics_plots(all_stats)
            
            self.logger.success(f"统计信息计算完成，处理了 {len(all_stats)} 个组装 | Statistics calculation completed for {len(all_stats)} assemblies")
            return True
            
        except Exception as e:
            self.logger.error(f"统计信息计算失败 | Statistics calculation failed: {e}")
            return False
    
    def _calculate_assembly_statistics(self, fasta_path: Path, assembly_type: str) -> Optional[Dict]:
        """计算单个组装的统计信息 | Calculate statistics for single assembly"""
        try:
            self.logger.info(f"计算 {assembly_type} 统计信息 | Calculating statistics for {assembly_type}")
            
            sequences = []
            sequence_info = []
            gc_counts = []
            n_counts = []
            current_seq = ""
            current_id = ""
            
            # 读取FASTA文件
            with open(fasta_path, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line.startswith('>'):
                        if current_seq:
                            seq_len = len(current_seq)
                            sequences.append(seq_len)
                            
                            # 计算GC和N含量
                            gc_count = current_seq.upper().count('G') + current_seq.upper().count('C')
                            n_count = current_seq.upper().count('N')
                            
                            sequence_info.append({
                                'id': current_id,
                                'length': seq_len,
                                'gc_content': (gc_count / seq_len * 100) if seq_len > 0 else 0,
                                'n_content': (n_count / seq_len * 100) if seq_len > 0 else 0
                            })
                            
                            gc_counts.append(gc_count)
                            n_counts.append(n_count)
                        
                        current_id = line[1:]
                        current_seq = ""
                    else:
                        current_seq += line.upper()
                
                # 处理最后一个序列
                if current_seq:
                    seq_len = len(current_seq)
                    sequences.append(seq_len)
                    
                    gc_count = current_seq.count('G') + current_seq.count('C')
                    n_count = current_seq.count('N')
                    
                    sequence_info.append({
                        'id': current_id,
                        'length': seq_len,
                        'gc_content': (gc_count / seq_len * 100) if seq_len > 0 else 0,
                        'n_content': (n_count / seq_len * 100) if seq_len > 0 else 0
                    })
                    
                    gc_counts.append(gc_count)
                    n_counts.append(n_count)
            
            if not sequences:
                self.logger.warning(f"未找到序列在文件中 | No sequences found in file: {fasta_path}")
                return None
            
            # 计算统计信息
            stats = self._compute_assembly_metrics(sequences, sequence_info, gc_counts, n_counts)
            stats['assembly_type'] = assembly_type
            stats['file_path'] = str(fasta_path)
            
            return stats
            
        except Exception as e:
            self.logger.error(f"计算统计信息失败 | Failed to calculate statistics: {e}")
            return None
    
    def _compute_assembly_metrics(self, sequences: List[int], sequence_info: List[Dict], 
                                  gc_counts: List[int], n_counts: List[int]) -> Dict:
        """计算组装指标 | Compute assembly metrics"""
        sequences_sorted = sorted(sequences, reverse=True)
        total_length = sum(sequences)
        num_sequences = len(sequences)
        
        # 基本统计
        stats = {
            'total_length': total_length,
            'num_contigs': num_sequences,
            'longest_contig': max(sequences),
            'shortest_contig': min(sequences),
            'mean_contig_length': total_length / num_sequences,
            'median_contig_length': sequences_sorted[num_sequences // 2]
        }
        
        # 计算Nx值
        nx_values = {}
        for x in [10, 20, 30, 40, 50, 60, 70, 80, 90]:
            nx_values[f'N{x}'] = self._calculate_nx(sequences_sorted, total_length, x)
        stats.update(nx_values)
        
        # 计算Lx值（达到Nx所需的contig数量）
        lx_values = {}
        for x in [50, 90]:
            lx_values[f'L{x}'] = self._calculate_lx(sequences_sorted, total_length, x)
        stats.update(lx_values)
        
        # 计算auN (area under Nx curve)
        aun = sum(self._calculate_nx(sequences_sorted, total_length, i) for i in range(1, 101))
        stats['auN'] = aun
        
        # GC含量统计
        total_gc = sum(gc_counts)
        total_bases = total_length
        stats['gc_content'] = (total_gc / total_bases * 100) if total_bases > 0 else 0
        
        # N含量统计
        total_n = sum(n_counts)
        stats['n_content'] = (total_n / total_bases * 100) if total_bases > 0 else 0
        
        # 长度分布统计
        length_thresholds = [1000, 5000, 10000, 50000, 100000, 500000, 1000000, 5000000, 10000000]
        for threshold in length_thresholds:
            count = len([s for s in sequences if s >= threshold])
            stats[f'contigs_>={threshold//1000}kb' if threshold < 1000000 else f'contigs_>={threshold//1000000}mb'] = count
        
        # 覆盖度分析（如果有长度信息）
        stats['gaps'] = len([info for info in sequence_info if info['n_content'] > 0])
        stats['gap_percentage'] = stats['n_content']
        
        return stats
    
    def _calculate_nx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        """计算Nx值 | Calculate Nx value"""
        target = total_length * (x / 100.0)
        cumulative = 0
        
        for length in sorted_lengths:
            cumulative += length
            if cumulative >= target:
                return length
        
        return 0
    
    def _calculate_lx(self, sorted_lengths: List[int], total_length: int, x: int) -> int:
        """计算Lx值 | Calculate Lx value"""
        target = total_length * (x / 100.0)
        cumulative = 0
        
        for i, length in enumerate(sorted_lengths):
            cumulative += length
            if cumulative >= target:
                return i + 1
        
        return len(sorted_lengths)
    
    def _generate_statistics_report(self, all_stats: Dict[str, Dict]):
        """生成统计报告 | Generate statistics report"""
        report_file = self.stats_dir / 'assembly_statistics_report.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("HiFiasm组装统计报告 | HiFiasm Assembly Statistics Report\n")
            f.write("="*80 + "\n")
            f.write(f"生成时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"样本: {self.config.prefix}\n\n")
            
            for assembly_type, stats in all_stats.items():
                f.write(f"\n{assembly_type.upper()} 组装统计 | {assembly_type.upper()} Assembly Statistics\n")
                f.write("-"*60 + "\n")
                
                # 基本信息
                f.write("基本信息 | Basic Information:\n")
                f.write(f"  总长度: {stats['total_length']:,} bp ({stats['total_length']/1e9:.2f} Gb)\n")
                f.write(f"  Contig数量: {stats['num_contigs']:,}\n")
                f.write(f"  最长contig: {stats['longest_contig']:,} bp ({stats['longest_contig']/1e6:.2f} Mb)\n")
                f.write(f"  最短contig: {stats['shortest_contig']:,} bp\n")
                f.write(f"  平均长度: {stats['mean_contig_length']:,.0f} bp\n")
                f.write(f"  中位数长度: {stats['median_contig_length']:,} bp\n\n")
                
                # Nx统计
                f.write("连续性统计 | Contiguity Statistics:\n")
                for nx in ['N10', 'N20', 'N30', 'N40', 'N50', 'N60', 'N70', 'N80', 'N90']:
                    f.write(f"  {nx}: {stats[nx]:,} bp ({stats[nx]/1e6:.2f} Mb)\n")
                f.write(f"  L50: {stats['L50']:,} contigs\n")
                f.write(f"  L90: {stats['L90']:,} contigs\n")
                f.write(f"  auN: {stats['auN']:,.0f} bp\n\n")
                
                # 质量统计
                f.write("质量统计 | Quality Statistics:\n")
                f.write(f"  GC含量: {stats['gc_content']:.2f}%\n")
                f.write(f"  N含量: {stats['n_content']:.2f}%\n")
                f.write(f"  包含gaps的contig数: {stats['gaps']:,}\n\n")
                
                # 长度分布
                f.write("长度分布 | Length Distribution:\n")
                length_categories = [
                    ('>=1kb', 'contigs_>=1kb'),
                    ('>=5kb', 'contigs_>=5kb'), 
                    ('>=10kb', 'contigs_>=10kb'),
                    ('>=50kb', 'contigs_>=50kb'),
                    ('>=100kb', 'contigs_>=100kb'),
                    ('>=500kb', 'contigs_>=500kb'),
                    ('>=1mb', 'contigs_>=1mb'),
                    ('>=5mb', 'contigs_>=5mb'),
                    ('>=10mb', 'contigs_>=10mb')
                ]
                
                for label, key in length_categories:
                    if key in stats:
                        f.write(f"  Contigs {label}: {stats[key]:,}\n")
                
                f.write("\n")
    
    def _generate_comparison_table(self, all_stats: Dict[str, Dict]):
        """生成比较表格 | Generate comparison table"""
        try:
            # 创建DataFrame
            comparison_data = []
            
            for assembly_type, stats in all_stats.items():
                comparison_data.append({
                    'Assembly': assembly_type,
                    'Total_Length_Gb': round(stats['total_length'] / 1e9, 2),
                    'Num_Contigs': stats['num_contigs'],
                    'N50_Mb': round(stats['N50'] / 1e6, 2),
                    'N90_Mb': round(stats['N90'] / 1e6, 2),
                    'L50': stats['L50'],
                    'Longest_Contig_Mb': round(stats['longest_contig'] / 1e6, 2),
                    'GC_Content_%': round(stats['gc_content'], 2),
                    'N_Content_%': round(stats['n_content'], 2),
                    'Contigs_>=1Mb': stats.get('contigs_>=1mb', 0),
                    'auN_Mb': round(stats['auN'] / 1e6, 2)
                })
            
            df = pd.DataFrame(comparison_data)
            
            # 保存CSV
            csv_file = self.stats_dir / 'assembly_statistics_comparison.csv'
            df.to_csv(csv_file, index=False)
            
            # 保存Excel（如果可用）
            try:
                excel_file = self.stats_dir / 'assembly_statistics_comparison.xlsx'
                df.to_excel(excel_file, index=False)
            except ImportError:
                self.logger.warning("pandas Excel支持不可用，跳过Excel文件生成 | pandas Excel support not available, skipping Excel file")
            
            self.logger.info(f"比较表格已生成 | Comparison table generated: {csv_file}")
            
        except Exception as e:
            self.logger.error(f"生成比较表格失败 | Failed to generate comparison table: {e}")
    
    def _generate_statistics_plots(self, all_stats: Dict[str, Dict]):
        """生成统计图表 | Generate statistics plots"""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            import numpy as np
            
            plt.style.use('seaborn-v0_8')
            
            # 创建多子图
            fig, axes = plt.subplots(2, 3, figsize=(18, 12))
            fig.suptitle('HiFiasm组装统计分析 | HiFiasm Assembly Statistics Analysis', 
                        fontsize=16, fontweight='bold')
            
            assemblies = list(all_stats.keys())
            colors = plt.cm.Set2(np.linspace(0, 1, len(assemblies)))
            
            # 1. 基因组大小比较
            ax1 = axes[0, 0]
            sizes = [all_stats[asm]['total_length'] / 1e9 for asm in assemblies]
            bars1 = ax1.bar(assemblies, sizes, color=colors, alpha=0.7)
            ax1.set_title('基因组大小比较 | Genome Size Comparison', fontweight='bold')
            ax1.set_ylabel('大小 (Gb) | Size (Gb)')
            ax1.tick_params(axis='x', rotation=45)
            
            for bar, size in zip(bars1, sizes):
                ax1.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.01,
                        f'{size:.2f}', ha='center', va='bottom')
            
            # 2. N50比较
            ax2 = axes[0, 1]
            n50s = [all_stats[asm]['N50'] / 1e6 for asm in assemblies]
            bars2 = ax2.bar(assemblies, n50s, color=colors, alpha=0.7)
            ax2.set_title('N50比较 | N50 Comparison', fontweight='bold')
            ax2.set_ylabel('N50 (Mb)')
            ax2.tick_params(axis='x', rotation=45)
            
            for bar, n50 in zip(bars2, n50s):
                ax2.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.1,
                        f'{n50:.1f}', ha='center', va='bottom')
            
            # 3. Contig数量比较
            ax3 = axes[0, 2]
            contig_nums = [all_stats[asm]['num_contigs'] for asm in assemblies]
            bars3 = ax3.bar(assemblies, contig_nums, color=colors, alpha=0.7)
            ax3.set_title('Contig数量比较 | Contig Number Comparison', fontweight='bold')
            ax3.set_ylabel('Contig数量 | Number of Contigs')
            ax3.tick_params(axis='x', rotation=45)
            
            for bar, num in zip(bars3, contig_nums):
                ax3.text(bar.get_x() + bar.get_width()/2., bar.get_height() + max(contig_nums)*0.01,
                        f'{num:,}', ha='center', va='bottom')
            
            # 4. GC含量比较
            ax4 = axes[1, 0]
            gc_contents = [all_stats[asm]['gc_content'] for asm in assemblies]
            bars4 = ax4.bar(assemblies, gc_contents, color=colors, alpha=0.7)
            ax4.set_title('GC含量比较 | GC Content Comparison', fontweight='bold')
            ax4.set_ylabel('GC含量 (%) | GC Content (%)')
            ax4.tick_params(axis='x', rotation=45)
            
            for bar, gc in zip(bars4, gc_contents):
                ax4.text(bar.get_x() + bar.get_width()/2., bar.get_height() + 0.1,
                        f'{gc:.1f}%', ha='center', va='bottom')
            
            # 5. 长contig数量比较
            ax5 = axes[1, 1]
            long_contigs = [all_stats[asm].get('contigs_>=1mb', 0) for asm in assemblies]
            bars5 = ax5.bar(assemblies, long_contigs, color=colors, alpha=0.7)
            ax5.set_title('长Contig数量 (≥1Mb) | Long Contigs (≥1Mb)', fontweight='bold')
            ax5.set_ylabel('数量 | Number')
            ax5.tick_params(axis='x', rotation=45)
            
            for bar, num in zip(bars5, long_contigs):
                ax5.text(bar.get_x() + bar.get_width()/2., bar.get_height() + max(long_contigs)*0.01,
                        f'{num}', ha='center', va='bottom')
            
            # 6. N含量比较
            ax6 = axes[1, 2]
            n_contents = [all_stats[asm]['n_content'] for asm in assemblies]
            bars6 = ax6.bar(assemblies, n_contents, color=colors, alpha=0.7)
            ax6.set_title('N含量比较 | N Content Comparison', fontweight='bold')
            ax6.set_ylabel('N含量 (%) | N Content (%)')
            ax6.tick_params(axis='x', rotation=45)
            
            for bar, n_content in zip(bars6, n_contents):
                ax6.text(bar.get_x() + bar.get_width()/2., bar.get_height() + max(n_contents)*0.01,
                        f'{n_content:.2f}%', ha='center', va='bottom')
            
            plt.tight_layout()
            plot_file = self.stats_dir / 'assembly_statistics_plots.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"统计图表已生成 | Statistics plots generated: {plot_file}")
            
        except ImportError:
            self.logger.warning("matplotlib/seaborn未安装，跳过统计图表生成 | matplotlib/seaborn not installed, skipping statistics plots")
        except Exception as e:
            self.logger.warning(f"生成统计图表失败 | Failed to generate statistics plots: {e}")

class FormatConverter:
    """格式转换器 | Format Converter"""
    
    def __init__(self, config, logger: logging.Logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.conversion_dir = Path(config.output_dir) / 'format_conversion'
        self.conversion_dir.mkdir(parents=True, exist_ok=True)
    
    def convert_all_formats(self) -> bool:
        """转换所有格式 | Convert all formats"""
        try:
            self.logger.info("开始格式转换 | Starting format conversion")
            
            success = True
            
            # 压缩输出文件（如果需要）
            if self.config.compress_output:
                if not self._compress_output_files():
                    success = False
            
            # 生成索引文件
            if not self._generate_index_files():
                success = False
            
            # 生成BED文件（如果需要）
            if not self._generate_bed_files():
                success = False
            
            if success:
                self.logger.success("格式转换完成 | Format conversion completed")
            else:
                self.logger.warning("部分格式转换失败 | Some format conversions failed")
            
            return success
            
        except Exception as e:
            self.logger.error(f"格式转换失败 | Format conversion failed: {e}")
            return False
    
    def _compress_output_files(self) -> bool:
        """压缩输出文件 | Compress output files"""
        try:
            self.logger.info("压缩输出文件 | Compressing output files")
            
            # 获取FASTA文件
            from .assembly import GFAConverter
            converter = GFAConverter(self.config, self.logger, self.cmd_runner)
            fasta_files = converter.get_converted_files()
            
            compressed_files = {}
            
            for assembly_type, fasta_path in fasta_files.items():
                compressed_path = self._compress_single_file(fasta_path)
                if compressed_path:
                    compressed_files[assembly_type] = compressed_path
            
            self.logger.success(f"压缩了 {len(compressed_files)} 个文件 | Compressed {len(compressed_files)} files")
            return True
            
        except Exception as e:
            self.logger.error(f"文件压缩失败 | File compression failed: {e}")
            return False
    
    def _compress_single_file(self, file_path: Path) -> Optional[Path]:
        """压缩单个文件 | Compress single file"""
        try:
            compressed_path = file_path.with_suffix(file_path.suffix + '.gz')
            
            cmd = ['gzip', '-c', str(file_path)]
            
            result = self.cmd_runner.run(
                cmd=cmd,
                description=f"压缩文件 | Compress file: {file_path.name}"
            )
            
            with open(compressed_path, 'wb') as f:
                f.write(result.stdout.encode())
            
            self.logger.info(f"文件已压缩 | File compressed: {compressed_path.name}")
            return compressed_path
            
        except Exception as e:
            self.logger.error(f"压缩文件失败 | Failed to compress file: {e}")
            return None
    
    def _generate_index_files(self) -> bool:
        """生成索引文件 | Generate index files"""
        try:
            self.logger.info("生成索引文件 | Generating index files")
            
            # 获取FASTA文件
            from .assembly import GFAConverter
            converter = GFAConverter(self.config, self.logger, self.cmd_runner)
            fasta_files = converter.get_converted_files()
            
            if not fasta_files:
                return True
            
            indexed_files = 0
            
            for assembly_type, fasta_path in fasta_files.items():
                if self._generate_single_index(fasta_path):
                    indexed_files += 1
            
            self.logger.success(f"生成了 {indexed_files} 个索引文件 | Generated {indexed_files} index files")
            return True
            
        except Exception as e:
            self.logger.error(f"生成索引文件失败 | Failed to generate index files: {e}")
            return False
    
    def _generate_single_index(self, fasta_path: Path) -> bool:
        """生成单个索引文件 | Generate single index file"""
        try:
            # 使用samtools生成fai索引
            cmd = [self.config.samtools_path, 'faidx', str(fasta_path)]
            
            result = self.cmd_runner.run(
                cmd=cmd,
                description=f"生成索引 | Generate index: {fasta_path.name}",
                check=False
            )
            
            if result.returncode == 0:
                index_file = fasta_path.with_suffix(fasta_path.suffix + '.fai')
                if index_file.exists():
                    self.logger.info(f"索引文件已生成 | Index file generated: {index_file.name}")
                    return True
            
            return False
            
        except Exception as e:
            self.logger.warning(f"生成索引失败 | Failed to generate index: {e}")
            return False
    
    def _generate_bed_files(self) -> bool:
        """生成BED文件 | Generate BED files"""
        try:
            self.logger.info("生成BED格式文件 | Generating BED format files")
            
            # 获取FASTA文件
            from .assembly import GFAConverter
            converter = GFAConverter(self.config, self.logger, self.cmd_runner)
            fasta_files = converter.get_converted_files()
            
            if not fasta_files:
                return True
            
            bed_files = 0
            
            for assembly_type, fasta_path in fasta_files.items():
                if self._generate_single_bed(fasta_path, assembly_type):
                    bed_files += 1
            
            self.logger.success(f"生成了 {bed_files} 个BED文件 | Generated {bed_files} BED files")
            return True
            
        except Exception as e:
            self.logger.error(f"生成BED文件失败 | Failed to generate BED files: {e}")
            return False
    
    def _generate_single_bed(self, fasta_path: Path, assembly_type: str) -> bool:
        """生成单个BED文件 | Generate single BED file"""
        try:
            bed_path = self.conversion_dir / f"{assembly_type}_contigs.bed"
            
            with open(fasta_path, 'r') as fasta_file, open(bed_path, 'w') as bed_file:
                current_seq_id = None
                current_seq_len = 0
                
                for line in fasta_file:
                    line = line.strip()
                    if line.startswith('>'):
                        # 写入上一个序列的BED记录
                        if current_seq_id and current_seq_len > 0:
                            bed_file.write(f"{current_seq_id}\t0\t{current_seq_len}\t{current_seq_id}\n")
                        
                        # 开始新序列
                        current_seq_id = line[1:].split()[0]  # 取ID的第一部分
                        current_seq_len = 0
                    else:
                        current_seq_len += len(line)
                
                # 写入最后一个序列
                if current_seq_id and current_seq_len > 0:
                    bed_file.write(f"{current_seq_id}\t0\t{current_seq_len}\t{current_seq_id}\n")
            
            self.logger.info(f"BED文件已生成 | BED file generated: {bed_path.name}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成BED文件失败 | Failed to generate BED file: {e}")
            return False