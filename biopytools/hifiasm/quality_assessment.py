"""
HiFiasm质量评估模块 | HiFiasm Quality Assessment Module
"""

import os
import logging
import pandas as pd
import subprocess
from pathlib import Path
from typing import Dict, List, Optional, Any, Tuple
import json
import re

class BUSCOAssessor:
    """BUSCO质量评估器 | BUSCO Quality Assessor"""
    
    def __init__(self, config, logger: logging.Logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.busco_dir = Path(config.output_dir) / 'quality_assessment' / 'busco'
        self.busco_dir.mkdir(parents=True, exist_ok=True)
    
    def run_busco_assessment(self) -> bool:
        """运行BUSCO质量评估 | Run BUSCO quality assessment"""
        try:
            self.logger.info("开始BUSCO质量评估 | Starting BUSCO quality assessment")
            
            # 获取FASTA文件
            from .assembly import GFAConverter
            converter = GFAConverter(self.config, self.logger, self.cmd_runner)
            fasta_files = converter.get_converted_files()
            
            if not fasta_files:
                self.logger.error("未找到FASTA文件进行BUSCO评估 | No FASTA files found for BUSCO assessment")
                return False
            
            # 准备BUSCO数据集
            lineage = self._prepare_busco_lineage()
            if not lineage:
                return False
            
            # 运行BUSCO评估
            busco_results = {}
            for assembly_type, fasta_path in fasta_files.items():
                result = self._run_single_busco(fasta_path, assembly_type, lineage)
                if result:
                    busco_results[assembly_type] = result
            
            if not busco_results:
                self.logger.error("所有BUSCO评估都失败了 | All BUSCO assessments failed")
                return False
            
            # 生成比较报告
            self._generate_busco_comparison(busco_results)
            
            self.logger.success(f"BUSCO评估完成，评估了 {len(busco_results)} 个组装 | BUSCO assessment completed for {len(busco_results)} assemblies")
            return True
            
        except Exception as e:
            self.logger.error(f"BUSCO评估失败 | BUSCO assessment failed: {e}")
            return False
    
    def _prepare_busco_lineage(self) -> Optional[str]:
        """准备BUSCO谱系数据集 | Prepare BUSCO lineage dataset"""
        lineage = self.config.get_busco_lineage()
        
        self.logger.info(f"使用BUSCO谱系 | Using BUSCO lineage: {lineage}")
        
        # 检查数据集是否可用
        try:
            cmd = [self.config.busco_path, '--list-datasets']
            result = self.cmd_runner.run(cmd, "检查BUSCO数据集 | Check BUSCO datasets")
            
            available_datasets = result.stdout
            if lineage in available_datasets:
                self.logger.info(f"BUSCO数据集可用 | BUSCO dataset available: {lineage}")
                return lineage
            else:
                self.logger.warning(f"数据集不可用，尝试下载 | Dataset not available, trying to download: {lineage}")
                return self._download_busco_dataset(lineage)
                
        except Exception as e:
            self.logger.error(f"检查BUSCO数据集失败 | Failed to check BUSCO datasets: {e}")
            return None
    
    def _download_busco_dataset(self, lineage: str) -> Optional[str]:
        """下载BUSCO数据集 | Download BUSCO dataset"""
        try:
            self.logger.info(f"下载BUSCO数据集 | Downloading BUSCO dataset: {lineage}")
            
            cmd = [self.config.busco_path, '--download', lineage]
            
            if self.config.busco_download_path:
                cmd.extend(['--download_path', self.config.busco_download_path])
            
            result = self.cmd_runner.run(
                cmd=cmd,
                description=f"下载BUSCO数据集 | Download BUSCO dataset: {lineage}",
                timeout=1800  # 30分钟超时
            )
            
            if result.returncode == 0:
                self.logger.success(f"BUSCO数据集下载成功 | BUSCO dataset downloaded successfully: {lineage}")
                return lineage
            else:
                self.logger.error(f"BUSCO数据集下载失败 | BUSCO dataset download failed: {lineage}")
                return None
                
        except Exception as e:
            self.logger.error(f"下载BUSCO数据集异常 | BUSCO dataset download exception: {e}")
            return None
    
    def _run_single_busco(self, fasta_path: Path, assembly_type: str, lineage: str) -> Optional[Dict]:
        """运行单个BUSCO评估 | Run single BUSCO assessment"""
        try:
            self.logger.info(f"评估 {assembly_type}: {fasta_path.name}")
            
            output_name = f"busco_{assembly_type}"
            output_path = self.busco_dir / output_name
            
            # 构建BUSCO命令
            cmd = [
                self.config.busco_path,
                '-i', str(fasta_path),
                '-l', lineage,
                '-m', self.config.busco_mode,
                '-c', str(self.config.threads),
                '-o', output_name,
                '--offline',
                '--quiet'
            ]
            
            # 设置工作目录
            original_working_dir = self.cmd_runner.working_dir
            self.cmd_runner.working_dir = self.busco_dir
            
            try:
                result = self.cmd_runner.run(
                    cmd=cmd,
                    description=f"BUSCO评估 | BUSCO assessment: {assembly_type}",
                    timeout=3600  # 1小时超时
                )
                
                if result.returncode == 0:
                    # 解析BUSCO结果
                    busco_result = self._parse_busco_result(output_path, assembly_type)
                    if busco_result:
                        self.logger.info(f"✓ {assembly_type} BUSCO评估完成 | BUSCO assessment completed")
                        return busco_result
                    
            finally:
                self.cmd_runner.working_dir = original_working_dir
                
        except Exception as e:
            self.logger.error(f"BUSCO评估异常 | BUSCO assessment exception: {e}")
        
        return None
    
    def _parse_busco_result(self, output_path: Path, assembly_type: str) -> Optional[Dict]:
        """解析BUSCO结果 | Parse BUSCO result"""
        try:
            # 查找summary文件
            summary_files = list(output_path.glob("short_summary.*.txt"))
            if not summary_files:
                self.logger.error(f"未找到BUSCO摘要文件 | BUSCO summary file not found: {output_path}")
                return None
            
            summary_file = summary_files[0]
            
            with open(summary_file, 'r') as f:
                content = f.read()
            
            # 解析关键统计信息
            pattern = r'C:(\d+\.?\d*)%\[S:(\d+\.?\d*)%,D:(\d+\.?\d*)%\],F:(\d+\.?\d*)%,M:(\d+\.?\d*)%,n:(\d+)'
            match = re.search(pattern, content)
            
            if not match:
                self.logger.error(f"无法解析BUSCO结果 | Cannot parse BUSCO result: {summary_file}")
                return None
            
            result = {
                'assembly_type': assembly_type,
                'complete': float(match.group(1)),
                'single': float(match.group(2)),
                'duplicated': float(match.group(3)),
                'fragmented': float(match.group(4)),
                'missing': float(match.group(5)),
                'total_genes': int(match.group(6)),
                'summary_file': str(summary_file),
                'lineage': self._extract_lineage_from_summary(content)
            }
            
            # 计算质量等级
            result['quality_grade'] = self._calculate_quality_grade(result['complete'])
            result['polyploidy_signal'] = self._analyze_polyploidy_signal(result['duplicated'])
            
            return result
            
        except Exception as e:
            self.logger.error(f"解析BUSCO结果失败 | Failed to parse BUSCO result: {e}")
            return None
    
    def _extract_lineage_from_summary(self, content: str) -> str:
        """从摘要中提取谱系信息 | Extract lineage info from summary"""
        match = re.search(r'The lineage dataset is: (\w+)', content)
        return match.group(1) if match else 'unknown'
    
    def _calculate_quality_grade(self, complete_percent: float) -> str:
        """计算质量等级 | Calculate quality grade"""
        if complete_percent >= 95:
            return "优秀 ⭐⭐⭐⭐⭐"
        elif complete_percent >= 90:
            return "良好 ⭐⭐⭐⭐"
        elif complete_percent >= 85:
            return "中等 ⭐⭐⭐"
        elif complete_percent >= 80:
            return "一般 ⭐⭐"
        else:
            return "较差 ⭐"
    
    def _analyze_polyploidy_signal(self, duplicated_percent: float) -> str:
        """分析古多倍体信号 | Analyze ancient polyploidy signal"""
        if duplicated_percent >= 15:
            return "强烈古多倍体信号"
        elif duplicated_percent >= 10:
            return "明显古多倍体信号"
        elif duplicated_percent >= 5:
            return "轻微古多倍体信号"
        else:
            return "无明显古多倍体信号"
    
    def _generate_busco_comparison(self, results: Dict[str, Dict]):
        """生成BUSCO比较报告 | Generate BUSCO comparison report"""
        try:
            self.logger.info("生成BUSCO比较报告 | Generating BUSCO comparison report")
            
            # 创建比较表格
            comparison_data = []
            for assembly_type, result in results.items():
                comparison_data.append({
                    'Assembly': result['assembly_type'],
                    'Complete(%)': result['complete'],
                    'Single(%)': result['single'],
                    'Duplicated(%)': result['duplicated'],
                    'Fragmented(%)': result['fragmented'],
                    'Missing(%)': result['missing'],
                    'Total_Genes': result['total_genes'],
                    'Quality_Grade': result['quality_grade'].split()[0],
                    'Polyploidy_Signal': result['polyploidy_signal']
                })
            
            # 保存CSV报告
            df = pd.DataFrame(comparison_data)
            csv_file = self.busco_dir / 'busco_comparison_report.csv'
            df.to_csv(csv_file, index=False)
            
            # 生成文本报告
            self._generate_busco_text_report(results, df)
            
            # 如果启用绘图，生成图表
            if self.config.generate_plots:
                self._generate_busco_plots(df)
            
            self.logger.success(f"BUSCO比较报告已生成 | BUSCO comparison report generated: {csv_file}")
            
        except Exception as e:
            self.logger.error(f"生成BUSCO报告失败 | Failed to generate BUSCO report: {e}")
    
    def _generate_busco_text_report(self, results: Dict, df: pd.DataFrame):
        """生成BUSCO文本报告 | Generate BUSCO text report"""
        report_file = self.busco_dir / 'busco_assessment_report.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("HiFiasm基因组BUSCO质量评估报告 | HiFiasm Genome BUSCO Quality Assessment Report\n")
            f.write("="*80 + "\n")
            f.write(f"评估时间: {pd.Timestamp.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"使用谱系: {list(results.values())[0]['lineage']}\n")
            f.write(f"总基因数: {list(results.values())[0]['total_genes']}\n\n")
            
            # 详细结果对比
            f.write("详细结果对比 | Detailed Results Comparison\n")
            f.write("-"*80 + "\n")
            f.write(f"{'组装名称':<15} {'完整性':<10} {'单拷贝':<10} {'重复':<10} {'片段':<8} {'缺失':<8} {'质量等级':<15}\n")
            f.write("-"*80 + "\n")
            
            for _, row in df.iterrows():
                f.write(f"{row['Assembly']:<15} {row['Complete(%)']:<10.1f}% {row['Single(%)']:<10.1f}% "
                       f"{row['Duplicated(%)']:<10.1f}% {row['Fragmented(%)']:<8.1f}% "
                       f"{row['Missing(%)']:<8.1f}% {row['Quality_Grade']:<15}\n")
            
            # 推荐和建议
            f.write("\n" + "="*80 + "\n")
            f.write("质量评估和推荐 | Quality Assessment and Recommendations\n")
            f.write("="*80 + "\n")
            
            best_assembly = df.loc[df['Complete(%)'].idxmax()]
            f.write(f"最佳组装: {best_assembly['Assembly']}\n")
            f.write(f"完整性: {best_assembly['Complete(%)']:.1f}%\n")
            f.write(f"质量等级: {best_assembly['Quality_Grade']}\n\n")
            
            # 古多倍体化分析
            f.write("古多倍体化特征分析 | Ancient Polyploidy Analysis:\n")
            for _, row in df.iterrows():
                f.write(f"{row['Assembly']}: {row['Duplicated(%)']:.1f}% 重复基因 - {row['Polyploidy_Signal']}\n")
    
    def _generate_busco_plots(self, df: pd.DataFrame):
        """生成BUSCO可视化图表 | Generate BUSCO visualization plots"""
        try:
            import matplotlib.pyplot as plt
            import seaborn as sns
            
            plt.style.use('seaborn-v0_8')
            fig, axes = plt.subplots(2, 2, figsize=(12, 10))
            fig.suptitle('HiFiasm基因组BUSCO质量评估结果', fontsize=16, fontweight='bold')
            
            # 1. 完整性比较
            ax1 = axes[0, 0]
            bars1 = ax1.bar(df['Assembly'], df['Complete(%)'], color='steelblue', alpha=0.7)
            ax1.set_title('基因组完整性 (Complete %)', fontweight='bold')
            ax1.set_ylabel('完整性 (%)')
            ax1.set_ylim(80, 100)
            
            for bar in bars1:
                height = bar.get_height()
                ax1.text(bar.get_x() + bar.get_width()/2., height + 0.3,
                        f'{height:.1f}%', ha='center', va='bottom')
            
            # 2. BUSCO分类详细对比
            ax2 = axes[0, 1]
            categories = ['Single(%)', 'Duplicated(%)', 'Fragmented(%)', 'Missing(%)']
            colors = ['lightblue', 'orange', 'lightcoral', 'lightgray']
            
            x = range(len(df))
            width = 0.2
            
            for i, category in enumerate(categories):
                ax2.bar([j + i*width for j in x], df[category], width, 
                        label=category.replace('(%)', ''), color=colors[i], alpha=0.7)
            
            ax2.set_title('BUSCO分类详细对比', fontweight='bold')
            ax2.set_ylabel('百分比 (%)')
            ax2.set_xticks([j + width*1.5 for j in x])
            ax2.set_xticklabels(df['Assembly'])
            ax2.legend()
            
            # 3. 重复基因（古多倍体信号）
            ax3 = axes[1, 0]
            bars3 = ax3.bar(df['Assembly'], df['Duplicated(%)'], color='darkorange', alpha=0.7)
            ax3.set_title('重复基因比例 (古多倍体信号)', fontweight='bold')
            ax3.set_ylabel('重复基因 (%)')
            ax3.axhline(y=10, color='red', linestyle='--', alpha=0.5, label='古多倍体阈值')
            ax3.legend()
            
            for bar in bars3:
                height = bar.get_height()
                ax3.text(bar.get_x() + bar.get_width()/2., height + 0.1,
                        f'{height:.1f}%', ha='center', va='bottom')
            
            # 4. 质量缺陷（片段化+缺失）
            ax4 = axes[1, 1]
            defects = df['Fragmented(%)'] + df['Missing(%)']
            bars4 = ax4.bar(df['Assembly'], defects, color='lightcoral', alpha=0.7)
            ax4.set_title('组装缺陷 (片段化 + 缺失)', fontweight='bold')
            ax4.set_ylabel('缺陷比例 (%)')
            ax4.axhline(y=5, color='red', linestyle='--', alpha=0.5, label='期望阈值(<5%)')
            ax4.legend()
            
            for bar in bars4:
                height = bar.get_height()
                ax4.text(bar.get_x() + bar.get_width()/2., height + 0.05,
                        f'{height:.1f}%', ha='center', va='bottom')
            
            plt.tight_layout()
            plot_file = self.busco_dir / 'busco_comparison_plot.png'
            plt.savefig(plot_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"BUSCO可视化图表已生成 | BUSCO plots generated: {plot_file}")
            
        except ImportError:
            self.logger.warning("matplotlib未安装，跳过绘图 | matplotlib not installed, skipping plots")
        except Exception as e:
            self.logger.warning(f"生成BUSCO图表失败 | Failed to generate BUSCO plots: {e}")

class QUASTAssessor:
    """QUAST质量评估器 | QUAST Quality Assessor"""
    
    def __init__(self, config, logger: logging.Logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.quast_dir = Path(config.output_dir) / 'quality_assessment' / 'quast'
        self.quast_dir.mkdir(parents=True, exist_ok=True)
    
    def run_quast_assessment(self) -> bool:
        """运行QUAST质量评估 | Run QUAST quality assessment"""
        try:
            self.logger.info("开始QUAST质量评估 | Starting QUAST quality assessment")
            
            # 获取FASTA文件
            from .assembly import GFAConverter
            converter = GFAConverter(self.config, self.logger, self.cmd_runner)
            fasta_files = converter.get_converted_files()
            
            if not fasta_files:
                self.logger.error("未找到FASTA文件进行QUAST评估 | No FASTA files found for QUAST assessment")
                return False
            
            # 运行QUAST
            if not self._run_quast(fasta_files):
                return False
            
            # 解析QUAST结果
            if not self._parse_quast_results():
                return False
            
            self.logger.success("QUAST评估完成 | QUAST assessment completed")
            return True
            
        except Exception as e:
            self.logger.error(f"QUAST评估失败 | QUAST assessment failed: {e}")
            return False
    
    def _run_quast(self, fasta_files: Dict[str, Path]) -> bool:
        """运行QUAST命令 | Run QUAST command"""
        try:
            # 构建QUAST命令
            cmd = [
                self.config.quast_path,
                '--threads', str(self.config.threads),
                '--output-dir', str(self.quast_dir),
                '--labels', ','.join(fasta_files.keys())
            ]
            
            # 添加参考基因组（如果有）
            if self.config.reference_genome:
                cmd.extend(['-r', self.config.reference_genome])
            
            # 添加输入文件
            cmd.extend([str(path) for path in fasta_files.values()])
            
            result = self.cmd_runner.run(
                cmd=cmd,
                description="QUAST组装质量评估 | QUAST assembly quality assessment",
                timeout=1800  # 30分钟超时
            )
            
            return result.returncode == 0
            
        except Exception as e:
            self.logger.error(f"运行QUAST失败 | Failed to run QUAST: {e}")
            return False
    
    def _parse_quast_results(self) -> bool:
        """解析QUAST结果 | Parse QUAST results"""
        try:
            # 查找QUAST报告文件
            report_file = self.quast_dir / 'report.txt'
            if not report_file.exists():
                self.logger.error(f"QUAST报告文件不存在 | QUAST report file not found: {report_file}")
                return False
            
            # 解析报告
            quast_data = self._parse_quast_report(report_file)
            
            # 生成总结
            self._generate_quast_summary(quast_data)
            
            return True
            
        except Exception as e:
            self.logger.error(f"解析QUAST结果失败 | Failed to parse QUAST results: {e}")
            return False
    
    def _parse_quast_report(self, report_file: Path) -> Dict:
        """解析QUAST报告文件 | Parse QUAST report file"""
        quast_data = {}
        
        with open(report_file, 'r') as f:
            lines = f.readlines()
        
        # 解析关键统计信息
        for line in lines:
            line = line.strip()
            if line.startswith('# contigs'):
                quast_data['num_contigs'] = self._extract_values_from_line(line)
            elif line.startswith('Total length'):
                quast_data['total_length'] = self._extract_values_from_line(line)
            elif line.startswith('N50'):
                quast_data['n50'] = self._extract_values_from_line(line)
            elif line.startswith('N90'):
                quast_data['n90'] = self._extract_values_from_line(line)
            elif line.startswith('L50'):
                quast_data['l50'] = self._extract_values_from_line(line)
            elif line.startswith('Largest contig'):
                quast_data['largest_contig'] = self._extract_values_from_line(line)
            elif line.startswith('GC (%)'):
                quast_data['gc_content'] = self._extract_values_from_line(line)
        
        return quast_data
    
    def _extract_values_from_line(self, line: str) -> List:
        """从行中提取数值 | Extract values from line"""
        parts = line.split('\t')
        return parts[1:] if len(parts) > 1 else []
    
    def _generate_quast_summary(self, quast_data: Dict):
        """生成QUAST总结 | Generate QUAST summary"""
        summary_file = self.quast_dir / 'quast_summary.txt'
        
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("QUAST质量评估总结 | QUAST Quality Assessment Summary\n")
            f.write("="*60 + "\n\n")
            
            for metric, values in quast_data.items():
                f.write(f"{metric}:\n")
                for i, value in enumerate(values):
                    f.write(f"  Assembly {i+1}: {value}\n")
                f.write("\n")

class HaplotypeAnalyzer:
    """单倍型分析器 | Haplotype Analyzer"""
    
    def __init__(self, config, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.haplotype_dir = Path(config.output_dir) / 'quality_assessment' / 'haplotype_analysis'
        self.haplotype_dir.mkdir(parents=True, exist_ok=True)
    
    def analyze_haplotypes(self) -> bool:
        """分析单倍型差异 | Analyze haplotype differences"""
        try:
            self.logger.info("开始单倍型差异分析 | Starting haplotype difference analysis")
            
            # 获取单倍型文件
            haplotype_files = self._get_haplotype_files()
            if len(haplotype_files) < 2:
                self.logger.warning("未找到足够的单倍型文件进行比较 | Not enough haplotype files for comparison")
                return True
            
            # 计算基本统计信息
            hap_stats = self._calculate_haplotype_statistics(haplotype_files)
            
            # 分析差异
            differences = self._analyze_haplotype_differences(hap_stats)
            
            # 生成分析报告
            self._generate_haplotype_report(hap_stats, differences)
            
            # 生成选择策略建议
            self._generate_selection_strategy(hap_stats, differences)
            
            self.logger.success("单倍型差异分析完成 | Haplotype difference analysis completed")
            return True
            
        except Exception as e:
            self.logger.error(f"单倍型分析失败 | Haplotype analysis failed: {e}")
            return False
    
    def _get_haplotype_files(self) -> Dict[str, Path]:
        """获取单倍型文件 | Get haplotype files"""
        from .assembly import GFAConverter
        converter = GFAConverter(self.config, self.logger, None)
        return converter.get_converted_files()
    
    def _calculate_haplotype_statistics(self, files: Dict[str, Path]) -> Dict:
        """计算单倍型统计信息 | Calculate haplotype statistics"""
        stats = {}
        
        for hap_type, file_path in files.items():
            try:
                hap_stat = self._analyze_single_assembly(file_path)
                stats[hap_type] = hap_stat
            except Exception as e:
                self.logger.warning(f"分析 {hap_type} 失败 | Failed to analyze {hap_type}: {e}")
        
        return stats
    
    def _analyze_single_assembly(self, fasta_path: Path) -> Dict:
        """分析单个组装文件 | Analyze single assembly file"""
        sequences = []
        gc_content = []
        n_content = []
        current_seq = ""
        
        with open(fasta_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    if current_seq:
                        sequences.append(len(current_seq))
                        gc_count = current_seq.upper().count('G') + current_seq.upper().count('C')
                        n_count = current_seq.upper().count('N')
                        total_bases = len(current_seq)
                        if total_bases > 0:
                            gc_content.append(gc_count / total_bases * 100)
                            n_content.append(n_count / total_bases * 100)
                    current_seq = ""
                else:
                    current_seq += line.upper()
        
        if current_seq:
            sequences.append(len(current_seq))
            gc_count = current_seq.count('G') + current_seq.count('C')
            n_count = current_seq.count('N')
            total_bases = len(current_seq)
            if total_bases > 0:
                gc_content.append(gc_count / total_bases * 100)
                n_content.append(n_count / total_bases * 100)
        
        if not sequences:
            return {}
        
        sequences.sort(reverse=True)
        total_length = sum(sequences)
        
        # 计算N50
        def calculate_nx(lengths, x):
            target = total_length * (x / 100.0)
            cumulative = 0
            for length in lengths:
                cumulative += length
                if cumulative >= target:
                    return length
            return 0
        
        return {
            'total_length': total_length,
            'num_contigs': len(sequences),
            'longest_contig': max(sequences),
            'n50': calculate_nx(sequences, 50),
            'n90': calculate_nx(sequences, 90),
            'mean_gc': sum(gc_content) / len(gc_content) if gc_content else 0,
            'mean_n_content': sum(n_content) / len(n_content) if n_content else 0,
            'contigs_1mb': len([s for s in sequences if s >= 1000000]),
            'contigs_5mb': len([s for s in sequences if s >= 5000000]),
            'contigs_10mb': len([s for s in sequences if s >= 10000000])
        }
    
    def _analyze_haplotype_differences(self, stats: Dict) -> Dict:
        """分析单倍型差异 | Analyze haplotype differences"""
        if 'haplotype1' not in stats or 'haplotype2' not in stats:
            return {}
        
        hap1 = stats['haplotype1']
        hap2 = stats['haplotype2']
        
        differences = {
            'size_diff_percent': abs(hap1['total_length'] - hap2['total_length']) / max(hap1['total_length'], hap2['total_length']) * 100,
            'n50_diff_percent': abs(hap1['n50'] - hap2['n50']) / max(hap1['n50'], hap2['n50']) * 100,
            'contig_diff_percent': abs(hap1['num_contigs'] - hap2['num_contigs']) / max(hap1['num_contigs'], hap2['num_contigs']) * 100,
            'gc_diff': abs(hap1['mean_gc'] - hap2['mean_gc']),
            'n_content_diff': abs(hap1['mean_n_content'] - hap2['mean_n_content'])
        }
        
        return differences
    
    def _generate_haplotype_report(self, stats: Dict, differences: Dict):
        """生成单倍型分析报告 | Generate haplotype analysis report"""
        report_file = self.haplotype_dir / 'haplotype_analysis_report.txt'
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("单倍型差异分析报告 | Haplotype Difference Analysis Report\n")
            f.write("="*80 + "\n\n")
            
            # 基本统计信息
            f.write("基本统计信息 | Basic Statistics\n")
            f.write("-"*40 + "\n")
            for hap_type, stat in stats.items():
                f.write(f"\n{hap_type}:\n")
                f.write(f"  总长度: {stat['total_length']/1e9:.2f} Gb\n")
                f.write(f"  Contig数量: {stat['num_contigs']:,}\n")
                f.write(f"  N50: {stat['n50']/1e6:.2f} Mb\n")
                f.write(f"  最长contig: {stat['longest_contig']/1e6:.2f} Mb\n")
                f.write(f"  >1Mb contigs: {stat['contigs_1mb']:,}\n")
                f.write(f"  平均GC含量: {stat['mean_gc']:.1f}%\n")
                f.write(f"  平均N含量: {stat['mean_n_content']:.2f}%\n")
            
            # 差异分析
            if differences:
                f.write(f"\n\n单倍型差异分析 | Haplotype Difference Analysis\n")
                f.write("-"*40 + "\n")
                f.write(f"基因组大小差异: {differences['size_diff_percent']:.1f}%\n")
                f.write(f"N50差异: {differences['n50_diff_percent']:.1f}%\n")
                f.write(f"Contig数量差异: {differences['contig_diff_percent']:.1f}%\n")
                f.write(f"GC含量差异: {differences['gc_diff']:.2f}%\n")
                f.write(f"N含量差异: {differences['n_content_diff']:.2f}%\n")
    
    def _generate_selection_strategy(self, stats: Dict, differences: Dict):
        """生成选择策略建议 | Generate selection strategy recommendations"""
        strategy_file = self.haplotype_dir / 'selection_strategy.txt'
        
        with open(strategy_file, 'w', encoding='utf-8') as f:
            f.write("组装选择策略建议 | Assembly Selection Strategy Recommendations\n")
            f.write("="*80 + "\n\n")
            
            # 基于质量的推荐
            if stats:
                best_assembly = max(stats.keys(), 
                                  key=lambda x: stats[x]['n50'] * (1 - stats[x]['mean_n_content']/100))
                
                f.write(f"推荐的最佳组装: {best_assembly}\n")
                f.write(f"推荐原因: 综合考虑N50和序列质量\n\n")
            
            # 针对不同用途的建议
            strategies = {
                "基因注释和预测": "使用Primary Assembly或质量最好的单倍型",
                "变异检测和群体遗传学": "保留两个单倍型进行比较分析",
                "比较基因组学": "使用质量最好的单倍型作为参考",
                "功能基因挖掘": "Primary Assembly + 单倍型验证"
            }
            
            f.write("针对不同分析目的的建议 | Recommendations for Different Analysis Purposes\n")
            f.write("-"*60 + "\n")
            
            for purpose, strategy in strategies.items():
                f.write(f"{purpose}: {strategy}\n")