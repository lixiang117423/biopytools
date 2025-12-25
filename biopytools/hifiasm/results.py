"""
HiFiasm结果处理模块 | HiFiasm Results Processing Module
"""

import os
import logging
import json
import shutil
from pathlib import Path
from typing import Dict, List, Optional, Any
import pandas as pd
from datetime import datetime

class ResultsProcessor:
    """结果处理器 | Results Processor"""
    
    def __init__(self, config, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.results_dir = Path(config.output_dir) / 'final_results'
        self.results_dir.mkdir(parents=True, exist_ok=True)
    
    def process_all_results(self) -> bool:
        """处理所有结果 | Process all results"""
        try:
            self.logger.info("开始结果处理 | Starting results processing")
            
            # 整理主要组装文件
            if not self._organize_assembly_files():
                return False
            
            # 整理质量评估结果
            if not self._organize_quality_results():
                return False
            
            # 生成结果清单
            if not self._generate_file_manifest():
                return False
            
            # 创建结果摘要
            if not self._create_results_summary():
                return False
            
            self.logger.success("结果处理完成 | Results processing completed")
            return True
            
        except Exception as e:
            self.logger.error(f"结果处理失败 | Results processing failed: {e}")
            return False
    
    def _organize_assembly_files(self) -> bool:
        """整理组装文件 | Organize assembly files"""
        try:
            self.logger.info("整理组装文件 | Organizing assembly files")
            
            # 创建组装结果目录
            assembly_results_dir = self.results_dir / 'assemblies'
            assembly_results_dir.mkdir(exist_ok=True)
            
            # 获取并复制主要组装文件
            from .assembly import GFAConverter
            converter = GFAConverter(self.config, self.logger, None)
            fasta_files = converter.get_converted_files()
            
            organized_files = {}
            
            for assembly_type, fasta_path in fasta_files.items():
                # 标准化文件名
                standard_name = f"{self.config.prefix}.{assembly_type}.fasta"
                target_path = assembly_results_dir / standard_name
                
                # 复制文件
                shutil.copy2(fasta_path, target_path)
                organized_files[assembly_type] = target_path
                
                self.logger.info(f"已整理组装文件 | Organized assembly file: {standard_name}")
            
            # 保存文件信息
            self._save_file_info(organized_files, 'assembly_files.json')
            
            return True
            
        except Exception as e:
            self.logger.error(f"整理组装文件失败 | Failed to organize assembly files: {e}")
            return False
    
    def _organize_quality_results(self) -> bool:
        """整理质量评估结果 | Organize quality assessment results"""
        try:
            self.logger.info("整理质量评估结果 | Organizing quality assessment results")
            
            # 创建质量评估结果目录
            quality_results_dir = self.results_dir / 'quality_assessment'
            quality_results_dir.mkdir(exist_ok=True)
            
            organized_results = {}
            
            # 整理BUSCO结果
            busco_results = self._organize_busco_results(quality_results_dir)
            if busco_results:
                organized_results['busco'] = busco_results
            
            # 整理QUAST结果
            quast_results = self._organize_quast_results(quality_results_dir)
            if quast_results:
                organized_results['quast'] = quast_results
            
            # 整理统计结果
            stats_results = self._organize_statistics_results(quality_results_dir)
            if stats_results:
                organized_results['statistics'] = stats_results
            
            # 整理单倍型分析结果
            haplotype_results = self._organize_haplotype_results(quality_results_dir)
            if haplotype_results:
                organized_results['haplotype_analysis'] = haplotype_results
            
            # 保存质量评估信息
            self._save_file_info(organized_results, 'quality_assessment_files.json')
            
            return True
            
        except Exception as e:
            self.logger.error(f"整理质量评估结果失败 | Failed to organize quality assessment results: {e}")
            return False
    
    def _organize_busco_results(self, quality_dir: Path) -> Optional[Dict]:
        """整理BUSCO结果 | Organize BUSCO results"""
        try:
            busco_source_dir = Path(self.config.output_dir) / 'quality_assessment' / 'busco'
            if not busco_source_dir.exists():
                return None
            
            busco_target_dir = quality_dir / 'busco'
            busco_target_dir.mkdir(exist_ok=True)
            
            organized_files = {}
            
            # 复制主要BUSCO文件
            for file_pattern in ['*.csv', '*.txt', '*.png']:
                for file_path in busco_source_dir.glob(file_pattern):
                    target_path = busco_target_dir / file_path.name
                    shutil.copy2(file_path, target_path)
                    organized_files[file_path.stem] = str(target_path)
            
            return organized_files
            
        except Exception as e:
            self.logger.warning(f"整理BUSCO结果失败 | Failed to organize BUSCO results: {e}")
            return None
    
    def _organize_quast_results(self, quality_dir: Path) -> Optional[Dict]:
        """整理QUAST结果 | Organize QUAST results"""
        try:
            quast_source_dir = Path(self.config.output_dir) / 'quality_assessment' / 'quast'
            if not quast_source_dir.exists():
                return None
            
            quast_target_dir = quality_dir / 'quast'
            quast_target_dir.mkdir(exist_ok=True)
            
            organized_files = {}
            
            # 复制主要QUAST文件
            important_files = ['report.txt', 'report.html', 'report.pdf']
            for filename in important_files:
                source_file = quast_source_dir / filename
                if source_file.exists():
                    target_file = quast_target_dir / filename
                    shutil.copy2(source_file, target_file)
                    organized_files[filename] = str(target_file)
            
            return organized_files
            
        except Exception as e:
            self.logger.warning(f"整理QUAST结果失败 | Failed to organize QUAST results: {e}")
            return None
    
    def _organize_statistics_results(self, quality_dir: Path) -> Optional[Dict]:
        """整理统计结果 | Organize statistics results"""
        try:
            stats_source_dir = Path(self.config.output_dir) / 'statistics'
            if not stats_source_dir.exists():
                return None
            
            stats_target_dir = quality_dir / 'statistics'
            stats_target_dir.mkdir(exist_ok=True)
            
            organized_files = {}
            
            # 复制统计文件
            for file_pattern in ['*.txt', '*.csv', '*.xlsx', '*.png']:
                for file_path in stats_source_dir.glob(file_pattern):
                    target_path = stats_target_dir / file_path.name
                    shutil.copy2(file_path, target_path)
                    organized_files[file_path.stem] = str(target_path)
            
            return organized_files
            
        except Exception as e:
            self.logger.warning(f"整理统计结果失败 | Failed to organize statistics results: {e}")
            return None
    
    def _organize_haplotype_results(self, quality_dir: Path) -> Optional[Dict]:
        """整理单倍型分析结果 | Organize haplotype analysis results"""
        try:
            hap_source_dir = Path(self.config.output_dir) / 'quality_assessment' / 'haplotype_analysis'
            if not hap_source_dir.exists():
                return None
            
            hap_target_dir = quality_dir / 'haplotype_analysis'
            hap_target_dir.mkdir(exist_ok=True)
            
            organized_files = {}
            
            # 复制单倍型分析文件
            for file_path in hap_source_dir.glob('*.txt'):
                target_path = hap_target_dir / file_path.name
                shutil.copy2(file_path, target_path)
                organized_files[file_path.stem] = str(target_path)
            
            return organized_files
            
        except Exception as e:
            self.logger.warning(f"整理单倍型分析结果失败 | Failed to organize haplotype analysis results: {e}")
            return None
    
    def _save_file_info(self, file_info: Dict, filename: str):
        """保存文件信息 | Save file information"""
        info_file = self.results_dir / filename
        with open(info_file, 'w', encoding='utf-8') as f:
            json.dump(file_info, f, indent=2, ensure_ascii=False)
    
    def _generate_file_manifest(self) -> bool:
        """生成文件清单 | Generate file manifest"""
        try:
            self.logger.info("生成文件清单 | Generating file manifest")
            
            manifest_file = self.results_dir / 'file_manifest.txt'
            
            with open(manifest_file, 'w', encoding='utf-8') as f:
                f.write("HiFiasm分析结果文件清单 | HiFiasm Analysis Results File Manifest\n")
                f.write("="*80 + "\n")
                f.write(f"生成时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"样本: {self.config.prefix}\n")
                f.write(f"输出目录: {self.config.output_dir}\n\n")
                
                # 遍历结果目录生成清单
                self._write_directory_contents(f, self.results_dir, "最终结果 | Final Results")
            
            self.logger.success(f"文件清单已生成 | File manifest generated: {manifest_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成文件清单失败 | Failed to generate file manifest: {e}")
            return False
    
    def _write_directory_contents(self, file_handle, directory: Path, title: str, level: int = 0):
        """写入目录内容 | Write directory contents"""
        indent = "  " * level
        file_handle.write(f"{indent}{title}\n")
        file_handle.write(f"{indent}{'-' * len(title)}\n")
        
        try:
            items = sorted(directory.iterdir())
            
            for item in items:
                if item.is_file():
                    size = item.stat().st_size
                    size_str = self._format_file_size(size)
                    file_handle.write(f"{indent}  📄 {item.name} ({size_str})\n")
                elif item.is_dir():
                    file_handle.write(f"{indent}  📁 {item.name}/\n")
                    self._write_directory_contents(file_handle, item, item.name, level + 2)
            
            file_handle.write("\n")
            
        except PermissionError:
            file_handle.write(f"{indent}  ❌ 权限不足 | Permission denied\n\n")
    
    def _format_file_size(self, size_bytes: int) -> str:
        """格式化文件大小 | Format file size"""
        for unit in ['B', 'KB', 'MB', 'GB']:
            if size_bytes < 1024.0:
                return f"{size_bytes:.1f} {unit}"
            size_bytes /= 1024.0
        return f"{size_bytes:.1f} TB"
    
    def _create_results_summary(self) -> bool:
        """创建结果摘要 | Create results summary"""
        try:
            self.logger.info("创建结果摘要 | Creating results summary")
            
            summary_file = self.results_dir / 'analysis_summary.txt'
            
            with open(summary_file, 'w', encoding='utf-8') as f:
                f.write("HiFiasm基因组组装分析摘要 | HiFiasm Genome Assembly Analysis Summary\n")
                f.write("="*80 + "\n\n")
                
                # 基本信息
                f.write("基本信息 | Basic Information\n")
                f.write("-"*40 + "\n")
                f.write(f"样本名称: {self.config.prefix}\n")
                f.write(f"输入文件: {self.config.input_reads}\n")
                f.write(f"基因组大小估计: {self.config.estimate_genome_size()}\n")
                f.write(f"组装类型: {self.config.assembly_type}\n")
                f.write(f"线程数: {self.config.threads}\n")
                f.write(f"分析时间: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n\n")
                
                # 主要结果文件
                f.write("主要结果文件 | Main Result Files\n")
                f.write("-"*40 + "\n")
                
                assemblies_dir = self.results_dir / 'assemblies'
                if assemblies_dir.exists():
                    for fasta_file in assemblies_dir.glob('*.fasta'):
                        f.write(f"  • {fasta_file.name}\n")
                
                f.write("\n")
                
                # 质量评估摘要
                self._write_quality_summary(f)
                
                # 分析建议
                self._write_analysis_recommendations(f)
            
            self.logger.success(f"结果摘要已创建 | Results summary created: {summary_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"创建结果摘要失败 | Failed to create results summary: {e}")
            return False
    
    def _write_quality_summary(self, file_handle):
        """写入质量评估摘要 | Write quality assessment summary"""
        file_handle.write("质量评估摘要 | Quality Assessment Summary\n")
        file_handle.write("-"*40 + "\n")
        
        # 尝试读取BUSCO结果
        try:
            busco_csv = Path(self.config.output_dir) / 'quality_assessment' / 'busco' / 'busco_comparison_report.csv'
            if busco_csv.exists():
                df = pd.read_csv(busco_csv)
                best_assembly = df.loc[df['Complete(%)'].idxmax()]
                file_handle.write(f"BUSCO评估结果:\n")
                file_handle.write(f"  最佳组装: {best_assembly['Assembly']}\n")
                file_handle.write(f"  完整性: {best_assembly['Complete(%)']:.1f}%\n")
                file_handle.write(f"  质量等级: {best_assembly['Quality_Grade']}\n\n")
        except:
            file_handle.write("BUSCO评估结果: 未找到或解析失败\n\n")
        
        # 尝试读取统计结果
        try:
            stats_csv = Path(self.config.output_dir) / 'statistics' / 'assembly_statistics_comparison.csv'
            if stats_csv.exists():
                df = pd.read_csv(stats_csv)
                best_stats = df.loc[df['N50_Mb'].idxmax()]
                file_handle.write(f"组装统计摘要:\n")
                file_handle.write(f"  最高N50组装: {best_stats['Assembly']}\n")
                file_handle.write(f"  N50: {best_stats['N50_Mb']:.1f} Mb\n")
                file_handle.write(f"  总长度: {best_stats['Total_Length_Gb']:.2f} Gb\n")
                file_handle.write(f"  Contig数量: {best_stats['Num_Contigs']:,}\n\n")
        except:
            file_handle.write("组装统计摘要: 未找到或解析失败\n\n")
    
    def _write_analysis_recommendations(self, file_handle):
        """写入分析建议 | Write analysis recommendations"""
        file_handle.write("分析建议 | Analysis Recommendations\n")
        file_handle.write("-"*40 + "\n")
        
        recommendations = [
            "1. 查看 quality_assessment/ 目录中的详细质量报告",
            "2. 根据研究目标选择合适的组装版本：",
            "   - 基因注释：使用primary组装",
            "   - 变异分析：使用单倍型组装",
            "   - 比较基因组学：选择质量最好的组装",
            "3. 检查BUSCO完整性评估结果",
            "4. 根据N50和contig数量评估组装连续性",
            "5. 如果进行功能性二倍体研究，关注单倍型差异分析",
            "6. 建议进行进一步的质量控制和验证"
        ]
        
        for rec in recommendations:
            file_handle.write(f"{rec}\n")
        
        file_handle.write("\n")
        file_handle.write("详细使用说明请参考各目录下的README文件和报告。\n")
        file_handle.write("For detailed usage instructions, please refer to README files and reports in each directory.\n")

class SummaryGenerator:
    """摘要生成器 | Summary Generator"""
    
    def __init__(self, config, logger: logging.Logger):
        self.config = config
        self.logger = logger
        self.output_dir = Path(config.output_dir)
    
    def generate_summary(self) -> bool:
        """生成分析摘要 | Generate analysis summary"""
        try:
            self.logger.info("生成最终分析摘要 | Generating final analysis summary")
            
            # 生成HTML报告
            if not self._generate_html_report():
                return False
            
            # 生成README文件
            if not self._generate_readme():
                return False
            
            # 生成使用指南
            if not self._generate_usage_guide():
                return False
            
            self.logger.success("最终摘要生成完成 | Final summary generation completed")
            return True
            
        except Exception as e:
            self.logger.error(f"摘要生成失败 | Summary generation failed: {e}")
            return False
    
    def _generate_html_report(self) -> bool:
        """生成HTML报告 | Generate HTML report"""
        try:
            self.logger.info("生成HTML报告 | Generating HTML report")
            
            html_file = self.output_dir / 'HiFiasm_Analysis_Report.html'
            
            # 收集所有分析数据
            analysis_data = self._collect_analysis_data()
            
            # 生成HTML内容
            html_content = self._create_html_content(analysis_data)
            
            with open(html_file, 'w', encoding='utf-8') as f:
                f.write(html_content)
            
            self.logger.success(f"HTML报告已生成 | HTML report generated: {html_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成HTML报告失败 | Failed to generate HTML report: {e}")
            return False
    
    def _collect_analysis_data(self) -> Dict:
        """收集分析数据 | Collect analysis data"""
        data = {
            'config': self.config.to_dict(),
            'timestamp': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'busco_results': None,
            'statistics': None,
            'haplotype_analysis': None
        }
        
        # 尝试读取BUSCO结果
        try:
            busco_csv = self.output_dir / 'quality_assessment' / 'busco' / 'busco_comparison_report.csv'
            if busco_csv.exists():
                data['busco_results'] = pd.read_csv(busco_csv).to_dict('records')
        except:
            pass
        
        # 尝试读取统计结果
        try:
            stats_csv = self.output_dir / 'statistics' / 'assembly_statistics_comparison.csv'
            if stats_csv.exists():
                data['statistics'] = pd.read_csv(stats_csv).to_dict('records')
        except:
            pass
        
        return data
    
    def _create_html_content(self, data: Dict) -> str:
        """创建HTML内容 | Create HTML content"""
        html_template = """
        <!DOCTYPE html>
        <html lang="zh-CN">
        <head>
            <meta charset="UTF-8">
            <meta name="viewport" content="width=device-width, initial-scale=1.0">
            <title>HiFiasm基因组组装分析报告 | HiFiasm Genome Assembly Analysis Report</title>
            <style>
                body {{ font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 20px; background-color: #f5f5f5; }}
                .container {{ max-width: 1200px; margin: 0 auto; background-color: white; padding: 30px; border-radius: 10px; box-shadow: 0 0 20px rgba(0,0,0,0.1); }}
                .header {{ text-align: center; margin-bottom: 40px; }}
                .header h1 {{ color: #2c3e50; margin-bottom: 10px; }}
                .header p {{ color: #7f8c8d; font-size: 18px; }}
                .section {{ margin-bottom: 30px; }}
                .section h2 {{ color: #34495e; border-bottom: 2px solid #3498db; padding-bottom: 10px; }}
                .info-grid {{ display: grid; grid-template-columns: repeat(auto-fit, minmax(300px, 1fr)); gap: 20px; }}
                .info-card {{ background-color: #ecf0f1; padding: 20px; border-radius: 8px; border-left: 4px solid #3498db; }}
                .info-card h3 {{ margin-top: 0; color: #2c3e50; }}
                .table {{ width: 100%; border-collapse: collapse; margin-top: 20px; }}
                .table th, .table td {{ padding: 12px; text-align: left; border-bottom: 1px solid #ddd; }}
                .table th {{ background-color: #3498db; color: white; }}
                .table tr:hover {{ background-color: #f5f5f5; }}
                .status-good {{ color: #27ae60; font-weight: bold; }}
                .status-warning {{ color: #f39c12; font-weight: bold; }}
                .status-error {{ color: #e74c3c; font-weight: bold; }}
                .footer {{ text-align: center; margin-top: 40px; color: #7f8c8d; }}
            </style>
        </head>
        <body>
            <div class="container">
                <div class="header">
                    <h1>HiFiasm基因组组装分析报告<br>HiFiasm Genome Assembly Analysis Report</h1>
                    <p>样本 | Sample: {sample_name}</p>
                    <p>生成时间 | Generated: {timestamp}</p>
                </div>
                
                <div class="section">
                    <h2>分析配置 | Analysis Configuration</h2>
                    <div class="info-grid">
                        <div class="info-card">
                            <h3>输入参数 | Input Parameters</h3>
                            <p><strong>输入文件:</strong> {input_file}</p>
                            <p><strong>基因组大小:</strong> {genome_size}</p>
                            <p><strong>线程数:</strong> {threads}</p>
                            <p><strong>组装类型:</strong> {assembly_type}</p>
                        </div>
                        <div class="info-card">
                            <h3>质量控制 | Quality Control</h3>
                            <p><strong>BUSCO评估:</strong> {busco_status}</p>
                            <p><strong>QUAST评估:</strong> {quast_status}</p>
                            <p><strong>单倍型分析:</strong> {haplotype_status}</p>
                        </div>
                    </div>
                </div>
                
                {busco_section}
                {statistics_section}
                
                <div class="section">
                    <h2>建议与下一步 | Recommendations & Next Steps</h2>
                    <div class="info-card">
                        <ul>
                            <li>查看详细的质量评估报告以了解组装质量</li>
                            <li>根据研究目标选择合适的组装版本</li>
                            <li>进行进一步的功能注释和分析</li>
                            <li>Check detailed quality assessment reports for assembly quality</li>
                            <li>Choose appropriate assembly version based on research objectives</li>
                            <li>Proceed with functional annotation and downstream analysis</li>
                        </ul>
                    </div>
                </div>
                
                <div class="footer">
                    <p>该报告由HiFiasm分析流水线自动生成 | This report was automatically generated by HiFiasm analysis pipeline</p>
                </div>
            </div>
        </body>
        </html>
        """
        
        # 格式化模板
        busco_section = self._create_busco_section(data.get('busco_results'))
        statistics_section = self._create_statistics_section(data.get('statistics'))
        
        return html_template.format(
            sample_name=self.config.prefix,
            timestamp=data['timestamp'],
            input_file=os.path.basename(self.config.input_reads),
            genome_size=self.config.estimate_genome_size(),
            threads=self.config.threads,
            assembly_type=self.config.assembly_type,
            busco_status="已启用" if not self.config.skip_busco else "已跳过",
            quast_status="已启用" if not self.config.skip_quast else "已跳过",
            haplotype_status="已启用" if self.config.analyze_haplotypes else "已跳过",
            busco_section=busco_section,
            statistics_section=statistics_section
        )
    
    def _create_busco_section(self, busco_data: Optional[List[Dict]]) -> str:
        """创建BUSCO部分 | Create BUSCO section"""
        if not busco_data:
            return '<div class="section"><h2>BUSCO评估结果 | BUSCO Assessment Results</h2><p>无可用数据 | No data available</p></div>'
        
        table_rows = ""
        for result in busco_data:
            table_rows += f"""
            <tr>
                <td>{result['Assembly']}</td>
                <td class="status-good">{result['Complete(%)']}%</td>
                <td>{result['Single(%)']}%</td>
                <td>{result['Duplicated(%)']}%</td>
                <td>{result['Fragmented(%)']}%</td>
                <td>{result['Missing(%)']}%</td>
                <td>{result['Quality_Grade']}</td>
            </tr>
            """
        
        return f"""
        <div class="section">
            <h2>BUSCO评估结果 | BUSCO Assessment Results</h2>
            <table class="table">
                <thead>
                    <tr>
                        <th>组装 | Assembly</th>
                        <th>完整性 | Complete (%)</th>
                        <th>单拷贝 | Single (%)</th>
                        <th>重复 | Duplicated (%)</th>
                        <th>片段 | Fragmented (%)</th>
                        <th>缺失 | Missing (%)</th>
                        <th>质量等级 | Quality Grade</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>
        """
    
    def _create_statistics_section(self, stats_data: Optional[List[Dict]]) -> str:
        """创建统计部分 | Create statistics section"""
        if not stats_data:
            return '<div class="section"><h2>组装统计 | Assembly Statistics</h2><p>无可用数据 | No data available</p></div>'
        
        table_rows = ""
        for result in stats_data:
            table_rows += f"""
            <tr>
                <td>{result['Assembly']}</td>
                <td>{result['Total_Length_Gb']} Gb</td>
                <td>{result['Num_Contigs']:,}</td>
                <td>{result['N50_Mb']} Mb</td>
                <td>{result['Longest_Contig_Mb']} Mb</td>
                <td>{result['GC_Content_%']}%</td>
                <td>{result['N_Content_%']}%</td>
            </tr>
            """
        
        return f"""
        <div class="section">
            <h2>组装统计 | Assembly Statistics</h2>
            <table class="table">
                <thead>
                    <tr>
                        <th>组装 | Assembly</th>
                        <th>总长度 | Total Length</th>
                        <th>Contig数 | Num Contigs</th>
                        <th>N50</th>
                        <th>最长Contig | Longest Contig</th>
                        <th>GC含量 | GC Content</th>
                        <th>N含量 | N Content</th>
                    </tr>
                </thead>
                <tbody>
                    {table_rows}
                </tbody>
            </table>
        </div>
        """
    
    def _generate_readme(self) -> bool:
        """生成README文件 | Generate README file"""
        try:
            readme_file = self.output_dir / 'README.md'
            
            readme_content = f"""# HiFiasm基因组组装分析结果 | HiFiasm Genome Assembly Analysis Results

## 基本信息 | Basic Information

- **样本名称 | Sample Name**: {self.config.prefix}
- **输入文件 | Input File**: {self.config.input_reads}
- **分析时间 | Analysis Time**: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}
- **基因组大小估计 | Estimated Genome Size**: {self.config.estimate_genome_size()}
- **组装类型 | Assembly Type**: {self.config.assembly_type}

## 目录结构 | Directory Structure

```
{self.config.output_dir}/
├── final_results/           # 最终结果文件 | Final result files
│   ├── assemblies/         # 组装序列文件 | Assembly sequence files
│   └── quality_assessment/ # 质量评估结果 | Quality assessment results
├── assembly/               # 原始组装输出 | Raw assembly output
├── logs/                   # 分析日志 | Analysis logs
├── statistics/             # 统计分析 | Statistical analysis
└── HiFiasm_Analysis_Report.html  # HTML报告 | HTML report
```

## 主要文件说明 | Main Files Description

### 组装文件 | Assembly Files
- `*.primary.fasta`: 主要组装序列 | Primary assembly sequences
- `*.haplotype1.fasta`: 单倍型1序列 | Haplotype 1 sequences  
- `*.haplotype2.fasta`: 单倍型2序列 | Haplotype 2 sequences

### 质量评估 | Quality Assessment
- `busco/`: BUSCO完整性评估结果 | BUSCO completeness assessment results
- `statistics/`: 组装统计信息 | Assembly statistics
- `haplotype_analysis/`: 单倍型差异分析 | Haplotype difference analysis

## 使用建议 | Usage Recommendations

1. **基因注释 | Gene Annotation**: 使用primary组装 | Use primary assembly
2. **变异分析 | Variant Analysis**: 使用单倍型组装 | Use haplotype assemblies
3. **比较基因组学 | Comparative Genomics**: 选择质量最好的组装 | Choose the best quality assembly

## 质量评估标准 | Quality Assessment Criteria

- **优秀 | Excellent**: BUSCO完整性 ≥95% | BUSCO completeness ≥95%
- **良好 | Good**: BUSCO完整性 90-95% | BUSCO completeness 90-95%
- **中等 | Fair**: BUSCO完整性 85-90% | BUSCO completeness 85-90%

## 技术支持 | Technical Support

如有问题请查看日志文件或联系分析团队。
For questions, please check log files or contact the analysis team.
"""
            
            with open(readme_file, 'w', encoding='utf-8') as f:
                f.write(readme_content)
            
            self.logger.success(f"README文件已生成 | README file generated: {readme_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成README失败 | Failed to generate README: {e}")
            return False
    
    def _generate_usage_guide(self) -> bool:
        """生成使用指南 | Generate usage guide"""
        try:
            guide_file = self.output_dir / 'USAGE_GUIDE.txt'
            
            guide_content = f"""HiFiasm分析结果使用指南 | HiFiasm Analysis Results Usage Guide
{'='*80}

1. 选择合适的组装文件 | Choose Appropriate Assembly File
   - 基因注释: {self.config.prefix}.primary.fasta
   - 变异检测: {self.config.prefix}.haplotype1.fasta, {self.config.prefix}.haplotype2.fasta
   - 比较基因组学: 根据质量评估选择最佳组装

2. 质量检查 | Quality Check
   - 查看 BUSCO 评估报告了解基因完整性
   - 检查组装统计了解连续性 (N50, contig数量)
   - 查看GC含量和N含量评估序列质量

3. 下游分析建议 | Downstream Analysis Recommendations
   - 使用 MAKER 或 Augustus 进行基因预测
   - 使用 RepeatMasker 进行重复序列注释
   - 进行功能注释 (InterProScan, eggNOG-mapper)

4. 数据发布准备 | Data Publishing Preparation
   - 检查序列命名规范
   - 准备元数据描述文件
   - 验证文件完整性

5. 常见问题 | Common Issues
   - 如果BUSCO完整性较低，检查输入数据质量
   - 如果N50较低，可能需要额外的脚手架构建
   - 单倍型差异大可能反映真实的生物学特征

详细技术文档请参考 HiFiasm 官方文档。
For detailed technical documentation, please refer to HiFiasm official documentation.
"""
            
            with open(guide_file, 'w', encoding='utf-8') as f:
                f.write(guide_content)
            
            self.logger.success(f"使用指南已生成 | Usage guide generated: {guide_file}")
            return True
            
        except Exception as e:
            self.logger.error(f"生成使用指南失败 | Failed to generate usage guide: {e}")
            return False