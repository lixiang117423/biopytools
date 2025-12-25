"""
基因组组装结果处理模块 | Genome Assembly Results Processing Module 📋
"""

import os
import json
import shutil
from pathlib import Path
from typing import Dict, List, Tuple, Optional, Any
from .utils import CommandRunner, format_time, format_size
from datetime import datetime

class ResultsManager:
    """结果管理器 | Results Manager"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.results_dir = Path(config.output_dir) / "results"
        self.results_dir.mkdir(parents=True, exist_ok=True)
        
        # 创建子目录
        (self.results_dir / "final_assemblies").mkdir(exist_ok=True)
        (self.results_dir / "scaffolds").mkdir(exist_ok=True)
        (self.results_dir / "for_juicebox").mkdir(exist_ok=True)
        (self.results_dir / "statistics").mkdir(exist_ok=True)
    
    def organize_final_results(self, assembly_files: Dict[str, List[str]], 
                             scaffold_results: Dict[str, Dict[str, str]] = None,
                             quality_results: Dict[str, Dict[str, Any]] = None) -> Dict[str, Any]:
        """整理最终结果 | Organize final results"""
        self.logger.info("📋 整理最终结果 | Organizing final results")
        
        final_results = {
            'project_info': self._get_project_info(),
            'assemblies': {},
            'scaffolds': {},
            'quality_summary': {},
            'juicebox_files': {},
            'file_manifest': {}
        }
        
        # 整理组装结果
        final_results['assemblies'] = self._organize_assemblies(assembly_files)
        
        # 整理挂载结果
        if scaffold_results:
            final_results['scaffolds'] = self._organize_scaffolds(scaffold_results)
            
            # 准备Juicebox文件
            final_results['juicebox_files'] = self._prepare_juicebox_files(scaffold_results)
        
        # 整理质量评估结果
        if quality_results:
            final_results['quality_summary'] = self._summarize_quality_results(quality_results)
        
        # 生成文件清单
        final_results['file_manifest'] = self._generate_file_manifest()
        
        # 生成最终报告
        self._generate_final_report(final_results)
        
        # 生成JSON格式的结果摘要
        self._save_results_json(final_results)
        
        # 创建README文件
        self._create_readme_file(final_results)
        
        self.logger.info("✅ 最终结果整理完成 | Final results organization completed")
        return final_results
    
    def _get_project_info(self) -> Dict[str, Any]:
        """获取项目信息 | Get project information"""
        return {
            'project_name': self.config.project_name,
            'assembly_strategy': self.config.assembly_strategy,
            'hic_strategy': self.config.hic_strategy,
            'detected_data_types': self.config.detected_data_types,
            'genome_size': self.config.genome_size,
            'species_type': self.config.species_type,
            'analysis_date': datetime.now().strftime('%Y-%m-%d %H:%M:%S'),
            'output_directory': self.config.output_dir
        }
    
    def _organize_assemblies(self, assembly_files: Dict[str, List[str]]) -> Dict[str, Any]:
        """整理组装结果 | Organize assembly results"""
        organized = {}
        
        for assembly_type, fasta_list in assembly_files.items():
            organized[assembly_type] = {}
            
            for fasta_file in fasta_list:
                file_path = Path(fasta_file)
                assembly_name = file_path.stem
                
                # 复制到最终结果目录
                final_path = self.results_dir / "final_assemblies" / file_path.name
                if not final_path.exists():
                    shutil.copy2(fasta_file, final_path)
                    self.logger.info(f"📁 复制组装文件: {file_path.name}")
                
                # 创建索引
                self._create_fasta_index(final_path)
                
                organized[assembly_type][assembly_name] = {
                    'original_path': fasta_file,
                    'final_path': str(final_path),
                    'index_path': str(final_path) + '.fai',
                    'file_size': format_size(final_path.stat().st_size),
                    'description': self._generate_assembly_description(assembly_type, assembly_name)
                }
        
        return organized
    
    def _organize_scaffolds(self, scaffold_results: Dict[str, Dict[str, str]]) -> Dict[str, Any]:
        """整理挂载结果 | Organize scaffolding results"""
        organized = {}
        
        for assembly_type, results in scaffold_results.items():
            organized[assembly_type] = {}
            
            for assembly_name, files in results.items():
                organized[assembly_type][assembly_name] = {}
                
                # 整理不同类型的输出文件
                for file_type, file_path in files.items():
                    if os.path.exists(file_path):
                        file_path_obj = Path(file_path)
                        
                        # 根据文件类型决定目标目录
                        if file_type in ['scaffolds_final', 'final_fasta', 'scaffolds_FINAL']:
                            target_dir = self.results_dir / "scaffolds"
                            new_name = f"{self.config.project_name}_{assembly_type}_{assembly_name}_scaffolded.fa"
                        elif file_type == 'hic_file':
                            target_dir = self.results_dir / "for_juicebox"
                            new_name = f"{self.config.project_name}_{assembly_type}_{assembly_name}.hic"
                        else:
                            target_dir = self.results_dir / "scaffolds"
                            new_name = f"{self.config.project_name}_{assembly_type}_{assembly_name}_{file_type}.{file_path_obj.suffix[1:]}"
                        
                        # 复制文件
                        final_path = target_dir / new_name
                        if not final_path.exists():
                            shutil.copy2(file_path, final_path)
                            self.logger.info(f"📁 复制挂载文件: {new_name}")
                        
                        # 为FASTA文件创建索引
                        if final_path.suffix.lower() in ['.fa', '.fasta']:
                            self._create_fasta_index(final_path)
                        
                        organized[assembly_type][assembly_name][file_type] = {
                            'original_path': file_path,
                            'final_path': str(final_path),
                            'file_size': format_size(final_path.stat().st_size),
                            'description': self._generate_file_description(file_type)
                        }
        
        return organized
    
    def _prepare_juicebox_files(self, scaffold_results: Dict[str, Dict[str, str]]) -> Dict[str, Any]:
        """准备Juicebox可视化文件 | Prepare Juicebox visualization files"""
        juicebox_files = {}
        juicebox_dir = self.results_dir / "for_juicebox"
        
        # 创建Juicebox使用说明
        instructions_file = juicebox_dir / "JUICEBOX_INSTRUCTIONS.txt"
        with open(instructions_file, 'w', encoding='utf-8') as f:
            f.write("🔗 Juicebox手动校正指南 | Juicebox Manual Correction Guide\n")
            f.write("=" * 60 + "\n\n")
            f.write("1. 下载并安装Juicebox Desktop版本\n")
            f.write("   Download and install Juicebox Desktop\n")
            f.write("   网址 | URL: https://github.com/aidenlab/Juicebox\n\n")
            f.write("2. 加载文件 | Load files:\n")
            f.write("   - .hic文件：Hi-C接触图 | Hi-C contact map\n")
            f.write("   - .assembly文件：组装信息 | Assembly information\n\n")
            f.write("3. 手动校正步骤 | Manual correction steps:\n")
            f.write("   a) 检查染色体边界清晰度\n")
            f.write("   b) 识别组装错误（倒位、易位等）\n")
            f.write("   c) 标记和修正错误区域\n")
            f.write("   d) 导出修正后的.review.assembly文件\n\n")
            f.write("4. 应用修正 | Apply corrections:\n")
            f.write("   使用3D-DNA的post-review脚本应用修正\n")
            f.write("   Use 3D-DNA post-review script to apply corrections\n\n")
        
        # 收集所有Juicebox相关文件
        for assembly_type, results in scaffold_results.items():
            juicebox_files[assembly_type] = {}
            
            for assembly_name, files in results.items():
                assembly_juicebox_files = {}
                
                # 收集.hic, .assembly, .review.assembly文件
                for file_type, file_path in files.items():
                    if file_type in ['hic_file', 'assembly_file', 'review_file'] and os.path.exists(file_path):
                        # 文件已经在_organize_scaffolds中复制了
                        file_name = Path(file_path).name
                        juicebox_path = juicebox_dir / f"{self.config.project_name}_{assembly_type}_{assembly_name}_{file_name}"
                        
                        if juicebox_path.exists() or file_type == 'hic_file':
                            assembly_juicebox_files[file_type] = str(juicebox_path)
                
                if assembly_juicebox_files:
                    juicebox_files[assembly_type][assembly_name] = assembly_juicebox_files
        
        return juicebox_files
    
    def _summarize_quality_results(self, quality_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """汇总质量评估结果 | Summarize quality assessment results"""
        summary = {
            'overall_grade': 'N/A',
            'key_metrics': {},
            'recommendations': [],
            'best_assembly': None,
            'best_scaffold': None
        }
        
        try:
            # 提取关键指标
            if 'assemblies' in quality_results:
                assembly_metrics = []
                for assembly_type, results in quality_results['assemblies'].items():
                    for name, data in results.items():
                        metrics = {
                            'name': f"{assembly_type}_{name}",
                            'type': 'assembly'
                        }
                        
                        # 提取N50
                        if 'continuity' in data:
                            metrics['n50'] = data['continuity'].get('N50', 0)
                        
                        # 提取大小信息
                        if 'basic_stats' in data:
                            stats = data['basic_stats']
                            if 'sum_len' in stats:
                                metrics['size'] = int(stats['sum_len'].replace(',', '')) if stats['sum_len'] != 'N/A' else 0
                            if 'num' in stats:
                                metrics['sequences'] = int(stats['num']) if stats['num'] != 'N/A' else 0
                        
                        assembly_metrics.append(metrics)
                
                # 找到最佳组装
                if assembly_metrics:
                    best_assembly = max(assembly_metrics, key=lambda x: x.get('n50', 0))
                    summary['best_assembly'] = best_assembly['name']
                    summary['key_metrics']['best_assembly_n50'] = best_assembly.get('n50', 0)
            
            # 处理挂载结果
            if 'scaffolds' in quality_results:
                scaffold_metrics = []
                for assembly_type, results in quality_results['scaffolds'].items():
                    for name, data in results.items():
                        metrics = {
                            'name': f"{assembly_type}_{name}",
                            'type': 'scaffold'
                        }
                        
                        # 提取挂载指标
                        if 'scaffolding_metrics' in data:
                            sm = data['scaffolding_metrics']
                            metrics['n50'] = sm.get('N50', 0)
                            
                            if 'chromosome_analysis' in sm:
                                ca = sm['chromosome_analysis']
                                metrics['large_scaffolds'] = ca.get('large_scaffolds', 0)
                                metrics['estimated_chromosomes'] = ca.get('estimated_chromosomes', 0)
                        
                        scaffold_metrics.append(metrics)
                
                # 找到最佳挂载
                if scaffold_metrics:
                    best_scaffold = max(scaffold_metrics, key=lambda x: x.get('n50', 0))
                    summary['best_scaffold'] = best_scaffold['name']
                    summary['key_metrics']['best_scaffold_n50'] = best_scaffold.get('n50', 0)
                    summary['key_metrics']['estimated_chromosomes'] = best_scaffold.get('estimated_chromosomes', 0)
            
            # 处理BUSCO结果
            if 'busco' in quality_results and quality_results['busco']:
                best_busco = max(quality_results['busco'].items(), 
                               key=lambda x: x[1].get('complete_percent', 0))
                summary['key_metrics']['best_busco_name'] = best_busco[0]
                summary['key_metrics']['best_busco_complete'] = best_busco[1].get('complete_percent', 0)
            
            # 生成总体评级
            summary['overall_grade'] = self._calculate_overall_grade(quality_results)
            
            # 生成建议
            summary['recommendations'] = self._generate_recommendations(quality_results, summary)
        
        except Exception as e:
            self.logger.warning(f"⚠️ 质量结果汇总失败: {e}")
        
        return summary
    
    def _calculate_overall_grade(self, quality_results: Dict[str, Dict[str, Any]]) -> str:
        """计算总体评级 | Calculate overall grade"""
        try:
            grades = quality_results.get('grades', {})
            if not grades:
                return "未评级"
            
            # 提取所有评级的分数
            grade_scores = []
            for grade_str in grades.values():
                if 'A+' in grade_str:
                    grade_scores.append(95)
                elif 'A' in grade_str:
                    grade_scores.append(90)
                elif 'B' in grade_str:
                    grade_scores.append(80)
                elif 'C' in grade_str:
                    grade_scores.append(70)
                elif 'D' in grade_str:
                    grade_scores.append(60)
                elif 'F' in grade_str:
                    grade_scores.append(50)
            
            if grade_scores:
                avg_score = sum(grade_scores) / len(grade_scores)
                if avg_score >= 90:
                    return "A (优秀)"
                elif avg_score >= 80:
                    return "B (良好)"
                elif avg_score >= 70:
                    return "C (中等)"
                elif avg_score >= 60:
                    return "D (及格)"
                else:
                    return "F (不及格)"
        except:
            pass
        
        return "未评级"
    
    def _generate_recommendations(self, quality_results: Dict[str, Dict[str, Any]], summary: Dict[str, Any]) -> List[str]:
        """生成改进建议 | Generate improvement recommendations"""
        recommendations = []
        
        try:
            # 基于BUSCO结果的建议
            if 'busco' in quality_results:
                best_complete = summary['key_metrics'].get('best_busco_complete', 0)
                if best_complete < 90:
                    recommendations.append("BUSCO完整性偏低，建议检查数据质量和组装参数")
                elif best_complete < 95:
                    recommendations.append("BUSCO完整性良好，可考虑进一步优化组装参数")
            
            # 基于连续性的建议
            best_n50 = summary['key_metrics'].get('best_assembly_n50', 0)
            if best_n50 < 100000:  # <100kb
                recommendations.append("N50偏低，建议增加长读长数据或优化组装策略")
            elif best_n50 < 1000000:  # <1Mb
                recommendations.append("N50中等，可考虑添加Hi-C数据进行染色体挂载")
            
            # 基于挂载结果的建议
            if 'Hi-C' in self.config.detected_data_types:
                estimated_chroms = summary['key_metrics'].get('estimated_chromosomes', 0)
                if estimated_chroms == 0:
                    recommendations.append("Hi-C挂载失败，建议检查Hi-C数据质量和比对参数")
                elif estimated_chroms > 100:
                    recommendations.append("染色体数量过多，建议进行手动校正或调整挂载参数")
            else:
                recommendations.append("建议补充Hi-C数据以获得染色体级别的组装")
            
            # 基于数据类型的建议
            if 'ONT' not in self.config.detected_data_types:
                recommendations.append("可考虑补充ONT数据以提高组装连续性")
            
            if not recommendations:
                recommendations.append("组装质量良好，建议进行下游分析")
        
        except Exception as e:
            self.logger.warning(f"⚠️ 建议生成失败: {e}")
            recommendations.append("请参考质量评估报告进行进一步分析")
        
        return recommendations
    
    def _create_fasta_index(self, fasta_file: Path):
        """创建FASTA索引 | Create FASTA index"""
        try:
            cmd = f"{self.config.samtools_path} faidx {fasta_file}"
            self.cmd_runner.run(cmd, f"创建索引 - {fasta_file.name}")
        except Exception as e:
            self.logger.warning(f"⚠️ 无法创建索引 {fasta_file.name}: {e}")
    
    def _generate_assembly_description(self, assembly_type: str, assembly_name: str) -> str:
        """生成组装描述 | Generate assembly description"""
        descriptions = {
            'primary': f"主要组装 - {assembly_name}",
            'alternate': f"备选组装 - {assembly_name}",
            'haplotype_1': f"单倍型1组装 - {assembly_name}",
            'haplotype_2': f"单倍型2组装 - {assembly_name}"
        }
        return descriptions.get(assembly_type, f"{assembly_type}组装 - {assembly_name}")
    
    def _generate_file_description(self, file_type: str) -> str:
        """生成文件描述 | Generate file description"""
        descriptions = {
            'hic_file': 'Hi-C接触图文件，用于Juicebox可视化',
            'assembly_file': '3D-DNA组装信息文件',
            'review_file': '待手动校正的组装文件',
            'final_fasta': '最终染色体级基因组序列',
            'scaffolds_final': 'SALSA2最终挂载序列',
            'scaffolds_FINAL': 'SALSA2最终挂载序列',
            'agp_file': 'AGP格式的挂载信息文件'
        }
        return descriptions.get(file_type, f'{file_type}文件')
    
    def _generate_file_manifest(self) -> Dict[str, List[str]]:
        """生成文件清单 | Generate file manifest"""
        manifest = {}
        
        # 扫描结果目录
        for subdir in self.results_dir.iterdir():
            if subdir.is_dir():
                manifest[subdir.name] = []
                for file_path in subdir.rglob('*'):
                    if file_path.is_file():
                        manifest[subdir.name].append(file_path.name)
        
        return manifest
    
    def _generate_final_report(self, final_results: Dict[str, Any]):
        """生成最终报告 | Generate final report"""
        report_file = self.results_dir / "FINAL_ASSEMBLY_REPORT.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🧬 基因组组装最终报告 | Final Genome Assembly Report\n")
            f.write("=" * 70 + "\n\n")
            
            # 项目信息
            project_info = final_results['project_info']
            f.write("📋 项目信息 | Project Information:\n")
            f.write(f"  项目名称 | Project Name: {project_info['project_name']}\n")
            f.write(f"  组装策略 | Assembly Strategy: {project_info['assembly_strategy']}\n")
            f.write(f"  Hi-C策略 | Hi-C Strategy: {project_info['hic_strategy']}\n")
            f.write(f"  数据类型 | Data Types: {', '.join(project_info['detected_data_types'])}\n")
            f.write(f"  预估基因组大小 | Estimated Genome Size: {project_info['genome_size']}\n")
            f.write(f"  分析日期 | Analysis Date: {project_info['analysis_date']}\n\n")
            
            # 质量摘要
            quality_summary = final_results.get('quality_summary', {})
            if quality_summary:
                f.write("🏆 质量评估摘要 | Quality Assessment Summary:\n")
                f.write(f"  总体评级 | Overall Grade: {quality_summary.get('overall_grade', 'N/A')}\n")
                
                key_metrics = quality_summary.get('key_metrics', {})
                if 'best_assembly_n50' in key_metrics:
                    f.write(f"  最佳组装N50 | Best Assembly N50: {key_metrics['best_assembly_n50']:,} bp\n")
                if 'best_scaffold_n50' in key_metrics:
                    f.write(f"  最佳挂载N50 | Best Scaffold N50: {key_metrics['best_scaffold_n50']:,} bp\n")
                if 'best_busco_complete' in key_metrics:
                    f.write(f"  最高BUSCO完整性 | Best BUSCO Completeness: {key_metrics['best_busco_complete']:.1f}%\n")
                if 'estimated_chromosomes' in key_metrics:
                    f.write(f"  预估染色体数 | Estimated Chromosomes: {key_metrics['estimated_chromosomes']}\n")
                
                f.write("\n")
            
            # 主要输出文件
            f.write("📁 主要输出文件 | Main Output Files:\n")
            f.write("\n组装文件 | Assembly Files:\n")
            for assembly_type, assemblies in final_results.get('assemblies', {}).items():
                for name, info in assemblies.items():
                    f.write(f"  📄 {info['description']}\n")
                    f.write(f"     文件 | File: {Path(info['final_path']).name}\n")
                    f.write(f"     大小 | Size: {info['file_size']}\n\n")
            
            if final_results.get('scaffolds'):
                f.write("挂载文件 | Scaffold Files:\n")
                for assembly_type, scaffolds in final_results['scaffolds'].items():
                    for name, files in scaffolds.items():
                        f.write(f"  🔗 {assembly_type}_{name} 挂载结果:\n")
                        for file_type, info in files.items():
                            if file_type in ['scaffolds_final', 'final_fasta', 'scaffolds_FINAL']:
                                f.write(f"     📄 {Path(info['final_path']).name} ({info['file_size']})\n")
                        f.write("\n")
            
            # Juicebox文件
            if final_results.get('juicebox_files'):
                f.write("🔗 Juicebox手动校正文件 | Juicebox Manual Correction Files:\n")
                f.write("位置 | Location: for_juicebox/\n")
                f.write("说明文件 | Instructions: JUICEBOX_INSTRUCTIONS.txt\n\n")
            
            # 建议和下一步
            if quality_summary.get('recommendations'):
                f.write("💡 改进建议 | Recommendations:\n")
                for i, rec in enumerate(quality_summary['recommendations'], 1):
                    f.write(f"  {i}. {rec}\n")
                f.write("\n")
            
            f.write("🎯 下一步分析 | Next Steps:\n")
            f.write("  1. 如有Juicebox文件，进行手动校正\n")
            f.write("  2. 基因组注释 (基因预测、重复序列注释等)\n")
            f.write("  3. 比较基因组学分析\n")
            f.write("  4. 变异检测和群体基因组学分析\n")
            f.write("  5. 表观基因组学分析 (如有Hi-C数据)\n\n")
            
            f.write("📞 技术支持 | Technical Support:\n")
            f.write("  如需技术支持，请提供完整的日志文件和错误信息\n")
            f.write("  For technical support, please provide complete log files and error messages\n")
        
        self.logger.info(f"📋 最终报告已生成: {report_file}")
    
    def _save_results_json(self, final_results: Dict[str, Any]):
        """保存JSON格式的结果摘要 | Save JSON format results summary"""
        json_file = self.results_dir / "assembly_results_summary.json"
        
        try:
            with open(json_file, 'w', encoding='utf-8') as f:
                json.dump(final_results, f, indent=2, ensure_ascii=False)
            self.logger.info(f"📋 JSON结果摘要已保存: {json_file}")
        except Exception as e:
            self.logger.warning(f"⚠️ JSON保存失败: {e}")
    
    def _create_readme_file(self, final_results: Dict[str, Any]):
        """创建README文件 | Create README file"""
        readme_file = self.results_dir / "README.md"
        
        with open(readme_file, 'w', encoding='utf-8') as f:
            f.write("# 基因组组装结果 | Genome Assembly Results\n\n")
            
            # 项目信息
            project_info = final_results['project_info']
            f.write("## 项目信息 | Project Information\n\n")
            f.write(f"- **项目名称 | Project Name**: {project_info['project_name']}\n")
            f.write(f"- **组装策略 | Assembly Strategy**: {project_info['assembly_strategy']}\n")
            f.write(f"- **Hi-C策略 | Hi-C Strategy**: {project_info['hic_strategy']}\n")
            f.write(f"- **数据类型 | Data Types**: {', '.join(project_info['detected_data_types'])}\n")
            f.write(f"- **分析日期 | Analysis Date**: {project_info['analysis_date']}\n\n")
            
            # 目录结构
            f.write("## 目录结构 | Directory Structure\n\n")
            f.write("```\n")
            f.write("results/\n")
            f.write("├── final_assemblies/     # 最终组装序列 | Final assembly sequences\n")
            f.write("├── scaffolds/           # 染色体挂载结果 | Chromosome scaffolding results\n") 
            f.write("├── for_juicebox/        # Juicebox手动校正文件 | Files for Juicebox manual correction\n")
            f.write("├── statistics/          # 统计和质量报告 | Statistics and quality reports\n")
            f.write("├── FINAL_ASSEMBLY_REPORT.txt  # 最终报告 | Final report\n")
            f.write("└── README.md           # 本文件 | This file\n")
            f.write("```\n\n")
            
            # 主要文件
            f.write("## 主要文件 | Main Files\n\n")
            f.write("### 组装序列 | Assembly Sequences\n\n")
            
            for assembly_type, assemblies in final_results.get('assemblies', {}).items():
                f.write(f"**{assembly_type.upper()}**:\n")
                for name, info in assemblies.items():
                    f.write(f"- `{Path(info['final_path']).name}` - {info['description']} ({info['file_size']})\n")
                f.write("\n")
            
            if final_results.get('scaffolds'):
                f.write("### 挂载序列 | Scaffolded Sequences\n\n")
                for assembly_type, scaffolds in final_results['scaffolds'].items():
                    f.write(f"**{assembly_type.upper()} Scaffolds**:\n")
                    for name, files in scaffolds.items():
                        for file_type, info in files.items():
                            if file_type in ['scaffolds_final', 'final_fasta', 'scaffolds_FINAL']:
                                f.write(f"- `{Path(info['final_path']).name}` ({info['file_size']})\n")
                    f.write("\n")
            
            # 质量评估
            quality_summary = final_results.get('quality_summary', {})
            if quality_summary:
                f.write("## 质量评估 | Quality Assessment\n\n")
                f.write(f"- **总体评级 | Overall Grade**: {quality_summary.get('overall_grade', 'N/A')}\n")
                
                key_metrics = quality_summary.get('key_metrics', {})
                if key_metrics:
                    f.write("- **关键指标 | Key Metrics**:\n")
                    if 'best_assembly_n50' in key_metrics:
                        f.write(f"  - 最佳组装N50 | Best Assembly N50: {key_metrics['best_assembly_n50']:,} bp\n")
                    if 'best_scaffold_n50' in key_metrics:
                        f.write(f"  - 最佳挂载N50 | Best Scaffold N50: {key_metrics['best_scaffold_n50']:,} bp\n")
                    if 'best_busco_complete' in key_metrics:
                        f.write(f"  - 最高BUSCO完整性 | Best BUSCO Completeness: {key_metrics['best_busco_complete']:.1f}%\n")
                f.write("\n")
            
            # 使用说明
            f.write("## 使用说明 | Usage Instructions\n\n")
            f.write("### 组装序列使用 | Using Assembly Sequences\n\n")
            f.write("组装序列文件为标准FASTA格式，可直接用于:\n\n")
            f.write("- 基因组注释 (基因预测、重复序列注释)\n")
            f.write("- 比较基因组学分析\n")
            f.write("- 变异检测\n")
            f.write("- 其他下游分析\n\n")
            
            if final_results.get('juicebox_files'):
                f.write("### Juicebox手动校正 | Juicebox Manual Correction\n\n")
                f.write("1. 安装Juicebox Desktop: https://github.com/aidenlab/Juicebox\n")
                f.write("2. 加载`for_juicebox/`目录中的.hic和.assembly文件\n")
                f.write("3. 进行手动校正\n")
                f.write("4. 导出.review.assembly文件\n")
                f.write("5. 使用3D-DNA post-review脚本应用修正\n\n")
                f.write("详细说明请参考: `for_juicebox/JUICEBOX_INSTRUCTIONS.txt`\n\n")
            
            # 引用信息
            f.write("## 引用信息 | Citation\n\n")
            f.write("如果使用了这些结果，请引用相关工具:\n\n")
            f.write("- **Hifiasm**: Cheng, H., et al. (2021). Haplotype-resolved de novo assembly using phased assembly graphs with hifiasm. Nature methods, 18(2), 170-175.\n")
            if 'Hi-C' in project_info['detected_data_types']:
                if project_info['hic_strategy'] == 'complete_juicer':
                    f.write("- **Juicer**: Durand, N.C., et al. (2016). Juicer provides a one-click system for analyzing loop-resolution Hi-C experiments. Cell systems, 3(1), 95-98.\n")
                    f.write("- **3D-DNA**: Dudchenko, O., et al. (2017). De novo assembly of the Aedes aegypti genome using Hi-C yields chromosome-length scaffolds. Science, 356(6333), 92-95.\n")
                elif project_info['hic_strategy'] == 'simplified_salsa2':
                    f.write("- **SALSA2**: Ghurye, J., et al. (2019). Integrating Hi-C links with assembly graphs for chromosome-scale assembly. PLoS computational biology, 15(8), e1007273.\n")
        
        self.logger.info(f"📋 README文件已创建: {readme_file}")
