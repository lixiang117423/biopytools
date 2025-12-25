"""
🌾 EDTA结果处理模块 | EDTA Results Processing Module
"""

import json
from pathlib import Path
from datetime import datetime
from typing import Dict, Any
from .utils import save_checkpoint

class ResultsGenerator:
    """结果生成器 | Results Generator"""
    
    def __init__(self, config, logger, directories: Dict[str, Path]):
        self.config = config
        self.logger = logger
        self.directories = directories
    
    def generate_all_reports(self, edta_results: Dict[str, Dict], 
                           processed_results: Dict[str, Any]):
        """生成所有报告 | Generate all reports"""
        self.logger.info("📋 生成综合分析报告 | Generating comprehensive analysis reports")
        
        try:
            # 生成主要摘要报告
            self._generate_main_summary_report(edta_results, processed_results)
            
            # 生成详细统计报告
            self._generate_detailed_statistics_report(processed_results)
            
            # 生成结果比较报告（如果有多个结果）
            successful_results = {k: v for k, v in processed_results.items() 
                                if v.get("status") == "processed"}
            if len(successful_results) > 1 and self.config.compare_results:
                self._generate_comparative_report(successful_results)
            
            # 保存检查点
            save_checkpoint(self.config.output_path, "reports_generated", 
                          {"timestamp": datetime.now().isoformat()}, self.logger)
            
            self.logger.info("✅ 所有报告生成完成 | All reports generation completed")
            
        except Exception as e:
            self.logger.error(f"❌ 报告生成失败 | Reports generation failed: {e}")
    
    def _generate_main_summary_report(self, edta_results: Dict[str, Dict], 
                                    processed_results: Dict[str, Any]):
        """生成主要摘要报告 | Generate main summary report"""
        self.logger.info("📄 生成主要摘要报告 | Generating main summary report")
        
        report_file = self.directories["reports"] / "edta_analysis_summary.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🌾 EDTA植物基因组TE注释分析摘要报告 | EDTA Plant Genome TE Annotation Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 基本信息
            f.write("📋 基本分析信息 | Basic Analysis Information:\n")
            f.write("-" * 40 + "\n")
            f.write(f"分析时间 | Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"输出目录 | Output directory: {self.config.output_dir}\n")
            
            if self.config.batch_mode:
                f.write(f"分析模式 | Analysis mode: 批量模式 | Batch mode\n")
                f.write(f"基因组数量 | Number of genomes: {len(self.config.genome_list)}\n")
            else:
                f.write(f"分析模式 | Analysis mode: 单基因组模式 | Single genome mode\n")
                f.write(f"基因组文件 | Genome file: {self.config.genome}\n")
            
            f.write(f"物种设置 | Species setting: {self.config.species}\n")
            f.write(f"分析步骤 | Analysis step: {self.config.step}\n")
            f.write(f"使用线程数 | Threads used: {self.config.threads}\n")
            f.write("\n")
            
            # 分析结果概览
            f.write("📊 分析结果概览 | Analysis Results Overview:\n")
            f.write("-" * 40 + "\n")
            
            total_analyses = len(edta_results)
            successful_analyses = len([r for r in edta_results.values() if r.get("status") == "success"])
            
            f.write(f"总分析数 | Total analyses: {total_analyses}\n")
            f.write(f"成功分析数 | Successful analyses: {successful_analyses}\n")
            f.write(f"失败分析数 | Failed analyses: {total_analyses - successful_analyses}\n")
            if total_analyses > 0:
                f.write(f"成功率 | Success rate: {successful_analyses/total_analyses*100:.1f}%\n")
            f.write("\n")
            
            # 详细结果
            for analysis_name, result_data in processed_results.items():
                if result_data.get("status") == "processed":
                    f.write(f"📋 分析详情 | Analysis Details: {analysis_name}\n")
                    f.write("-" * 30 + "\n")
                    
                    annotation_stats = result_data.get("annotation_stats", {})
                    if annotation_stats and "error" not in annotation_stats:
                        f.write(f"  TE元件总数 | Total TE elements: {annotation_stats.get('total_elements', 0):,}\n")
                        f.write(f"  TE总长度 | Total TE length: {annotation_stats.get('total_length', 0):,} bp\n")
                        
                        # 染色体分布
                        chr_dist = annotation_stats.get("chromosome_distribution", {})
                        if chr_dist:
                            f.write(f"  涉及染色体数 | Chromosomes involved: {len(chr_dist)}\n")
                        
                        # 特征分布
                        feature_dist = annotation_stats.get("feature_distribution", {})
                        if feature_dist:
                            f.write("  主要TE特征类型 | Main TE feature types:\n")
                            for feature, count in sorted(feature_dist.items(), key=lambda x: x[1], reverse=True)[:3]:
                                percentage = count/annotation_stats.get('total_elements', 1)*100
                                f.write(f"    - {feature}: {count:,} ({percentage:.1f}%)\n")
                    else:
                        f.write("  ⚠️ 注释统计数据不可用 | Annotation statistics not available\n")
                    
                    f.write("\n")
            
            # 使用建议
            f.write("💡 结果使用建议 | Recommendations for Using Results:\n")
            f.write("-" * 40 + "\n")
            f.write("1. 主要注释文件 | Main annotation files:\n")
            f.write("   - *.TEanno.gff3: 用于基因组注释和TE分析 | For genome annotation and TE analysis\n")
            f.write("   - *.TElib.fa: TE库文件，可用于其他基因组 | TE library for other genomes\n")
            f.write("   - *.MAKER.masked: 屏蔽基因组，用于基因预测 | Masked genome for gene prediction\n")
            f.write("2. 统计图表 | Statistical plots:\n")
            f.write("   - 查看04_visualization/目录中的图表 | Check plots in 04_visualization/ directory\n")
            f.write("\n")
        
        self.logger.info(f"✅ 主要摘要报告已保存 | Main summary report saved: {report_file}")
    
    def _generate_detailed_statistics_report(self, processed_results: Dict[str, Any]):
        """生成详细统计报告 | Generate detailed statistics report"""
        self.logger.info("📊 生成详细统计报告 | Generating detailed statistics report")
        
        report_file = self.directories["reports"] / "detailed_statistics_report.json"
        
        # 准备详细统计数据
        detailed_stats = {
            "analysis_timestamp": datetime.now().isoformat(),
            "configuration": {
                "species": self.config.species,
                "step": self.config.step,
                "threads": self.config.threads,
                "sensitive_mode": bool(self.config.sensitive),
                "batch_mode": self.config.batch_mode
            },
            "detailed_results": processed_results
        }
        
        # 保存为JSON格式
        with open(report_file, 'w', encoding='utf-8') as f:
            json.dump(detailed_stats, f, indent=2, ensure_ascii=False)
        
        self.logger.info(f"✅ 详细统计报告已保存 | Detailed statistics report saved: {report_file}")
    
    def _generate_comparative_report(self, successful_results: Dict[str, Any]):
        """生成比较分析报告 | Generate comparative analysis report"""
        self.logger.info("🔍 生成比较分析报告 | Generating comparative analysis report")
        
        report_file = self.directories["reports"] / "comparative_analysis_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🔍 EDTA分析结果比较报告 | EDTA Analysis Results Comparison Report\n")
            f.write("=" * 70 + "\n\n")
            
            f.write("📋 比较概览 | Comparison Overview:\n")
            f.write("-" * 30 + "\n")
            f.write(f"比较的分析数量 | Number of analyses compared: {len(successful_results)}\n\n")
            
            # TE数量比较
            f.write("📊 TE元件数量比较 | TE Element Count Comparison:\n")
            f.write("-" * 40 + "\n")
            element_counts = {}
            for name, result in successful_results.items():
                count = result.get("annotation_stats", {}).get("total_elements", 0)
                element_counts[name] = count
                f.write(f"  {name}: {count:,} elements\n")
            
            if element_counts:
                max_analysis = max(element_counts.items(), key=lambda x: x[1])
                min_analysis = min(element_counts.items(), key=lambda x: x[1])
                f.write(f"\n  最多TE元件 | Most TE elements: {max_analysis[0]} ({max_analysis[1]:,})\n")
                f.write(f"  最少TE元件 | Fewest TE elements: {min_analysis[0]} ({min_analysis[1]:,})\n")
            f.write("\n")
        
        self.logger.info(f"✅ 比较分析报告已保存 | Comparative analysis report saved: {report_file}")
