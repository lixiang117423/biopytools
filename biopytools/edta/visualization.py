"""
🌾 EDTA可视化模块 | EDTA Visualization Module
"""

import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
from pathlib import Path
from typing import Dict, List, Any
from collections import Counter

class EDTAVisualizer:
    """EDTA结果可视化器 | EDTA Results Visualizer"""
    
    def __init__(self, config, logger, directories: Dict[str, Path]):
        self.config = config
        self.logger = logger
        self.directories = directories
        
        # 设置绘图样式
        plt.style.use('seaborn-v0_8')
        sns.set_palette("husl")
        
        # 设置中文字体支持
        plt.rcParams['font.sans-serif'] = ['DejaVu Sans', 'SimHei', 'Arial Unicode MS']
        plt.rcParams['axes.unicode_minus'] = False
    
    def generate_all_plots(self, processed_results: Dict[str, Any]):
        """生成所有可视化图表 | Generate all visualization plots"""
        if not self.config.generate_plots:
            self.logger.info("📊 跳过图表生成 | Skipping plot generation")
            return
        
        self.logger.info("🎨 生成可视化图表 | Generating visualization plots")
        
        try:
            # 为每个成功的分析生成图表
            for analysis_name, result_data in processed_results.items():
                if result_data.get("status") == "processed":
                    self._generate_single_analysis_plots(analysis_name, result_data)
            
            # 如果有多个结果，生成比较图表
            successful_results = {k: v for k, v in processed_results.items() 
                                if v.get("status") == "processed"}
            
            if len(successful_results) > 1 and self.config.compare_results:
                self._generate_comparative_plots(successful_results)
            
            self.logger.info("✅ 可视化图表生成完成 | Visualization plots generation completed")
            
        except Exception as e:
            self.logger.error(f"❌ 可视化图表生成失败 | Visualization plots generation failed: {e}")
    
    def _generate_single_analysis_plots(self, analysis_name: str, result_data: Dict[str, Any]):
        """为单个分析生成图表 | Generate plots for single analysis"""
        self.logger.info(f"📊 生成图表 | Generating plots for: {analysis_name}")
        
        # 创建分析专用的可视化目录
        vis_dir = self.directories["visualization"] / analysis_name
        vis_dir.mkdir(parents=True, exist_ok=True)
        
        annotation_stats = result_data.get("annotation_stats", {})
        
        if annotation_stats and "error" not in annotation_stats:
            # 1. TE元件数量分布饼图
            self._plot_te_feature_distribution(annotation_stats, vis_dir, analysis_name)
            
            # 2. TE分类分布条形图
            self._plot_te_classification_distribution(annotation_stats, vis_dir, analysis_name)
            
            # 3. 染色体分布柱状图
            self._plot_chromosome_distribution(annotation_stats, vis_dir, analysis_name)
            
            # 4. TE链方向分布
            self._plot_strand_distribution(annotation_stats, vis_dir, analysis_name)
        else:
            self.logger.warning(f"⚠️ 无有效注释统计数据，跳过图表生成 | No valid annotation stats, skipping plots: {analysis_name}")
    
    def _plot_te_feature_distribution(self, stats: Dict, vis_dir: Path, analysis_name: str):
        """绘制TE特征分布饼图 | Plot TE feature distribution pie chart"""
        try:
            feature_dist = stats.get("feature_distribution", {})
            if not feature_dist:
                return
            
            plt.figure(figsize=(10, 8))
            
            labels = list(feature_dist.keys())
            sizes = list(feature_dist.values())
            colors = sns.color_palette("husl", len(labels))
            
            plt.pie(sizes, labels=labels, colors=colors, autopct='%1.1f%%', startangle=90)
            plt.title(f'🧬 TE特征类型分布 | TE Feature Type Distribution\n{analysis_name}', 
                     fontsize=14, fontweight='bold')
            plt.axis('equal')
            
            plt.tight_layout()
            output_file = vis_dir / f"{analysis_name}_te_feature_distribution.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 TE特征分布图已保存 | TE feature distribution plot saved: {output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ TE特征分布图生成失败 | TE feature distribution plot failed: {e}")
    
    def _plot_te_classification_distribution(self, stats: Dict, vis_dir: Path, analysis_name: str):
        """绘制TE分类分布条形图 | Plot TE classification distribution bar chart"""
        try:
            class_dist = stats.get("classification_distribution", {})
            if not class_dist or len(class_dist) == 0:
                self.logger.warning(f"⚠️ 无TE分类数据 | No TE classification data: {analysis_name}")
                return
            
            plt.figure(figsize=(12, 8))
            
            # 准备数据
            classifications = list(class_dist.keys())
            counts = list(class_dist.values())
            
            # 如果分类太多，只显示前15个
            if len(classifications) > 15:
                sorted_items = sorted(class_dist.items(), key=lambda x: x[1], reverse=True)
                top_items = sorted_items[:14]
                others_count = sum(item[1] for item in sorted_items[14:])
                
                classifications = [item[0] for item in top_items] + ["Others"]
                counts = [item[1] for item in top_items] + [others_count]
            
            # 创建条形图
            bars = plt.bar(range(len(classifications)), counts, 
                          color=sns.color_palette("husl", len(classifications)))
            
            plt.title(f'🏷️ TE分类分布 | TE Classification Distribution\n{analysis_name}', 
                     fontsize=14, fontweight='bold')
            plt.xlabel('TE分类 | TE Classification', fontsize=12)
            plt.ylabel('数量 | Count', fontsize=12)
            
            # 设置x轴标签
            plt.xticks(range(len(classifications)), classifications, rotation=45, ha='right')
            
            # 在条形图上添加数值标签
            for bar, count in zip(bars, counts):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(counts)*0.01, 
                        str(count), ha='center', va='bottom', fontsize=9)
            
            plt.tight_layout()
            output_file = vis_dir / f"{analysis_name}_te_classification_distribution.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 TE分类分布图已保存 | TE classification distribution plot saved: {output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ TE分类分布图生成失败 | TE classification distribution plot failed: {e}")
    
    def _plot_chromosome_distribution(self, stats: Dict, vis_dir: Path, analysis_name: str):
        """绘制染色体分布柱状图 | Plot chromosome distribution bar chart"""
        try:
            chr_dist = stats.get("chromosome_distribution", {})
            if not chr_dist:
                return
            
            plt.figure(figsize=(14, 8))
            
            chromosomes = list(chr_dist.keys())
            counts = list(chr_dist.values())
            
            bars = plt.bar(chromosomes, counts, color=sns.color_palette("viridis", len(chromosomes)))
            
            plt.title(f'🧬 TE染色体分布 | TE Chromosome Distribution\n{analysis_name}', 
                     fontsize=14, fontweight='bold')
            plt.xlabel('染色体 | Chromosome', fontsize=12)
            plt.ylabel('TE数量 | TE Count', fontsize=12)
            
            # 旋转x轴标签以避免重叠
            plt.xticks(rotation=45, ha='right')
            
            plt.tight_layout()
            output_file = vis_dir / f"{analysis_name}_chromosome_distribution.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 染色体分布图已保存 | Chromosome distribution plot saved: {output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 染色体分布图生成失败 | Chromosome distribution plot failed: {e}")
    
    def _plot_strand_distribution(self, stats: Dict, vis_dir: Path, analysis_name: str):
        """绘制TE链方向分布 | Plot TE strand distribution"""
        try:
            strand_dist = stats.get("strand_distribution", {})
            if not strand_dist:
                return
            
            plt.figure(figsize=(8, 8))
            
            strands = list(strand_dist.keys())
            counts = list(strand_dist.values())
            colors = ['#ff9999', '#66b3ff', '#99ff99', '#ffcc99'][:len(strands)]
            
            plt.pie(counts, labels=strands, colors=colors, autopct='%1.1f%%', startangle=90)
            plt.title(f'⬅️➡️ TE链方向分布 | TE Strand Distribution\n{analysis_name}', 
                     fontsize=14, fontweight='bold')
            plt.axis('equal')
            
            plt.tight_layout()
            output_file = vis_dir / f"{analysis_name}_strand_distribution.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 链方向分布图已保存 | Strand distribution plot saved: {output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 链方向分布图生成失败 | Strand distribution plot failed: {e}")
    
    def _generate_comparative_plots(self, successful_results: Dict[str, Any]):
        """生成比较图表 | Generate comparative plots"""
        self.logger.info("📊 生成比较图表 | Generating comparative plots")
        
        try:
            # 创建比较图表目录
            comp_dir = self.directories["visualization"] / "comparative_analysis"
            comp_dir.mkdir(parents=True, exist_ok=True)
            
            # 比较TE总数
            self._plot_comparative_te_counts(successful_results, comp_dir)
            
            self.logger.info("✅ 比较图表生成完成 | Comparative plots generation completed")
            
        except Exception as e:
            self.logger.error(f"❌ 比较图表生成失败 | Comparative plots generation failed: {e}")
    
    def _plot_comparative_te_counts(self, results: Dict[str, Any], comp_dir: Path):
        """绘制TE数量比较图 | Plot comparative TE counts"""
        try:
            analysis_names = []
            te_counts = []
            
            for analysis_name, result_data in results.items():
                annotation_stats = result_data.get("annotation_stats", {})
                if annotation_stats and "error" not in annotation_stats:
                    analysis_names.append(analysis_name)
                    te_counts.append(annotation_stats.get("total_elements", 0))
            
            if len(analysis_names) < 2:
                return
            
            plt.figure(figsize=(10, 6))
            
            bars = plt.bar(analysis_names, te_counts, 
                          color=sns.color_palette("viridis", len(analysis_names)))
            
            plt.title('🔍 TE数量比较 | TE Count Comparison', fontsize=14, fontweight='bold')
            plt.xlabel('分析 | Analysis', fontsize=12)
            plt.ylabel('TE数量 | TE Count', fontsize=12)
            
            # 旋转x轴标签
            plt.xticks(rotation=45, ha='right')
            
            # 添加数值标签
            for bar, count in zip(bars, te_counts):
                plt.text(bar.get_x() + bar.get_width()/2, bar.get_height() + max(te_counts)*0.01, 
                        str(count), ha='center', va='bottom', fontsize=10)
            
            plt.tight_layout()
            output_file = comp_dir / "comparative_te_counts.png"
            plt.savefig(output_file, dpi=300, bbox_inches='tight')
            plt.close()
            
            self.logger.info(f"📊 TE数量比较图已保存 | TE count comparison plot saved: {output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ TE数量比较图生成失败 | TE count comparison plot failed: {e}")
