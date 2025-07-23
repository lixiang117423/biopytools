"""
主分析器模块 | Main Analyzer Module
"""

import sys
import numpy as np
from .vcf_reader import VCFReader
from .data_filter import DataFilter
from .ld_calculator import LDCalculator
from .plot_generator import PlotGenerator
from .matrix_exporter import MatrixExporter

class LDHeatmapAnalyzer:
    """LD热图分析器 | LD Heatmap Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        
        # 初始化各个组件 | Initialize components
        self.vcf_reader = VCFReader(config, logger)
        self.data_filter = DataFilter(config, logger)
        self.ld_calculator = LDCalculator(config, logger)
        self.plot_generator = PlotGenerator(config, logger)
        self.matrix_exporter = MatrixExporter(config, logger)
    
    def run_analysis(self):
        """运行LD分析流程 | Run LD analysis pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始LD热图分析 | Starting LD heatmap analysis")
            self.logger.info("=" * 60)
            
            # 步骤1: 读取VCF数据 | Step 1: Read VCF data
            self.logger.info("\n步骤1: 读取VCF数据 | Step 1: Reading VCF data")
            gt, positions, chrom, samples = self.vcf_reader.read_vcf_data()
            
            if gt.shape[0] == 0:
                raise ValueError("错误: 没有找到有效的变异位点 | Error: No valid variants found")
            
            # 步骤2: 过滤变异位点 | Step 2: Filter variants
            self.logger.info("\n步骤2: 过滤变异位点 | Step 2: Filtering variants")
            gt_filtered, positions_filtered = self.data_filter.filter_variants(gt, positions)
            
            if gt_filtered.shape[0] < 2:
                raise ValueError("错误: 过滤后的变异位点数量不足 | Error: Insufficient variants after filtering")
            
            # 步骤3: 计算LD矩阵 | Step 3: Calculate LD matrix
            self.logger.info("\n步骤3: 计算LD矩阵 | Step 3: Calculating LD matrix")
            ld_matrix = self.ld_calculator.calculate_ld_matrix(gt_filtered)
            
            # 步骤4: 生成热图 | Step 4: Generate heatmap
            self.logger.info("\n步骤4: 生成热图 | Step 4: Generating heatmap")
            self.plot_generator.create_ld_heatmap(ld_matrix, positions_filtered)
            
            # 步骤5: 保存LD矩阵 | Step 5: Save LD matrix
            if self.config.save_matrix:
                self.logger.info("\n步骤5: 保存LD矩阵 | Step 5: Saving LD matrix")
                self.matrix_exporter.save_ld_matrix(ld_matrix, positions_filtered)
            
            # 完成信息 | Completion information
            self.logger.info("\n" + "=" * 60)
            self.logger.info("LD热图分析完成！| LD heatmap analysis completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"最终包含 {len(positions_filtered)} 个SNP | Final {len(positions_filtered)} SNPs included")
            self.logger.info(f"平均r² = {np.nanmean(ld_matrix):.3f} | Average r² = {np.nanmean(ld_matrix):.3f}")
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_file}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)
