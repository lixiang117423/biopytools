# ===== FILE: orthofinder_pangenome/main.py =====
"""
OrthoFinder泛基因组分析主程序模块 | OrthoFinder Pangenome Analysis Main Module
"""

import argparse
import sys
from pathlib import Path
import pandas as pd
from .config import PangenomeConfig
from .utils import PangenomeLogger, CommandRunner, check_dependencies, validate_fasta_files, count_fasta_files
from .orthofinder_runner import OrthoFinderRunner
from .pangenome_classifier import PangenomeClassifier
from .rarefaction_analysis import RarefactionAnalyzer
from .visualization import PangenomeVisualizer
from .results import ResultsProcessor
from .sequence_extractor import SequenceExtractor

class OrthoFinderPangenomeAnalyzer:
    """OrthoFinder泛基因组分析主类 | Main OrthoFinder Pangenome Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PangenomeConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PangenomeLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        # self.orthofinder_runner = OrthoFinderRunner(self.config, self.logger, self.cmd_runner)
        # self.pangenome_classifier = PangenomeClassifier(self.config, self.logger)
        # self.sequence_extractor = SequenceExtractor(self.config, self.logger)
        # self.rarefaction_analyzer = RarefactionAnalyzer(self.config, self.logger)
        # self.visualizer = PangenomeVisualizer(self.config, self.logger)
        # self.results_processor = ResultsProcessor(self.config, self.logger)

        # 初始化各个处理器 | Initialize processors
        self.orthofinder_runner = OrthoFinderRunner(self.config, self.logger, self.cmd_runner)
        self.pangenome_classifier = PangenomeClassifier(self.config, self.logger)
        self.rarefaction_analyzer = RarefactionAnalyzer(self.config, self.logger)
        self.visualizer = PangenomeVisualizer(self.config, self.logger)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.sequence_extractor = SequenceExtractor(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)

    # def run_analysis(self):
    #     """运行完整的泛基因组分析流程 | Run complete pangenome analysis pipeline"""
    #     try:
    #         self.logger.info("开始OrthoFinder泛基因组分析流程 | Starting OrthoFinder pangenome analysis pipeline")
    #         self.logger.info(f"输入目录 | Input directory: {self.config.input_dir}")
    #         self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")

    #         # 检查依赖 | Check dependencies
    #         # self.check_dependencies()  # 临时跳过依赖检查
    #         self.logger.info("跳过依赖检查 | Skipping dependency check")

    #         # 验证输入文件 | Validate input files
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info("步骤1: 验证输入数据 | Step 1: Validate input data")
    #         self.logger.info(f"{'=' * 60}")

    #         if not validate_fasta_files(self.config.input_dir, self.logger):
    #             raise RuntimeError("输入文件验证失败 | Input file validation failed")

    #         # 获取基因组数量信息 | Get genome count information
    #         genome_count, sample_names = count_fasta_files(self.config.input_dir)
    #         self.logger.info(f"检测到 {genome_count} 个基因组 | Detected {genome_count} genomes")

    #         # 运行OrthoFinder | Run OrthoFinder
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info("步骤2: 运行OrthoFinder同源基因群分析 | Step 2: Run OrthoFinder orthogroups analysis")
    #         self.logger.info(f"{'=' * 60}")

    #         if not self.orthofinder_runner.run_orthofinder():
    #             raise RuntimeError("OrthoFinder分析失败 | OrthoFinder analysis failed")

    #         # 加载OrthoFinder结果 | Load OrthoFinder results
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info("步骤3: 加载和解析OrthoFinder结果 | Step 3: Load and parse OrthoFinder results")
    #         self.logger.info(f"{'=' * 60}")

    #         orthogroups_file = self.orthofinder_runner.get_orthogroups_file()
    #         gene_count_file = self.orthofinder_runner.get_gene_count_file()

    #         self.pangenome_classifier.load_orthofinder_results(orthogroups_file, gene_count_file)

    #         # 进行泛基因组分类 | Perform pangenome classification
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info("步骤4: 泛基因组分类分析 | Step 4: Pangenome classification analysis")
    #         self.logger.info(f"{'=' * 60}")

    #         classification_results = self.pangenome_classifier.classify_pangenome()

    #         # 单拷贝基因分析 | Single copy gene analysis
    #         if self.config.enable_single_copy_analysis:
    #             step_num += 1
    #             self.logger.info(f"{'=' * 60}")
    #             self.logger.info(f"步骤{step_num}: 单拷贝基因分析和序列提取 | Step {step_num}: Single copy gene analysis and sequence extraction")
    #             self.logger.info(f"{'=' * 60}")
                
    #             single_copy_count = len(classification_results['single_copy'])
    #             self.logger.info(f"识别到 {single_copy_count} 个单拷贝同源群 | Identified {single_copy_count} single copy orthogroups")
                
    #             if self.config.extract_sequences:
    #                 # 加载基因组序列 | Load genome sequences
    #                 self.sequence_extractor.load_genome_sequences()
                    
    #                 # 提取单拷贝基因序列 | Extract single copy gene sequences
    #                 self.sequence_extractor.extract_single_copy_sequences(
    #                     classification_results['single_copy'], 
    #                     self.config.output_path
    #                 )

    #         # 步骤5: 保存分类结果 (移动到稀释分析之前)
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info("步骤5: 保存分类结果 | Step 5: Save classification results")
    #         self.logger.info(f"{'=' * 60}")
            
    #         # 先保存主要结果文件，这样稀释分析才能使用它
    #         main_output_file = self.pangenome_classifier.save_classification_results(
    #             classification_results, self.config.output_path)
            
    #         # 初始化 rarefaction_results
    #         rarefaction_results = None

    #         # 稀释曲线分析 | Rarefaction curve analysis
    #         if self.config.enable_rarefaction:
    #             self.logger.info(f"{'=' * 60}")
    #             self.logger.info("步骤6: 泛基因组稀释曲线分析 | Step 6: Pangenome rarefaction curve analysis")
    #             self.logger.info(f"{'=' * 60}")
                
    #             # 加载分类结果到稀释分析器 (现在 main_output_file 已定义)
    #             self.rarefaction_analyzer.load_classification_results(main_output_file)
                
    #             # 运行稀释分析 | Run rarefaction analysis
    #             rarefaction_results = self.rarefaction_analyzer.run_rarefaction_analysis()
                
    #             # 保存稀释分析结果 | Save rarefaction results
    #             self.rarefaction_analyzer.save_rarefaction_results(rarefaction_results, self.config.output_path)
            
    #         # 获取其他结果 | Get other results
    #         step_num = 7 if self.config.enable_rarefaction else 6
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info(f"步骤{step_num}: 准备最终报告数据 | Step {step_num}: Prepare final report data")
    #         self.logger.info(f"{'=' * 60}")

    #         # 获取基因组统计信息 | Get genome statistics
    #         genome_stats = self.pangenome_classifier.get_genome_gene_counts(classification_results)
            
    #         # 获取频率分布 | Get frequency distribution
    #         frequency_distributions = self.pangenome_classifier.calculate_frequency_distribution(classification_results)
            
    #         # 保存详细结果表格 | Save detailed results table
    #         detailed_table_file = self.results_processor.save_gene_families_table(
    #             classification_results, self.config.output_path)
            
    #         # 生成可视化图表 | Generate visualization plots
    #         if self.config.generate_plots:
    #             step_num += 1
    #             self.logger.info(f"{'=' * 60}")
    #             self.logger.info(f"步骤{step_num}: 生成可视化图表 | Step {step_num}: Generate visualization plots")
    #             self.logger.info(f"{'=' * 60}")
                
    #             # 创建分类饼图 | Create classification pie chart
    #             self.visualizer.create_classification_pie_chart(classification_results, self.config.output_path)
                
    #             # 创建频率分布图 | Create frequency distribution plot
    #             self.visualizer.create_frequency_distribution_plot(frequency_distributions, self.config.output_path)
                
    #             # 创建稀释曲线图 | Create rarefaction curve plot
    #             if self.config.enable_rarefaction and rarefaction_results is not None:
    #                 # 注意: rarefaction_results['pan_counts'] 是详细列表，绘图需要的是统计摘要
    #                 # 从文件中读取摘要数据
    #                 summary_file = self.config.output_path / "rarefaction_curve_summary.tsv"
    #                 if summary_file.exists():
    #                     rarefaction_summary_data = pd.read_csv(summary_file, sep='\t')
    #                     self.visualizer.create_rarefaction_plot(rarefaction_summary_data, self.config.output_path)
    #                 else:
    #                     self.logger.warning("稀释曲线摘要文件未找到，无法绘图 | Rarefaction summary file not found, skipping plot.")

    #         # 生成综合报告 | Generate comprehensive report
    #         step_num += 1
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info(f"步骤{step_num}: 生成分析报告 | Step {step_num}: Generate analysis report")
    #         self.logger.info(f"{'=' * 60}")
            
    #         report_file = self.results_processor.generate_comprehensive_report(
    #             classification_results, genome_stats, frequency_distributions,
    #             rarefaction_results, # 传递原始字典给报告
    #             self.config.output_path)
            
    #         # 分析完成总结 | Analysis completion summary
    #         # self.logger.info(f"{'=' * 60}")
    #         # self.logger.info("泛基因组分析完成！| Pangenome analysis completed!")
    #         # self.logger.info(f"{'=' * 60}")
    #         # self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")

    #         # 分析完成总结 | Analysis completion summary
    #         self.logger.info(f"{'=' * 60}")
    #         self.logger.info("泛基因组分析完成！| Pangenome analysis completed!")
    #         self.logger.info(f"{'=' * 60}")

    #         # 拷贝OrthoFinder原始结果到输出目录 | Copy OrthoFinder original results to output directory
    #         self._copy_orthofinder_results_to_output()

    #         self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")
            
    #         # 输出关键统计信息 | Output key statistics
    #         total_orthogroups = sum(len(results) for results in classification_results.values())
    #         total_genes = sum(sum(len(gene_info[1]) for gene_info in results) 
    #                         for results in classification_results.values())
            
    #         self.logger.info(f"关键统计 | Key Statistics:")
    #         self.logger.info(f"   总基因组数 | Total genomes: {genome_count}")
    #         self.logger.info(f"   同源基因群数 | Orthogroups: {total_orthogroups:,}")
    #         self.logger.info(f"   总基因数 | Total genes: {total_genes:,}")
    #         self.logger.info(f"   核心基因群 | Core gene families: {len(classification_results['core']):,}")
    #         self.logger.info(f"   软核心基因群 | Softcore gene families: {len(classification_results['softcore']):,}")
    #         self.logger.info(f"   非必需基因群 | Dispensable gene families: {len(classification_results['dispensable']):,}")
    #         self.logger.info(f"   私有基因群 | Private gene families: {len(classification_results['private']):,}")
            
    #         self.logger.info(f"\n主要输出文件 | Main output files:")
    #         self.logger.info(f"   基因家族信息 | Gene families: {main_output_file.name}")
    #         self.logger.info(f"   详细数据表格 | Detailed table: {detailed_table_file.name}")
    #         self.logger.info(f"   综合报告 | Comprehensive report: {report_file.name}")
            
    #     except Exception as e:
    #         self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}", exc_info=True) # 添加 exc_info=True 获取更详细的回溯
    #         sys.exit(1)

    def run_analysis(self):
        """运行完整的泛基因组分析流程 | Run complete pangenome analysis pipeline"""
        try:
            self.logger.info("开始OrthoFinder泛基因组分析流程 | Starting OrthoFinder pangenome analysis pipeline")
            self.logger.info(f"输入目录 | Input directory: {self.config.input_dir}")
            self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
            
            # 检查依赖 | Check dependencies
            # self.check_dependencies()  # 临时跳过依赖检查
            self.logger.info("跳过依赖检查 | Skipping dependency check")
            
            # 验证输入文件 | Validate input files
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤1: 验证输入数据 | Step 1: Validate input data")
            self.logger.info(f"{'=' * 60}")
            
            if not validate_fasta_files(self.config.input_dir, self.logger):
                raise RuntimeError("输入文件验证失败 | Input file validation failed")
            
            # 获取基因组数量信息 | Get genome count information
            genome_count, sample_names = count_fasta_files(self.config.input_dir)
            self.logger.info(f"检测到 {genome_count} 个基因组 | Detected {genome_count} genomes")
            
            # 运行OrthoFinder | Run OrthoFinder
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤2: 运行OrthoFinder同源基因群分析 | Step 2: Run OrthoFinder orthogroups analysis")
            self.logger.info(f"{'=' * 60}")
            
            if not self.orthofinder_runner.run_orthofinder():
                raise RuntimeError("OrthoFinder分析失败 | OrthoFinder analysis failed")
            
            # 加载OrthoFinder结果 | Load OrthoFinder results
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤3: 加载和解析OrthoFinder结果 | Step 3: Load and parse OrthoFinder results")
            self.logger.info(f"{'=' * 60}")
            
            orthogroups_file = self.orthofinder_runner.get_orthogroups_file()
            gene_count_file = self.orthofinder_runner.get_gene_count_file()
            
            self.pangenome_classifier.load_orthofinder_results(orthogroups_file, gene_count_file)
            
            # 进行泛基因组分类 | Perform pangenome classification
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤4: 泛基因组分类分析 | Step 4: Pangenome classification analysis")
            self.logger.info(f"{'=' * 60}")
            
            classification_results = self.pangenome_classifier.classify_pangenome()
            
            # 初始化步骤计数器 | Initialize step counter
            current_step = 5
            
            # 稀释曲线分析 | Rarefaction curve analysis
            if self.config.enable_rarefaction:
                self.logger.info(f"{'=' * 60}")
                self.logger.info(f"步骤{current_step}: 泛基因组稀释曲线分析 | Step {current_step}: Pangenome rarefaction curve analysis")
                self.logger.info(f"{'=' * 60}")
                
                # 加载分类结果到稀释分析器 | Load classification results to rarefaction analyzer
                main_output_file = self.pangenome_classifier.save_classification_results(
                    classification_results, self.config.output_path)
                self.rarefaction_analyzer.load_classification_results(main_output_file)
                
                # 运行稀释分析 | Run rarefaction analysis
                rarefaction_results = self.rarefaction_analyzer.run_rarefaction_analysis()
                
                # 保存稀释分析结果 | Save rarefaction results
                self.rarefaction_analyzer.save_rarefaction_results(rarefaction_results, self.config.output_path)
                current_step += 1
            else:
                rarefaction_results = None
            
            # # 保存分类结果 | Save classification results
            # self.logger.info(f"{'=' * 60}")
            # self.logger.info(f"步骤{current_step}: 保存分析结果 | Step {current_step}: Save analysis results")
            # self.logger.info(f"{'=' * 60}")
            
            # # 如果没有进行稀释分析，需要先保存分类结果 | Save classification results if not done in rarefaction step
            # if not self.config.enable_rarefaction:
            #     main_output_file = self.pangenome_classifier.save_classification_results(
            #         classification_results, self.config.output_path)
            
            # # 获取基因组统计信息 | Get genome statistics
            # genome_stats = self.pangenome_classifier.get_genome_gene_counts(classification_results)
            
            # # 获取频率分布 | Get frequency distribution
            # frequency_distributions = self.pangenome_classifier.calculate_frequency_distribution(classification_results)
            
            # # 保存详细结果表格 | Save detailed results table
            # detailed_table_file = self.results_processor.save_gene_families_table(
            #     classification_results, self.config.output_path)
            
            # current_step += 1
            # 保存分类结果 | Save classification results
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"步骤{current_step}: 保存分析结果 | Step {current_step}: Save analysis results")
            self.logger.info(f"{'=' * 60}")

            # 如果没有进行稀释分析，需要先保存分类结果
            if not self.config.enable_rarefaction:
                self.logger.info("正在保存分类结果...")
                main_output_file = self.pangenome_classifier.save_classification_results(
                    classification_results, self.config.output_path)
                self.logger.info("分类结果保存完成")

            # 获取基因组统计信息
            self.logger.info("正在获取基因组统计信息...")
            genome_stats = self.pangenome_classifier.get_genome_gene_counts(classification_results)
            self.logger.info("基因组统计信息获取完成")

            # 获取频率分布
            self.logger.info("正在获取频率分布...")
            frequency_distributions = self.pangenome_classifier.calculate_frequency_distribution(classification_results)
            self.logger.info("频率分布获取完成")

            # 保存详细结果表格
            self.logger.info("正在保存详细结果表格...")
            detailed_table_file = self.results_processor.save_gene_families_table(
                classification_results, self.config.output_path)
            self.logger.info("详细结果表格保存完成")

            current_step += 1
            
            # 单拷贝基因分析 | Single copy gene analysis
            if self.config.enable_single_copy_analysis:
                self.logger.info(f"{'=' * 60}")
                self.logger.info(f"步骤{current_step}: 单拷贝基因分析和序列提取 | Step {current_step}: Single copy gene analysis and sequence extraction")
                self.logger.info(f"{'=' * 60}")
                
                single_copy_count = len(classification_results['single_copy'])
                self.logger.info(f"识别到 {single_copy_count} 个单拷贝同源群 | Identified {single_copy_count} single copy orthogroups")
                
                if self.config.extract_sequences:
                    # 加载基因组序列 | Load genome sequences
                    self.sequence_extractor.load_genome_sequences()
                    
                    # 提取单拷贝基因序列 | Extract single copy gene sequences
                    self.sequence_extractor.extract_single_copy_sequences(
                        classification_results['single_copy'], 
                        self.config.output_path
                    )
                current_step += 1
            
            # 生成可视化图表 | Generate visualization plots
            if self.config.generate_plots:
                self.logger.info(f"{'=' * 60}")
                self.logger.info(f"步骤{current_step}: 生成可视化图表 | Step {current_step}: Generate visualization plots")
                self.logger.info(f"{'=' * 60}")
                
                # 创建分类饼图 | Create classification pie chart
                self.visualizer.create_classification_pie_chart(classification_results, self.config.output_path)
                
                # 创建频率分布图 | Create frequency distribution plot
                self.visualizer.create_frequency_distribution_plot(frequency_distributions, self.config.output_path)
                
                # 创建稀释曲线图 | Create rarefaction curve plot
                if self.config.enable_rarefaction:
                    self.visualizer.create_rarefaction_plot(rarefaction_results, self.config.output_path)
                
                current_step += 1
            
            # 生成综合报告 | Generate comprehensive report
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"步骤{current_step}: 生成分析报告 | Step {current_step}: Generate analysis report")
            self.logger.info(f"{'=' * 60}")
            
            report_file = self.results_processor.generate_comprehensive_report(
                classification_results, genome_stats, frequency_distributions,
                rarefaction_results if self.config.enable_rarefaction else None,
                self.config.output_path)
            
            # 分析完成总结 | Analysis completion summary
            self.logger.info(f"{'=' * 60}")
            self.logger.info("泛基因组分析完成！| Pangenome analysis completed!")
            self.logger.info(f"{'=' * 60}")
            
            # 拷贝OrthoFinder原始结果到输出目录 | Copy OrthoFinder original results to output directory
            self._copy_orthofinder_results_to_output()
            
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")
            
            # 输出关键统计信息 | Output key statistics
            total_orthogroups = sum(len(results) for results in classification_results.values())
            total_genes = sum(sum(len(gene_info[1]) for gene_info in results) 
                            for results in classification_results.values())
            
            self.logger.info(f"关键统计 | Key Statistics:")
            self.logger.info(f"   总基因组数 | Total genomes: {genome_count}")
            self.logger.info(f"   同源基因群数 | Orthogroups: {total_orthogroups:,}")
            self.logger.info(f"   总基因数 | Total genes: {total_genes:,}")
            self.logger.info(f"   核心基因群 | Core gene families: {len(classification_results['core']):,}")
            self.logger.info(f"   软核心基因群 | Softcore gene families: {len(classification_results['softcore']):,}")
            self.logger.info(f"   非必需基因群 | Dispensable gene families: {len(classification_results['dispensable']):,}")
            self.logger.info(f"   私有基因群 | Private gene families: {len(classification_results['private']):,}")
            self.logger.info(f"   单拷贝基因群 | Single copy gene families: {len(classification_results['single_copy']):,}")
            
            self.logger.info(f"\n主要输出文件 | Main output files:")
            self.logger.info(f"   基因家族信息 | Gene families: {main_output_file.name}")
            self.logger.info(f"   详细数据表格 | Detailed table: {detailed_table_file.name}")
            self.logger.info(f"   综合报告 | Comprehensive report: {report_file.name}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

    def _copy_orthofinder_results_to_output(self):
        """拷贝OrthoFinder原始结果到输出目录 | Copy OrthoFinder original results to output directory"""
        if not self.orthofinder_runner.orthofinder_output_dir:
            self.logger.warning("未找到OrthoFinder结果目录，跳过拷贝 | OrthoFinder results directory not found, skipping copy")
            return
        
        import shutil
        
        source_dir = Path(self.orthofinder_runner.orthofinder_output_dir)
        target_dir = Path(self.config.output_dir) / "OrthoFinder_Results"
        
        # 可能存在的重复目录列表
        duplicate_dirs = [
            Path(self.config.output_dir) / "Results_OrthoFinder",
            Path(self.config.output_dir) / "OrthoFinder_WorkingDirectory"
        ]
        
        try:
            self.logger.info("拷贝OrthoFinder原始结果到输出目录 | Copying OrthoFinder original results to output directory")
            
            # 如果目标目录已存在，先删除
            if target_dir.exists():
                shutil.rmtree(target_dir)
            
            # 拷贝整个结果目录
            shutil.copytree(str(source_dir), str(target_dir))
            self.logger.info(f"OrthoFinder结果已拷贝到 | OrthoFinder results copied to: {target_dir}")
            
            # 同时拷贝工作目录（如果存在且还在原位置）
            input_path = Path(self.config.input_dir)
            orthofinder_work_dir = input_path / "OrthoFinder"
            
            if orthofinder_work_dir.exists():
                target_work_dir = Path(self.config.output_dir) / "OrthoFinder_WorkingDirectory"
                
                if target_work_dir.exists():
                    shutil.rmtree(target_work_dir)
                
                shutil.copytree(str(orthofinder_work_dir), str(target_work_dir))
                self.logger.info(f"OrthoFinder工作目录已拷贝到 | OrthoFinder working directory copied to: {target_work_dir}")
            
            # 清理可能存在的重复目录
            for duplicate_dir in duplicate_dirs:
                if duplicate_dir.exists() and duplicate_dir != target_dir:
                    self.logger.info(f"删除重复目录 | Removing duplicate directory: {duplicate_dir}")
                    shutil.rmtree(duplicate_dir)
            
            self.logger.info("OrthoFinder结果整理完成 | OrthoFinder results organization completed")
            
        except Exception as e:
            self.logger.error(f"拷贝OrthoFinder结果失败 | Failed to copy OrthoFinder results: {e}")

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='OrthoFinder泛基因组分析脚本 | OrthoFinder Pangenome Analysis Script',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -i protein_sequences/ -o pangenome_results
  %(prog)s --input data/proteomes/ --output results/ --softcore-threshold 1
  %(prog)s -i sequences/ -o output/ -s diamond -t 32 --dispensable-threshold 3
  %(prog)s -i proteomes/ -o results/ --force --generate-trees
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True, 
                       help='输入蛋白质序列文件目录 | Input protein sequence files directory')
    parser.add_argument('-o', '--output', required=True,
                       help='输出目录 | Output directory')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-n', '--project-name', 
                       help='项目名称 | Project name')
    
    # 泛基因组分类参数 | Pangenome classification parameters
    parser.add_argument('--softcore-threshold', type=int, default=2,
                       help='Softcore基因缺失阈值 | Softcore gene missing threshold (default: <=2)')
    parser.add_argument('--dispensable-threshold', type=int, default=2, 
                       help='Dispensable基因缺失阈值 | Dispensable gene missing threshold (default: >2)')
    
    # OrthoFinder参数 | OrthoFinder parameters
    parser.add_argument('-t', '--threads', type=int, default=88, 
                       help='线程数 | Number of threads')
    parser.add_argument('-s', '--search', '--search-program', default='blast', 
                       choices=['blast', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl'],
                       help='序列搜索程序 | Sequence search program')
    parser.add_argument('--mcl-inflation', type=float, default=1.2,
                       help='MCL inflation参数 | MCL inflation parameter')

    # 单拷贝基因分析参数 | Single copy gene analysis parameters
    parser.add_argument('--enable-single-copy', action='store_true', default=True,
                    help='启用单拷贝基因分析 | Enable single copy gene analysis')
    parser.add_argument('--extract-sequences', action='store_true', default=True,
                    help='提取单拷贝基因序列 | Extract single copy gene sequences')
    parser.add_argument('--single-copy-format', default='both',
                    choices=['by_orthogroup', 'by_genome', 'both'],
                    help='单拷贝基因序列输出格式 | Single copy gene sequence output format')
    
    # 稀释分析参数 | Rarefaction analysis parameters
    parser.add_argument('--enable-rarefaction', action='store_true', default=True,
                       help='启用稀释曲线分析 | Enable rarefaction curve analysis')
    parser.add_argument('--rarefaction-iterations', type=int, default=100,
                       help='稀释分析迭代次数 | Rarefaction analysis iterations')
    parser.add_argument('--disable-rarefaction', action='store_true',
                       help='禁用稀释曲线分析 | Disable rarefaction curve analysis')
    
    # 高级参数 | Advanced parameters
    parser.add_argument('-d', '--dna', action='store_true',
                       help='输入序列为DNA序列 | Input sequences are DNA')
    parser.add_argument('--basic-only', action='store_true',
                       help='仅进行基础分析 | Perform basic analysis only')
    parser.add_argument('--generate-trees', action='store_true',
                       help='生成系统发育树 | Generate phylogenetic trees')
    parser.add_argument('--msa-program', default='mafft', 
                       choices=['mafft', 'muscle'],
                       help='多序列比对程序 | Multiple sequence alignment program')
    parser.add_argument('--tree-program', default='fasttree',
                       choices=['fasttree', 'fasttree_fastest', 'raxml', 'raxml-ng', 'iqtree'], 
                       help='系统发育树构建程序 | Phylogenetic tree inference program')
    
    # 断点续跑参数 | Resume parameters
    parser.add_argument('--resume', action='store_true', default=True,
                       help='使用已有结果续跑 | Resume from existing results (default: enabled)')
    parser.add_argument('--force', action='store_true',
                       help='强制重新分析覆盖已有结果 | Force reanalysis overwriting existing results')
    parser.add_argument('--skip-orthofinder', action='store_true',
                       help='跳过OrthoFinder步骤直接进行分类 | Skip OrthoFinder step and go directly to classification')
    
    # 可视化参数 | Visualization parameters
    parser.add_argument('--no-plots', action='store_true',
                       help='不生成图表 | Do not generate plots')
    parser.add_argument('--plot-format', default='png', choices=['png', 'pdf', 'svg'],
                       help='图表格式 | Plot format')
    
    # 工具路径 | Tool paths
    parser.add_argument('--orthofinder-path', default='orthofinder',
                       help='OrthoFinder程序路径 | OrthoFinder program path')
    
    args = parser.parse_args()
    
    # 处理分析模式逻辑 | Handle analysis mode logic
    if args.generate_trees:
        basic_analysis_only = False
        generate_trees = True
    elif args.basic_only:
        basic_analysis_only = True
        generate_trees = False
    else:
        basic_analysis_only = True
        generate_trees = False
    
    # 处理断点续跑逻辑 | Handle resume logic
    resume_from_existing = args.resume and not args.force
    skip_orthofinder = args.skip_orthofinder
    force_overwrite = args.force
    
    # 处理稀释分析逻辑 | Handle rarefaction analysis logic
    enable_rarefaction = args.enable_rarefaction and not args.disable_rarefaction
    
    # 处理可视化逻辑 | Handle visualization logic
    generate_plots = not args.no_plots
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = OrthoFinderPangenomeAnalyzer(
        input_dir=args.input,
        output_dir=args.output,
        project_name=args.project_name,
        softcore_missing_threshold=args.softcore_threshold,
        dispensable_missing_threshold=args.dispensable_threshold,
        enable_single_copy_analysis=args.enable_single_copy,
        extract_sequences=args.extract_sequences,
        single_copy_output_format=args.single_copy_format,
        threads=args.threads,
        search_program=args.search,
        mcl_inflation=args.mcl_inflation,
        sequence_type='dna' if args.dna else 'protein',
        basic_analysis_only=basic_analysis_only,
        generate_trees=generate_trees,
        msa_program=args.msa_program,
        tree_program=args.tree_program,
        orthofinder_path=args.orthofinder_path,
        resume_from_existing=resume_from_existing,
        skip_orthofinder=skip_orthofinder,
        force_overwrite=force_overwrite,
        enable_rarefaction=enable_rarefaction,
        rarefaction_iterations=args.rarefaction_iterations,
        generate_plots=generate_plots,
        plot_format=args.plot_format
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()

# ===== END FILE =====