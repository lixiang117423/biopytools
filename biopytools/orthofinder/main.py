"""
OrthoFinder泛基因组分析主程序模块|OrthoFinder Pangenome Analysis Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import PangenomeConfig
from .utils import PangenomeLogger, CommandRunner, check_dependencies, validate_fasta_files, clean_fasta_sequences, count_fasta_files, format_number
from .orthofinder_runner import OrthoFinderRunner
from .pangenome_classifier import PangenomeClassifier
from .rarefaction_analysis import RarefactionAnalyzer
from .visualization import PangenomeVisualizer
from .results import ResultsProcessor
from .sequence_extractor import SequenceExtractor
from .pipeline_info import PipelineInfoGenerator


class OrthoFinderPangenomeAnalyzer:
    """OrthoFinder泛基因组分析主类|Main OrthoFinder Pangenome Analyzer Class"""

    def __init__(self, **kwargs):
        self.config = PangenomeConfig(**kwargs)
        self.config.validate()

        self.logger_manager = PangenomeLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()

        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        self.orthofinder_runner = OrthoFinderRunner(self.config, self.logger, self.cmd_runner)
        self.pangenome_classifier = PangenomeClassifier(self.config, self.logger)
        self.rarefaction_analyzer = RarefactionAnalyzer(self.config, self.logger)
        self.visualizer = PangenomeVisualizer(self.config, self.logger)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.sequence_extractor = SequenceExtractor(self.config, self.logger)
        self.pipeline_info_generator = PipelineInfoGenerator(self.logger, self.config)

    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)

    def run_analysis(self):
        """运行完整的泛基因组分析流程|Run complete pangenome analysis pipeline"""
        try:
            self.logger.info("开始OrthoFinder泛基因组分析流程|Starting OrthoFinder pangenome analysis pipeline")
            self.logger.info(f"输入目录|Input directory: {self.config.input_dir}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

            self.logger.info("生成流程信息|Generating pipeline information")
            self.pipeline_info_generator.generate_pipeline_info(self.config.output_dir)

            self.check_dependencies()

            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤1: 验证输入数据|Step 1: Validate input data")
            self.logger.info(f"{'=' * 60}")

            if not validate_fasta_files(self.config.input_dir, self.logger):
                self.logger.info("尝试清洗非法字符后重新验证|Attempting to clean invalid characters and revalidate")
                clean_fasta_sequences(self.config.input_dir, self.logger)
                if not validate_fasta_files(self.config.input_dir, self.logger):
                    raise RuntimeError("输入文件验证失败|Input file validation failed")

            genome_count, sample_names = count_fasta_files(self.config.input_dir)
            self.logger.info(f"检测到 {genome_count} 个基因组|Detected {genome_count} genomes")

            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤2: 运行OrthoFinder同源基因群分析|Step 2: Run OrthoFinder orthogroups analysis")
            self.logger.info(f"{'=' * 60}")

            if not self.orthofinder_runner.run_orthofinder():
                raise RuntimeError("OrthoFinder分析失败|OrthoFinder analysis failed")

            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤3: 加载和解析OrthoFinder结果|Step 3: Load and parse OrthoFinder results")
            self.logger.info(f"{'=' * 60}")

            orthogroups_file = self.orthofinder_runner.get_orthogroups_file()
            gene_count_file = self.orthofinder_runner.get_gene_count_file()

            self.pangenome_classifier.load_orthofinder_results(orthogroups_file, gene_count_file)

            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤4: 泛基因组分类分析|Step 4: Pangenome classification analysis")
            self.logger.info(f"{'=' * 60}")

            classification_results = self.pangenome_classifier.classify_pangenome()

            current_step = 5

            if self.config.enable_rarefaction:
                self.logger.info(f"{'=' * 60}")
                self.logger.info(f"步骤{current_step}: 泛基因组稀释曲线分析|Step {current_step}: Pangenome rarefaction curve analysis")
                self.logger.info(f"{'=' * 60}")

                main_output_file = self.pangenome_classifier.save_classification_results(
                    classification_results, self.config.output_path)
                self.rarefaction_analyzer.load_classification_results(main_output_file)

                rarefaction_results = self.rarefaction_analyzer.run_rarefaction_analysis()
                self.rarefaction_analyzer.save_rarefaction_results(rarefaction_results, self.config.output_path)
                current_step += 1
            else:
                rarefaction_results = None

            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"步骤{current_step}: 保存分析结果|Step {current_step}: Save analysis results")
            self.logger.info(f"{'=' * 60}")

            if not self.config.enable_rarefaction:
                self.logger.info("正在保存分类结果|Saving classification results")
                main_output_file = self.pangenome_classifier.save_classification_results(
                    classification_results, self.config.output_path)
                self.logger.info("分类结果保存完成|Classification results saved")

            self.logger.info("正在获取基因组统计信息|Getting genome statistics")
            genome_stats = self.pangenome_classifier.get_genome_gene_counts(classification_results)
            self.logger.info("基因组统计信息获取完成|Genome statistics retrieved")

            self.logger.info("正在获取频率分布|Getting frequency distribution")
            frequency_distributions = self.pangenome_classifier.calculate_frequency_distribution(classification_results)
            self.logger.info("频率分布获取完成|Frequency distribution retrieved")

            self.logger.info("正在保存详细结果表格|Saving detailed results table")
            detailed_table_file = self.results_processor.save_gene_families_table(
                classification_results, self.config.output_path)
            self.logger.info("详细结果表格保存完成|Detailed results table saved")

            current_step += 1

            if self.config.enable_single_copy_analysis:
                self.logger.info(f"{'=' * 60}")
                self.logger.info(f"步骤{current_step}: 单拷贝基因分析和序列提取|Step {current_step}: Single copy gene analysis and sequence extraction")
                self.logger.info(f"{'=' * 60}")

                single_copy_count = len(classification_results['single_copy'])
                self.logger.info(f"识别到 {single_copy_count} 个单拷贝同源群|Identified {single_copy_count} single copy orthogroups")

                if self.config.extract_sequences:
                    self.sequence_extractor.load_genome_sequences()
                    self.sequence_extractor.extract_single_copy_sequences(
                        classification_results['single_copy'],
                        self.config.output_path
                    )
                current_step += 1

            if self.config.generate_plots:
                self.logger.info(f"{'=' * 60}")
                self.logger.info(f"步骤{current_step}: 生成可视化图表|Step {current_step}: Generate visualization plots")
                self.logger.info(f"{'=' * 60}")

                self.visualizer.create_classification_pie_chart(classification_results, self.config.output_path)
                self.visualizer.create_frequency_distribution_plot(frequency_distributions, self.config.output_path)

                if self.config.enable_rarefaction:
                    self.visualizer.create_rarefaction_plot(rarefaction_results, self.config.output_path)

                current_step += 1

            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"步骤{current_step}: 生成分析报告|Step {current_step}: Generate analysis report")
            self.logger.info(f"{'=' * 60}")

            report_file = self.results_processor.generate_comprehensive_report(
                classification_results, genome_stats, frequency_distributions,
                rarefaction_results if self.config.enable_rarefaction else None,
                self.config.output_path)

            self.logger.info(f"{'=' * 60}")
            self.logger.info("泛基因组分析完成|Pangenome analysis completed!")
            self.logger.info(f"{'=' * 60}")

            self._copy_orthofinder_results_to_output()

            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")

            total_orthogroups = sum(len(classification_results[c]) for c in ['core', 'softcore', 'dispensable', 'private'])
            total_genes = sum(sum(len(gene_info[1]) for gene_info in classification_results[c])
                            for c in ['core', 'softcore', 'dispensable', 'private'])

            self.logger.info(f"关键统计|Key Statistics:")
            self.logger.info(f"  总基因组数|Total genomes: {genome_count}")
            self.logger.info(f"  同源基因群数|Orthogroups: {format_number(total_orthogroups)}")
            self.logger.info(f"  总基因数|Total genes: {format_number(total_genes)}")
            self.logger.info(f"  核心基因群|Core gene families: {format_number(len(classification_results['core']))}")
            self.logger.info(f"    其中单拷贝|  of which single copy: {format_number(len(classification_results['single_copy']))}")
            self.logger.info(f"  软核心基因群|Softcore gene families: {format_number(len(classification_results['softcore']))}")
            self.logger.info(f"  非必需基因群|Dispensable gene families: {format_number(len(classification_results['dispensable']))}")
            self.logger.info(f"  私有基因群|Private gene families: {format_number(len(classification_results['private']))}")

            self.logger.info(f"主要输出文件|Main output files:")
            self.logger.info(f"  基因家族信息|Gene families: {main_output_file.name}")
            self.logger.info(f"  详细数据表格|Detailed table: {detailed_table_file.name}")
            self.logger.info(f"  综合报告|Comprehensive report: {report_file.name}")

            self.pipeline_info_generator.finalize_pipeline_info(self.config.output_dir)

        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

    def _copy_orthofinder_results_to_output(self):
        """清理可能的旧目录|Clean up possible old directories"""
        import shutil

        old_dirs = [
            Path(self.config.output_dir) / "Results_OrthoFinder",
            Path(self.config.output_dir) / "OrthoFinder_Results",
            Path(self.config.output_dir) / "OrthoFinder_WorkingDirectory"
        ]

        try:
            for old_dir in old_dirs:
                if old_dir.exists():
                    self.logger.info(f"删除旧目录|Removing old directory: {old_dir}")
                    shutil.rmtree(old_dir)

        except Exception as e:
            self.logger.error(f"清理旧目录失败|Failed to clean up old directories: {e}")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='OrthoFinder泛基因组分析脚本|OrthoFinder Pangenome Analysis Script',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i protein_sequences/ -o pangenome_results
        """
    )

    parser.add_argument('-i', '--input', required=True,
                       help='输入蛋白质序列文件目录|Input protein sequence files directory')
    parser.add_argument('-o', '--output', required=True,
                       help='输出目录|Output directory')

    parser.add_argument('-n', '--project-name',
                       help='项目名称|Project name')

    parser.add_argument('--softcore-threshold', type=int, default=1,
                       help='Softcore基因缺失阈值|Softcore gene missing threshold')
    parser.add_argument('--dispensable-threshold', type=int, default=1,
                       help='Dispensable基因缺失阈值|Dispensable gene missing threshold')

    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')
    parser.add_argument('-s', '--search', '--search-program', default='blast',
                       choices=['blast', 'diamond', 'diamond_ultra_sens', 'mmseqs', 'blast_nucl'],
                       help='序列搜索程序|Sequence search program')
    parser.add_argument('--mcl-inflation', type=float, default=1.2,
                       help='MCL inflation参数|MCL inflation parameter')

    parser.add_argument('--disable-single-copy', action='store_true',
                       help='禁用单拷贝基因分析|Disable single copy gene analysis')
    parser.add_argument('--extract-sequences', action='store_true', default=True,
                       help='提取单拷贝基因序列|Extract single copy gene sequences')
    parser.add_argument('--single-copy-format', default='both',
                       choices=['by_orthogroup', 'by_genome', 'both'],
                       help='单拷贝基因序列输出格式|Single copy gene sequence output format')

    parser.add_argument('--disable-rarefaction', action='store_true',
                       help='禁用稀释曲线分析|Disable rarefaction curve analysis')
    parser.add_argument('--rarefaction-iterations', type=int, default=100,
                       help='稀释分析迭代次数|Rarefaction analysis iterations')

    parser.add_argument('-d', '--dna', action='store_true',
                       help='输入序列为DNA序列|Input sequences are DNA')
    parser.add_argument('--basic-only', action='store_true',
                       help='仅进行基础分析|Perform basic analysis only')
    parser.add_argument('--generate-trees', action='store_true',
                       help='生成系统发育树|Generate phylogenetic trees')
    parser.add_argument('--msa-program', default='mafft',
                       choices=['mafft', 'muscle'],
                       help='多序列比对程序|Multiple sequence alignment program')
    parser.add_argument('--tree-program', default='fasttree',
                       choices=['fasttree', 'fasttree_fastest', 'raxml', 'raxml-ng', 'iqtree'],
                       help='系统发育树构建程序|Phylogenetic tree inference program')

    parser.add_argument('--force', action='store_true',
                       help='强制重新分析覆盖已有结果|Force reanalysis overwriting existing results')
    parser.add_argument('--skip-orthofinder', action='store_true',
                       help='跳过OrthoFinder步骤直接进行分类|Skip OrthoFinder step and go directly to classification')

    parser.add_argument('--no-plots', action='store_true',
                       help='不生成图表|Do not generate plots')
    parser.add_argument('--plot-format', default='png', choices=['png', 'pdf', 'svg'],
                       help='图表格式|Plot format')

    parser.add_argument('--orthofinder-path',
                       default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/orthofinder_v.3.1.5/bin/orthofinder',
                       help='OrthoFinder程序路径|OrthoFinder program path')

    args = parser.parse_args()

    if args.generate_trees:
        basic_analysis_only = False
        generate_trees = True
    elif args.basic_only:
        basic_analysis_only = True
        generate_trees = False
    else:
        basic_analysis_only = True
        generate_trees = False

    analyzer = OrthoFinderPangenomeAnalyzer(
        input_dir=args.input,
        output_dir=args.output,
        project_name=args.project_name,
        softcore_missing_threshold=args.softcore_threshold,
        dispensable_missing_threshold=args.dispensable_threshold,
        enable_single_copy_analysis=not args.disable_single_copy,
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
        resume_from_existing=True,
        skip_orthofinder=args.skip_orthofinder,
        force_overwrite=args.force,
        enable_rarefaction=not args.disable_rarefaction,
        rarefaction_iterations=args.rarefaction_iterations,
        generate_plots=not args.no_plots,
        plot_format=args.plot_format
    )

    analyzer.run_analysis()


if __name__ == "__main__":
    main()
