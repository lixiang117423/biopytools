"""
GCTB主程序模块|GCTB Main Program Module
"""

import argparse
import sys
import os
from pathlib import Path
from .config import GCTBConfig
from .utils import GCTBLogger, CommandRunner
from .data_processing import DataProcessor
from .analyzer import GCTBAnalyzer


class GCTBRunner:
    """GCTB运行器|GCTB Runner"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = GCTBConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = GCTBLogger(self.config.logs_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_dir)

        # 初始化数据处理器|Initialize data processor
        self.data_processor = DataProcessor(
            self.config,
            self.logger,
            self.cmd_runner
        )

        # 初始化分析器|Initialize analyzer
        self.analyzer = GCTBAnalyzer(
            self.config,
            self.logger,
            self.cmd_runner
        )

    def _is_step_completed(self, output_file: str) -> bool:
        """检查步骤是否已完成|Check if step is completed"""
        return os.path.exists(output_file)

    def step_convert_vcf(self) -> bool:
        """步骤1: VCF转换为PLINK格式|Step 1: Convert VCF to PLINK format"""
        output_file = f"{self.data_processor.plink_prefix}.bed"

        if self._is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: VCF转换|VCF conversion")
            return True

        self.logger.info("执行步骤1: VCF转换为PLINK格式|Executing step 1: Convert VCF to PLINK")
        return self.data_processor.vcf_to_plink()

    def step_quality_control(self) -> bool:
        """步骤2: 质量控制|Step 2: Quality control"""
        output_file = f"{self.data_processor.plink_qc_prefix}.bed"

        if self._is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 质量控制|Quality control")
            return True

        self.logger.info("执行步骤2: 质量控制|Executing step 2: Quality control")
        return self.data_processor.quality_control()

    def step_calculate_freq(self) -> bool:
        """步骤3: 计算等位基因频率|Step 3: Calculate allele frequencies"""
        output_file = f"{self.data_processor.plink_qc_prefix}.frq"

        if self._is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 计算频率|Calculate frequencies")
            return True

        self.logger.info("执行步骤3: 计算等位基因频率|Executing step 3: Calculate allele frequencies")
        return self.data_processor.calculate_freq()

    def step_build_ld_matrix(self) -> bool:
        """步骤4: 构建LD矩阵|Step 4: Build LD matrix"""
        ld_prefix = str(self.config.ld_dir / "ld_matrix")

        if self.config.ld_matrix_type == "sparse":
            output_file = f"{ld_prefix}.ldm.sparse.info"
            build_func = self.data_processor.build_sparse_ld_matrix
        elif self.config.ld_matrix_type == "block":
            output_file = f"{ld_prefix}.ldm.block.info"
            build_func = self.data_processor.build_block_ld_matrix
        elif self.config.ld_matrix_type == "eigen":
            output_file = f"{ld_prefix}.ldm.eigen.info"
            build_func = self.data_processor.eigen_decompose_ld_matrix
        else:
            self.logger.error(f"未知的LD矩阵类型|Unknown LD matrix type: {self.config.ld_matrix_type}")
            return False

        if self._is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: 构建LD矩阵|Build LD matrix")
            return True

        self.logger.info(f"执行步骤4: 构建{self.config.ld_matrix_type}LD矩阵|Executing step 4: Build {self.config.ld_matrix_type} LD matrix")
        return build_func()

    def step_run_analysis(self) -> bool:
        """步骤5: 运行GCTB分析|Step 5: Run GCTB analysis"""
        if self.config.analysis_mode == "individual":
            output_file = f"{self.config.analysis_dir}/bayes{self.config.bayes_type.lower()}.snpRes"
            run_func = self.analyzer.run_individual_analysis
        else:
            output_file = f"{self.config.analysis_dir}/sbayes{self.config.bayes_type.lower()}.snpRes"
            run_func = self.analyzer.run_summary_analysis

        if self._is_step_completed(output_file):
            self.logger.info(f"跳过已完成步骤|Skipping completed step: GCTB分析|GCTB analysis")
            return True

        self.logger.info(f"执行步骤5: 运行{self.config.analysis_mode}分析|Executing step 5: Run {self.config.analysis_mode} analysis")
        return run_func()

    def run_single_step(self, step: str) -> bool:
        """
        运行单个步骤|Run single step

        Args:
            step: 步骤名称|Step name (convert/qc/freq/ld/analysis)

        Returns:
            bool: 执行是否成功|Whether execution succeeded
        """
        step_functions = {
            'convert': (self.step_convert_vcf, "VCF转换|VCF conversion"),
            'qc': (self.step_quality_control, "质量控制|Quality control"),
            'freq': (self.step_calculate_freq, "计算频率|Calculate frequencies"),
            'ld': (self.step_build_ld_matrix, "构建LD矩阵|Build LD matrix"),
            'analysis': (self.step_run_analysis, "运行分析|Run analysis")
        }

        if step not in step_functions:
            self.logger.error(f"无效的步骤名称|Invalid step name: {step}")
            return False

        step_func, step_name = step_functions[step]
        self.logger.info(f"执行步骤|Executing step: {step_name}")

        success = step_func()
        if success:
            self.logger.info(f"步骤完成|Step completed: {step_name}")
        else:
            self.logger.error(f"步骤失败|Step failed: {step_name}")

        return success

    def run_batch_analysis(self) -> bool:
        """
        运行批量表型分析|Run batch phenotype analysis

        自动转换表型文件并为每个表型运行独立的分析
        Automatically convert phenotype file and run independent analysis for each phenotype

        Returns:
            bool: 执行是否成功|Whether execution succeeded
        """
        self.logger.info("="*60)
        self.logger.info("开始批量GCTB分析|Starting batch GCTB analysis")
        self.logger.info("="*60)

        # 步骤0: 转换表型文件|Step 0: Convert phenotype file
        self.logger.info("步骤 0: 转换表型文件|Step 0: Convert phenotype file")
        self.logger.info("-"*60)

        converted_pheno_files = self.data_processor.convert_and_split_phenotype()

        if not converted_pheno_files:
            self.logger.error("表型文件转换失败|Phenotype file conversion failed")
            return False

        self.logger.info(f"✓ 生成了 {len(converted_pheno_files)} 个表型文件|Generated {len(converted_pheno_files)} phenotype files")
        self.logger.info("")

        # 步骤1-3: VCF转换、质控、计算频率（只执行一次）|Steps 1-3: VCF conversion, QC, calculate frequencies (run once)
        prep_steps = [
            (self.step_convert_vcf, "VCF转换|VCF conversion"),
            (self.step_quality_control, "质量控制|Quality control"),
            (self.step_calculate_freq, "计算频率|Calculate frequencies"),
        ]

        for i, (step_func, step_name) in enumerate(prep_steps, 1):
            self.logger.info(f"步骤 {i}/{len(prep_steps)}: {step_name}|Step {i}/{len(prep_steps)}: {step_name}")
            self.logger.info("-"*60)

            if not step_func():
                self.logger.error(f"\n步骤 {i} 失败，终止流程|Step {i} failed, terminating pipeline")
                return False

            self.logger.info(f"✓ 步骤 {i} 完成|Step {i} completed")
            self.logger.info("")

        # 为每个表型运行分析|Run analysis for each phenotype
        all_success = True
        summary_results = []

        # 创建共享的分析结果目录|Create shared analysis results directory
        shared_analysis_dir = self.config.output_path / "04_gctb_analysis"
        shared_analysis_dir.mkdir(parents=True, exist_ok=True)

        for idx, pheno_file in enumerate(converted_pheno_files, 1):
            pheno_name = Path(pheno_file).stem
            self.logger.info("="*60)
            self.logger.info(f"处理表型 {idx}/{len(converted_pheno_files)}: {pheno_name}|Processing phenotype {idx}/{len(converted_pheno_files)}: {pheno_name}")
            self.logger.info("="*60)

            # 创建单个表型的结果子目录|Create single phenotype result subdirectory
            trait_result_dir = shared_analysis_dir / pheno_name
            trait_result_dir.mkdir(parents=True, exist_ok=True)

            # 更新配置|Update configuration
            self.config.pheno_file = pheno_file
            self.config.analysis_dir = trait_result_dir
            self.config.analysis_dir.mkdir(parents=True, exist_ok=True)

            # 运行分析|Run analysis
            success = self.step_run_analysis()

            if success:
                self.logger.info(f"✓ 表型 {pheno_name} 分析完成|Phenotype {pheno_name} analysis completed")
                summary_results.append({
                    'phenotype': pheno_name,
                    'status': 'success',
                    'output_dir': str(trait_result_dir)
                })
            else:
                self.logger.error(f"✗ 表型 {pheno_name} 分析失败|Phenotype {pheno_name} analysis failed")
                summary_results.append({
                    'phenotype': pheno_name,
                    'status': 'failed',
                    'output_dir': str(trait_result_dir)
                })
                all_success = False

            self.logger.info("")

        # 打印总结|Print summary
        self.logger.info("="*60)
        self.logger.info("批量分析完成|Batch analysis completed")
        self.logger.info("="*60)
        self.logger.info(f"总计|Total: {len(converted_pheno_files)} 个表型|phenotypes")
        self.logger.info(f"成功|Success: {sum(1 for r in summary_results if r['status'] == 'success')}")
        self.logger.info(f"失败|Failed: {sum(1 for r in summary_results if r['status'] == 'failed')}")
        self.logger.info("")

        for result in summary_results:
            status_symbol = "✓" if result['status'] == 'success' else "✗"
            self.logger.info(f"{status_symbol} {result['phenotype']}: {result['output_dir']}")

        # 步骤5: 汇总所有表型的结果|Step 5: Summarize results for all phenotypes
        if all_success:
            self.logger.info("")
            self.logger.info("="*60)
            self.logger.info("步骤 5: 汇总分析结果|Step 5: Summarize analysis results")
            self.logger.info("="*60)

            try:
                from .results_summary import ResultsSummarizer

                summarizer = ResultsSummarizer(self.logger, self.config.output_dir)

                # 收集成功的表型及其分析目录|Collect successful traits and their analysis directories
                trait_names = [r['phenotype'] for r in summary_results if r['status'] == 'success']
                analysis_dirs = {r['phenotype']: Path(r['output_dir']) for r in summary_results if r['status'] == 'success'}

                # 生成汇总结果|Generate summary results
                summarizer.summarize_batch_results(
                    trait_names=trait_names,
                    analysis_dirs=analysis_dirs,
                    bayes_type=self.config.bayes_type
                )

            except Exception as e:
                self.logger.warning(f"结果汇总失败，但分析已完成|Results summarization failed, but analysis completed: {e}")

        return all_success

    def run_full_pipeline(self) -> bool:
        """
        运行完整流程|Run full pipeline

        Returns:
            bool: 执行是否成功|Whether execution succeeded
        """
        self.logger.info("="*60)
        self.logger.info("开始GCTB分析流程|Starting GCTB analysis pipeline")
        self.logger.info("="*60)

        steps = [
            (self.step_convert_vcf, "VCF转换|VCF conversion"),
            (self.step_quality_control, "质量控制|Quality control"),
            (self.step_calculate_freq, "计算频率|Calculate frequencies"),
        ]

        # 如果是汇总水平分析，需要构建LD矩阵|If summary-level analysis, need to build LD matrix
        if self.config.analysis_mode == "summary":
            steps.append((self.step_build_ld_matrix, f"构建{self.config.ld_matrix_type}LD矩阵|Build {self.config.ld_matrix_type} LD matrix"))

        steps.append((self.step_run_analysis, f"运行{self.config.analysis_mode}分析|Run {self.config.analysis_mode} analysis"))

        for i, (step_func, step_name) in enumerate(steps, 1):
            self.logger.info(f"步骤 {i}/{len(steps)}: {step_name}|Step {i}/{len(steps)}: {step_name}")
            self.logger.info("-"*60)

            if not step_func():
                self.logger.error(f"\n步骤 {i} 失败，终止流程|Step {i} failed, terminating pipeline")
                return False

            self.logger.info(f"✓ 步骤 {i} 完成|Step {i} completed")

        self.logger.info("="*60)
        self.logger.info("GCTB分析流程全部完成|GCTB analysis pipeline completed successfully!")
        self.logger.info("="*60)

        # 输出结果文件位置|Output result file locations
        self._print_output_files()

        return True

    def _print_output_files(self):
        """打印输出文件位置|Print output file locations"""
        self.logger.info("\n输出文件|Output files:")

        if self.config.analysis_mode == "individual":
            prefix = f"{self.config.analysis_dir}/bayes{self.config.bayes_type.lower()}"
        else:
            prefix = f"{self.config.analysis_dir}/sbayes{self.config.bayes_type.lower()}"

        self.logger.info(f"  SNP结果|SNP results: {prefix}.snpRes")
        self.logger.info(f"  参数结果|Parameter results: {prefix}.parRes")
        self.logger.info(f"  日志文件|Log file: {prefix}.log")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='GCTB全基因组复杂性状贝叶斯分析自动化工具|GCTB Genome-wide Complex Trait Bayesian Analysis Automation Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--vcf-file', required=True,
                       help='VCF变异文件路径|VCF variant file path')
    parser.add_argument('-p', '--pheno-file', required=True,
                       help='表型文件路径（个体水平分析）或GWAS汇总统计文件（汇总水平分析）|'
                            'Phenotype file path (individual-level) or GWAS summary statistics file (summary-level)')

    # 输出配置|Output configuration
    parser.add_argument('-o', '--output-dir', default='./gctb_output',
                       help='输出目录|Output directory')

    # 软件路径配置|Software path configuration
    parser.add_argument('--gctb-path',
                       default='~/miniforge3/envs/gctb/bin/gctb',
                       help='GCTB软件路径|GCTB software path')
    parser.add_argument('--plink-path',
                       default='~/miniforge3/envs/Population_genetics/bin/plink',
                       help='PLINK软件路径|PLINK software path')

    # 质量控制参数|Quality control parameters
    parser.add_argument('--maf-threshold', type=float, default=0.01,
                       help='MAF阈值|MAF threshold')
    parser.add_argument('--miss-threshold', type=float, default=0.1,
                       help='缺失率阈值|Missing rate threshold')
    parser.add_argument('--hwe-p', type=float, default=1e-6,
                       help='HWE p值阈值|HWE p-value threshold')

    # 分析参数|Analysis parameters
    parser.add_argument('--bayes-type', choices=['S', 'R', 'C'], default='S',
                       help='贝叶斯模型类型|Bayesian model type (S/R/C)')
    parser.add_argument('--analysis-mode', choices=['individual', 'summary'], default='individual',
                       help='分析模式|Analysis mode')
    parser.add_argument('--ld-matrix-type', choices=['sparse', 'block', 'eigen'], default='sparse',
                       help='LD矩阵类型|LD matrix type (sparse/block/eigen)')

    # GCTB 高级参数|GCTB advanced parameters
    parser.add_argument('--threads', type=int, default=12,
                       help='线程数|Number of threads')
    parser.add_argument('--seed', type=int,
                       help='随机种子|Random seed')
    parser.add_argument('--pi', type=float,
                       help='polygenicity参数|Polygenicity parameter')
    parser.add_argument('--sigma-g', type=float,
                       help='遗传方差|Genetic variance')
    parser.add_argument('--rho',
                       help='SNP效应与MAF关系参数|Parameter for effect-MAF relationship')

    # 步骤控制|Step control
    parser.add_argument('--step', choices=['convert', 'qc', 'freq', 'ld', 'analysis'],
                       help='只运行指定步骤|Run only specified step')
    parser.add_argument('--batch', action='store_true', default=True,
                       help='批量处理多个表型（默认开启，使用--no-batch禁用）|Batch process multiple phenotypes (enabled by default, use --no-batch to disable)')
    parser.add_argument('--no-batch', dest='batch', action='store_false',
                       help='禁用批量处理模式|Disable batch processing mode')
    parser.add_argument('--keep-intermediate', action='store_true',
                       help='保留中间文件|Keep intermediate files')

    args = parser.parse_args()

    try:
        # 创建运行器并运行|Create runner and run
        runner = GCTBRunner(
            vcf_file=args.vcf_file,
            pheno_file=args.pheno_file,
            output_dir=args.output_dir,
            gctb_path=args.gctb_path,
            plink_path=args.plink_path,
            maf_threshold=args.maf_threshold,
            miss_threshold=args.miss_threshold,
            hwe_p=args.hwe_p,
            bayes_type=args.bayes_type,
            analysis_mode=args.analysis_mode,
            ld_matrix_type=args.ld_matrix_type,
            threads=args.threads,
            seed=args.seed,
            pi=args.pi,
            sigma_g=args.sigma_g,
            rho=args.rho,
            step=args.step,
            batch=args.batch,
            keep_intermediate=args.keep_intermediate
        )

        if args.step:
            success = runner.run_single_step(args.step)
        elif args.batch:
            success = runner.run_batch_analysis()
        else:
            success = runner.run_full_pipeline()

        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"错误|Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
