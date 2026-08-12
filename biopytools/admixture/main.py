"""
ADMIXTURE分析主程序模块|ADMIXTURE Analysis Main Module
"""

import argparse
import sys
import os
import time
import subprocess
from pathlib import Path
from .config import AdmixtureConfig
from .utils import AdmixtureLogger, CommandRunner, SoftwareChecker, build_conda_command
from .data_processing import VCFProcessor, PlinkProcessor
from .analysis import AdmixtureAnalyzer as CoreAnalyzer, ResultsProcessor
from .results import CovariateGenerator, PlotGenerator, SummaryGenerator

# 版本信息|Version information
VERSION = "1.0.0"


class AdmixtureAnalyzer:
    """ADMIXTURE分析主类|Main ADMIXTURE Analyzer Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = AdmixtureConfig(**kwargs)
        self.config.validate()

        # 创建分层输出目录(§12.2)|Create layered output directories
        self._create_output_dirs()

        # 初始化日志(落99_logs/)|Initialize logging (into 99_logs/)
        self.logger_manager = AdmixtureLogger(
            Path(self.config.logs_dir),
            log_name="admixture_analysis.log",
            log_level=self.config.log_level,
            quiet=self.config.quiet
        )
        self.logger = self.logger_manager.get_logger()

        # 检查软件环境(非阻断,仅warning)|Check software env (non-blocking, warning only)
        self.software_checker = SoftwareChecker(self.logger, self.config)
        self.software_checker.check_dependencies()

        # 初始化命令执行器(透传dry_run)|Initialize command runner (forward dry_run)
        self.cmd_runner = CommandRunner(
            self.logger, self.config.output_path, dry_run=self.config.dry_run)

        # 初始化各个处理器|Initialize processors
        self.vcf_processor = VCFProcessor(self.config, self.logger, self.cmd_runner)
        self.plink_processor = PlinkProcessor(self.config, self.logger, self.cmd_runner)
        self.admixture_analyzer = CoreAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.covariate_generator = CovariateGenerator(self.config, self.logger)
        self.plot_generator = PlotGenerator(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)

    def _create_output_dirs(self):
        """创建所有分层输出目录(§12.2)|Create all layered output directories"""
        for d in (self.config.pipeline_info_dir, self.config.preprocessing_dir,
                  self.config.plink_dir, self.config.admixture_dir,
                  self.config.results_dir, self.config.logs_dir):
            os.makedirs(d, exist_ok=True)

    def _step(self, num, desc: str):
        """打印步骤横幅|Print step banner"""
        self.logger.info("=" * 60)
        self.logger.info(f"STEP {num}: {desc}")
        self.logger.info("=" * 60)

    def _generate_software_versions(self):
        """生成software_versions.yml到00_pipeline_info(§12.5)|Generate software_versions.yml"""
        import yaml
        tools = {
            'plink': self.config.plink_path,
            'bcftools': self.config.bcftools_path,
        }
        if self.config.method == 'adamixture':
            tools['adamixture'] = self.config.adamixture_path
        else:
            tools['admixture'] = self.config.admixture_path

        versions = {}
        for name, path in tools.items():
            try:
                cmd = build_conda_command(path, ['--version'])
                self.logger.info(f"命令|Command: {' '.join(cmd)}")
                r = subprocess.run(cmd, capture_output=True, text=True, timeout=30)
                ver = (r.stdout or r.stderr or '').strip().split('\n')[0]
                versions[name] = {'version': ver or 'unknown', 'path': path}
            except Exception as e:
                self.logger.warning(f"获取{name}版本失败|Failed to get {name} version: {e}")
                versions[name] = {'version': 'unknown', 'path': path}

        info = {
            'pipeline': {'name': 'biopytools admixture', 'version': VERSION},
            'tools': versions,
            'parameters': {
                'method': self.config.method,
                'min_k': self.config.min_k, 'max_k': self.config.max_k,
                'cv_folds': self.config.cv_folds, 'threads': self.config.threads,
                'maf': self.config.maf, 'missing_rate': self.config.missing_rate,
                'hwe_pvalue': self.config.hwe_pvalue,
                'ld_prune': self.config.ld_prune, 'ld_window': self.config.ld_window,
                'ld_step': self.config.ld_step, 'ld_r2': self.config.ld_r2,
                'skip_preprocessing': self.config.skip_preprocessing,
            },
        }

        out_file = os.path.join(self.config.pipeline_info_dir, 'software_versions.yml')
        try:
            with open(out_file, 'w', encoding='utf-8') as f:
                yaml.dump(info, f, default_flow_style=False, allow_unicode=True)
            self.logger.info(f"软件版本信息已保存|Software versions saved: {out_file}")
        except Exception as e:
            self.logger.warning(f"保存软件版本信息失败|Failed to save software versions: {e}")

    def _dry_run_preview(self):
        """模拟运行:打印计划步骤与代表性命令,不执行|Dry run: print planned steps and commands, no execution"""
        self.logger.info("=" * 60)
        self.logger.info("模拟运行模式|DRY RUN MODE - 以下命令不会实际执行|commands below will NOT be executed")
        self.logger.info("=" * 60)

        vcf = self.config.vcf_file
        biallelic = os.path.join(self.config.preprocessing_dir, "biallelic.vcf.gz")
        raw = os.path.join(self.config.plink_dir, "raw_data")
        qc_prefix = os.path.join(self.config.plink_dir, self.config.base_name)

        # bcftools 双等位过滤|bcftools biallelic filter
        cmd_bcf = build_conda_command(self.config.bcftools_path,
            ['view', '-m2', '-M2', '-v', 'snps', vcf, '-Oz', '-o', biallelic])
        self.logger.info(f"[STEP 1 预览|preview] 命令|Command: {' '.join(cmd_bcf)}")

        # PLINK VCF转换|PLINK VCF->BED
        cmd_plink = build_conda_command(self.config.plink_path,
            ['--vcf', vcf, '--make-bed', '--out', raw, '--allow-extra-chr',
             '--double-id', '--threads', str(self.config.threads)])
        self.logger.info(f"[STEP 2 预览|preview] 命令|Command: {' '.join(cmd_plink)}")

        # PLINK 质控|PLINK QC
        if not self.config.skip_preprocessing:
            cmd_qc = build_conda_command(self.config.plink_path,
                ['--bfile', raw, '--maf', str(self.config.maf), '--hwe', str(self.config.hwe_pvalue),
                 '--geno', str(self.config.missing_rate), '--mind', str(self.config.missing_rate),
                 '--make-bed', '--out', qc_prefix, '--allow-extra-chr'])
            self.logger.info(f"[STEP 3 预览|preview] 命令|Command: {' '.join(cmd_qc)}")

        # LD剪枝|LD pruning
        if self.config.ld_prune:
            self.logger.info(f"[STEP LD 预览|preview] plink --indep-pairwise {self.config.ld_window} {self.config.ld_step} {self.config.ld_r2} (基于质控后数据|on QC'd data)")

        # ADMIXTURE 各K|ADMIXTURE per K
        bed = os.path.join(self.config.plink_dir, "admixture_chr_fixed.bed")
        for k in range(self.config.min_k, self.config.max_k + 1):
            if self.config.method == 'adamixture':
                cmd_k = build_conda_command(self.config.adamixture_path,
                    ['--k', str(k), '--data_path', bed, '--save_dir', self.config.admixture_dir])
            else:
                cmd_k = build_conda_command(self.config.admixture_path,
                    ['--cv', str(self.config.cv_folds), '-j' + str(self.config.threads), bed, str(k)])
            self.logger.info(f"[STEP 4 预览|preview] K={k} 命令|Command: {' '.join(cmd_k)}")

        self.logger.info("模拟运行结束(实际未执行任何命令)|Dry run finished (no commands actually executed)")

    def run_analysis(self):
        """运行完整的ADMIXTURE分析流程|Run complete ADMIXTURE analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始ADMIXTURE群体结构分析|Starting ADMIXTURE Population Structure Analysis")
            self.logger.info("=" * 80)
            self.logger.info(f"输入VCF文件|Input VCF file: {self.config.vcf_file}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"方法|Method: {self.config.method} | K值范围|K range: {self.config.min_k}-{self.config.max_k}")
            self.logger.info(f"线程数|Threads: {self.config.threads} | LD剪枝|LD pruning: {self.config.ld_prune}")
            self.logger.info("=" * 80)

            # 写软件版本信息(§12.5)|Write software versions up front
            self._generate_software_versions()

            # dry_run:打印计划步骤与命令即返回,不执行链(避免Python步骤读不存在的中间文件)
            # |Dry run: print planned steps/commands and return without executing the chain
            if self.config.dry_run:
                self._dry_run_preview()
                return True

            # STEP 1: VCF预处理|Step 1: VCF preprocessing
            self._step(1, "VCF文件预处理|VCF file preprocessing")
            processed_vcf = self.vcf_processor.preprocess_vcf()
            if processed_vcf is None:
                processed_vcf = self.config.vcf_file

            # STEP 2: 转换为PLINK格式|Step 2: Convert to PLINK format
            self._step(2, "转换为PLINK格式|Convert to PLINK format")
            raw_prefix = self.plink_processor.convert_vcf_to_plink(processed_vcf)

            # STEP 3: 质量控制(可跳过)|Step 3: Quality control (skippable)
            if self.config.skip_preprocessing:
                self.logger.info("跳过质量控制(--skip-preprocessing)|Skipping QC (--skip-preprocessing)")
                qc_prefix = raw_prefix
            else:
                self._step(3, "质量控制|Quality control")
                qc_prefix = self.plink_processor.quality_control(raw_prefix)

            # STEP: LD剪枝|LD pruning
            self._step("LD", "LD剪枝|LD pruning")
            pruned_prefix = self.plink_processor.ld_prune(qc_prefix)

            # STEP: 修复染色体编号|Fix chromosome codes
            self._step("CHR", "修复染色体编号|Fix chromosome codes")
            fixed_prefix = self.plink_processor.fix_chromosome_codes(pruned_prefix)

            # 更新base_name以匹配ADMIXTURE实际输入文件名|Update base_name to match actual input
            self.config.base_name = os.path.basename(fixed_prefix)
            self.logger.info(f"更新分析文件基础名称|Update base name: {self.config.base_name}")

            # STEP 4: ADMIXTURE分析|Step 4: ADMIXTURE analysis
            self._step(4, "ADMIXTURE分析|ADMIXTURE analysis")
            best_k = self.admixture_analyzer.run_admixture_analysis(fixed_prefix)

            if best_k is None:
                self.logger.warning("未确定最优K值,跳过结果处理|Best K undetermined, skipping results")
                return False

            # STEP 5: 结果处理|Step 5: Results processing
            self._step(5, "结果处理|Results processing")
            q_data, p_data, stats = self.results_processor.process_results(best_k)

            # STEP 6: GWAS协变量|Step 6: GWAS covariates
            self._step(6, "生成GWAS协变量|Generate GWAS covariates")
            self.covariate_generator.generate_gwas_covariates(best_k)

            # STEP: 结果可视化|Visualization
            self._step("PLOT", "结果可视化|Visualization")
            self.plot_generator.generate_plots(q_data, best_k)

            # STEP 7: 总结报告|Step 7: Summary report
            self._step(7, "生成总结报告|Generate summary report")
            self.summary_generator.generate_summary(best_k, stats)

            # 完成信息|Completion information
            self.logger.info("=" * 80)
            self.logger.info("ADMIXTURE分析圆满完成|ADMIXTURE analysis completed successfully!")
            self.logger.info("=" * 80)
            self.logger.info(f"最优K值|Best K value: {best_k}")
            self.logger.info(f"结果目录|Results directory: {self.config.output_dir}")
            self.logger.info("主要输出文件|Main output files:")
            self.logger.info("  - 04_results/admixture_proportions.csv: 个体祖先成分|Individual ancestry proportions")
            self.logger.info("  - 04_results/gwas_covariates.txt: GWAS协变量文件|GWAS covariate file")
            self.logger.info("  - 03_admixture/cv_results.csv: 交叉验证结果|Cross-validation results")
            self.logger.info("  - 04_results/*.pdf: 可视化图表|Visualization plots")
            self.logger.info("  - 04_results/analysis_summary.txt: 分析总结报告|Analysis summary report")
            return True

        except Exception as e:
            self.logger.error(f"分析过程中发生严重错误|A critical error occurred during analysis: {e}", exc_info=True)
            sys.exit(1)


def main():
    """主函数入口|Main function entry point"""
    start_time = time.time()

    parser = argparse.ArgumentParser(
        description="ADMIXTURE群体结构分析工具(模块化版本)|ADMIXTURE Population Structure Analysis Tool (Modular Version)",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i input.vcf -o results
        '''
    )

    # 必需参数|Required arguments
    parser.add_argument("-i", "--vcf", required=True,
                       help="输入VCF文件路径|Input VCF file path")

    # 可选参数|Optional arguments
    parser.add_argument("-o", "--output", default="admixture_results",
                       help="输出目录|Output directory")
    parser.add_argument("--method", default="admixture", choices=["admixture", "adamixture"],
                       help="分析方法|Analysis method (admixture or adamixture)")
    parser.add_argument("-k", "--min-k", type=int, default=2,
                       help="最小K值|Minimum K value")
    parser.add_argument("-K", "--max-k", type=int, default=10,
                       help="最大K值|Maximum K value")
    parser.add_argument("-c", "--cv-folds", type=int, default=5,
                       help="交叉验证折数|Cross-validation folds")
    parser.add_argument("-t", "--threads", type=int, default=12,
                       help="线程数|Number of threads")

    # ADAMIXTURE 参数|ADAMIXTURE parameters
    adamixture_group = parser.add_argument_group('ADAMIXTURE参数|ADAMIXTURE parameters')
    adamixture_group.add_argument("--adamixture-path",
                                  default="~/miniforge3/envs/adamixture_v.1.0.2/bin/adamixture",
                                  help="ADAMIXTURE可执行文件路径|ADAMIXTURE executable path")
    adamixture_group.add_argument("--adamixture-lr", type=float, default=0.005,
                                  help="ADAMIXTURE学习率|ADAMIXTURE learning rate")
    adamixture_group.add_argument("--adamixture-beta1", type=float, default=0.80,
                                  help="ADAMIXTURE beta1参数|ADAMIXTURE beta1 parameter")
    adamixture_group.add_argument("--adamixture-beta2", type=float, default=0.88,
                                  help="ADAMIXTURE beta2参数|ADAMIXTURE beta2 parameter")
    adamixture_group.add_argument("--adamixture-max-iter", type=int, default=1500,
                                  help="ADAMIXTURE最大迭代次数|ADAMIXTURE maximum iterations")
    adamixture_group.add_argument("--adamixture-seed", type=int, default=42,
                                  help="ADAMIXTURE随机种子|ADAMIXTURE random seed")

    # 质控参数|Quality control parameters
    parser.add_argument("-m", "--maf", type=float, default=0.05,
                       help="MAF阈值|MAF threshold")
    parser.add_argument("-M", "--missing", type=float, default=0.1,
                       help="缺失率阈值|Missing rate threshold")
    parser.add_argument("-H", "--hwe", type=float, default=1e-6,
                       help="HWE p值阈值|HWE p-value threshold")

    # LD剪枝参数|LD pruning parameters
    parser.add_argument("--no-ld-prune", action="store_true",
                       help="关闭LD剪枝(默认开启)|Disable LD pruning (on by default)")
    parser.add_argument("--ld-window", default="3000kb",
                       help="LD剪枝窗口(kb或SNP数)|LD pruning window (kb or SNP count)")
    parser.add_argument("--ld-step", type=int, default=1,
                       help="LD剪枝步长|LD pruning step size")
    parser.add_argument("--ld-r2", type=float, default=0.2,
                       help="LD剪枝r2阈值|LD pruning r2 threshold")

    # 处理选项|Processing options
    parser.add_argument("-s", "--skip-preprocessing", action="store_true",
                       help="跳过VCF预处理和质控|Skip VCF preprocessing and QC")
    parser.add_argument("--keep-intermediate", action="store_true",
                       help="保留中间文件|Keep intermediate files")

    # 日志参数|Logging parameters
    log_group = parser.add_argument_group('日志选项|Logging options')
    log_group.add_argument("--verbose", action="count", default=0,
                          help="详细输出模式(-v: INFO, -vv: DEBUG)|Verbose mode (-v: INFO, -vv: DEBUG)")
    log_group.add_argument("--quiet", action="store_true",
                          help="静默模式(只输出ERROR)|Quiet mode (ERROR only)")
    log_group.add_argument("--log-level",
                          help="日志级别(DEBUG/INFO/WARNING/ERROR/CRITICAL)|Log level (default: INFO)")
    log_group.add_argument("--log-file",
                          help="日志文件路径|Log file path")

    # 执行控制|Execution control
    exec_group = parser.add_argument_group('执行选项|Execution options')
    exec_group.add_argument("-f", "--force", action="store_true",
                           help="强制覆盖已存在文件|Force overwrite existing files")
    exec_group.add_argument("--dry-run", action="store_true",
                           help="模拟运行(不实际执行)|Dry run without execution")

    # 版本信息|Version information
    parser.add_argument("-V", "--version", action="version",
                       version=f'%(prog)s {VERSION}')

    args = parser.parse_args()

    # 确定日志级别|Determine log level
    if args.log_level:
        log_level = args.log_level
    elif args.verbose >= 2:
        log_level = "DEBUG"
    elif args.verbose == 1:
        log_level = "INFO"
    elif args.quiet:
        log_level = "ERROR"
    else:
        log_level = "INFO"

    # 创建分析器并运行|Create analyzer and run
    analyzer = AdmixtureAnalyzer(
        vcf_file=args.vcf,
        output_dir=args.output,
        method=args.method,
        min_k=args.min_k,
        max_k=args.max_k,
        cv_folds=args.cv_folds,
        threads=args.threads,
        adamixture_path=args.adamixture_path,
        adamixture_lr=args.adamixture_lr,
        adamixture_beta1=args.adamixture_beta1,
        adamixture_beta2=args.adamixture_beta2,
        adamixture_max_iter=args.adamixture_max_iter,
        adamixture_seed=args.adamixture_seed,
        maf=args.maf,
        missing_rate=args.missing,
        hwe_pvalue=args.hwe,
        ld_prune=not args.no_ld_prune,
        ld_window=args.ld_window,
        ld_step=args.ld_step,
        ld_r2=args.ld_r2,
        skip_preprocessing=args.skip_preprocessing,
        keep_intermediate=args.keep_intermediate,
        log_level=log_level,
        quiet=args.quiet,
        verbose=args.verbose,
        force=args.force,
        dry_run=args.dry_run
    )

    try:
        # 输出程序信息|Output program information
        analyzer.logger.info("=" * 80)
        analyzer.logger.info("程序|Program: ADMIXTURE Population Structure Analysis")
        analyzer.logger.info(f"版本|Version: {VERSION}")
        analyzer.logger.info("=" * 80)

        if args.dry_run:
            analyzer.logger.info("模拟运行模式-不会实际执行命令|DRY RUN mode - commands will not be executed")

        # 执行分析|Run analysis
        analyzer.run_analysis()

        # 输出总结信息|Output summary
        elapsed_time = time.time() - start_time
        analyzer.logger.info("=" * 80)
        analyzer.logger.info("流程总结|Pipeline Summary")
        analyzer.logger.info("=" * 80)
        analyzer.logger.info(f"总运行时间|Total runtime: {elapsed_time:.2f} seconds")
        analyzer.logger.info(f"输出目录|Output directory: {args.output}")
        analyzer.logger.info("流程完成|Pipeline completed successfully")

    except KeyboardInterrupt:
        analyzer.logger.warning("用户中断程序执行|Process interrupted by user")
        sys.exit(130)
    except Exception as e:
        analyzer.logger.critical(f"Pipeline failed: {str(e)}", exc_info=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
