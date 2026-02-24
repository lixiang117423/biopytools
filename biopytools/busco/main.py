"""
BUSCO质量评估分析主程序模块|BUSCO Quality Assessment Analysis Main Module
"""

import argparse
import sys
from typing import Tuple
from .config import BUSCOConfig
from .utils import BUSCOLogger, CommandRunner, FileManager, check_dependencies
from .analyzer import BUSCORunner
from .results import ResultsProcessor, SummaryGenerator

class BUSCOAnalyzer:
    """BUSCO质量评估分析主类|Main BUSCO Quality Assessment Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = BUSCOConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志|Initialize logging
        self.logger_manager = BUSCOLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器|Initialize processors
        self.file_manager = FileManager(self.config, self.logger)
        self.busco_runner = BUSCORunner(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件|Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的BUSCO分析流程|Run complete BUSCO analysis pipeline"""
        try:
            self.logger.info("开始BUSCO质量评估分析流程|Starting BUSCO quality assessment analysis pipeline")
            self.logger.info(f"输入路径|Input path: {self.config.input_path}")
            self.logger.info(f"数据库|Database: {self.config.lineage}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

            # 检查依赖|Check dependencies
            self.check_dependencies()

            # 获取输入文件列表|Get input file list
            self.logger.info("="*80)
            self.logger.info("步骤1: 获取输入文件|Step 1: Get input files")
            self.logger.info("="*80)

            input_files = self.file_manager.get_input_files()

            if not input_files:
                raise RuntimeError("未找到有效的输入文件|No valid input files found")

            # 批量运行BUSCO分析|Batch run BUSCO analyses
            self.logger.info("="*80)
            self.logger.info("步骤2: 运行BUSCO分析|Step 2: Run BUSCO analyses")
            self.logger.info("="*80)

            results = []
            failed_samples = []

            for i, (file_path, sample_name) in enumerate(input_files, 1):
                self.logger.info(f"处理样本 {i}/{len(input_files)}|Processing sample {i}/{len(input_files)}: {sample_name}")

                success, result_data = self.busco_runner.run_busco_analysis(file_path, sample_name)

                if success and result_data:
                    results.append(result_data)
                    self.logger.info(f"样本分析成功|Sample analysis successful: {sample_name}")
                else:
                    failed_samples.append(sample_name)
                    self.logger.error(f"样本分析失败|Sample analysis failed: {sample_name}")

            # 处理和保存结果|Process and save results
            self.logger.info("="*80)
            self.logger.info("步骤3: 处理和保存结果|Step 3: Process and save results")
            self.logger.info("="*80)

            # 编译结果表格|Compile results table
            results_df = self.results_processor.compile_results(results)

            # 添加失败样本|Add failed samples
            if failed_samples:
                results_df = self.results_processor.add_failed_samples(results_df, failed_samples)

            # 保存结果|Save results
            self.results_processor.save_results(results_df)

            # 生成总结报告|Generate summary report
            self.summary_generator.generate_summary_report(results_df)

            self.logger.info("="*80)
            self.logger.info("BUSCO质量评估分析完成|BUSCO quality assessment analysis completed")
            self.logger.info("="*80)
            self.logger.info(f"总样本数|Total samples: {len(input_files)}")
            self.logger.info(f"成功样本数|Successful samples: {len(results)}")
            self.logger.info(f"失败样本数|Failed samples: {len(failed_samples)}")
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")

        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止|Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='BUSCO质量评估分析脚本|BUSCO Quality Assessment Analysis Script',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog='''
示例|Examples:
  %(prog)s -i genome.fa -l eukaryota_odb12
        '''
    )

    # 必需参数|Required arguments
    parser.add_argument('-i', '--input', required=True,
                       help='输入文件或目录路径|Input file or directory path')
    parser.add_argument('-l', '--lineage', required=True,
                       help='BUSCO数据库/谱系名称|BUSCO database lineage name')

    # 可选参数|Optional arguments
    parser.add_argument('-o', '--output-dir', default='./busco_output',
                       help='输出目录|Output directory')
    parser.add_argument('-m', '--mode', default='genome',
                       choices=['genome', 'geno', 'transcriptome', 'tran', 'proteins', 'prot'],
                       help='BUSCO分析模式|BUSCO analysis mode')
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='CPU线程数|Number of CPU threads')
    parser.add_argument('--sample-suffix', default='*.fa',
                       help='样本名称提取后缀模式|Sample name extraction suffix pattern')

    # 输出格式|Output format
    parser.add_argument('--output-format', default='txt',
                       choices=['txt', 'csv', 'xlsx'],
                       help='输出文件格式|Output file format')

    # BUSCO高级参数|BUSCO advanced parameters
    parser.add_argument('-f', '--force', action='store_true',
                       help='强制重写现有文件|Force rewriting of existing files')
    parser.add_argument('--augustus', action='store_true',
                       help='使用Augustus基因预测器|Use Augustus gene predictor')
    parser.add_argument('--augustus-parameters',
                       help='Augustus附加参数|Additional Augustus parameters')
    parser.add_argument('--augustus-species',
                       help='Augustus物种名称|Augustus species name')
    parser.add_argument('--auto-lineage', action='store_true',
                       help='自动选择谱系|Automatically select lineage')
    parser.add_argument('--auto-lineage-euk', action='store_true',
                       help='自动选择真核生物谱系|Automatically select eukaryote lineage')
    parser.add_argument('--auto-lineage-prok', action='store_true',
                       help='自动选择原核生物谱系|Automatically select prokaryote lineage')
    parser.add_argument('--contig-break', type=int, default=10,
                       help='Contig断点N数量|Number of Ns for contig break')
    parser.add_argument('--datasets-version', default='odb12',
                       help='数据集版本|Dataset version')
    parser.add_argument('--download-path',
                       help='数据集下载路径|Dataset download path')
    parser.add_argument('-e', '--evalue', type=float, default=1e-3,
                       help='BLAST E值阈值|BLAST E-value threshold')
    parser.add_argument('--limit', type=int, default=3,
                       help='候选区域数量限制|Candidate region limit')
    parser.add_argument('--long', action='store_true',
                       help='启用Augustus长模式优化|Enable Augustus long mode optimization')
    parser.add_argument('--metaeuk', action='store_true',
                       help='使用Metaeuk基因预测器|Use Metaeuk gene predictor')
    parser.add_argument('--metaeuk-parameters',
                       help='Metaeuk附加参数|Additional Metaeuk parameters')
    parser.add_argument('--metaeuk-rerun-parameters',
                       help='Metaeuk重运行参数|Metaeuk rerun parameters')
    parser.add_argument('--miniprot', action='store_true',
                       help='使用Miniprot基因预测器|Use Miniprot gene predictor')
    parser.add_argument('--skip-bbtools', action='store_true',
                       help='跳过BBTools统计|Skip BBTools statistics')
    parser.add_argument('--offline', action='store_true',
                       help='离线模式|Offline mode')
    parser.add_argument('-r', '--restart', action='store_true',
                       help='重启未完成的分析|Restart incomplete analysis')
    parser.add_argument('-q', '--quiet', action='store_true',
                       help='静默模式|Quiet mode')
    parser.add_argument('--scaffold-composition', action='store_true',
                       help='生成scaffold组成文件|Generate scaffold composition file')
    parser.add_argument('--tar', action='store_true',
                       help='压缩子目录|Compress subdirectories')

    # 工具路径|Tool paths
    parser.add_argument('--busco-path', default='busco',
                       help='BUSCO软件路径|BUSCO software path')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    analyzer = BUSCOAnalyzer(
        input_path=args.input,
        lineage=args.lineage,
        output_dir=args.output_dir,
        mode=args.mode,
        threads=args.threads,
        sample_suffix=args.sample_suffix,
        output_format=args.output_format,
        force=args.force,
        augustus=args.augustus,
        augustus_parameters=args.augustus_parameters,
        augustus_species=args.augustus_species,
        auto_lineage=args.auto_lineage,
        auto_lineage_euk=args.auto_lineage_euk,
        auto_lineage_prok=args.auto_lineage_prok,
        contig_break=args.contig_break,
        datasets_version=args.datasets_version,
        download_path=args.download_path,
        evalue=args.evalue,
        limit=args.limit,
        long=args.long,
        metaeuk=args.metaeuk,
        metaeuk_parameters=args.metaeuk_parameters,
        metaeuk_rerun_parameters=args.metaeuk_rerun_parameters,
        miniprot=args.miniprot,
        skip_bbtools=args.skip_bbtools,
        offline=args.offline,
        restart=args.restart,
        quiet=args.quiet,
        scaffold_composition=args.scaffold_composition,
        tar=args.tar,
        busco_path=args.busco_path
    )

    analyzer.run_analysis()

if __name__ == "__main__":
    main()
