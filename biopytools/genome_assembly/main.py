"""
基因组组装分析主程序模块 | Genome Assembly Analysis Main Module
"""

import argparse
import sys
from .config import AssemblyConfig
from .utils import AssemblyLogger, CommandRunner, check_dependencies
from .data_processing import ReadProcessor
from .assembly import VerkkoAssembler, HifiasmAssembler, GraphasingProcessor
from .quality_control import ContaminationScreener, AssemblyAnnotator
from .quality_assessment import QualityAssessor
from .alignment import AlignmentAnalyzer, TrioAnalyzer
from .results import ResultsProcessor

class GenomeAssemblyAnalyzer:
    """基因组组装分析主类 | Main Genome Assembly Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = AssemblyConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = AssemblyLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.read_processor = ReadProcessor(self.config, self.logger, self.cmd_runner)
        self.verkko_assembler = VerkkoAssembler(self.config, self.logger, self.cmd_runner)
        self.hifiasm_assembler = HifiasmAssembler(self.config, self.logger, self.cmd_runner)
        self.graphasing_processor = GraphasingProcessor(self.config, self.logger, self.cmd_runner)
        self.contamination_screener = ContaminationScreener(self.config, self.logger, self.cmd_runner)
        self.assembly_annotator = AssemblyAnnotator(self.config, self.logger, self.cmd_runner)
        self.quality_assessor = QualityAssessor(self.config, self.logger, self.cmd_runner)
        self.alignment_analyzer = AlignmentAnalyzer(self.config, self.logger, self.cmd_runner)
        self.trio_analyzer = TrioAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的基因组组装分析流程 | Run complete genome assembly analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始基因组组装分析流程 | Starting genome assembly analysis pipeline")
            self.logger.info("=" * 80)
            
            # 步骤1: 检查依赖软件 | Step 1: Check dependencies
            self.logger.info("\n步骤1: 检查依赖软件 | Step 1: Check dependencies")
            self.check_dependencies()
            
            # 步骤2: 验证输入数据 | Step 2: Validate input data
            self.logger.info("\n步骤2: 验证输入数据 | Step 2: Validate input data")
            if not self.read_processor.validate_reads():
                raise RuntimeError("输入数据验证失败 | Input data validation failed")
            
            assembly_results = {}
            
            # 步骤3: Verkko主要组装 | Step 3: Primary assembly with Verkko
            self.logger.info("\n步骤3: Verkko主要组装 | Step 3: Primary assembly with Verkko")
            verkko_assembly = self.verkko_assembler.run_verkko_assembly()
            if verkko_assembly:
                assembly_results['verkko_assembly'] = verkko_assembly
            
            # 步骤4: hifiasm补充组装 | Step 4: Complementary assembly with hifiasm
            self.logger.info("\n步骤4: hifiasm补充组装 | Step 4: Complementary assembly with hifiasm")
            hifiasm_assembly = self.hifiasm_assembler.run_hifiasm_assembly()
            if hifiasm_assembly:
                assembly_results['hifiasm_assembly'] = hifiasm_assembly
            
            # 选择主要组装文件 | Choose primary assembly file
            primary_assembly = verkko_assembly if verkko_assembly else hifiasm_assembly
            if not primary_assembly:
                raise RuntimeError("所有组装都失败了 | All assemblies failed")
            
            # 步骤5: Graphasing分期分析 | Step 5: Graphasing phasing analysis
            self.logger.info("\n步骤5: Graphasing分期分析 | Step 5: Graphasing phasing analysis")
            phased_assembly = self.graphasing_processor.run_graphasing(primary_assembly)
            assembly_results['phased_assembly'] = phased_assembly
            
            # 步骤6: 外源污染筛查 | Step 6: Contamination screening
            self.logger.info("\n步骤6: 外源污染筛查 | Step 6: Contamination screening")
            clean_assembly = self.contamination_screener.screen_contamination(phased_assembly)
            assembly_results['clean_assembly'] = clean_assembly
            
            # 步骤7: 组装注释 | Step 7: Assembly annotation
            self.logger.info("\n步骤7: 组装注释 | Step 7: Assembly annotation")
            flagger_results = self.assembly_annotator.run_flagger(clean_assembly)
            nucfreq_results = self.assembly_annotator.run_nucfreq(clean_assembly)
            assembly_results.update(flagger_results)
            assembly_results.update(nucfreq_results)
            
            # 步骤8: 质量评估 | Step 8: Quality assessment
            self.logger.info("\n步骤8: 质量评估 | Step 8: Quality assessment")
            merqury_results = self.quality_assessor.run_merqury(clean_assembly)
            inspector_results = self.quality_assessor.run_inspector(clean_assembly)
            deepvariant_results = self.quality_assessor.run_deepvariant_qv(clean_assembly)
            compleasm_results = self.quality_assessor.run_compleasm(clean_assembly)
            
            assembly_results.update(merqury_results)
            assembly_results.update(inspector_results)
            assembly_results.update(deepvariant_results)
            assembly_results.update(compleasm_results)
            
            # 步骤9: 参考基因组比对分析 | Step 9: Reference genome alignment analysis
            self.logger.info("\n步骤9: 参考基因组比对分析 | Step 9: Reference genome alignment analysis")
            alignment_results = self.alignment_analyzer.align_to_reference(clean_assembly)
            t2t_results = self.alignment_analyzer.analyze_t2t_status(alignment_results)
            assembly_results.update(alignment_results)
            assembly_results.update(t2t_results)
            
            # 步骤10: 家系分析(如果适用) | Step 10: Family analysis (if applicable)
            if self.config.trio_mode:
                self.logger.info("\n步骤10: 家系三元组分析 | Step 10: Family trio analysis")
                # 这里需要为每个家系成员运行完整的组装流程
                # For demonstration, we'll assume assemblies are already available
                trio_results = self.trio_analyzer.analyze_parental_support(
                    clean_assembly, clean_assembly, clean_assembly
                )
                assembly_results.update(trio_results)
            
            # 步骤11: 生成总结报告 | Step 11: Generate summary report
            self.logger.info("\n步骤11: 生成总结报告 | Step 11: Generate summary report")
            summary_report = self.results_processor.generate_summary_report(assembly_results)
            assembly_results['summary_report'] = summary_report
            
            # 完成信息 | Completion information
            self.logger.info("\n" + "=" * 80)
            self.logger.info("基因组组装分析完成！| Genome assembly analysis completed!")
            self.logger.info("=" * 80)
            self.logger.info(f"主要组装文件 | Primary assembly file: {clean_assembly}")
            self.logger.info(f"结果目录 | Results directory: {self.config.output_dir}")
            self.logger.info("主要输出文件 | Main output files:")
            self.logger.info("  - verkko_assembly/: Verkko组装结果 | Verkko assembly results")
            self.logger.info("  - hifiasm_assembly/: hifiasm组装结果 | hifiasm assembly results")
            self.logger.info("  - graphasing_output/: 分期分析结果 | Phasing analysis results")
            self.logger.info("  - contamination_screen/: 污染筛查结果 | Contamination screening results")
            self.logger.info("  - merqury_output/: 质量评估结果 | Quality assessment results")
            self.logger.info("  - alignment_output/: 比对分析结果 | Alignment analysis results")
            if self.config.trio_mode:
                self.logger.info("  - trio_analysis/: 家系分析结果 | Family analysis results")
            self.logger.info("  - assembly_summary_report.txt: 总结报告 | Summary report")
            self.logger.info("  - assembly_analysis.log: 分析日志 | Analysis log")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='基因组组装分析工具 (模块化版本) | Genome Assembly Analysis Tool (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
基于Verkko和hifiasm的杂合全基因组组装流程 | Hybrid Genome Assembly Pipeline using Verkko and hifiasm

主要功能 | Main Features:
  - Verkko主要组装器 | Primary assembly with Verkko
  - hifiasm超长读长组装 | Ultra-long read assembly with hifiasm
  - Graphasing分期信号 | Phasing signal with Graphasing
  - 质量控制和注释 | Quality control and annotation
  - 组装质量评估 | Assembly quality assessment
  - 参考基因组比对 | Reference genome alignment
  - 家系三元组分析 | Family trio analysis

示例 | Examples:
  %(prog)s --hifi-reads hifi.fastq.gz -o assembly_results
  %(prog)s --hifi-reads hifi.fq.gz --ont-reads ont.fq.gz --genome-size 3.2g
  %(prog)s --hifi-reads child_hifi.fq.gz --trio-mode --parent1-reads p1.fq.gz --parent2-reads p2.fq.gz --child-reads child.fq.gz
  %(prog)s --hifi-reads hifi.fastq.gz --reference-genome T2T-CHM13.fa --threads 64 --memory 256
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('--hifi-reads', required=True,
                       help='HiFi测序reads文件路径 | HiFi sequencing reads file path')
    
    # 可选输入参数 | Optional input arguments
    parser.add_argument('--ont-reads',
                       help='ONT测序reads文件路径 | ONT sequencing reads file path')
    parser.add_argument('--reference-genome',
                       help='参考基因组文件路径 (用于比对分析) | Reference genome file path (for alignment analysis)')
    parser.add_argument('-o', '--output-dir', default='./assembly_output',
                       help='输出目录 | Output directory')
    
    # 组装参数 | Assembly parameters
    parser.add_argument('--genome-size', default='3.2g',
                       help='预期基因组大小 | Expected genome size')
    parser.add_argument('--ploidy', type=int, choices=[1, 2], default=2,
                       help='基因组倍性 | Genome ploidy')
    
    # Verkko参数 | Verkko parameters
    parser.add_argument('--verkko-version', default='1.4.1',
                       help='Verkko版本 | Verkko version')
    parser.add_argument('--verkko-grid', action='store_true',
                       help='使用网格计算 | Use grid computing')
    parser.add_argument('--verkko-local-cpus', type=int, default=32,
                       help='Verkko本地CPU数量 | Verkko local CPU count')
    parser.add_argument('--verkko-local-memory', type=int, default=128,
                       help='Verkko本地内存(GB) | Verkko local memory (GB)')
    parser.add_argument('--no-verkko-screen', action='store_false', dest='verkko_screen',
                       help='禁用Verkko筛选 | Disable Verkko screening')
    
    # hifiasm参数 | hifiasm parameters
    parser.add_argument('--hifiasm-version', default='0.19.6',
                       help='hifiasm版本 | hifiasm version')
    parser.add_argument('--no-hifiasm-ultra-long', action='store_false', dest='hifiasm_ultra_long',
                       help='禁用hifiasm超长模式 | Disable hifiasm ultra-long mode')
    parser.add_argument('--hifiasm-purge-level', type=int, default=3, choices=[0, 1, 2, 3],
                       help='hifiasm purge级别 | hifiasm purge level')
    parser.add_argument('--hifiasm-similarity', type=float, default=0.8,
                       help='hifiasm相似性阈值 | hifiasm similarity threshold')
    
    # Graphasing参数 | Graphasing parameters
    parser.add_argument('--graphasing-version', default='0.3.1-alpha',
                       help='Graphasing版本 | Graphasing version')
    parser.add_argument('--graphasing-kmer-size', type=int, default=21,
                       help='Graphasing k-mer大小 | Graphasing k-mer size')
    
    # 质控参数 | Quality control parameters
    parser.add_argument('--skip-contamination-screen', action='store_false', dest='run_contamination_screen',
                       help='跳过外源污染筛查 | Skip contamination screening')
    parser.add_argument('--fcs-version', default='0.4.0',
                       help='NCBI FCS版本 | NCBI FCS version')
    parser.add_argument('--skip-flagger', action='store_false', dest='run_flagger',
                       help='跳过Flagger错误注释 | Skip Flagger error annotation')
    parser.add_argument('--flagger-version', default='0.3.3',
                       help='Flagger版本 | Flagger version')
    
    # 质量评估参数 | Quality assessment parameters
    parser.add_argument('--skip-merqury', action='store_false', dest='run_merqury',
                       help='跳过Merqury质量评估 | Skip Merqury quality assessment')
    parser.add_argument('--merqury-version', default='1.0',
                       help='Merqury版本 | Merqury version')
    parser.add_argument('--skip-inspector', action='store_false', dest='run_inspector',
                       help='跳过Inspector组装检查 | Skip Inspector assembly inspection')
    parser.add_argument('--inspector-version', default='1.2',
                       help='Inspector版本 | Inspector version')
    parser.add_argument('--skip-deepvariant', action='store_false', dest='run_deepvariant',
                       help='跳过DeepVariant质量值估计 | Skip DeepVariant quality value estimation')
    parser.add_argument('--deepvariant-version', default='1.6',
                       help='DeepVariant版本 | DeepVariant version')
    parser.add_argument('--skip-compleasm', action='store_false', dest='run_compleasm',
                       help='跳过compleasm基因完整性评估 | Skip compleasm gene completeness assessment')
    parser.add_argument('--compleasm-version', default='0.2.5',
                       help='compleasm版本 | compleasm version')
    parser.add_argument('--orthodb-version', default='10',
                       help='OrthoDB版本 | OrthoDB version')
    
    # 比对参数 | Alignment parameters
    parser.add_argument('--minimap2-version', default='2.26',
                       help='minimap2版本 | minimap2 version')
    parser.add_argument('--mashmap-version', default='3.1.3',
                       help='mashmap版本 | mashmap version')
    
    # 系统参数 | System parameters
    parser.add_argument('--threads', type=int, default=32,
                       help='线程数 | Number of threads')
    parser.add_argument('--memory', type=int, default=128,
                       help='内存大小(GB) | Memory size (GB)')
    parser.add_argument('--tmp-dir', default='/tmp',
                       help='临时目录 | Temporary directory')
    parser.add_argument('--keep-temp', action='store_true',
                       help='保留临时文件 | Keep temporary files')
    
    # 工具路径 | Tool paths
    parser.add_argument('--verkko-path', default='verkko',
                       help='Verkko软件路径 | Verkko software path')
    parser.add_argument('--hifiasm-path', default='hifiasm',
                       help='hifiasm软件路径 | hifiasm software path')
    parser.add_argument('--graphasing-path', default='graphasing',
                       help='Graphasing软件路径 | Graphasing software path')
    parser.add_argument('--fcs-path', default='fcs.py',
                       help='NCBI FCS软件路径 | NCBI FCS software path')
    parser.add_argument('--flagger-path', default='flagger',
                       help='Flagger软件路径 | Flagger software path')
    parser.add_argument('--merqury-path', default='merqury.sh',
                       help='Merqury软件路径 | Merqury software path')
    parser.add_argument('--nucfreq-path', default='nucfreq',
                       help='NucFreq软件路径 | NucFreq software path')
    parser.add_argument('--inspector-path', default='inspector.py',
                       help='Inspector软件路径 | Inspector software path')
    parser.add_argument('--deepvariant-path', default='run_deepvariant',
                       help='DeepVariant软件路径 | DeepVariant software path')
    parser.add_argument('--compleasm-path', default='compleasm',
                       help='compleasm软件路径 | compleasm software path')
    parser.add_argument('--minimap2-path', default='minimap2',
                       help='minimap2软件路径 | minimap2 software path')
    parser.add_argument('--mashmap-path', default='mashmap',
                       help='mashmap软件路径 | mashmap software path')
    
    # 家系分析参数 | Family analysis parameters
    parser.add_argument('--trio-mode', action='store_true',
                       help='启用家系三元组分析模式 | Enable family trio analysis mode')
    parser.add_argument('--parent1-reads',
                       help='亲本1的reads文件路径 | Parent1 reads file path')
    parser.add_argument('--parent2-reads',
                       help='亲本2的reads文件路径 | Parent2 reads file path')
    parser.add_argument('--child-reads',
                       help='子代的reads文件路径 | Child reads file path')
    
    args = parser.parse_args()
    
    # 验证家系模式参数 | Validate trio mode parameters
    if args.trio_mode:
        if not all([args.parent1_reads, args.parent2_reads, args.child_reads]):
            parser.error("家系模式需要提供 --parent1-reads, --parent2-reads, 和 --child-reads 参数 | "
                        "Trio mode requires --parent1-reads, --parent2-reads, and --child-reads arguments")
    
    # 转换参数为字典 | Convert arguments to dictionary
    config_dict = vars(args)
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = GenomeAssemblyAnalyzer(**config_dict)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
