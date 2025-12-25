#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
HiFiasm分析主程序模块 | HiFiasm Analysis Main Module
"""

# ==============================================================================
# 导入模块 | Import Modules (合并自两个文件)
# ==============================================================================
import os
import sys
import time
import argparse
from pathlib import Path
from typing import Dict, Any, Optional

from .config import HifiasmConfig
from .utils import (
    HifiasmLogger, CommandRunner, check_dependencies,
    estimate_resources, setup_signal_handlers, format_duration,
    create_directory_structure
)
from .assembly import HifiasmAssembler, GFAConverter
from .quality_assessment import BUSCOAssessor, QUASTAssessor, HaplotypeAnalyzer
from .data_processing import StatisticsCalculator, FormatConverter
from .results import ResultsProcessor, SummaryGenerator


# ==============================================================================
# 核心分析类 | Core Analyzer Class (来自原始 main.py)
# ==============================================================================

class HifiasmAnalyzer:
    """HiFiasm分析主类 | Main HiFiasm Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = HifiasmConfig(**kwargs)
        self.config.validate()
        
        # 保存配置文件 | Save configuration file
        self.config_file = self.config.save_config()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = HifiasmLogger(
            output_dir=self.config.output_dir,
            log_level=self.config.log_level,
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(
            logger=self.logger,
            working_dir=self.config.output_dir,
            timeout=self.config.max_runtime * 3600  # 转换为秒
        )
        
        # 设置信号处理器 | Setup signal handlers
        setup_signal_handlers(self.cmd_runner, self.logger)
        
        # 创建目录结构 | Create directory structure
        self.directories = create_directory_structure(self.config.output_dir)
        
        # 初始化各个处理器 | Initialize processors
        self._initialize_processors()
        
        # 记录初始化信息 | Log initialization info
        self._log_initialization_info()
    
    def _initialize_processors(self):
        """初始化处理器 | Initialize processors"""
        # 组装相关 | Assembly related
        self.hifiasm_assembler = HifiasmAssembler(self.config, self.logger, self.cmd_runner)
        self.gfa_converter = GFAConverter(self.config, self.logger, self.cmd_runner)
        
        # 质量评估相关 | Quality assessment related
        self.busco_assessor = BUSCOAssessor(self.config, self.logger, self.cmd_runner)
        self.quast_assessor = QUASTAssessor(self.config, self.logger, self.cmd_runner)
        self.haplotype_analyzer = HaplotypeAnalyzer(self.config, self.logger)
        
        # 数据处理相关 | Data processing related
        self.statistics_calculator = StatisticsCalculator(self.config, self.logger)
        self.format_converter = FormatConverter(self.config, self.logger, self.cmd_runner)
        
        # 结果处理相关 | Results processing related
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def _log_initialization_info(self):
        """记录初始化信息 | Log initialization information"""
        self.logger.info("="*80)
        self.logger.info("HiFiasm基因组组装分析流水线 | HiFiasm Genome Assembly Analysis Pipeline")
        self.logger.info("="*80)
        self.logger.info(f"输入文件 | Input file: {self.config.input_reads}")
        self.logger.info(f"输出目录 | Output directory: {self.config.output_dir}")
        self.logger.info(f"样本前缀 | Sample prefix: {self.config.prefix}")
        self.logger.info(f"线程数 | Threads: {self.config.threads}")
        self.logger.info(f"基因组大小 | Genome size: {self.config.estimate_genome_size()}")
        self.logger.info(f"组装类型 | Assembly type: {self.config.assembly_type}")
        
        if self.config.ont_reads:
            self.logger.info(f"ONT数据 | ONT data: {self.config.ont_reads}")
        if self.config.hi_c_1 and self.config.hi_c_2:
            self.logger.info(f"Hi-C数据 | Hi-C data: {self.config.hi_c_1}, {self.config.hi_c_2}")
        
        self.logger.info(f"配置文件 | Config file: {self.config_file}")
        self.logger.info("="*80)
    
    def run_analysis(self) -> bool:
        """运行完整的分析流程 | Run complete analysis pipeline"""
        start_time = time.time()
        
        try:
            self.logger.start("HiFiasm基因组组装分析流程 | HiFiasm genome assembly analysis pipeline")
            
            # 步骤1: 检查依赖 | Step 1: Check dependencies
            if not self._check_dependencies():
                return False
            
            # 步骤2: 估计资源需求 | Step 2: Estimate resource requirements
            self._estimate_resources()
            
            # 步骤3: 预处理检查 | Step 3: Preprocessing checks
            if not self._preprocessing_checks():
                return False
            
            # 步骤4: HiFiasm组装 | Step 4: HiFiasm assembly
            if not self._run_hifiasm_assembly():
                return False
            
            # 步骤5: 格式转换 | Step 5: Format conversion
            if not self._convert_formats():
                return False
            
            # 步骤6: 基本统计 | Step 6: Basic statistics
            if not self._calculate_statistics():
                return False
            
            # 步骤7: 质量评估 | Step 7: Quality assessment
            if not self._run_quality_assessment():
                return False
            
            # 步骤8: 单倍型分析 | Step 8: Haplotype analysis
            if self.config.analyze_haplotypes:
                if not self._analyze_haplotypes():
                    return False
            
            # 步骤9: 生成最终结果 | Step 9: Generate final results
            if not self._generate_final_results():
                return False
            
            # 步骤10: 清理临时文件 | Step 10: Cleanup temporary files
            if not self.config.keep_intermediate:
                self._cleanup_temporary_files()
            
            end_time = time.time()
            total_time = end_time - start_time
            
            self.logger.success("HiFiasm分析流程完成 | HiFiasm analysis pipeline completed")
            self.logger.info(f"总耗时 | Total time: {format_duration(total_time)}")
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")
            
            return True
            
        except KeyboardInterrupt:
            self.logger.warning("分析被用户中断 | Analysis interrupted by user")
            return False
        except Exception as e:
            self.logger.error(f"分析流程失败 | Analysis pipeline failed: {e}")
            if self.config.debug:
                import traceback
                self.logger.debug(traceback.format_exc())
            return False
        finally:
            # 清理命令执行器 | Cleanup command runner
            self.cmd_runner.cleanup()
    
    def _check_dependencies(self) -> bool:
        """检查依赖 | Check dependencies"""
        self.logger.start("检查依赖软件 | Checking dependencies")
        
        if self.config.dry_run:
            self.logger.info("试运行模式：跳过依赖检查 | Dry run mode: skipping dependency check")
            return True
        
        return check_dependencies(self.config, self.logger)
    
    def _estimate_resources(self):
        """估计资源需求 | Estimate resource requirements"""
        self.logger.start("估计资源需求 | Estimating resource requirements")
        
        resources = estimate_resources(self.config.input_reads, self.logger)
        
        if resources['recommended_memory_gb'] > self.config.memory:
            self.logger.warning(
                f"推荐内存({resources['recommended_memory_gb']}GB) 大于配置内存({self.config.memory}GB) | "
                f"Recommended memory({resources['recommended_memory_gb']}GB) > configured memory({self.config.memory}GB)"
            )
        
        if resources['recommended_threads'] > self.config.threads:
            self.logger.warning(
                f"推荐线程数({resources['recommended_threads']}) 大于配置线程数({self.config.threads}) | "
                f"Recommended threads({resources['recommended_threads']}) > configured threads({self.config.threads})"
            )
    
    def _preprocessing_checks(self) -> bool:
        """预处理检查 | Preprocessing checks"""
        self.logger.start("预处理检查 | Preprocessing checks")
        
        try:
            if not self._validate_input_format():
                return False
            if not self._check_disk_space():
                return False
            if self.config.resume:
                if not self._check_resume_status():
                    self.logger.warning("无法恢复，将重新开始分析 | Cannot resume, will start fresh analysis")
                    self.config.resume = False
            return True
        except Exception as e:
            self.logger.error(f"预处理检查失败 | Preprocessing checks failed: {e}")
            return False
    
    def _validate_input_format(self) -> bool:
        """验证输入文件格式 | Validate input file format"""
        valid_extensions = ['.fq', '.fastq', '.fq.gz', '.fastq.gz', '.fa', '.fasta', '.fa.gz', '.fasta.gz']
        file_path = Path(self.config.input_reads)
        if not any(str(file_path).endswith(ext) for ext in valid_extensions):
            self.logger.error(f"不支持的文件格式 | Unsupported file format: {file_path.suffix}")
            return False
        self.logger.info(f"输入文件格式验证通过 | Input file format validation passed: {file_path.suffix}")
        return True
    
    def _check_disk_space(self) -> bool:
        """检查磁盘空间 | Check disk space"""
        import shutil
        input_size = os.path.getsize(self.config.input_reads)
        required_space = input_size * 5
        free_space = shutil.disk_usage(self.config.output_dir).free
        if free_space < required_space:
            self.logger.error(
                f"磁盘空间不足 | Insufficient disk space: "
                f"需要 {required_space/1e9:.1f}GB, 可用 {free_space/1e9:.1f}GB"
            )
            return False
        self.logger.info(f"磁盘空间充足 | Sufficient disk space: 可用 {free_space/1e9:.1f}GB")
        return True
    
    def _check_resume_status(self) -> bool:
        """检查恢复状态 | Check resume status"""
        self.logger.info("恢复分析功能开发中 | Resume analysis functionality under development")
        return False
    
    def _run_hifiasm_assembly(self) -> bool:
        """运行HiFiasm组装 | Run HiFiasm assembly"""
        self.logger.start("HiFiasm基因组组装 | HiFiasm genome assembly")
        if self.config.dry_run:
            self.logger.info("试运行模式：模拟HiFiasm组装")
            return True
        return self.hifiasm_assembler.run_assembly()
    
    def _convert_formats(self) -> bool:
        """转换格式 | Convert formats"""
        self.logger.start("格式转换 | Format conversion")
        if self.config.dry_run:
            self.logger.info("试运行模式：模拟格式转换")
            return True
        if not self.gfa_converter.convert_gfa_to_fasta():
            return False
        return self.format_converter.convert_all_formats()
    
    def _calculate_statistics(self) -> bool:
        """计算统计信息 | Calculate statistics"""
        self.logger.start("计算组装统计信息 | Calculating assembly statistics")
        if self.config.dry_run:
            self.logger.info("试运行模式：模拟统计计算")
            return True
        return self.statistics_calculator.calculate_all_statistics()
    
    def _run_quality_assessment(self) -> bool:
        """运行质量评估 | Run quality assessment"""
        self.logger.start("质量评估 | Quality assessment")
        if self.config.dry_run:
            self.logger.info("试运行模式：模拟质量评估")
            return True
        success = True
        if not self.config.skip_busco:
            if not self.busco_assessor.run_busco_assessment():
                success = False
        if not self.config.skip_quast:
            if not self.quast_assessor.run_quast_assessment():
                success = False
        return success
    
    def _analyze_haplotypes(self) -> bool:
        """分析单倍型 | Analyze haplotypes"""
        self.logger.start("单倍型差异分析 | Haplotype difference analysis")
        if self.config.dry_run:
            self.logger.info("试运行模式：模拟单倍型分析")
            return True
        return self.haplotype_analyzer.analyze_haplotypes()
    
    def _generate_final_results(self) -> bool:
        """生成最终结果 | Generate final results"""
        self.logger.start("生成最终结果 | Generating final results")
        if self.config.dry_run:
            self.logger.info("试运行模式：模拟结果生成")
            return True
        if not self.results_processor.process_all_results():
            return False
        if not self.summary_generator.generate_summary():
            return False
        return True
    
    def _cleanup_temporary_files(self):
        """清理临时文件 | Cleanup temporary files"""
        self.logger.start("清理临时文件 | Cleaning up temporary files")
        try:
            import shutil
            tmp_dir = Path(self.directories['tmp'])
            if tmp_dir.exists():
                shutil.rmtree(tmp_dir)
                self.logger.info(f"已删除临时目录 | Deleted temporary directory: {tmp_dir}")
        except Exception as e:
            self.logger.warning(f"清理临时文件失败 | Failed to cleanup temporary files: {e}")
    
    def get_analysis_status(self) -> Dict[str, Any]:
        """获取分析状态 | Get analysis status"""
        return {'config': self.config.to_dict(), 'directories': self.directories, 'steps_completed': [], 'current_step': None, 'total_steps': 10}
    
    def resume_analysis(self) -> bool:
        """恢复分析 | Resume analysis"""
        self.logger.info("恢复分析功能开发中 | Resume analysis functionality under development")
        return self.run_analysis()


# ==============================================================================
# 命令行接口 | Command-Line Interface (来自原始 run_hifiasm.py)
# ==============================================================================

def create_argument_parser():
    """创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='HiFiasm基因组组装完整流水线 | HiFiasm Genome Assembly Complete Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 | Usage Examples:
  # 基本用法 | Basic usage
  %(prog)s -i sample.hifi.fq.gz -o hifiasm_results -p sample_prefix
  
  # 二倍体组装 | Diploid assembly
  %(prog)s -i reads.hifi.fq.gz -o output --hg-size 1.4g --purge-max 65 -s 0.75
  
  # 包含ONT数据的组装 | Assembly with ONT data
  %(prog)s -i hifi.fq.gz --ont-reads ont.fq.gz -o output -p sample
        """
    )
    
    required = parser.add_argument_group('必需参数 | Required arguments')
    required.add_argument('-i', '--input-reads', required=True, help='输入HiFi测序数据文件 | Input HiFi sequencing data file')
    
    basic = parser.add_argument_group('基本参数 | Basic arguments')
    basic.add_argument('-o', '--output-dir', default='./hifiasm_output', help='输出目录 | Output directory')
    basic.add_argument('-p', '--prefix', default='sample', help='输出文件前缀 | Output file prefix')
    basic.add_argument('-t', '--threads', type=int, default=32, help='线程数 | Number of threads')
    
    hifiasm = parser.add_argument_group('HiFiasm组装参数 | HiFiasm assembly parameters')
    hifiasm.add_argument('--hg-size', default='auto', help='基因组大小估计 (如 1.4g, 2.1g) | Genome size estimation (e.g., 1.4g, 2.1g)')
    hifiasm.add_argument('-l', '--purge-level', type=int, default=3, help='purge级别 (0-3) | Purge level (0-3)')
    hifiasm.add_argument('--purge-max', type=int, default=65, help='最大purge覆盖度 | Maximum purge coverage')
    hifiasm.add_argument('-s', '--similarity-threshold', type=float, default=0.75, help='相似性阈值 | Similarity threshold')
    hifiasm.add_argument('--ont-reads', help='ONT长读长数据文件 | ONT long-read data file')
    hifiasm.add_argument('--hi-c-1', help='Hi-C第一端数据文件 | Hi-C first-end data file')
    hifiasm.add_argument('--hi-c-2', help='Hi-C第二端数据文件 | Hi-C second-end data file')
    hifiasm.add_argument('--extra-hifiasm-args', default='', help='额外的HiFiasm参数 | Additional HiFiasm arguments')
    
    quality = parser.add_argument_group('质量评估参数 | Quality assessment parameters')
    quality.add_argument('--skip-busco', action='store_true', help='跳过BUSCO质量评估 | Skip BUSCO quality assessment')
    quality.add_argument('--busco-lineage', default='auto', help='BUSCO谱系数据集 (如 embryophyta_odb10) | BUSCO lineage dataset')
    quality.add_argument('--busco-mode', default='genome', choices=['genome', 'proteins', 'transcriptome'], help='BUSCO评估模式 | BUSCO assessment mode')
    quality.add_argument('--skip-quast', action='store_true', help='跳过QUAST质量评估 | Skip QUAST quality assessment')
    quality.add_argument('--reference-genome', help='参考基因组文件 (用于QUAST) | Reference genome file (for QUAST)')
    
    analysis = parser.add_argument_group('分析参数 | Analysis parameters')
    analysis.add_argument('--analyze-haplotypes', action='store_true', help='分析单倍型差异 | Analyze haplotype differences')
    analysis.add_argument('--min-contig-length', type=int, default=1000, help='最小contig长度过滤 | Minimum contig length filter')
    analysis.add_argument('--generate-plots', action='store_true', help='生成可视化图表 | Generate visualization plots')
    analysis.add_argument('--assembly-type', default='auto', choices=['auto', 'diploid', 'triploid', 'polyploid'], help='组装类型 | Assembly type')
    
    output = parser.add_argument_group('输出控制参数 | Output control parameters')
    output.add_argument('--keep-intermediate', action='store_true', help='保留中间文件 | Keep intermediate files')
    output.add_argument('--compress-output', action='store_true', help='压缩输出文件 | Compress output files')
    output.add_argument('--output-formats', nargs='+', choices=['fasta', 'gfa', 'both'], default=['both'], help='输出格式选择 | Output format selection')
    
    system = parser.add_argument_group('系统参数 | System parameters')
    system.add_argument('--memory', type=int, default=64, help='内存大小(GB) | Memory size (GB)')
    system.add_argument('--tmp-dir', default='/tmp', help='临时目录 | Temporary directory')
    system.add_argument('--max-runtime', type=int, default=48, help='最大运行时间(小时) | Maximum runtime (hours)')
    system.add_argument('--resume', action='store_true', help='恢复中断的分析 | Resume interrupted analysis')
    
    tools = parser.add_argument_group('工具路径参数 | Tool paths parameters')
    tools.add_argument('--hifiasm-path', default='hifiasm', help='HiFiasm软件路径 | HiFiasm software path')
    tools.add_argument('--busco-path', default='busco', help='BUSCO软件路径 | BUSCO software path')
    tools.add_argument('--quast-path', default='quast', help='QUAST软件路径 | QUAST software path')
    tools.add_argument('--python-path', default='python3', help='Python解释器路径 | Python interpreter path')
    tools.add_argument('--samtools-path', default='samtools', help='Samtools软件路径 | Samtools software path')
    
    databases = parser.add_argument_group('数据库路径参数 | Database paths parameters')
    databases.add_argument('--busco-db-path', help='BUSCO数据库路径 | BUSCO database path')
    databases.add_argument('--busco-download-path', help='BUSCO数据集下载路径 | BUSCO dataset download path')
    
    advanced = parser.add_argument_group('高级参数 | Advanced parameters')
    advanced.add_argument('--debug', action='store_true', help='启用调试模式 | Enable debug mode')
    advanced.add_argument('--verbose', '-v', action='count', default=0, help='详细输出模式 (-v, -vv, -vvv) | Verbose output mode')
    advanced.add_argument('--log-level', choices=['DEBUG', 'INFO', 'WARNING', 'ERROR'], default='INFO', help='日志级别 | Log level')
    advanced.add_argument('--config-file', help='配置文件路径 | Configuration file path')
    advanced.add_argument('--dry-run', action='store_true', help='试运行模式 (不执行实际命令) | Dry run mode (do not execute actual commands)')
    
    return parser

def validate_arguments(args):
    """验证命令行参数 | Validate command line arguments"""
    errors = []
    
    if not os.path.exists(args.input_reads):
        errors.append(f"输入文件不存在 | Input file does not exist: {args.input_reads}")
    
    if args.threads <= 0:
        errors.append("线程数必须大于0 | Number of threads must be greater than 0")
    
    if args.memory <= 0:
        errors.append("内存大小必须大于0 | Memory size must be greater than 0")
    
    if args.purge_level not in range(0, 4):
        errors.append("purge级别必须在0-3之间 | Purge level must be between 0-3")
    
    if args.similarity_threshold <= 0 or args.similarity_threshold > 1:
        errors.append("相似性阈值必须在0-1之间 | Similarity threshold must be between 0-1")
    
    if (args.hi_c_1 and not args.hi_c_2) or (args.hi_c_2 and not args.hi_c_1):
        errors.append("Hi-C数据需要同时提供两端数据 | Hi-C data requires both end files")
    
    if args.reference_genome and not os.path.exists(args.reference_genome):
        errors.append(f"参考基因组文件不存在 | Reference genome file does not exist: {args.reference_genome}")
    
    if errors:
        for error in errors:
            print(f"错误 | Error: {error}", file=sys.stderr)
        sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    validate_arguments(args)
    
    os.makedirs(args.output_dir, exist_ok=True)
    
    try:
        analyzer = HifiasmAnalyzer(**vars(args))
        success = analyzer.run_analysis()
        
        if success:
            print("\n✅ HiFiasm分析流程成功完成！| HiFiasm analysis pipeline completed successfully!")
            print(f"📁 结果保存在 | Results saved in: {args.output_dir}")
            sys.exit(0)
        else:
            print("\n❌ HiFiasm分析流程执行失败 | HiFiasm analysis pipeline failed")
            sys.exit(1)
            
    except KeyboardInterrupt:
        print("\n⚠️  分析被用户中断 | Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 意外错误 | Unexpected error: {e}")
        if 'args' in locals() and args.debug:
            import traceback
            traceback.print_exc()
        sys.exit(1)

# ==============================================================================
# 执行入口 | Execution Entrypoint
# ==============================================================================

if __name__ == "__main__":
    main()