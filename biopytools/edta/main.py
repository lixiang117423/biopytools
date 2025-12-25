"""
🌾 EDTA植物基因组TE注释工具主程序模块 | EDTA Plant Genome TE Annotation Tool Main Module
"""

import argparse
import sys
import signal
from pathlib import Path

from .config import EDTAConfig
from .utils import (EDTALogger, CommandRunner, check_edta_dependencies, 
                   setup_output_directories, check_resume_capability, save_checkpoint)
from .data_processing import EDTAProcessor
from .analysis import EDTAAnalysisProcessor
from .visualization import EDTAVisualizer
from .results import ResultsGenerator

class EDTAAnalyzer:
    """EDTA分析主类 | Main EDTA Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = EDTAConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = EDTALogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 设置信号处理 | Setup signal handling
        self.cmd_runner = None
        signal.signal(signal.SIGINT, self._signal_handler)
        signal.signal(signal.SIGTERM, self._signal_handler)
        
        # 初始化输出目录 | Initialize output directories
        self.directories = setup_output_directories(self.config.output_path, self.logger)
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.data_processor = EDTAProcessor(self.config, self.logger, self.cmd_runner, self.directories)
        self.analysis_processor = EDTAAnalysisProcessor(self.config, self.logger, self.directories)
        self.visualizer = EDTAVisualizer(self.config, self.logger, self.directories)
        self.results_generator = ResultsGenerator(self.config, self.logger, self.directories)
    
    def _signal_handler(self, signum, frame):
        """处理中断信号 | Handle interrupt signals"""
        self.logger.warning(f"🛑 接收到中断信号 | Received interrupt signal: {signum}")
        
        if self.cmd_runner:
            self.cmd_runner.terminate_process()
        
        self.logger.info("🔚 程序被用户中断 | Program interrupted by user")
        sys.exit(1)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_edta_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的EDTA分析流程 | Run complete EDTA analysis pipeline"""
        try:
            self.logger.info("🚀 开始EDTA植物基因组TE注释分析流程 | Starting EDTA plant genome TE annotation analysis pipeline")
            # self.logger.info("🌾 =" * 60)
            
            # 步骤1: 检查依赖
            self.logger.info("📋 步骤1: 检查依赖软件 | Step 1: Checking dependencies")
            self.check_dependencies()
            
            # 步骤2: 检查断点续跑
            checkpoint = None
            if self.config.resume:
                self.logger.info("🔄 步骤2: 检查断点续跑状态 | Step 2: Checking resume status")
                checkpoint = check_resume_capability(self.config.output_path, self.logger)
            
            # 步骤3: 准备输入文件
            if not checkpoint or "input_files_prepared" not in checkpoint.get("data", {}):
                self.logger.info("📂 步骤3: 准备输入文件 | Step 3: Preparing input files")
                input_files = self.data_processor.prepare_input_files()
                save_checkpoint(self.config.output_path, "input_files_prepared", 
                              {"input_files": input_files}, self.logger)
            else:
                self.logger.info("📂 步骤3: 从检查点恢复输入文件 | Step 3: Resuming input files from checkpoint")
                input_files = checkpoint["data"]["input_files"]
            
            # 步骤4: 运行EDTA分析
            if not checkpoint or "edta_completed" not in checkpoint.get("data", {}):
                self.logger.info("🧬 步骤4: 运行EDTA分析 | Step 4: Running EDTA analysis")
                edta_results = self.data_processor.run_edta_analysis(input_files)
            else:
                self.logger.info("🧬 步骤4: 从检查点恢复EDTA结果 | Step 4: Resuming EDTA results from checkpoint")
                edta_results = checkpoint["data"].get("edta_results", {})
            
            # 步骤5: 处理分析结果
            self.logger.info("📊 步骤5: 处理分析结果 | Step 5: Processing analysis results")
            processed_results = self.analysis_processor.process_results(edta_results)
            
            # 步骤6: 生成可视化图表
            self.logger.info("🎨 步骤6: 生成可视化图表 | Step 6: Generating visualization plots")
            self.visualizer.generate_all_plots(processed_results)
            
            # 步骤7: 生成综合报告
            self.logger.info("📋 步骤7: 生成综合报告 | Step 7: Generating comprehensive reports")
            self.results_generator.generate_all_reports(edta_results, processed_results)
            
            # 完成分析
            # self.logger.info("🌾 =" * 20)
            self.logger.info("✅ EDTA分析流程成功完成！| EDTA analysis pipeline completed successfully!")
            self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            # self.logger.info("🌾 =" * 20)
            
            # 清理检查点文件
            resume_file = self.config.output_path / "logs" / "resume_checkpoint.json"
            if resume_file.exists():
                resume_file.unlink()
                self.logger.info("🔄 清理检查点文件 | Cleaned up checkpoint file")
            
        except KeyboardInterrupt:
            self.logger.warning("🛑 分析被用户中断 | Analysis interrupted by user")
            sys.exit(1)
        except Exception as e:
            self.logger.error(f"❌ 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            raise

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🌾 EDTA植物基因组TE注释工具 | EDTA Plant Genome TE Annotation Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s --genome plant.fa --output edta_results
  %(prog)s -g genome.fa -o results --species rice --threads 88
  %(prog)s --genome plant.fa --cds proteins.fa --sensitive 1
        """
    )
    
    # 必需参数
    parser.add_argument('-g', '--genome', required=True,
                       help='输入基因组FASTA文件路径 | Input genome FASTA file path')
    
    # 基本参数
    parser.add_argument('-o', '--output-dir', default='./edta_output', 
                       help='输出目录 | Output directory')
    parser.add_argument('-s', '--species', default='others', choices=['Rice', 'Maize', 'others'],
                       help='物种类型 | Species for TIR candidate identification')
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='使用的线程数 | Number of threads to use')
    
    # EDTA核心参数
    parser.add_argument('--step', default='all', choices=['all', 'filter', 'final', 'anno'],
                       help='指定运行的步骤 | Specify which steps to run')
    parser.add_argument('--sensitive', type=int, default=0, choices=[0, 1],
                       help='使用RepeatModeler增强敏感性 | Use RepeatModeler to enhance sensitivity')
    parser.add_argument('--anno', type=int, default=1, choices=[0, 1],
                       help='执行全基因组TE注释 | Perform whole-genome TE annotation')
    parser.add_argument('--cds', 
                       help='提供CDS序列文件 | Provide CDS sequences file')
    
    # 批量处理参数
    parser.add_argument('--genome-list',
                       help='包含多个基因组文件路径的列表文件 | File containing list of genome file paths')
    parser.add_argument('--batch-mode', action='store_true',
                       help='启用批量处理模式 | Enable batch processing mode')
    
    # 其他参数
    parser.add_argument('--resume', action='store_true',
                       help='从上次中断的地方继续分析 | Resume analysis from last interruption')
    parser.add_argument('--generate-plots', action='store_true', default=True,
                       help='生成可视化图表 | Generate visualization plots')
    parser.add_argument('--compare-results', action='store_true',
                       help='启用结果比较 | Enable results comparison')
    parser.add_argument('--check-dependencies', action='store_true',
                       help='检查依赖软件后退出 | Check dependencies and exit')
    
    args = parser.parse_args()
    
    # 处理批量模式
    genome_list = None
    if args.genome_list:
        try:
            with open(args.genome_list, 'r') as f:
                genome_list = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            args.batch_mode = True
        except FileNotFoundError:
            print(f"❌ 基因组列表文件不存在 | Genome list file not found: {args.genome_list}")
            sys.exit(1)
    
    # 创建分析器并运行
    try:
        analyzer = EDTAAnalyzer(
            genome=args.genome,
            genome_list=genome_list,
            output_dir=args.output_dir,
            species=args.species,
            threads=args.threads,
            step=args.step,
            sensitive=args.sensitive,
            anno=args.anno,
            cds=args.cds,
            batch_mode=args.batch_mode,
            resume=args.resume,
            generate_plots=args.generate_plots,
            compare_results=args.compare_results,
            check_dependencies=args.check_dependencies
        )
        
        # 如果只是检查依赖
        if args.check_dependencies:
            analyzer.check_dependencies()
            print("✅ 所有依赖检查通过 | All dependencies check passed")
            sys.exit(0)
        
        # 运行完整分析
        analyzer.run_analysis()
        
    except Exception as e:
        print(f"❌ 分析失败 | Analysis failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
