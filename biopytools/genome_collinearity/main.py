"""
基因组共线性分析主程序模块 | Genome Collinearity Analysis Main Module
"""

import argparse
import sys
from .config import CollinearityConfig
from .utils import CollinearityLogger, CommandRunner, check_dependencies
from .alignment import GenomeAligner
from .syri_analysis import SyRIAnalyzer
from .visualization import CollinearityVisualizer
from .results import SummaryGenerator

class GenomeCollinearityAnalyzer:
    """基因组共线性分析主类 | Main Genome Collinearity Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = CollinearityConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = CollinearityLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.aligner = GenomeAligner(self.config, self.logger, self.cmd_runner)
        self.syri_analyzer = SyRIAnalyzer(self.config, self.logger, self.cmd_runner)
        self.visualizer = CollinearityVisualizer(self.config, self.logger, self.cmd_runner)
        self.summary_generator = SummaryGenerator(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    # def run_analysis(self):
    #     """运行完整的共线性分析流程 | Run complete collinearity analysis pipeline"""
    #     try:
    #         self.logger.info("🚀 开始基因组共线性分析流程 | Starting genome collinearity analysis pipeline")
    #         self.logger.info(f"📊 样本数量 | Number of samples: {len(self.config.sample_list)}")
    #         self.logger.info(f"🔄 样本顺序 | Sample order: {' -> '.join(self.config.sample_list)}")
            
    #         if self.config.chromosome:
    #             self.logger.info(f"🧬 指定染色体 | Specified chromosome: {self.config.chromosome}")
    #         else:
    #             self.logger.info("🌍 分析范围 | Analysis scope: 全基因组 | Whole genome")
            
    #         # 检查依赖 | Check dependencies
    #         self.check_dependencies()
            
    #         # Step 1: 基因组比对 | Genome alignment
    #         self.logger.info("\n" + "="*20 + " 🧬 步骤1: 基因组比对 | Step 1: Genome Alignment " + "="*20)
    #         alignment_files = self.aligner.run_alignments()
    #         if not alignment_files:
    #             raise RuntimeError("❌ 基因组比对失败 | Genome alignment failed")
            
    #         # Step 2: SyRI结构变异分析 | SyRI structural variation analysis
    #         self.logger.info("\n" + "="*20 + " 🔬 步骤2: SyRI结构变异分析 | Step 2: SyRI Structural Variation Analysis " + "="*20)
    #         syri_files = self.syri_analyzer.run_syri_analysis(alignment_files)
    #         if not syri_files:
    #             raise RuntimeError("❌ SyRI结构变异分析失败 | SyRI structural variation analysis failed")
            
    #         # Step 3: PlotSR可视化 | PlotSR visualization
    #         self.logger.info("\n" + "="*20 + " 🎨 步骤3: PlotSR可视化 | Step 3: PlotSR Visualization " + "="*20)
    #         visualization_success = self.visualizer.create_plotsr_visualization(syri_files)
    #         if not visualization_success:
    #             raise RuntimeError("❌ PlotSR可视化失败 | PlotSR visualization failed")
            
    #         # Step 4: 生成总结报告 | Generate summary report
    #         self.logger.info("\n" + "="*20 + " 📋 步骤4: 生成总结报告 | Step 4: Generate Summary Report " + "="*20)
    #         self.summary_generator.generate_summary_report(alignment_files, syri_files)
            
    #         # 分析完成 | Analysis completed
    #         self.logger.info("\n" + "="*20 + " 🎉 基因组共线性分析完成 | Genome Collinearity Analysis Completed " + "="*20)
    #         self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
    #     except Exception as e:
    #         self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
    #         sys.exit(1)

    def run_analysis(self):
        """运行完整的共线性分析流程 | Run complete collinearity analysis pipeline"""
        try:
            self.logger.info("🚀 开始基因组共线性分析流程 | Starting genome collinearity analysis pipeline")
            self.logger.info(f"📊 样本数量 | Number of samples: {len(self.config.sample_list)}")
            self.logger.info(f"🔄 样本顺序 | Sample order: {' -> '.join(self.config.sample_list)}")
            
            if self.config.chromosome:
                self.logger.info(f"🧬 指定染色体 | Specified chromosome: {self.config.chromosome}")
            else:
                self.logger.info("🌍 分析范围 | Analysis scope: 全基因组 | Whole genome")
            
            # 检查依赖 | Check dependencies
            self.check_dependencies()
            
            # Step 1: 基因组比对 | Genome alignment
            self.logger.info("\n" + "="*20 + " 🧬 步骤1: 基因组比对 | Step 1: Genome Alignment " + "="*20)
            
            if self.config.skip_completed and not self.config.force_restart:
                alignment_files = self.aligner.check_alignments_completed()
                if len(alignment_files) == len(self.config.sample_list) - 1:
                    self.logger.info("⏭️ 跳过基因组比对步骤（已完成）| Skipping genome alignment step (completed)")
                else:
                    self.logger.info("🔄 需要重新运行部分或全部比对 | Need to re-run partial or all alignments")
                    alignment_files = self.aligner.run_alignments()
            else:
                alignment_files = self.aligner.run_alignments()
            
            if not alignment_files:
                raise RuntimeError("❌ 基因组比对失败 | Genome alignment failed")
            
            # Step 2: SyRI结构变异分析 | SyRI structural variation analysis
            self.logger.info("\n" + "="*20 + " 🔬 步骤2: SyRI结构变异分析 | Step 2: SyRI Structural Variation Analysis " + "="*20)
            
            if self.config.skip_completed and not self.config.force_restart:
                syri_files = self.syri_analyzer.check_syri_completed()
                if len(syri_files) == len(self.config.sample_list) - 1:
                    self.logger.info("⏭️ 跳过SyRI分析步骤（已完成）| Skipping SyRI analysis step (completed)")
                else:
                    self.logger.info("🔄 需要重新运行部分或全部SyRI分析 | Need to re-run partial or all SyRI analyses")
                    syri_files = self.syri_analyzer.run_syri_analysis(alignment_files)
            else:
                syri_files = self.syri_analyzer.run_syri_analysis(alignment_files)
            
            if not syri_files:
                raise RuntimeError("❌ SyRI结构变异分析失败 | SyRI structural variation analysis failed")
            
            # Step 3: PlotSR可视化 | PlotSR visualization
            self.logger.info("\n" + "="*20 + " 🎨 步骤3: PlotSR可视化 | Step 3: PlotSR Visualization " + "="*20)
            
            if self.config.skip_completed and not self.config.force_restart:
                if self.visualizer.check_visualization_completed():
                    self.logger.info("⏭️ 跳过可视化步骤（已完成）| Skipping visualization step (completed)")
                    visualization_success = True
                else:
                    self.logger.info("🔄 需要重新运行可视化 | Need to re-run visualization")
                    visualization_success = self.visualizer.create_plotsr_visualization(syri_files)
            else:
                visualization_success = self.visualizer.create_plotsr_visualization(syri_files)
            
            if not visualization_success:
                raise RuntimeError("❌ PlotSR可视化失败 | PlotSR visualization failed")
            
            # Step 4: 生成总结报告 | Generate summary report
            self.logger.info("\n" + "="*20 + " 📋 步骤4: 生成总结报告 | Step 4: Generate Summary Report " + "="*20)
            self.summary_generator.generate_summary_report(alignment_files, syri_files)
            
            # 分析完成 | Analysis completed
            self.logger.info("\n" + "="*20 + " 🎉 基因组共线性分析完成 | Genome Collinearity Analysis Completed " + "="*20)
            self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 基因组共线性可视化分析脚本 (模块化版本) | Genome Collinearity Visualization Analysis Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🎯 示例 | Examples:
  %(prog)s --sample-order samples.txt -o collinearity_results
  %(prog)s --sample-order samples.txt -c chr1 -o results
  %(prog)s --sample-order samples.txt -t 8 --plotsr-format pdf
  %(prog)s --sample-order samples.txt -c scaffold_1 --figure-size 16 10
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('--sample-order', required=True, 
                       help='📋 样本顺序文件路径 (制表符或空格分隔: 基因组路径 样本名称) | Sample order file path (tab or space separated: genome_path sample_name)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./collinearity_output', 
                       help='📁 输出目录 | Output directory')
    parser.add_argument('-c', '--chromosome', 
                       help='🧬 指定展示的染色体名称 (基于第一个基因组的染色体名称) | Specify chromosome name to display (based on first genome chromosome names)')
    parser.add_argument('-t', '--threads', type=int, default=4, 
                       help='⚡ 线程数 | Number of threads')
    
    # 工具路径参数 | Tool path parameters
    parser.add_argument('--minimap2-path', default='minimap2', 
                       help='🔧 minimap2软件路径 | minimap2 software path')
    parser.add_argument('--samtools-path', default='samtools', 
                       help='🔧 samtools软件路径 | samtools software path')
    parser.add_argument('--syri-path', default='syri', 
                       help='🔧 syri软件路径 | syri software path')
    parser.add_argument('--plotsr-path', default='plotsr', 
                       help='🔧 plotsr软件路径 | plotsr software path')
    
    # 分析参数 | Analysis parameters
    parser.add_argument('--minimap2-preset', default='asm5', 
                       help='⚙️ minimap2比对预设 | minimap2 alignment preset')
    parser.add_argument('--min-alignment-length', type=int, default=1000, 
                       help='📏 最小比对长度 | Minimum alignment length')
    
    # 可视化参数 | Visualization parameters
    parser.add_argument('--plotsr-format', default='png', choices=['png', 'pdf', 'svg'], 
                       help='🖼️ 输出图像格式 | Output image format')
    parser.add_argument('--figure-size', nargs=2, type=int, default=[12, 8], metavar=('WIDTH', 'HEIGHT'),
                       help='📐 图像尺寸 (宽度 高度) | Figure size (width height)')
    parser.add_argument('--line-width', type=float, default=1.5, 
                       help='📏 线条宽度 | Line width')
    parser.add_argument('--skip-synteny', action='store_true',
                       help='⏭️ 跳过共线性区域的绘制 | Skip plotting syntenic regions')
    # 流程控制参数 | Pipeline control parameters
    parser.add_argument('--skip-completed', action='store_true', default=True,
                    help='⏭️ 自动跳过已完成的步骤 | Automatically skip completed steps')
    parser.add_argument('--force-restart', action='store_true',
                    help='🔄 强制重新运行所有步骤 | Force restart all steps')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = GenomeCollinearityAnalyzer(
        sample_order_file=args.sample_order,
        output_dir=args.output,
        chromosome=args.chromosome,
        threads=args.threads,
        minimap2_path=args.minimap2_path,
        samtools_path=args.samtools_path,
        syri_path=args.syri_path,
        plotsr_path=args.plotsr_path,
        minimap2_preset=args.minimap2_preset,
        plotsr_format=args.plotsr_format,
        figure_width=args.figure_size[0],
        figure_height=args.figure_size[1],
        line_width=args.line_width,
        skip_synteny=args.skip_synteny,
        min_alignment_length=args.min_alignment_length,
        skip_completed=args.skip_completed,
        force_restart=args.force_restart
    )
    
    # 运行分析 | Run analysis
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
