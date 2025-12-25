"""
重复序列分析主程序模块 | Repeat Sequence Analysis Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import RepeatConfig
from .utils import RepeatLogger, CommandRunner, ProgressTracker, check_dependencies, get_genome_stats
from .tools import RepeatModelerRunner, LTRFinderRunner, LTRHarvestRunner, LTRRetrieverRunner, RepeatMaskerRunner, TESorterRunner

class RepeatAnalyzer:
    """🧬 重复序列分析主类 | Main Repeat Sequence Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = RepeatConfig(**kwargs)
        self.config.validate()
        
        # 自动检测工具路径 | Auto-detect tool paths
        self.config.detect_tool_paths()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = RepeatLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化进度跟踪 | Initialize progress tracking
        total_steps = self._calculate_total_steps()
        self.progress = ProgressTracker(self.logger, total_steps)
        
        # 初始化各个工具运行器 | Initialize tool runners
        self.repeatmodeler = RepeatModelerRunner(self.config, self.logger, self.cmd_runner)
        self.ltr_finder = LTRFinderRunner(self.config, self.logger, self.cmd_runner)
        self.ltr_harvest = LTRHarvestRunner(self.config, self.logger, self.cmd_runner)
        self.ltr_retriever = LTRRetrieverRunner(self.config, self.logger, self.cmd_runner)
        self.repeatmasker = RepeatMaskerRunner(self.config, self.logger, self.cmd_runner)
        self.tesorter = TESorterRunner(self.config, self.logger, self.cmd_runner)
        
        # 存储分析结果 | Store analysis results
        self.results = {
            'genome_stats': {},
            'repeat_libraries': [],
            'repeatmasker_files': {},
            'tesorter_files': {},
            'combined_library': None
        }
    
    def _calculate_total_steps(self) -> int:
        """📊 计算总步骤数 | Calculate total steps"""
        steps = 1  # 依赖检查 | Dependency check
        
        if not self.config.skip_modeler:
            steps += 1  # RepeatModeler
        
        if not self.config.skip_ltr:
            steps += 1  # LTR分析
        
        steps += 3  # 库合并、RepeatMasker、TEsorter
        return steps
    
    def check_dependencies(self):
        """🔍 检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def analyze_genome(self):
        """📊 分析基因组统计信息 | Analyze genome statistics"""
        self.results['genome_stats'] = get_genome_stats(self.config.genome_file, self.logger)
    
    def build_denovo_library(self) -> Path:
        """🏗️ 构建从头重复库 | Build de novo repeat library"""
        if self.config.skip_modeler:
            self.logger.info("⏭️ 跳过RepeatModeler步骤 | Skipping RepeatModeler step")
            return None
        
        self.progress.start_step()
        
        # 构建数据库并运行RepeatModeler | Build database and run RepeatModeler
        if self.repeatmodeler.build_database():
            if self.repeatmodeler.run_repeatmodeler():
                library_path = self.repeatmodeler.get_library_path()
                if library_path:
                    self.results['repeat_libraries'].append(('RepeatModeler', library_path))
                    return library_path
        
        self.logger.warning("⚠️ RepeatModeler库构建失败，继续其他步骤 | RepeatModeler library construction failed, continuing with other steps")
        return None
    
    def run_ltr_analysis(self) -> Path:
        """🧬 运行LTR序列分析 | Run LTR sequence analysis"""
        if self.config.skip_ltr:
            self.logger.info("⏭️ 跳过LTR分析步骤 | Skipping LTR analysis step")
            return None
        
        self.progress.start_step()
        
        ltr_finder_out = None
        ltrharvest_out = None
        
        # 运行LTR_FINDER | Run LTR_FINDER
        if self.ltr_finder.run_ltr_finder():
            ltr_finder_out = self.ltr_finder.get_output_path()
        
        # 运行LTRharvest | Run LTRharvest  
        if self.ltr_harvest.run_ltrharvest():
            ltrharvest_out = self.ltr_harvest.get_output_path()
        
        # 运行LTR_retriever整合结果 | Run LTR_retriever to integrate results
        if ltr_finder_out or ltrharvest_out:
            if self.ltr_retriever.run_ltr_retriever(ltr_finder_out, ltrharvest_out):
                library_path = self.ltr_retriever.get_library_path()
                if library_path:
                    self.results['repeat_libraries'].append(('LTR_retriever', library_path))
                    return library_path
        
        self.logger.warning("⚠️ LTR分析失败，继续其他步骤 | LTR analysis failed, continuing with other steps")
        return None
    
    def combine_libraries(self, denovo_lib: Path, ltr_lib: Path) -> Path:
        """📚 合并重复库 | Combine repeat libraries"""
        self.progress.start_step()
        
        available_libs = []
        if denovo_lib and denovo_lib.exists():
            available_libs.append(denovo_lib)
        if ltr_lib and ltr_lib.exists():
            available_libs.append(ltr_lib)
        
        if not available_libs:
            self.logger.error("❌ 没有可用的重复库文件 | No available repeat library files")
            return None
        
        # 如果只有一个库，直接使用 | If only one library, use it directly
        if len(available_libs) == 1:
            combined_lib = available_libs[0]
            self.logger.info(f"📚 使用单个重复库 | Using single repeat library: {combined_lib}")
        else:
            # 合并多个库 | Combine multiple libraries
            combined_lib = self.config.output_path / f"{self.config.base_name}_combined.fa"
            
            self.logger.info("📚 合并重复库文件 | Combining repeat library files")
            with open(combined_lib, 'w') as outfile:
                for i, lib_file in enumerate(available_libs):
                    self.logger.info(f"  📖 添加库文件 | Adding library file: {lib_file}")
                    with open(lib_file, 'r') as infile:
                        for line in infile:
                            if line.startswith('>'):
                                # 添加库来源标识 | Add library source identifier
                                header = line.strip()
                                source = ['RepeatModeler', 'LTR'][i] if i < 2 else f'Lib{i+1}'
                                line = f"{header}#{source}\n"
                            outfile.write(line)
        
        self.results['combined_library'] = combined_lib
        self.logger.info(f"✅ 合并库文件创建完成 | Combined library file created: {combined_lib}")
        return combined_lib
    
    def run_repeatmasker(self, library_path: Path) -> bool:
        """🎭 运行RepeatMasker | Run RepeatMasker"""
        self.progress.start_step()
        
        if self.repeatmasker.run_repeatmasker(library_path):
            self.results['repeatmasker_files'] = self.repeatmasker.get_output_files()
            return True
        
        return False
    
    def run_tesorter(self, library_path: Path) -> bool:
        """🏷️ 运行TEsorter | Run TEsorter"""
        self.progress.start_step()
        
        if self.tesorter.run_tesorter(library_path):
            self.results['tesorter_files'] = self.tesorter.get_output_files()
            return True
        
        return False
    
    def generate_summary(self):
        """📋 生成分析总结 | Generate analysis summary"""
        summary_file = self.config.output_path / "repeat_analysis_summary.txt"
        
        with open(summary_file, 'w') as f:
            f.write("🧬 重复序列分析总结报告 | Repeat Sequence Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 基因组信息 | Genome information
            f.write("📊 基因组统计信息 | Genome Statistics:\n")
            stats = self.results['genome_stats']
            f.write(f"  序列数 | Sequences: {stats.get('sequences', 'N/A'):,}\n")
            f.write(f"  总长度 | Total length: {stats.get('total_length', 'N/A'):,} bp\n")
            f.write(f"  N50: {stats.get('n50', 'N/A'):,} bp\n")
            f.write(f"  GC含量 | GC content: {stats.get('gc_content', 'N/A'):.2f}%\n\n")
            
            # 重复库信息 | Repeat library information
            f.write("📚 重复库文件 | Repeat Library Files:\n")
            for lib_name, lib_path in self.results['repeat_libraries']:
                f.write(f"  - {lib_name}: {lib_path}\n")
            if self.results['combined_library']:
                f.write(f"  - 合并库文件 | Combined library: {self.results['combined_library']}\n")
            f.write("\n")
            
            # RepeatMasker结果 | RepeatMasker results
            f.write("🎭 RepeatMasker输出文件 | RepeatMasker Output Files:\n")
            for file_type, file_path in self.results['repeatmasker_files'].items():
                f.write(f"  - {file_type}: {file_path}\n")
            f.write("\n")
            
            # TEsorter结果 | TEsorter results
            f.write("🏷️ TEsorter输出文件 | TEsorter Output Files:\n")
            for file_type, file_path in self.results['tesorter_files'].items():
                f.write(f"  - {file_type}: {file_path}\n")
            f.write("\n")
            
            f.write(f"📁 输出目录 | Output directory: {self.config.output_dir}\n")
        
        self.logger.info(f"📋 分析总结报告已生成 | Analysis summary report generated: {summary_file}")
    
    def run_analysis(self):
        """🚀 运行完整的重复序列分析流程 | Run complete repeat sequence analysis pipeline"""
        try:
            self.logger.info("🧬 开始重复序列分析流程 | Starting repeat sequence analysis pipeline")
            self.logger.info(f"📁 输入基因组 | Input genome: {self.config.genome_file}")
            self.logger.info(f"📂 输出目录 | Output directory: {self.config.output_dir}")
            self.logger.info(f"🧵 线程数 | Threads: {self.config.threads}")
            
            # 分析基因组统计信息 | Analyze genome statistics
            self.analyze_genome()
            
            # 检查依赖软件 | Check dependencies
            self.progress.start_step()
            self.check_dependencies()
            
            # 构建从头重复库 | Build de novo repeat library
            denovo_lib = self.build_denovo_library()
            
            # 运行LTR分析 | Run LTR analysis
            ltr_lib = self.run_ltr_analysis()
            
            # 合并重复库 | Combine repeat libraries
            combined_lib = self.combine_libraries(denovo_lib, ltr_lib)
            
            if not combined_lib:
                raise RuntimeError("❌ 无法创建重复库，分析终止 | Cannot create repeat library, analysis terminated")
            
            # 运行RepeatMasker | Run RepeatMasker
            if not self.run_repeatmasker(combined_lib):
                self.logger.error("❌ RepeatMasker运行失败 | RepeatMasker execution failed")
            
            # 运行TEsorter | Run TEsorter  
            if not self.run_tesorter(combined_lib):
                self.logger.error("❌ TEsorter运行失败 | TEsorter execution failed")
            
            # 生成分析总结 | Generate analysis summary
            self.generate_summary()
            
            # 完成分析 | Complete analysis
            self.progress.complete()
            self.logger.info("🎉 重复序列分析流程完成！| Repeat sequence analysis pipeline completed!")
            self.logger.info(f"📁 结果保存在 | Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"💥 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """🚀 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 重复序列分析脚本 (模块化版本) | Repeat Sequence Analysis Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🌟 示例 | Examples:
  %(prog)s -i genome.fasta -o repeat_results
  %(prog)s -i large_genome.fa -o results -t 64
  %(prog)s -i genome.fasta -o repeat_analysis --skip-modeler
  %(prog)s -i plant_genome.fa -o results --skip-ltr -t 32
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True, dest='genome_file',
                       help='🧬 输入基因组FASTA文件路径 | Input genome FASTA file path')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./repeat_output', dest='output_dir',
                       help='📂 输出目录 | Output directory')
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🧵 线程数 | Number of threads')
    
    # 流程控制参数 | Pipeline control parameters
    parser.add_argument('--skip-modeler', action='store_true',
                       help='⏭️ 跳过RepeatModeler步骤 | Skip RepeatModeler step')
    parser.add_argument('--skip-ltr', action='store_true',
                       help='⏭️ 跳过LTR分析步骤 | Skip LTR analysis step')
    
    # 工具路径参数 | Tool path parameters
    parser.add_argument('--repeatmodeler-path', default='RepeatModeler',
                       help='🔧 RepeatModeler程序路径 | RepeatModeler program path')
    parser.add_argument('--ltr-finder-path', default='ltr_finder',
                       help='🔧 LTR_FINDER程序路径 | LTR_FINDER program path')  
    parser.add_argument('--ltrharvest-path', default='gt ltrharvest',
                       help='🔧 LTRharvest程序路径 | LTRharvest program path')
    parser.add_argument('--ltr-retriever-path', default='LTR_retriever',
                       help='🔧 LTR_retriever程序路径 | LTR_retriever program path')
    parser.add_argument('--repeatmasker-path', default='RepeatMasker',
                       help='🔧 RepeatMasker程序路径 | RepeatMasker program path')
    parser.add_argument('--tesorter-path', default='TEsorter',
                       help='🔧 TEsorter程序路径 | TEsorter program path')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = RepeatAnalyzer(
        genome_file=args.genome_file,
        output_dir=args.output_dir,
        threads=args.threads,
        skip_modeler=args.skip_modeler,
        skip_ltr=args.skip_ltr,
        repeatmodeler_path=args.repeatmodeler_path,
        ltr_finder_path=args.ltr_finder_path,
        ltrharvest_path=args.ltrharvest_path,
        ltr_retriever_path=args.ltr_retriever_path,
        repeatmasker_path=args.repeatmasker_path,
        tesorter_path=args.tesorter_path
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()
