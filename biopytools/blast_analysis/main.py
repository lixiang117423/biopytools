# """
# 🧬 BLAST分析主程序模块
# """

# import argparse
# import sys
# from pathlib import Path
# from .config import BLASTConfig
# from .utils import BLASTLogger, CommandRunner, check_dependencies, SampleMapGenerator
# from .database import DatabaseManager
# from .blast_runner import BLASTRunner
# from .results_processor import ResultsProcessor
# from .statistics import StatisticsGenerator

# class BLASTAnalyzer:
#     """🧬 BLAST分析主类"""
    
#     def __init__(self, **kwargs):
#         # 初始化配置
#         self.config = BLASTConfig(**kwargs)
#         self.config.validate()
        
#         # 初始化日志
#         self.logger_manager = BLASTLogger(self.config.output_path)
#         self.logger = self.logger_manager.get_logger()
        
#         # 初始化命令执行器
#         self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
#         # 初始化各个处理器
#         self.sample_map_generator = SampleMapGenerator(self.config, self.logger)
#         self.database_manager = DatabaseManager(self.config, self.logger, self.cmd_runner)
#         self.blast_runner = BLASTRunner(self.config, self.logger, self.cmd_runner)
#         self.results_processor = ResultsProcessor(self.config, self.logger)
#         self.statistics_generator = StatisticsGenerator(self.config, self.logger)
    
#     def check_dependencies(self):
#         """检查依赖软件"""
#         return check_dependencies(self.config, self.logger)
    
#     def run_analysis(self):
#         """运行完整的BLAST分析流程"""
#         try:
#             self.logger.info("🧬 开始BLAST分析流程")
#             self.logger.info(f"{'=' * 80}")
            
#             # 步骤0: 生成或验证样品映射文件
#             self.logger.info("📋 步骤0: 准备样品映射文件")
#             sample_map_file = self.sample_map_generator.generate_sample_map()
            
#             # 步骤1: 检查依赖软件
#             self.logger.info("🔍 步骤1: 检查依赖软件")
#             self.check_dependencies()
            
#             # 步骤2: 创建BLAST数据库
#             self.logger.info("🗄️  步骤2: 创建BLAST数据库")
#             target_db_path = self.database_manager.create_target_database()
            
#             # 步骤3: 运行BLAST比对
#             self.logger.info("🔬 步骤3: 运行BLAST比对")
#             blast_results = self.blast_runner.run_blast_analysis(target_db_path)
            
#             # 步骤4: 处理结果
#             self.logger.info("📊 步骤4: 处理和汇总结果")
#             summary_file = self.results_processor.process_blast_results(blast_results)
            
#             # 步骤5: 排序结果
#             self.logger.info("📊 步骤5: 按覆盖度排序结果")
#             sorted_file = self.results_processor.sort_results_by_coverage(summary_file)
            
#             # 步骤6: 创建高质量结果
#             self.logger.info("🌟 步骤6: 筛选高质量结果")
#             high_quality_file = self.results_processor.create_high_quality_results(sorted_file)
            
#             # 步骤7: 生成统计报告
#             self.logger.info("📈 步骤7: 生成统计报告")
#             stats_file = self.statistics_generator.generate_statistics_report(sorted_file)
            
#             # 完成分析
#             self._log_completion_summary(sample_map_file, sorted_file, high_quality_file, stats_file)
            
#         except Exception as e:
#             self.logger.error(f"❌ 分析流程在执行过程中意外终止: {e}")
#             sys.exit(1)
    
#     def _log_completion_summary(self, sample_map_file: str, sorted_file: str, high_quality_file: str, stats_file: str):
#         """记录完成摘要"""
#         self.logger.info(f"{'=' * 80}")
#         self.logger.info("🎉 BLAST分析完成！")
#         self.logger.info(f"{'=' * 80}")
#         self.logger.info("📁 输出文件:")
        
#         if self.config.auto_generated_map:
#             self.logger.info(f"  📋 自动生成的样品映射文件: {sample_map_file}")
#         else:
#             self.logger.info(f"  📋 使用的样品映射文件: {sample_map_file}")
        
#         self.logger.info(f"  📊 汇总结果（已排序）: {sorted_file}")
        
#         if high_quality_file:
#             self.logger.info(f"  🌟 高质量结果: {high_quality_file}")
        
#         if stats_file:
#             self.logger.info(f"  📈 统计报告: {stats_file}")
        
#         self.logger.info(f"  📁 详细结果目录: {self.config.output_dir}")
#         self.logger.info(f"  📝 日志文件: {self.logger_manager.log_file}")
#         self.logger.info(f"{'=' * 80}")
#         self.logger.info(f"✅ 结果保存在: {self.config.output_dir}")
        
#         if self.config.auto_generated_map:
#             self.logger.info(f"💡 提示：如需修改样品名称，请编辑 {sample_map_file} 文件后重新运行")

# def main():
#     """主函数"""
#     parser = argparse.ArgumentParser(
#         description='🧬 BLAST比对分析脚本 (模块化版本)',
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         epilog="""
# 🧬 示例:
#   # 使用目录（自动生成样品映射文件）
#   %(prog)s -i sequences/ -t nlr_genes.fa -o blast_results
  
#   # 使用单个文件
#   %(prog)s -i single_sequence.fa -t nlr_genes.fa -o blast_results
  
#   # 使用现有样品映射文件
#   %(prog)s -s sample_map.txt -t nlr_genes.fa -o blast_results
  
#   # BLASTP分析
#   %(prog)s -i proteins/ -t targets.fa -o results --blast-type blastp --threads 32
  
#   # 高质量筛选
#   %(prog)s -i sequences/ -t nlr.fa -o results --min-identity 80 --min-coverage 70

# 📋 样品映射文件格式:
#   /path/to/sequence1.fa<TAB>Sample1
#   /path/to/sequence2.fa<TAB>Sample2
#   # 注释行以#开头

# 💡 工作流程:
#   1. 如果只提供-i：自动扫描文件，生成样品映射文件
#   2. 如果只提供-s：直接使用现有样品映射文件
#   3. 如果同时提供-i和-s：优先使用-s，忽略-i
#         """
#     )
    
#     # 输入数据参数
#     parser.add_argument('-i', '--input', 
#                        help='📁 输入文件或目录路径（支持单个文件或包含多个文件的目录）')
#     parser.add_argument('-s', '--sample-map-file', 
#                        help='🧪 样品映射文件，格式：文件路径<TAB>样品名称')
    
#     # 必需参数
#     parser.add_argument('-t', '--target-file', required=True, 
#                        help='🎯 目标基因序列文件')
    
#     # 基本参数
#     parser.add_argument('-o', '--output-dir', default='./blast_output', 
#                        help='📁 输出目录')
#     parser.add_argument('--blast-type', choices=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'], 
#                        default='blastn', help='🔬 BLAST程序类型')
    
#     # BLAST参数
#     parser.add_argument('-e', '--evalue', type=float, default=1e-5, 
#                        help='🔢 E-value阈值')
#     parser.add_argument('--max-target-seqs', type=int, default=10, 
#                        help='🎯 最大目标序列数')
#     parser.add_argument('--word-size', type=int, default=11, 
#                        help='📏 词大小 (仅适用于blastn/tblastx)')
#     parser.add_argument('--threads', '-j', type=int, default=88, 
#                        help='⚡ 线程数')
    
#     # 文件格式参数
#     parser.add_argument('--input-suffix', default='*.fa', 
#                        help='📁 输入文件后缀模式（当输入为目录时使用）')
#     parser.add_argument('--target-db-type', choices=['nucl', 'prot'], default='nucl', 
#                        help='🗄️ 目标数据库类型')
    
#     # 过滤参数
#     parser.add_argument('--min-identity', type=float, default=70.0, 
#                        help='📊 最小序列相似度 (%%)')
#     parser.add_argument('--min-coverage', type=float, default=50.0, 
#                        help='📐 最小覆盖度 (%%)')
#     parser.add_argument('--high-quality-evalue', type=float, default=1e-10, 
#                        help='🌟 高质量比对E-value阈值')
    
#     # 样品信息参数
#     parser.add_argument('--auto-detect-samples', action='store_true', default=True, 
#                        help='🔍 自动检测样品名称（仅当使用-i时有效）')
#     parser.add_argument('--sample-name-pattern', default=r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$', 
#                        help='🔍 样品名称提取正则表达式')
    
#     # 工具路径参数
#     parser.add_argument('--makeblastdb-path', default='makeblastdb', 
#                        help='🗄️ makeblastdb程序路径')
#     parser.add_argument('--blastn-path', default='blastn', 
#                        help='🔬 blastn程序路径')
#     parser.add_argument('--blastp-path', default='blastp', 
#                        help='🔬 blastp程序路径')
#     parser.add_argument('--blastx-path', default='blastx', 
#                        help='🔬 blastx程序路径')
#     parser.add_argument('--tblastn-path', default='tblastn', 
#                        help='🔬 tblastn程序路径')
#     parser.add_argument('--tblastx-path', default='tblastx', 
#                        help='🔬 tblastx程序路径')
    
#     args = parser.parse_args()
    
#     # 验证输入参数
#     if not args.input and not args.sample_map_file:
#         parser.error("必须指定输入路径(-i)或样品映射文件(-s)中的一个")
    
#     # 创建分析器并运行
#     analyzer = BLASTAnalyzer(
#         input_path=args.input,
#         target_file=args.target_file,
#         output_dir=args.output_dir,
#         sample_map_file=args.sample_map_file,
#         blast_type=args.blast_type,
#         evalue=args.evalue,
#         max_target_seqs=args.max_target_seqs,
#         word_size=args.word_size,
#         threads=args.threads,
#         input_suffix=args.input_suffix,
#         target_db_type=args.target_db_type,
#         min_identity=args.min_identity,
#         min_coverage=args.min_coverage,
#         high_quality_evalue=args.high_quality_evalue,
#         auto_detect_samples=args.auto_detect_samples,
#         sample_name_pattern=args.sample_name_pattern,
#         makeblastdb_path=args.makeblastdb_path,
#         blastn_path=args.blastn_path,
#         blastp_path=args.blastp_path,
#         blastx_path=args.blastx_path,
#         tblastn_path=args.tblastn_path,
#         tblastx_path=args.tblastx_path
#     )
    
#     analyzer.run_analysis()

# if __name__ == "__main__":
#     main()


"""
🧬 BLAST分析主程序模块
"""

import argparse
import sys
from pathlib import Path
from .config import BLASTConfig
from .utils import BLASTLogger, CommandRunner, check_dependencies, SampleMapGenerator
from .database import DatabaseManager
from .blast_runner import BLASTRunner
from .results_processor import ResultsProcessor
from .statistics import StatisticsGenerator
from .alignment_visualizer import AlignmentVisualizer

class BLASTAnalyzer:
    """🧬 BLAST分析主类"""
    
    def __init__(self, **kwargs):
        # 初始化配置
        self.config = BLASTConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志
        self.logger_manager = BLASTLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器
        self.sample_map_generator = SampleMapGenerator(self.config, self.logger)
        self.database_manager = DatabaseManager(self.config, self.logger, self.cmd_runner)
        self.blast_runner = BLASTRunner(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
        self.statistics_generator = StatisticsGenerator(self.config, self.logger)
        self.alignment_visualizer = AlignmentVisualizer(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的BLAST分析流程"""
        try:
            self.logger.info("🧬 开始BLAST分析流程")
            self.logger.info(f"{'=' * 80}")
            
            # 步骤0: 生成或验证样品映射文件
            self.logger.info("📋 步骤0: 准备样品映射文件")
            sample_map_file = self.sample_map_generator.generate_sample_map()
            
            # 步骤1: 检查依赖软件
            self.logger.info("🔍 步骤1: 检查依赖软件")
            self.check_dependencies()
            
            # 步骤2: 创建BLAST数据库
            self.logger.info("🗄️  步骤2: 创建BLAST数据库")
            target_db_path = self.database_manager.create_target_database()
            
            # 步骤3: 运行BLAST比对
            self.logger.info("🔬 步骤3: 运行BLAST比对")
            blast_results = self.blast_runner.run_blast_analysis(target_db_path)
            
            # 步骤4: 处理结果
            self.logger.info("📊 步骤4: 处理和汇总结果")
            summary_file = self.results_processor.process_blast_results(blast_results)
            
            # 步骤5: 排序结果
            self.logger.info("📊 步骤5: 按覆盖度排序结果")
            sorted_file = self.results_processor.sort_results_by_coverage(summary_file)
            
            # 步骤6: 创建高质量结果
            self.logger.info("🌟 步骤6: 筛选高质量结果")
            high_quality_file = self.results_processor.create_high_quality_results(sorted_file)
            
            # 步骤7: 生成统计报告
            self.logger.info("📈 步骤7: 生成统计报告")
            stats_file = self.statistics_generator.generate_statistics_report(sorted_file)
            
            # 步骤8: 生成序列比对可视化（新增，可选）
            alignment_files = None
            if self.config.alignment_output != 'none':
                self.logger.info("🧬 步骤8: 生成序列比对可视化")
                alignment_files = self.alignment_visualizer.generate_visualizations(blast_results)
            
            # 完成分析
            self._log_completion_summary(
                sample_map_file, sorted_file, high_quality_file, 
                stats_file, alignment_files
            )
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程在执行过程中意外终止: {e}")
            sys.exit(1)
    
    def _log_completion_summary(self, sample_map_file: str, sorted_file: str, 
                               high_quality_file: str, stats_file: str, 
                               alignment_files: dict):
        """记录完成摘要"""
        self.logger.info(f"{'=' * 80}")
        self.logger.info("🎉 BLAST分析完成！")
        self.logger.info(f"{'=' * 80}")
        self.logger.info("📁 输出文件:")
        
        if self.config.auto_generated_map:
            self.logger.info(f"  📋 自动生成的样品映射文件: {sample_map_file}")
        else:
            self.logger.info(f"  📋 使用的样品映射文件: {sample_map_file}")
        
        self.logger.info(f"  📊 汇总结果（已排序）: {sorted_file}")
        
        if high_quality_file:
            self.logger.info(f"  🌟 高质量结果: {high_quality_file}")
        
        if stats_file:
            self.logger.info(f"  📈 统计报告: {stats_file}")
        
        # 显示比对可视化文件
        if alignment_files:
            if 'text' in alignment_files:
                self.logger.info(f"  📝 文本比对文件: {self.config.output_path / self.config.alignment_output_dir / 'text'}")
            if 'html' in alignment_files:
                html_index = alignment_files['html']['index']
                self.logger.info(f"  🌐 HTML比对文件: {html_index}")
                self.logger.info(f"     💡 在浏览器中打开上述文件以查看交互式比对")
        
        self.logger.info(f"  📁 详细结果目录: {self.config.output_dir}")
        self.logger.info(f"  📝 日志文件: {self.logger_manager.log_file}")
        self.logger.info(f"{'=' * 80}")
        self.logger.info(f"✅ 结果保存在: {self.config.output_dir}")
        
        if self.config.auto_generated_map:
            self.logger.info(f"💡 提示：如需修改样品名称，请编辑 {sample_map_file} 文件后重新运行")

def main():
    """主函数"""
    parser = argparse.ArgumentParser(
        description='🧬 BLAST比对分析脚本 (模块化版本 v1.0 - 支持比对可视化)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🧬 示例:
  # 基本用法（使用目录）
  %(prog)s -i sequences/ -t nlr_genes.fa -o blast_results
  
  # 生成文本格式比对可视化
  %(prog)s -i sequences/ -t nlr_genes.fa -o results --alignment-output text
  
  # 生成HTML格式比对可视化
  %(prog)s -i sequences/ -t nlr_genes.fa -o results --alignment-output html
  
  # 同时生成文本和HTML格式
  %(prog)s -i sequences/ -t nlr_genes.fa -o results --alignment-output both
  
  # 带筛选条件的比对可视化
  %(prog)s -i sequences/ -t nlr.fa -o results \\
      --alignment-output both \\
      --alignment-min-identity 80 \\
      --alignment-min-coverage 70
  
  # 使用现有样品映射文件
  %(prog)s -s sample_map.txt -t nlr_genes.fa -o blast_results --alignment-output html

📋 样品映射文件格式:
  /path/to/sequence1.fa<TAB>Sample1
  /path/to/sequence2.fa<TAB>Sample2
  # 注释行以#开头

🎨 比对可视化说明:
  --alignment-output 参数:
    none: 不生成比对可视化（默认）
    text: 生成纯文本格式比对
    html: 生成交互式HTML格式比对
    both: 同时生成文本和HTML格式

💡 工作流程:
  1. 如果只提供-i：自动扫描文件，生成样品映射文件
  2. 如果只提供-s：直接使用现有样品映射文件
  3. 如果同时提供-i和-s：优先使用-s，忽略-i
        """
    )
    
    # 输入数据参数
    parser.add_argument('-i', '--input', 
                       help='📁 输入文件或目录路径（支持单个文件或包含多个文件的目录）')
    parser.add_argument('-s', '--sample-map-file', 
                       help='🧪 样品映射文件，格式：文件路径<TAB>样品名称')
    
    # 必需参数
    parser.add_argument('-t', '--target-file', required=True, 
                       help='🎯 目标基因序列文件')
    
    # 基本参数
    parser.add_argument('-o', '--output-dir', default='./blast_output', 
                       help='📁 输出目录')
    parser.add_argument('--blast-type', choices=['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx'], 
                       default='blastn', help='🔬 BLAST程序类型')
    
    # BLAST参数
    parser.add_argument('-e', '--evalue', type=float, default=1e-5, 
                       help='🔢 E-value阈值')
    parser.add_argument('--max-target-seqs', type=int, default=10, 
                       help='🎯 最大目标序列数')
    parser.add_argument('--word-size', type=int, default=11, 
                       help='📏 词大小 (仅适用于blastn/tblastx)')
    parser.add_argument('--threads', '-j', type=int, default=88, 
                       help='⚡ 线程数')
    
    # 文件格式参数
    parser.add_argument('--input-suffix', default='*.fa', 
                       help='📁 输入文件后缀模式（当输入为目录时使用）')
    parser.add_argument('--target-db-type', choices=['nucl', 'prot'], default='nucl', 
                       help='🗄️ 目标数据库类型')
    
    # 过滤参数
    parser.add_argument('--min-identity', type=float, default=70.0, 
                       help='📊 最小序列相似度 (%%)')
    parser.add_argument('--min-coverage', type=float, default=50.0, 
                       help='📐 最小覆盖度 (%%)')
    parser.add_argument('--high-quality-evalue', type=float, default=1e-10, 
                       help='🌟 高质量比对E-value阈值')
    
    # 比对可视化参数（新增）
    parser.add_argument('--alignment-output', choices=['none', 'text', 'html', 'both'], 
                       default='both', 
                       help='🧬 比对可视化输出格式：none=不生成, text=纯文本, html=交互式HTML, both=两种都生成')
    parser.add_argument('--alignment-width', type=int, default=80, 
                       help='📏 比对每行显示的字符数')
    parser.add_argument('--alignment-min-identity', type=float, default=0.0, 
                       help='🔍 比对可视化最小相似度过滤 (0表示不限制)')
    parser.add_argument('--alignment-min-coverage', type=float, default=0.0, 
                       help='🔍 比对可视化最小覆盖度过滤 (0表示不限制)')
    parser.add_argument('--alignment-max-per-sample', type=int, default=100, 
                       help='📊 每个样品最多显示的比对数')
    parser.add_argument('--html-theme', choices=['modern', 'classic', 'dark'], 
                       default='modern', 
                       help='🎨 HTML主题样式')
    
    # 样品信息参数
    parser.add_argument('--auto-detect-samples', action='store_true', default=True, 
                       help='🔍 自动检测样品名称（仅当使用-i时有效）')
    parser.add_argument('--sample-name-pattern', default=r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$', 
                       help='🔍 样品名称提取正则表达式')
    
    # 工具路径参数
    parser.add_argument('--makeblastdb-path', default='makeblastdb', 
                       help='🗄️ makeblastdb程序路径')
    parser.add_argument('--blastn-path', default='blastn', 
                       help='🔬 blastn程序路径')
    parser.add_argument('--blastp-path', default='blastp', 
                       help='🔬 blastp程序路径')
    parser.add_argument('--blastx-path', default='blastx', 
                       help='🔬 blastx程序路径')
    parser.add_argument('--tblastn-path', default='tblastn', 
                       help='🔬 tblastn程序路径')
    parser.add_argument('--tblastx-path', default='tblastx', 
                       help='🔬 tblastx程序路径')
    
    args = parser.parse_args()
    
    # 验证输入参数
    if not args.input and not args.sample_map_file:
        parser.error("必须指定输入路径(-i)或样品映射文件(-s)中的一个")
    
    # 创建分析器并运行
    analyzer = BLASTAnalyzer(
        input_path=args.input,
        target_file=args.target_file,
        output_dir=args.output_dir,
        sample_map_file=args.sample_map_file,
        blast_type=args.blast_type,
        evalue=args.evalue,
        max_target_seqs=args.max_target_seqs,
        word_size=args.word_size,
        threads=args.threads,
        input_suffix=args.input_suffix,
        target_db_type=args.target_db_type,
        min_identity=args.min_identity,
        min_coverage=args.min_coverage,
        high_quality_evalue=args.high_quality_evalue,
        alignment_output=args.alignment_output,
        alignment_width=args.alignment_width,
        alignment_min_identity=args.alignment_min_identity,
        alignment_min_coverage=args.alignment_min_coverage,
        alignment_max_per_sample=args.alignment_max_per_sample,
        html_theme=args.html_theme,
        auto_detect_samples=args.auto_detect_samples,
        sample_name_pattern=args.sample_name_pattern,
        makeblastdb_path=args.makeblastdb_path,
        blastn_path=args.blastn_path,
        blastp_path=args.blastp_path,
        blastx_path=args.blastx_path,
        tblastn_path=args.tblastn_path,
        tblastx_path=args.tblastx_path
    )
    
    analyzer.run_analysis()

if __name__ == "__main__":
    main()