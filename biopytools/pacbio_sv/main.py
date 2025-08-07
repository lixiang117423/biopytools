"""
PacBio HiFi结构变异检测主程序 | PacBio HiFi Structural Variant Detection Main Program
"""

import sys
import argparse
from pathlib import Path

from .config import PacBioSVConfig
from .logger import PacBioSVLogger
from .utils import SVCommandRunner, check_dependencies
from .data_processing import InputValidator, SVCaller
from .analysis import SVFilter, SVMerger, LargeSVAnalyzer
from .reporter import SVReporter


class PacBioSVAnalyzer:
    """PacBio HiFi结构变异检测分析器 | PacBio HiFi SV Detection Analyzer"""
    
    def __init__(self, bam_file: str, ref_genome: str, sample_name: str, 
                 output_dir: str = "./SV_analysis", **kwargs):
        """初始化分析器 | Initialize analyzer"""
        
        self.config = PacBioSVConfig(
            bam_file=bam_file,
            ref_genome=ref_genome,
            sample_name=sample_name,
            output_dir=output_dir,
            **kwargs
        )
        
        # 设置日志 | Setup logging
        logger_manager = PacBioSVLogger(
            output_dir=self.config.logs_dir,
            log_name=f"{sample_name}_analysis.log"
        )
        self.logger = logger_manager.get_logger()
        
        # 初始化组件 | Initialize components
        self.runner = SVCommandRunner(str(self.config.output_dir), self.logger)
        self.validator = InputValidator(self.config, self.runner)
        self.caller = SVCaller(self.config, self.runner)
        self.filter = SVFilter(self.config, self.runner)
        self.merger = SVMerger(self.config, self.runner)
        self.large_analyzer = LargeSVAnalyzer(self.config, self.runner)
        self.reporter = SVReporter(self.config, self.runner)
    
    def run_analysis(self):
        """运行完整分析 | Run complete analysis"""
        try:
            self.logger.info(f"{'=' * 60}")
            self.logger.info("开始PacBio HiFi结构变异检测分析 | Starting PacBio HiFi SV detection analysis")
            self.logger.info(f"BAM文件 | BAM file: {self.config.bam_file}")
            self.logger.info(f"参考基因组 | Reference genome: {self.config.ref_genome}")
            self.logger.info(f"样本名称 | Sample name: {self.config.sample_name}")
            self.logger.info(f"{'=' * 60}")
            
            # 检查依赖 | Check dependencies
            check_dependencies(self.config, self.logger)
            
            # 执行分析流程 | Execute analysis pipeline
            steps = [
                (self.validator.validate_files, "输入文件验证"),
                (self.caller.run_pbsv, "pbsv检测"),
                (self.caller.run_sniffles2, "Sniffles2检测"),
                (self.caller.run_cutesv, "cuteSV检测"),
                (self.filter.filter_vcf_files, "结果过滤"),
                (self.merger.merge_results, "结果合并"),
                (self.large_analyzer.analyze_large_svs, "大片段SV分析"),
                (self.reporter.generate_report, "报告生成"),
                (self.reporter.generate_visualization_data, "可视化数据生成")
            ]
            
            for step_func, step_name in steps:
                self.logger.info(f"执行 | Executing: {step_name}")
                if not step_func():
                    raise Exception(f"{step_name}执行失败 | {step_name} execution failed")
            
            self.logger.info(f"{'=' * 60}")
            self.logger.info("PacBio HiFi结构变异检测分析完成! | PacBio HiFi SV detection analysis completed!")
            self.logger.info(f"结果保存在 | Results saved in: {self.config.output_dir}")
            self.logger.info(f"查看报告 | View report: {self.config.results_dir}/{self.config.sample_name}_SV_summary_report.txt")
            self.logger.info(f"{'=' * 60}")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='PacBio HiFi结构变异检测分析脚本 | PacBio HiFi Structural Variant Detection Script',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -b sample.bam -r ref.fa -s sample01 -o sv_results
  %(prog)s -b data.bam -r genome.fa -s sample02 -t 32 --min-sv-length 100
  %(prog)s -b aligned.bam -r reference.fa -s test --large-sv-threshold 5000
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-b', '--bam', required=True,
                       help='输入BAM文件路径 | Input BAM file path')
    parser.add_argument('-r', '--reference', required=True,
                       help='参考基因组文件路径 | Reference genome file path')
    parser.add_argument('-s', '--sample', required=True,
                       help='样本名称 | Sample name')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./SV_analysis',
                       help='输出目录 | Output directory')
    
    # 分析参数 | Analysis parameters
    parser.add_argument('-t', '--threads', type=int, default=16,
                       help='线程数 | Number of threads')
    parser.add_argument('--min-sv-length', type=int, default=50,
                       help='最小SV长度 | Minimum SV length')
    parser.add_argument('--min-support', type=int, default=3,
                       help='最小支持读数 | Minimum support reads')
    parser.add_argument('-q', '--quality', type=float, default=20.0,
                       help='质量阈值 | Quality threshold')
    
    # 大片段SV参数 | Large SV parameters
    parser.add_argument('--large-sv-threshold', type=int, default=1000,
                       help='大片段SV阈值 | Large SV threshold')
    parser.add_argument('--very-large-sv-threshold', type=int, default=10000,
                       help='超大片段SV阈值 | Very large SV threshold')
    
    # SURVIVOR参数 | SURVIVOR parameters
    parser.add_argument('--survivor-distance', type=int, default=1000,
                       help='SURVIVOR合并距离 | SURVIVOR merge distance')
    parser.add_argument('--survivor-min-callers', type=int, default=2,
                       help='SURVIVOR最小调用者数量 | SURVIVOR minimum callers')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = PacBioSVAnalyzer(
        bam_file=args.bam,
        ref_genome=args.reference,
        sample_name=args.sample,
        output_dir=args.output,
        threads=args.threads,
        min_sv_length=args.min_sv_length,
        min_support=args.min_support,
        quality_threshold=args.quality,
        large_sv_threshold=args.large_sv_threshold,
        very_large_sv_threshold=args.very_large_sv_threshold,
        survivor_distance=args.survivor_distance,
        survivor_min_callers=args.survivor_min_callers
    )
    
    analyzer.run_analysis()


if __name__ == "__main__":
    main()
