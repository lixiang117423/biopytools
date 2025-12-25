"""
🏭 Ka/Ks Calculator主分析流水线
功能: 协调各个模块完成完整的Ka/Ks分析流程 | Main analysis pipeline coordinating all modules
"""

import os
import sys
import argparse
import tempfile
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Tuple

from .config import KaKsConfig
from .logger import Logger
from .validator import SequenceValidator
from .calculator import KaKsCalculator
from .processor import ResultProcessor

class KaKsAnalyzer:
    """🏭 Ka/Ks主分析器 | Main Ka/Ks analyzer"""
    
    def __init__(self, fasta1: str, fasta2: str, pairs: str, output_dir: str,
                 method: str = None, kaks_path: str = "KaKs_Calculator",
                 verbose: bool = False, temp_dir: str = None, keep_temp: bool = False):
        """
        🏗️ 初始化分析器 | Initialize analyzer
        
        Args:
            fasta1: 第一个FASTA文件路径 | First FASTA file path
            fasta2: 第二个FASTA文件路径 | Second FASTA file path
            pairs: 配对文件路径 | Pairs file path
            output_dir: 输出目录 | Output directory
            method: 计算方法 | Calculation method
            kaks_path: KaKs_Calculator路径 | KaKs_Calculator path
            verbose: 详细模式 | Verbose mode
            temp_dir: 临时目录 | Temporary directory
            keep_temp: 保留临时文件 | Keep temporary files
        """
        self.fasta1 = fasta1
        self.fasta2 = fasta2
        self.pairs = pairs
        self.output_dir = output_dir
        self.method = method or KaKsConfig.DEFAULT_METHOD
        self.kaks_path = kaks_path
        self.verbose = verbose
        self.temp_dir = temp_dir
        self.keep_temp = keep_temp
        
        # 📂 创建输出目录 | Create output directory
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        
        # 📊 初始化日志器 | Initialize logger
        log_file = os.path.join(output_dir, KaKsConfig.OUTPUT_FILES['log'])
        self.logger = Logger(log_file, verbose, 'KaKsAnalyzer')
        
        # 🔧 初始化组件 | Initialize components
        self.config = KaKsConfig()
        self.validator = SequenceValidator(self.logger)
        self.calculator = KaKsCalculator(self.logger, kaks_path)
        self.processor = ResultProcessor(self.logger)
        
        # 🗂️ 临时目录管理 | Temporary directory management
        self._temp_dir_created = None
    
    def run_analysis(self):
        """🚀 执行完整的Ka/Ks分析流程 | Execute complete Ka/Ks analysis pipeline"""
        try:
            start_time = datetime.now()
            self.logger.separator("Ka/Ks Calculator 2.0 Analysis Pipeline")
            self.logger.info("开始Ka/Ks分析流程 | Starting Ka/Ks analysis pipeline", "🚀")
            
            # 🗂️ 设置临时目录 | Setup temporary directory
            self._setup_temp_directory()
            
            # ✅ 验证输入文件 | Validate input files
            seq1_dict, seq2_dict, pairs_list = self._validate_inputs()
            
            # 🧮 运行Ka/Ks计算 | Run Ka/Ks calculation
            results_df = self._run_calculation(seq1_dict, seq2_dict, pairs_list)
            
            # 📈 处理和保存结果 | Process and save results
            self._process_and_save_results(results_df)
            
            # 🧹 清理临时文件 | Cleanup temporary files
            self._cleanup()
            
            # ⏱️ 计算运行时间 | Calculate runtime
            runtime = datetime.now() - start_time
            self.logger.separator()
            self.logger.success(f"分析完成！用时: {runtime} | Analysis completed! Runtime: {runtime}")
            
        except Exception as e:
            self.logger.error(f"分析流程失败 | Analysis pipeline failed: {e}")
            self._cleanup()
            raise
    
    def _setup_temp_directory(self):
        """🗂️ 设置临时工作目录 | Setup temporary working directory"""
        if self.temp_dir:
            self._temp_dir_created = self.temp_dir
            Path(self._temp_dir_created).mkdir(parents=True, exist_ok=True)
        else:
            self._temp_dir_created = tempfile.mkdtemp(prefix="kakscalc_")
        
        self.logger.info(f"使用临时目录 | Using temporary directory: {self._temp_dir_created}", "🗂️")
    
    def _validate_inputs(self) -> Tuple[Dict[str, str], Dict[str, str], List[Tuple[str, str, str]]]:
        """✅ 验证所有输入文件 | Validate all input files"""
        self.logger.separator("输入文件验证 | Input File Validation")
        
        # 🧬 验证FASTA文件 | Validate FASTA files
        self.logger.info(f"验证第一个FASTA文件 | Validating first FASTA file: {self.fasta1}")
        valid1, seq1_dict = self.validator.validate_fasta_file(self.fasta1)
        if not valid1:
            raise ValueError(f"第一个FASTA文件无效 | Invalid first FASTA file: {self.fasta1}")
        
        self.logger.info(f"验证第二个FASTA文件 | Validating second FASTA file: {self.fasta2}")
        valid2, seq2_dict = self.validator.validate_fasta_file(self.fasta2)
        if not valid2:
            raise ValueError(f"第二个FASTA文件无效 | Invalid second FASTA file: {self.fasta2}")
        
        # 🔗 验证配对文件 | Validate pair file
        self.logger.info(f"验证配对文件 | Validating pair file: {self.pairs}")
        valid_pairs, pairs_list = self.validator.validate_pair_file(
            self.pairs, set(seq1_dict.keys()), set(seq2_dict.keys())
        )
        if not valid_pairs or not pairs_list:
            raise ValueError(f"配对文件无效或为空 | Invalid or empty pair file: {self.pairs}")
        
        # 📊 检查批处理大小 | Check batch size
        if len(pairs_list) > self.config.MAX_BATCH_SIZE:
            self.logger.warning(f"配对数量 ({len(pairs_list)}) 超过推荐最大值 ({self.config.MAX_BATCH_SIZE})")
            self.logger.warning("大批量分析可能需要较长时间 | Large batch analysis may take significant time")
        
        self.logger.success("所有输入文件验证通过 | All input files validated successfully")
        return seq1_dict, seq2_dict, pairs_list
    
    def _run_calculation(self, seq1_dict: Dict[str, str], seq2_dict: Dict[str, str], 
                        pairs_list: List[Tuple[str, str, str]]):
        """🧮 执行Ka/Ks计算 | Execute Ka/Ks calculation"""
        self.logger.separator("Ka/Ks计算 | Ka/Ks Calculation")
        
        # 📝 准备输入文件 | Prepare input file
        input_file = self.calculator.prepare_input_file(
            seq1_dict, seq2_dict, pairs_list, self._temp_dir_created
        )
        
        # 🚀 运行计算 | Run calculation
        output_file = self.calculator.run_calculation(
            input_file, self.method, self._temp_dir_created
        )
        
        # 📊 解析结果 | Parse results
        results_df = self.processor.parse_results(output_file)
        
        return results_df
    
    def _process_and_save_results(self, results_df):
        """📈 处理并保存结果 | Process and save results"""
        self.logger.separator("结果处理 | Result Processing")
        
        # 📊 生成汇总统计 | Generate summary statistics
        stats = self.processor.generate_summary_stats(results_df)
        
        # 💾 保存结果 | Save results
        self.processor.save_results(results_df, stats, self.output_dir)
        
        # 📋 打印汇总到控制台 | Print summary to console
        self._print_analysis_summary(stats)
    
    def _print_analysis_summary(self, stats: Dict):
        """📋 打印分析汇总到控制台 | Print analysis summary to console"""
        self.logger.separator("分析汇总 | Analysis Summary")
        
        analysis_info = stats.get('analysis_info', {})
        omega_stats = stats.get('omega_statistics', {})
        selection_dist = stats.get('selection_distribution', {})
        
        print(f"  📊 总序列对数 | Total sequence pairs: {analysis_info.get('total_pairs', 'N/A')}")
        print(f"  ✅ 成功计算 | Successful calculations: {analysis_info.get('successful_calculations', 'N/A')}")
        print(f"  ❌ 失败计算 | Failed calculations: {analysis_info.get('failed_calculations', 'N/A')}")
        print(f"  📈 成功率 | Success rate: {analysis_info.get('success_rate', 0):.1%}")
        
        if omega_stats.get('count', 0) > 0:
            print(f"  🎯 平均ω (Ka/Ks) | Mean ω: {omega_stats.get('mean', 0):.4f}")
            print(f"  📊 中位数ω | Median ω: {omega_stats.get('median', 0):.4f}")
            print(f"  📏 ω范围 | ω range: {omega_stats.get('min', 0):.4f} - {omega_stats.get('max', 0):.4f}")
        
        if selection_dist:
            print("  🔍 选择压力分布 | Selection pressure distribution:")
            total_successful = analysis_info.get('successful_calculations', 1)
            for sel_type, count in selection_dist.items():
                percentage = (count / total_successful) * 100
                print(f"    {sel_type}: {count} ({percentage:.1f}%)")
        
        # 输出文件信息 | Output file information
        print("  📁 生成的输出文件 | Generated output files:")
        print(f"    📊 kaks_detailed.csv - CSV详细结果 | CSV detailed results")
        print(f"    📋 kaks_detailed.tsv - TSV详细结果 | TSV detailed results")
        print(f"    📈 kaks_summary.xlsx - Excel汇总文件 | Excel summary file")
        print(f"    💾 summary_stats.json - JSON统计数据 | JSON statistics")
        
        # 生物学解释 | Biological interpretation
        bio_interp = stats.get('biological_interpretation', {})
        if bio_interp:
            print("  🧠 生物学解释 | Biological interpretation:")
            for aspect, interpretation in bio_interp.items():
                print(f"    {interpretation}")
    
    def _cleanup(self):
        """🧹 清理临时文件 | Clean up temporary files"""
        if self._temp_dir_created and not self.keep_temp:
            try:
                shutil.rmtree(self._temp_dir_created)
                self.logger.info("临时文件已清理 | Temporary files cleaned up", "🧹")
            except Exception as e:
                self.logger.warning(f"清理临时文件失败 | Failed to clean up temporary files: {e}")
        elif self.keep_temp:
            self.logger.info(f"保留临时文件 | Keeping temporary files: {self._temp_dir_created}", "🗂️")

def create_argument_parser() -> argparse.ArgumentParser:
    """🔧 创建命令行参数解析器 | Create command line argument parser"""
    
    parser = argparse.ArgumentParser(
        description="🧬 Ka/Ks Calculator 2.0 Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
📋 使用示例 | Examples:
  %(prog)s -1 species1.fasta -2 species2.fasta -p pairs.txt -o results/
  %(prog)s --fasta1 human.fa --fasta2 mouse.fa --pairs orthologs.tsv --output analysis/ --method gamma-MYN
  %(prog)s -1 seqs1.fa -2 seqs2.fa -p pairs.csv -o out/ -m YN --verbose
  %(prog)s --fasta1 cds1.fa --fasta2 cds2.fa --pairs pairs.txt --output results/ --kaks-path /usr/local/bin/KaKs_Calculator

🧮 支持的计算方法 | Supported methods:
  GMYN (推荐), MYN, YN, NG, LWL, LPB, MLWL, MLPB, GY, MS, MA, GNG, GLWL, GLPB, GMLWL, GMLPB, GYN

📄 输入文件格式 | Input file formats:
  - FASTA文件: 标准CDS序列格式 | Standard CDS sequence format
  - 配对文件: TSV/CSV格式，包含序列ID对应关系 | TSV/CSV format with sequence ID pairs

📄 输出文件格式 | Output file formats:
  - Excel汇总文件: 多sheet Excel格式 | Multi-sheet Excel format
  - CSV详细结果: 逗号分隔值格式 | Comma-separated values format  
  - TSV详细结果: 制表符分隔值格式 | Tab-separated values format
  - JSON统计文件: 统计数据JSON格式 | Statistical data in JSON format
        """
    )
    
    # 📁 必需输入文件 | Required input files
    parser.add_argument('-1', '--fasta1', '--species1',
                       required=True, type=str, metavar='FILE',
                       help='🧬 第一个FASTA文件 (物种1 CDS序列) | First FASTA file (species 1 CDS sequences)')
    
    parser.add_argument('-2', '--fasta2', '--species2', 
                       required=True, type=str, metavar='FILE',
                       help='🧬 第二个FASTA文件 (物种2 CDS序列) | Second FASTA file (species 2 CDS sequences)')
    
    parser.add_argument('-p', '--pairs', '--pair-file',
                       required=True, type=str, metavar='FILE', 
                       help='🔗 序列配对文件 (TSV/CSV格式) | Sequence pair file (TSV/CSV format)')
    
    # 📂 输出目录 | Output directory
    parser.add_argument('-o', '--output', '--output-dir',
                       required=True, type=str, metavar='DIR',
                       help='📂 输出目录 | Output directory for results')
    
    # 🔧 分析参数 | Analysis parameters
    parser.add_argument('-m', '--method', '--calc-method',
                       default=KaKsConfig.DEFAULT_METHOD, 
                       choices=KaKsConfig.SUPPORTED_METHODS,
                       help=f'🧮 计算方法 (默认: {KaKsConfig.DEFAULT_METHOD}) | Calculation method (default: {KaKsConfig.DEFAULT_METHOD})')
    
    parser.add_argument('--kaks-path', '--calculator-path',
                       default='KaKs_Calculator', type=str, metavar='PATH',
                       help='🛠️ KaKs_Calculator可执行文件路径 | Path to KaKs_Calculator executable')
    
    # 🎛️ 运行选项 | Runtime options
    parser.add_argument('-v', '--verbose', '--debug',
                       action='store_true',
                       help='🔍 启用详细日志记录 | Enable verbose logging')
    
    parser.add_argument('--temp-dir', '--tmp-dir',
                       type=str, metavar='DIR',
                       help='🗂️ 自定义临时目录 | Custom temporary directory')
    
    parser.add_argument('--keep-temp', '--no-cleanup',
                       action='store_true',
                       help='🗑️ 保留临时文件 (用于调试) | Keep temporary files (for debugging)')
    
    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')
    
    return parser

def main():
    """🏁 主入口点 | Main entry point"""
    try:
        # 🔧 解析命令行参数 | Parse command line arguments
        parser = create_argument_parser()
        args = parser.parse_args()
        
        # 🏭 初始化并运行分析器 | Initialize and run analyzer
        analyzer = KaKsAnalyzer(
            fasta1=args.fasta1,
            fasta2=args.fasta2,
            pairs=args.pairs,
            output_dir=args.output,
            method=args.method,
            kaks_path=args.kaks_path,
            verbose=args.verbose,
            temp_dir=args.temp_dir,
            keep_temp=args.keep_temp
        )
        
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        print("\n❌ 分析被用户中断 | Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"❌ 致命错误 | Fatal error: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
