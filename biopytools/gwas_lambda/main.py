"""
GWAS Lambda计算主程序模块 | GWAS Lambda Calculator Main Module
"""

import argparse
import sys
import os
import glob
from typing import List, Dict

from .config import GWASLambdaConfig
from .analyzer import GWASResultAnalyzer


class GWASLambdaCalculator:
    """GWAS Lambda计算主类 | Main GWAS Lambda Calculator Class"""

    def __init__(self, **kwargs):
        """
        初始化计算器 | Initialize calculator

        Args:
            **kwargs: 配置参数 | Configuration parameters
        """
        # 初始化配置 | Initialize configuration
        self.config = GWASLambdaConfig(**kwargs)
        self.config.validate()

        # 初始化分析器 | Initialize analyzer
        self.analyzer = GWASResultAnalyzer(
            significance_threshold=self.config.significance_threshold,
            p_value_column=self.config.p_value_column,
            min_p_value=self.config.min_p_value,
            max_p_value=self.config.max_p_value,
            expected_median=self.config.expected_median
        )

        # 结果存储 | Results storage
        self.results = []

    def find_result_files(self) -> List[str]:
        """
        查找GWAS结果文件 | Find GWAS result files

        Returns:
            List[str]: 文件路径列表 | List of file paths
        """
        files = glob.glob(self.config.search_pattern)

        if not files:
            print(f"❌ 未找到匹配文件 | No files found matching pattern: {self.config.search_pattern}")
            sys.exit(1)

        print(f"📁 找到 {len(files)} 个文件，开始分析... | Found {len(files)} files, starting analysis...")
        return files

    def analyze_single_file(self, file_path: str) -> Dict:
        """
        分析单个文件 | Analyze a single file

        Args:
            file_path: 文件路径 | File path

        Returns:
            Dict: 分析结果 | Analysis result
        """
        folder_name = os.path.dirname(file_path)

        # 调用分析器 | Call analyzer
        lambda_val, total_snps, sig_count, error_msg = self.analyzer.analyze_file(file_path)

        # 生成结果字典 | Generate result dictionary
        result = {
            'folder': folder_name,
            'lambda_gc': lambda_val,
            'sig_count': sig_count,
            'total_snps': total_snps,
            'error_msg': error_msg
        }

        # 获取状态 | Get status
        if error_msg == "无显著位点 | No significant variants found":
            result['status'] = self.analyzer.get_status(None, 0)
            result['lambda_display'] = 'NA'
        elif lambda_val is not None:
            result['status'] = self.analyzer.get_status(lambda_val, sig_count)
            result['lambda_display'] = f"{lambda_val:.4f}"

            # 高亮显示需要注意的结果 | Highlight results that need attention
            if self.analyzer.should_highlight(result['status']):
                print(f"⚠️  {folder_name}: Lambda={lambda_val:.4f}, SigHits={sig_count} ({result['status']})")
        else:
            result['status'] = f"Error: {error_msg}"
            result['lambda_display'] = 'NA'
            print(f"❌ {folder_name}: Error - {error_msg}")

        return result

    def run_batch_analysis(self) -> bool:
        """
        运行批量分析 | Run batch analysis

        Returns:
            bool: 是否成功 | Whether successful
        """
        try:
            # 查找文件 | Find files
            files = self.find_result_files()

            # 分析每个文件 | Analyze each file
            for file_path in files:
                result = self.analyze_single_file(file_path)
                self.results.append(result)

            # 写入结果 | Write results
            self.write_results()

            # 打印总结 | Print summary
            self.print_summary()

            return True

        except Exception as e:
            print(f"❌ 分析过程中出现错误 | Error during analysis: {e}")
            return False

    def write_results(self):
        """写入结果到文件 | Write results to file"""
        output_path = self.config.get_output_path()

        with open(output_path, 'w', encoding='utf-8') as out:
            # 写入表头 | Write header
            out.write("Folder\tLambda_GC\tSig_SNPs(<1e-5)\tTotal_SNPs\tStatus\n")

            # 写入结果 | Write results
            for result in self.results:
                line = f"{result['folder']}\t{result['lambda_display']}\t{result['sig_count']}\t{result['total_snps']}\t{result['status']}\n"
                out.write(line)

        print(f"\n📝 结果已保存至 | Results saved to: {output_path}")

    def print_summary(self):
        """打印分析总结 | Print analysis summary"""
        total_files = len(self.results)
        no_signal_count = sum(1 for r in self.results if "No Signals" in r['status'])
        ideal_count = sum(1 for r in self.results if "Ideal" in r['status'])
        acceptable_count = sum(1 for r in self.results if "Acceptable" in r['status'])
        inflated_count = sum(1 for r in self.results if "Inflated" in r['status'])
        deflated_count = sum(1 for r in self.results if "Deflated" in r['status'])
        error_count = sum(1 for r in self.results if "Error" in r['status'])

        print(f"\n📊 分析总结 | Analysis Summary:")
        print(f"   总文件数 | Total files: {total_files}")
        print(f"   无显著信号 | No signals: {no_signal_count}")
        print(f"   理想结果 | Ideal: {ideal_count}")
        print(f"   可接受结果 | Acceptable: {acceptable_count}")
        print(f"   膨胀结果 | Inflated: {inflated_count}")
        print(f"   压缩结果 | Deflated: {deflated_count}")
        print(f"   错误 | Errors: {error_count}")

        print(f"\n💡 结果说明 | Result Interpretation:")
        print(f"1. Lambda_GC 显示 NA 且 Status 为 'No Signals': 说明该性状没有检测到 P<1e-5 的位点，无需关注。")
        print(f"2. 关注 Status 为 'Ideal' 或 'Acceptable' 且 Sig_SNPs 数量合理的性状。")
        print(f"3. 'Inflated' 状态表明可能存在群体分层或假阳性。")


def main():
    """命令行入口函数 | Command line entry function"""
    parser = argparse.ArgumentParser(
        description="GWAS Lambda GC计算工具 | GWAS Lambda GC Calculator",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例 | Examples:

  # 使用默认设置 | Use default settings
  python -m biopytools.gwas_lambda

  # 指定搜索模式和输出文件 | Specify search pattern and output file
  python -m biopytools.gwas_lambda --pattern "results/*/gwas.txt" --output "my_lambda_assessment.txt"

  # 调整显著性阈值 | Adjust significance threshold
  python -m biopytools.gwas_lambda --threshold 1e-6

  # 指定P值列 | Specify P-value column
  python -m biopytools.gwas_lambda --p-column 2
        """
    )

    parser.add_argument('--pattern', '-p',
                       default="feture_*/GWAS_Result.mlm.manht_input",
                       help='文件搜索模式 | File search pattern (default: "feture_*/GWAS_Result.mlm.manht_input")')

    parser.add_argument('--output', '-o',
                       default="Batch_Lambda_Assessment.txt",
                       help='输出文件名 | Output filename (default: "Batch_Lambda_Assessment.txt")')

    parser.add_argument('--threshold', '-t',
                       type=float,
                       default=1e-5,
                       help='显著性阈值 | Significance threshold (default: 1e-5)')

    parser.add_argument('--p-column', '-c',
                       type=int,
                       default=3,
                       help='P值所在列索引 | P-value column index (0-based, default: 3)')

    parser.add_argument('--output-dir', '-d',
                       default="./gwas_lambda_output",
                       help='输出目录 | Output directory (default: "./gwas_lambda_output")')

    args = parser.parse_args()

    # 创建计算器并运行分析 | Create calculator and run analysis
    calculator = GWASLambdaCalculator(
        search_pattern=args.pattern,
        output_file=args.output,
        significance_threshold=args.threshold,
        p_value_column=args.p_column,
        output_dir=args.output_dir
    )

    success = calculator.run_batch_analysis()
    sys.exit(0 if success else 1)


if __name__ == "__main__":
    main()