"""
RxLR扫描主程序模块|RxLR Scanner Main Module
"""

import argparse
import os
import sys
import time
import pandas as pd
from typing import List, Dict

from .config import RxLRConfig
from .utils import RxLRLogger, FastaParser
from .scanner import RxLRMotifScanner, ScanResult


class RxLRScanner:
    """RxLR扫描主类|Main RxLR Scanner Class"""

    def __init__(self, **kwargs):
        """初始化扫描器|Initialize scanner"""
        # 初始化配置|Initialize configuration
        self.config = RxLRConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = RxLRLogger(
            log_file=self.config.log_file_path,
            verbose=self.config.verbose
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化扫描器|Initialize motif scanner
        self.motif_scanner = RxLRMotifScanner(
            window_start=self.config.window_start,
            window_end=self.config.window_end,
            logger=self.logger
        )

        # 结果存储|Results storage
        self.scan_results: List[ScanResult] = []

        # 记录配置信息|Log configuration
        self.logger.info("RxLR效应蛋白扫描器已初始化|RxLR Effector Protein Scanner initialized")
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"搜索窗口|Search window: {self.config.window_start + 1}-{self.config.window_end}")
        self.logger.info(f"最小序列长度|Minimum sequence length: {self.config.min_length}")

    def run_scan(self) -> bool:
        """
        运行扫描流程|Run scanning pipeline

        Returns:
            是否成功|Success status
        """
        try:
            start_time = time.time()

            # 读取和解析FASTA文件|Read and parse FASTA file
            self.logger.info("正在读取FASTA文件|Reading FASTA file...")
            sequences = list(FastaParser.parse_fasta(self.config.input_file))

            if not sequences:
                self.logger.error("FASTA文件为空或格式错误|FASTA file is empty or malformed")
                return False

            self.logger.info(f"成功读取|Successfully read {len(sequences)} 条序列|sequences")

            # 扫描每条序列|Scan each sequence
            self.logger.info("开始扫描序列|Starting sequence scanning...")
            self.scan_results = []

            valid_count = 0
            invalid_count = 0
            candidate_count = 0

            for seq_id, sequence in sequences:
                result = self.motif_scanner.scan_sequence(
                    seq_id,
                    sequence,
                    self.config.min_length
                )

                self.scan_results.append(result)

                # 统计|Statistics
                if result.is_valid_length:
                    valid_count += 1
                else:
                    invalid_count += 1

                if result.is_candidate:
                    candidate_count += 1

            # 记录统计信息|Log statistics
            elapsed_time = time.time() - start_time
            self.logger.info(f"扫描完成|Scanning completed in {elapsed_time:.2f} 秒|seconds")
            self.logger.info(f"总序列数|Total sequences: {len(sequences)}")
            self.logger.info(f"有效序列数(≥{self.config.min_length}aa)|Valid sequences: {valid_count}")
            self.logger.info(f"无效序列数(<{self.config.min_length}aa)|Invalid sequences: {invalid_count}")
            self.logger.info(f"候选RxLR蛋白|Candidate RxLR proteins: {candidate_count}")
            self.logger.info(f"候选率|Candidate rate: {candidate_count/len(sequences)*100:.2f}%")

            # 保存结果|Save results
            self._save_results()

            return True

        except Exception as e:
            self.logger.error(f"扫描过程中发生错误|Error during scanning: {e}")
            return False

    def _save_results(self):
        """保存结果到文件|Save results to files"""
        self.logger.info("正在保存结果|Saving results...")

        # 格式化结果用于输出|Format results for output
        formatted_results = []
        for result in self.scan_results:
            formatted_results.append(
                self.motif_scanner.format_result_for_output(result)
            )

        # 创建DataFrame|Create DataFrame
        df = pd.DataFrame(formatted_results)

        # 保存Excel文件|Save Excel file
        if self.config.excel_output and self.config.excel_file:
            try:
                df.to_excel(self.config.excel_file, index=False, engine='openpyxl')
                self.logger.info(f"Excel结果已保存|Excel results saved: {self.config.excel_file}")
            except Exception as e:
                self.logger.error(f"保存Excel文件失败|Failed to save Excel file: {e}")

        # 保存TSV文件|Save TSV file
        if self.config.tsv_output and self.config.tsv_file:
            try:
                df.to_csv(self.config.tsv_file, sep='\t', index=False)
                self.logger.info(f"TSV结果已保存|TSV results saved: {self.config.tsv_file}")
            except Exception as e:
                self.logger.error(f"保存TSV文件失败|Failed to save TSV file: {e}")

        # 保存候选蛋白列表|Save candidate protein list
        candidate_df = df[df['RxLR_Candidate'] == 'Yes']
        if len(candidate_df) > 0:
            candidate_file = os.path.join(
                self.config.output_dir,
                f"{self.config.output_prefix}_candidates_only.tsv"
            )
            try:
                candidate_df.to_csv(candidate_file, sep='\t', index=False)
                self.logger.info(f"候选蛋白列表已保存|Candidate list saved: {candidate_file}")
            except Exception as e:
                self.logger.error(f"保存候选列表失败|Failed to save candidate list: {e}")

    def get_results(self) -> List[ScanResult]:
        """获取扫描结果|Get scan results"""
        return self.scan_results

    def get_summary(self) -> Dict:
        """
        获取结果摘要|Get results summary

        Returns:
            摘要字典|Summary dictionary
        """
        total = len(self.scan_results)
        if total == 0:
            return {}

        valid = sum(1 for r in self.scan_results if r.is_valid_length)
        invalid = total - valid
        candidates = sum(1 for r in self.scan_results if r.is_candidate)

        # 统计基序类型|Count motif types
        motif_type_counts = {}
        for result in self.scan_results:
            for match in result.motif_matches:
                motif_type_counts[match.motif_type] = motif_type_counts.get(match.motif_type, 0) + 1

        return {
            'total_sequences': total,
            'valid_sequences': valid,
            'invalid_sequences': invalid,
            'candidate_proteins': candidates,
            'candidate_rate_percent': candidates / total * 100 if total > 0 else 0,
            'motif_type_counts': motif_type_counts
        }


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description="RxLR效应蛋白扫描器|RxLR Effector Protein Scanner",
        formatter_class=argparse.RawTextHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-i', '--input',
                       required=True,
                       help='输入FASTA文件|Input FASTA file path')
    parser.add_argument('-o', '--output-prefix',
                       required=True,
                       help='输出文件前缀|Output file prefix')

    # 搜索窗口参数|Search window parameters
    parser.add_argument('--window-start',
                       type=int,
                       default=20,
                       help='窗口起始位置(默认20对应第21位)|Window start position (default: 20 for position 21)')
    parser.add_argument('--window-end',
                       type=int,
                       default=120,
                       help='窗口结束位置(默认120)|Window end position (default: 120)')
    parser.add_argument('--min-length',
                       type=int,
                       default=120,
                       help='最小序列长度(默认120)|Minimum sequence length (default: 120)')

    # 输出选项|Output options
    parser.add_argument('--output-dir',
                       default='./rxlr_scanner_output',
                       help='输出目录|Output directory (default: ./rxlr_scanner_output)')
    parser.add_argument('--no-excel',
                       action='store_true',
                       help='不生成Excel输出|Do not generate Excel output')
    parser.add_argument('--no-tsv',
                       action='store_true',
                       help='不生成TSV输出|Do not generate TSV output')

    # 日志选项|Logging options
    parser.add_argument('-v', '--verbose',
                       action='store_true',
                       help='详细输出模式|Verbose output mode')
    parser.add_argument('--log-file',
                       action='store_true',
                       help='生成日志文件|Generate log file')

    args = parser.parse_args()

    try:
        # 创建扫描器|Create scanner
        scanner = RxLRScanner(
            input_file=args.input,
            output_prefix=args.output_prefix,
            window_start=args.window_start,
            window_end=args.window_end,
            min_length=args.min_length,
            output_dir=args.output_dir,
            excel_output=not args.no_excel,
            tsv_output=not args.no_tsv,
            verbose=args.verbose,
            log_file=f"{args.output_prefix}.log" if args.log_file else None
        )

        # 运行扫描|Run scanning
        success = scanner.run_scan()

        # 退出|Exit
        sys.exit(0 if success else 1)

    except Exception as e:
        print(f"错误|Error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
