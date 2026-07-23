"""
Ka/Ks Calculator主分析流水线|Ka/Ks Calculator Main Analysis Pipeline
功能: 协调各个模块完成完整的Ka/Ks分析流程|Coordinate all modules for complete Ka/Ks analysis
"""

import os
import sys
import argparse
import subprocess
import tempfile
import shutil
from pathlib import Path
from datetime import datetime
from typing import Dict, List, Set, Tuple

from .config import KaKsConfig
from .utils import KaKsLogger
from .validator import SequenceValidator
from .calculator import KaKsCalculator
from .processor import ResultProcessor


class KaKsAnalyzer:
    """Ka/Ks主分析器|Main Ka/Ks analyzer"""

    def __init__(self, fasta1: str, fasta2: str, pairs: str, output_dir: str,
                 method: str = None, kaks_path: str = "KaKs_Calculator",
                 threads: int = 12, verbose: bool = False,
                 temp_dir: str = None, keep_temp: bool = False):
        """
        初始化分析器|Initialize analyzer

        Args:
            fasta1: 第一个FASTA文件路径|First FASTA file path
            fasta2: 第二个FASTA文件路径|Second FASTA file path
            pairs: 配对文件路径|Pairs file path
            output_dir: 输出目录|Output directory
            method: 计算方法|Calculation method
            kaks_path: KaKs_Calculator路径|KaKs_Calculator path
            threads: 线程数|Thread count
            verbose: 详细模式|Verbose mode
            temp_dir: 临时目录|Temporary directory
            keep_temp: 保留临时文件|Keep temporary files
        """
        self.fasta1 = fasta1
        self.fasta2 = fasta2
        self.pairs = pairs
        self.output_dir = output_dir
        self.method = method or KaKsConfig().method
        self.kaks_path = kaks_path
        self.threads = threads
        self.verbose = verbose
        self.temp_dir = temp_dir
        self.keep_temp = keep_temp

        output_path = Path(output_dir)
        output_path.mkdir(parents=True, exist_ok=True)

        # 创建标准输出目录结构|Create standard output directory structure
        (output_path / "00_pipeline_info").mkdir(parents=True, exist_ok=True)

        # 初始化日志器(自动放到99_logs/)|Initialize logger (auto-placed in 99_logs/)
        self.logger = KaKsLogger(output_path, verbose=verbose)

        # 初始化组件|Initialize components
        self.config = KaKsConfig()
        self.validator = SequenceValidator(self.logger)
        self.calculator = KaKsCalculator(self.logger, kaks_path)
        self.processor = ResultProcessor(self.logger)

        self._temp_dir_created = None

    def run_analysis(self):
        """执行完整的Ka/Ks分析流程|Execute complete Ka/Ks analysis pipeline"""
        try:
            start_time = datetime.now()
            self.logger.separator("Ka/Ks Calculator 2.0 Analysis Pipeline")
            self.logger.info("开始Ka/Ks分析流程|Starting Ka/Ks analysis pipeline")

            self._setup_temp_directory()

            seq1_dict, seq2_dict, pairs_list = self._validate_inputs()

            results_df = self._run_calculation(seq1_dict, seq2_dict, pairs_list)

            self._process_and_save_results(results_df)

            self._cleanup()

            runtime = datetime.now() - start_time
            self.logger.separator()
            self.logger.success(f"分析完成|Analysis completed. 用时|Runtime: {runtime}")

        except Exception as e:
            self.logger.error(f"分析流程失败|Analysis pipeline failed: {e}")
            self._cleanup()
            raise

    def _extract_ids_from_pairs(self) -> Tuple[Set[str], Set[str], List[Tuple[str, str, str]]]:
        """从配对文件提取需要的序列ID|Extract needed sequence IDs from pairs file"""
        self.logger.info(f"解析配对文件|Parsing pair file: {self.pairs}")
        pairs_list = self.validator.parse_pair_file(self.pairs)
        if not pairs_list:
            raise ValueError(f"配对文件为空|Empty pair file: {self.pairs}")

        seq1_ids = {p[0] for p in pairs_list}
        seq2_ids = {p[1] for p in pairs_list}
        self.logger.info(f"配对文件包含|Pair file contains: {len(pairs_list)} 对, "
                        f"涉及 {len(seq1_ids)} 个物种1序列|seq1 IDs, "
                        f"{len(seq2_ids)} 个物种2序列|seq2 IDs")
        return seq1_ids, seq2_ids, pairs_list

    def _extract_subset_fasta(self, fasta_file: str, ids: Set[str], label: str) -> str:
        """用seqkit grep提取需要的序列子集|Extract sequence subset using seqkit grep"""
        id_file = os.path.join(self._temp_dir_created, f"{label}_ids.txt")
        output_fasta = os.path.join(self._temp_dir_created, f"{label}_subset.fasta")

        with open(id_file, 'w') as f:
            for sid in sorted(ids):
                f.write(f"{sid}\n")

        seqkit_path = "~/miniforge3/envs/BioinfTools/bin/seqkit"
        from ..common.paths import expand_path
        seqkit_path = expand_path(seqkit_path)

        cmd = [seqkit_path, 'grep', '-f', id_file, fasta_file, '-o', output_fasta]
        self.logger.info(f"命令|Command: {' '.join(cmd)}")
        self.logger.info(f"提取{label}序列子集|Extracting {label} sequence subset: {len(ids)} sequences")

        result = subprocess.run(cmd, capture_output=True, text=True)
        if result.returncode != 0:
            self.logger.error(f"seqkit grep失败|seqkit grep failed: {result.stderr}")
            raise RuntimeError(f"seqkit grep failed: {result.stderr}")

        return output_fasta

    def _setup_temp_directory(self):
        """设置临时工作目录|Setup temporary working directory"""
        if self.temp_dir:
            self._temp_dir_created = self.temp_dir
            Path(self._temp_dir_created).mkdir(parents=True, exist_ok=True)
        else:
            # 临时目录落到 <output>/tmp 下,避免系统 /tmp 爆满|
            # Temp dir under <output>/tmp to avoid system /tmp overflow
            tmp_root = Path(self.output_dir) / "tmp"
            tmp_root.mkdir(parents=True, exist_ok=True)
            self._temp_dir_created = tempfile.mkdtemp(prefix="kakscalc_", dir=str(tmp_root))

        self.logger.info(f"使用临时目录|Using temporary directory: {self._temp_dir_created}")

    def _validate_inputs(self) -> Tuple[Dict[str, str], Dict[str, str], List[Tuple[str, str, str]]]:
        """验证输入文件|Validate input files"""
        self.logger.separator("输入文件验证|Input File Validation")

        # 1. 先解析配对文件，获取需要的序列ID|Parse pairs file first to get needed IDs
        seq1_ids, seq2_ids, pairs_list = self._extract_ids_from_pairs()

        # 2. 用seqkit提取需要的序列子集|Extract sequence subset with seqkit
        subset_fasta1 = self._extract_subset_fasta(self.fasta1, seq1_ids, "seq1")
        subset_fasta2 = self._extract_subset_fasta(self.fasta2, seq2_ids, "seq2")

        # 3. 只验证子集序列|Validate only subset sequences
        self.logger.info(f"验证物种1序列子集|Validating species1 sequence subset: {subset_fasta1}")
        valid1, seq1_dict = self.validator.validate_fasta_file(subset_fasta1)
        if not valid1:
            raise ValueError(f"第一个FASTA文件无效|Invalid first FASTA file: {self.fasta1}")

        self.logger.info(f"验证物种2序列子集|Validating species2 sequence subset: {subset_fasta2}")
        valid2, seq2_dict = self.validator.validate_fasta_file(subset_fasta2)
        if not valid2:
            raise ValueError(f"第二个FASTA文件无效|Invalid second FASTA file: {self.fasta2}")

        # 4. 过滤掉含无效序列的配对|Filter out pairs with invalid sequences
        valid_seq1_ids = set(seq1_dict.keys())
        valid_seq2_ids = set(seq2_dict.keys())
        original_count = len(pairs_list)
        pairs_list = [(s1, s2, name) for s1, s2, name in pairs_list
                      if s1 in valid_seq1_ids and s2 in valid_seq2_ids]
        filtered_count = original_count - len(pairs_list)
        if filtered_count > 0:
            self.logger.warning(f"因序列验证失败，过滤 {filtered_count} 个配对|Filtered {filtered_count} pairs due to sequence validation failure")
            self.logger.warning(f"剩余有效配对|Remaining valid pairs: {len(pairs_list)}")

        if len(pairs_list) > self.config.max_batch_size:
            self.logger.warning(f"配对数量 ({len(pairs_list)}) 超过推荐最大值 ({self.config.max_batch_size})|Pair count exceeds recommended max")
            self.logger.warning("大批量分析可能需要较长时间|Large batch analysis may take significant time")

        self.logger.success(f"输入文件验证通过|All input files validated: {len(seq1_dict)} + {len(seq2_dict)} 条序列, {len(pairs_list)} 对|sequences, pairs")
        return seq1_dict, seq2_dict, pairs_list

    def _run_calculation(self, seq1_dict: Dict[str, str], seq2_dict: Dict[str, str],
                        pairs_list: List[Tuple[str, str, str]]):
        """执行Ka/Ks计算|Execute Ka/Ks calculation"""
        self.logger.separator("Ka/Ks计算|Ka/Ks Calculation")

        input_file = self.calculator.prepare_input_file(
            seq1_dict, seq2_dict, pairs_list, self._temp_dir_created
        )

        output_file = self.calculator.run_calculation(
            input_file, self.method, self._temp_dir_created
        )

        results_df = self.processor.parse_results(output_file)

        # 添加seq1_id和seq2_id列|Add seq1_id and seq2_id columns
        pair_map = {pair_name: (seq1_id, seq2_id) for seq1_id, seq2_id, pair_name in pairs_list}
        id_col = 'Pair_ID' if 'Pair_ID' in results_df.columns else 'Sequence'
        results_df['seq1_id'] = results_df[id_col].map(lambda x: pair_map.get(x, ('', ''))[0])
        results_df['seq2_id'] = results_df[id_col].map(lambda x: pair_map.get(x, ('', ''))[1])

        return results_df

    def _process_and_save_results(self, results_df):
        """处理并保存结果|Process and save results"""
        self.logger.separator("结果处理|Result Processing")

        stats = self.processor.generate_summary_stats(results_df)

        self.processor.save_results(results_df, stats, self.output_dir)

        self._print_analysis_summary(stats)

    def _print_analysis_summary(self, stats: Dict):
        """打印分析汇总到控制台|Print analysis summary to console"""
        self.logger.separator("分析汇总|Analysis Summary")

        analysis_info = stats.get('analysis_info', {})
        omega_stats = stats.get('omega_statistics', {})
        selection_dist = stats.get('selection_distribution', {})

        self.logger.info(f"  总序列对数|Total sequence pairs: {analysis_info.get('total_pairs', 'N/A')}")
        self.logger.info(f"  成功计算|Successful calculations: {analysis_info.get('successful_calculations', 'N/A')}")
        self.logger.info(f"  失败计算|Failed calculations: {analysis_info.get('failed_calculations', 'N/A')}")
        self.logger.info(f"  成功率|Success rate: {analysis_info.get('success_rate', 0):.1%}")

        if omega_stats.get('count', 0) > 0:
            self.logger.info(f"  平均omega (Ka/Ks)|Mean omega: {omega_stats.get('mean', 0):.4f}")
            self.logger.info(f"  中位数omega|Median omega: {omega_stats.get('median', 0):.4f}")
            self.logger.info(f"  omega范围|Omega range: {omega_stats.get('min', 0):.4f} - {omega_stats.get('max', 0):.4f}")

        if selection_dist:
            self.logger.info("  选择压力分布|Selection pressure distribution:")
            total_successful = analysis_info.get('successful_calculations', 1)
            for sel_type, count in selection_dist.items():
                percentage = (count / total_successful) * 100
                self.logger.info(f"    {sel_type}: {count} ({percentage:.1f}%)")

        self.logger.info("  生成的输出文件|Generated output files:")
        self.logger.info(f"    kaks_detailed.csv - CSV详细结果|CSV detailed results")
        self.logger.info(f"    kaks_detailed.tsv - TSV详细结果|TSV detailed results")
        self.logger.info(f"    kaks_summary.xlsx - Excel汇总文件|Excel summary file")
        self.logger.info(f"    summary_stats.json - JSON统计数据|JSON statistics")

        bio_interp = stats.get('biological_interpretation', {})
        if bio_interp:
            self.logger.info("  生物学解释|Biological interpretation:")
            for aspect, interpretation in bio_interp.items():
                self.logger.info(f"    {interpretation}")

    def _cleanup(self):
        """清理临时文件|Clean up temporary files"""
        if self._temp_dir_created and not self.keep_temp:
            try:
                shutil.rmtree(self._temp_dir_created)
                self.logger.info("临时文件已清理|Temporary files cleaned up")
            except Exception as e:
                self.logger.warning(f"清理临时文件失败|Failed to clean up temporary files: {e}")
        elif self.keep_temp:
            self.logger.info(f"保留临时文件|Keeping temporary files: {self._temp_dir_created}")


def create_argument_parser() -> argparse.ArgumentParser:
    """创建命令行参数解析器|Create command line argument parser"""

    parser = argparse.ArgumentParser(
        description="Ka/Ks Calculator 2.0 Analysis Pipeline",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -1 species1.fasta -2 species2.fasta -p pairs.txt -o results/
        """
    )

    # 必需输入文件|Required input files
    parser.add_argument('-1', '--fasta1', '--species1',
                       required=True, type=str, metavar='FILE',
                       help='第一个FASTA文件 (物种1 CDS序列)|First FASTA file (species 1 CDS sequences)')

    parser.add_argument('-2', '--fasta2', '--species2',
                       required=True, type=str, metavar='FILE',
                       help='第二个FASTA文件 (物种2 CDS序列)|Second FASTA file (species 2 CDS sequences)')

    parser.add_argument('-p', '--pairs', '--pair-file',
                       required=True, type=str, metavar='FILE',
                       help='序列配对文件 (TSV/CSV格式)|Sequence pair file (TSV/CSV format)')

    # 输出目录|Output directory
    parser.add_argument('-o', '--output', '--output-dir',
                       required=True, type=str, metavar='DIR',
                       help='输出目录|Output directory for results')

    # 分析参数|Analysis parameters
    parser.add_argument('-m', '--method', '--calc-method',
                       default='GMYN',
                       choices=['GMYN', 'MYN', 'YN', 'NG', 'LWL', 'LPB',
                                'MLWL', 'MLPB', 'GY', 'MS', 'MA', 'GNG',
                                'GLWL', 'GLPB', 'GMLWL', 'GMLPB', 'GYN'],
                       help='计算方法 (默认: GMYN)|Calculation method (default: GMYN)')

    parser.add_argument('-t', '--threads',
                       default=12, type=int, metavar='INT',
                       help='线程数 (默认: 12)|Thread count (default: 12)')

    parser.add_argument('--kaks-path', '--calculator-path',
                       default='KaKs_Calculator', type=str, metavar='PATH',
                       help='KaKs_Calculator可执行文件路径|Path to KaKs_Calculator executable')

    # 运行选项|Runtime options
    parser.add_argument('-v', '--verbose', '--debug',
                       action='store_true',
                       help='启用详细日志记录|Enable verbose logging')

    parser.add_argument('--temp-dir', '--tmp-dir',
                       type=str, metavar='DIR',
                       help='自定义临时目录|Custom temporary directory')

    parser.add_argument('--keep-temp', '--no-cleanup',
                       action='store_true',
                       help='保留临时文件 (用于调试)|Keep temporary files (for debugging)')

    parser.add_argument('--version', action='version', version='%(prog)s 1.0.0')

    return parser


def main():
    """主入口点|Main entry point"""
    try:
        parser = create_argument_parser()
        args = parser.parse_args()

        analyzer = KaKsAnalyzer(
            fasta1=args.fasta1,
            fasta2=args.fasta2,
            pairs=args.pairs,
            output_dir=args.output,
            method=args.method,
            kaks_path=args.kaks_path,
            threads=args.threads,
            verbose=args.verbose,
            temp_dir=args.temp_dir,
            keep_temp=args.keep_temp
        )

        analyzer.run_analysis()

    except KeyboardInterrupt:
        print("分析被用户中断|Analysis interrupted by user")
        sys.exit(1)
    except Exception as e:
        print(f"致命错误|Fatal error: {e}")
        sys.exit(1)


if __name__ == "__main__":
    main()
