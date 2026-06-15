"""
GWAS GEC分析主程序模块（简化版）|GWAS GEC Analysis Main Module (Simplified)
功能: 使用GEC计算GWAS显著性阈值，支持PLINK binary和VCF格式|
Features: Calculate GWAS significance thresholds using GEC, supports PLINK binary and VCF formats
"""

import argparse
import gzip
import os
import subprocess
import sys
from .config import GECConfig
from .utils import GECLogger, CommandRunner, validate_chromosome_format
from .calculator import GECThresholdCalculator


class GECAnalyzer:
    """GWAS GEC分析器|GWAS GEC Analyzer"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = GECConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        self.logger_manager = GECLogger(
            log_file=self.config.log_file,
            log_level="INFO"
        )
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)

        # 初始化计算器|Initialize calculator
        self.calculator = GECThresholdCalculator(self.logger)

        # 自动转换染色体格式（如果需要）|Auto-convert chromosome format (if needed)
        self.temp_pfile = None
        if self.config.convert_chrom_format:
            self._convert_chromosome_format()

    def _convert_chromosome_format(self) -> bool:
        """
        自动转换P值文件和VCF文件的染色体格式|Auto-convert chromosome format in P-value and VCF files
        KGGSee不支持Chr01格式，需要转换为数字格式|KGGSee does not support Chr01 format, need to convert to numeric

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("\n" + "=" * 80)
        self.logger.info("检测并转换染色体格式|Detecting and converting chromosome format")
        self.logger.info("=" * 80)

        import shutil
        import gzip

        # 检查P值文件是否需要转换|Check if P-value file needs conversion
        needs_conversion = False
        sample_chrom = None

        # 读取P值文件前几行检查格式|Read first few lines of P-value file to check format
        open_func = gzip.open if self.config.pfile.endswith('.gz') else open
        try:
            with open_func(self.config.pfile, 'rt') as f:
                # 读取第一行（可能是标题）|Read first line (might be header)
                header = f.readline()
                # 读取第二行（数据）|Read second line (data)
                line = f.readline()

                if line:
                    first_chrom = line.strip().split('\t')[0]
                    self.logger.debug(f"检查第一个染色体|Checking first chromosome: {first_chrom}")

                    # 检查是否需要转换|Check if conversion is needed
                    if first_chrom.startswith('Chr') or first_chrom.startswith('chr'):
                        # Chr01 或 chr01 需要转换|Chr01 or chr01 needs conversion
                        if first_chrom[3:].isdigit() or first_chrom[3:] in ['X', 'Y', 'M', 'MT']:
                            needs_conversion = True
                            sample_chrom = first_chrom
        except Exception as e:
            self.logger.warning(f"染色体格式检测失败|Chromosome format detection failed: {str(e)}")

        if not needs_conversion:
            self.logger.info("P值文件染色体格式无需转换|P-value file chromosome format does not need conversion")
            return True

        self.logger.info(f"检测到需要转换的格式|Detected format needing conversion: {sample_chrom}")
        self.logger.info("开始转换P值文件染色体格式|Converting P-value file chromosome format")

        # 创建转换后的P值文件|Create converted P-value file
        self.temp_pfile = os.path.join(self.config.output_dir, 'gwas_domain_p_file_converted.txt')

        # 读取并转换|Read and convert
        with open(self.config.pfile, 'r') as fin, open(self.temp_pfile, 'w') as fout:
            for line in fin:
                if line.startswith('#'):
                    fout.write(line)
                    continue

                parts = line.strip().split('\t')
                if len(parts) >= 3:
                    # 转换染色体格式|Convert chromosome format
                    chrom = parts[0]
                    if chrom.startswith('Chr'):
                        # Chr01 -> 1, ChrX -> X
                        chrom = chrom[3:]
                        # 移除前导零|Remove leading zeros
                        if chrom.isdigit():
                            chrom = str(int(chrom))
                    elif chrom.startswith('chr'):
                        # chr01 -> 1, chrX -> X
                        chrom = chrom[3:]
                        if chrom.isdigit():
                            chrom = str(int(chrom))

                    parts[0] = chrom
                    fout.write('\t'.join(parts) + '\n')
                else:
                    fout.write(line)

        self.logger.info(f"P值文件转换完成|P-value file conversion completed: {self.temp_pfile}")

        # 更新配置使用转换后的文件|Update config to use converted file
        self.config.pfile = self.temp_pfile

        return True

    def build_java_command(self) -> str:
        """
        构建Java命令|Build Java command

        Returns:
            完整的Java命令字符串|Complete Java command string
        """
        cmd_parts = []

        # Java基础命令|Java base command
        cmd_parts.append(f"java -Xmx{self.config.memory}")
        cmd_parts.append(f"-jar {self.config.kggsee_jar}")

        # GEC算法选项|GEC algorithm option
        cmd_parts.append("--var-gec")

        # P值文件（可能是转换后的）|P-value file (possibly converted)
        cmd_parts.append(f"--pfile {self.config.pfile}")

        # 参考文件（VCF格式）|Reference file (VCF format)
        # KGGSee只支持VCF格式|KGGSee only supports VCF format
        cmd_parts.append(f"--vcf-ref {self.config.reference}")

        # 输出前缀|Output prefix
        output_prefix = os.path.join(self.config.output_dir, 'gwas_gec')
        cmd_parts.append(f"--out {output_prefix}")

        # 线程数|Thread number
        cmd_parts.append(f"--nt {self.config.threads}")

        # MAF过滤|MAF filter
        cmd_parts.append(f"--filter-maf-le {self.config.filter_maf_le}")

        # P值阈值|P-value threshold
        cmd_parts.append(f"--p-value-cutoff {self.config.p_value_cutoff}")

        # 列名配置|Column names configuration
        cmd_parts.append(f"--chrom-col {self.config.chrom_col}")
        cmd_parts.append(f"--pos-col {self.config.pos_col}")
        cmd_parts.append(f"--p-col {self.config.p_col}")

        # 保留参考文件缓存|Keep reference file cache
        if self.config.keep_ref:
            cmd_parts.append("--keep-ref")

        # 染色体范围|Chromosome range
        if self.config.chromosomes:
            cmd_parts.append(f"--chrom {self.config.chromosomes}")

        return " ".join(cmd_parts)

    def run_gec(self) -> bool:
        """
        运行GEC算法|Run GEC algorithm

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 80)
        self.logger.info("GWAS GEC分析|GWAS GEC Analysis")
        self.logger.info("=" * 80)
        self.logger.info(f"P值文件|P-value file: {self.config.pfile}")
        self.logger.info(f"参考文件|Reference file: {self.config.reference}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"线程数|Thread number: {self.config.threads}")
        self.logger.info(f"内存分配|Memory allocation: {self.config.memory}")
        self.logger.info("=" * 80)

        # 检测P值文件染色体格式（仅用于提示）|Detect P-value chromosome format (for hint only)
        self.logger.info("检测P值文件染色体格式|Detecting P-value file chromosome format")
        is_valid, format_type, examples = validate_chromosome_format(
            self.config.pfile,
            self.config.chrom_col,
            self.logger
        )

        if is_valid:
            self.logger.info(f"  检测到的格式|Detected format: {format_type}")
            self.logger.info(f"  示例染色体|Example chromosomes: {', '.join(examples)}")

        # 构建并执行命令|Build and execute command
        cmd = self.build_java_command()
        self.logger.info(f"执行GEC算法|Running GEC algorithm")
        self.logger.debug(f"完整命令|Full command:\n{cmd}\n")

        success = self.cmd_runner.run(
            cmd,
            description="GEC算法计算|GEC algorithm calculation"
        )

        if success:
            self.logger.info("GEC算法执行成功|GEC algorithm executed successfully")
        else:
            self.logger.error("GEC算法执行失败|GEC algorithm execution failed")

        return success

    def calculate_threshold(self) -> dict:
        """
        计算显著性阈值|Calculate significance threshold

        Returns:
            包含分析结果的字典|Dictionary containing analysis results
        """
        self.logger.info("\n计算显著性阈值|Calculating significance threshold")

        # 查找结果文件|Find result file
        result_file = os.path.join(self.config.output_dir, 'gwas_gec.effective.size.txt.gz')

        if not os.path.exists(result_file):
            self.logger.error(f"GEC输出文件未找到|GEC output file not found: {result_file}")
            raise FileNotFoundError(f"GEC输出文件未找到|GEC output file not found")

        self.logger.info(f"找到GEC输出文件|Found GEC output file: {result_file}")

        # 运行分析|Run analysis
        output_prefix = os.path.join(self.config.output_dir, 'gwas_gec')
        results = self.calculator.run_analysis(
            result_file=result_file,
            alpha=self.config.alpha,
            output_prefix=output_prefix
        )

        return results

    def run_analysis(self) -> dict:
        """
        运行完整的分析流程|Run complete analysis pipeline

        Returns:
            包含分析结果的字典|Dictionary containing analysis results
        """
        # 步骤1: 运行GEC算法|Step 1: Run GEC algorithm
        if not self.run_gec():
            self.logger.error("GEC算法执行失败，退出|GEC algorithm execution failed, exiting")
            self._cleanup_temp_files()
            sys.exit(1)

        # 步骤2: 计算显著性阈值|Step 2: Calculate significance threshold
        try:
            results = self.calculate_threshold()

            self.logger.info("\n" + "=" * 80)
            self.logger.info("分析完成|Analysis Completed")
            self.logger.info("=" * 80)
            self.logger.info(f"总有效检验数|Total effective tests: {results['total_effective_tests']:.2f}")
            self.logger.info(f"校正后显著性阈值|Adjusted significance threshold: {results['significance_threshold']:.2e}")
            output_prefix = os.path.join(self.config.output_dir, 'gwas_gec')
            self.logger.info(f"汇总报告|Summary report: {output_prefix}_summary.txt")
            self.logger.info("=" * 80)

            # 清理临时文件|Cleanup temporary files
            self._cleanup_temp_files()

            return results

        except Exception as e:
            self.logger.error(f"阈值计算失败|Threshold calculation failed: {str(e)}")
            self._cleanup_temp_files()
            sys.exit(1)

    def _cleanup_temp_files(self):
        """清理临时文件|Cleanup temporary files"""
        import os

        files_to_cleanup = []

        # 清理转换后的P值文件|Cleanup converted P-value file
        if self.temp_pfile and os.path.exists(self.temp_pfile):
            files_to_cleanup.append(self.temp_pfile)

        if files_to_cleanup:
            self.logger.info("\n清理临时文件|Cleaning up temporary files")
            for temp_file in files_to_cleanup:
                try:
                    os.remove(temp_file)
                    self.logger.debug(f"已删除|Deleted: {temp_file}")
                except Exception as e:
                    self.logger.warning(f"删除文件失败|Failed to delete file {temp_file}: {str(e)}")
            self.logger.info("临时文件清理完成|Temporary files cleanup completed")


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='GWAS基因组范围多重检验校正(GEC) - 简化版|GWAS Genome-wide Error Correction (GEC) - Simplified',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-p', '--pfile', required=True,
                       help='GWAS P值汇总统计文件|GWAS P-value summary statistics file')
    parser.add_argument('-r', '--reference', required=True,
                       help='参考文件（VCF文件或PLINK binary前缀）|Reference file (VCF or PLINK binary prefix)')

    # 可选参数|Optional parameters
    parser.add_argument('-o', '--output-dir',
                       default='./gec_output',
                       help='输出目录|Output directory')
    parser.add_argument('-t', '--threads',
                       type=int, default=16,
                       help='线程数|Number of threads')
    parser.add_argument('-m', '--memory',
                       default='8g',
                       help='Java内存分配|Java memory allocation')

    # 处理参数|Processing parameters
    parser.add_argument('--maf-filter',
                       type=float, default=0.05,
                       dest='filter_maf_le',
                       help='MAF过滤阈值|MAF filter threshold')
    parser.add_argument('--p-cutoff',
                       type=float, default=0.05,
                       dest='p_value_cutoff',
                       help='P值阈值|P-value threshold')
    parser.add_argument('--no-keep-ref',
                       action='store_false',
                       dest='keep_ref',
                       help='不保留参考文件缓存|Do not keep reference file cache')

    # 染色体范围|Chromosome range
    parser.add_argument('--chrom',
                       help='染色体范围|Chromosome range')

    # 染色体格式转换|Chromosome format conversion
    parser.add_argument('--no-convert-chrom',
                       action='store_false',
                       dest='convert_chrom_format',
                       help='禁用自动染色体格式转换|Disable automatic chromosome format conversion')

    # P值文件列名配置|P-value file column names
    parser.add_argument('--chrom-col',
                       default='CHR',
                       help='染色体列名|Chromosome column name')
    parser.add_argument('--pos-col',
                       default='BP',
                       help='位置列名|Position column name')
    parser.add_argument('--p-col',
                       default='P',
                       help='P值列名|P-value column name')

    # Alpha水平|Alpha level
    parser.add_argument('--alpha',
                       type=float, default=0.05,
                       help='显著性水平(FWER)|Significance level (FWER)')

    # KGGSee路径|KGGSee path
    parser.add_argument('--kggsee-jar',
                       default='~/software/kmmsee/kggsee.jar',
                       help='KGGSee JAR文件路径|KGGSee JAR file path')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    try:
        analyzer = GECAnalyzer(
            pfile=args.pfile,
            reference=args.reference,
            output_dir=args.output_dir,
            threads=args.threads,
            memory=args.memory,
            filter_maf_le=args.filter_maf_le,
            p_value_cutoff=args.p_value_cutoff,
            keep_ref=args.keep_ref,
            chromosomes=args.chrom,
            chrom_col=args.chrom_col,
            pos_col=args.pos_col,
            p_col=args.p_col,
            alpha=args.alpha,
            kggsee_jar=args.kggsee_jar,
            convert_chrom_format=args.convert_chrom_format
        )

        results = analyzer.run_analysis()

        # 打印最终结果摘要|Print final result summary
        print("\n" + "=" * 80)
        print("分析结果摘要|Analysis Result Summary")
        print("=" * 80)
        print(f"总有效检验数|Total effective tests: {results['total_effective_tests']:.2f}")
        print(f"校正后显著性阈值|Adjusted significance threshold: {results['significance_threshold']:.2e}")
        print("=" * 80)

    except Exception as e:
        print(f"错误|Error: {str(e)}")
        sys.exit(1)


if __name__ == "__main__":
    main()
