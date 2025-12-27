"""
PopLD主程序模块 | PopLD Main Module
用于连锁不平衡(LD)衰减分析 | Linkage Disequilibrium Decay Analysis
"""

import os
import sys
import argparse
import time
from collections import defaultdict
from .config import PopLDdecayConfig
from .utils import PopLDdecayLogger, CommandRunner


class PopLDRunner:
    """PopLD运行器类 | PopLD Runner Class"""

    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PopLDdecayConfig(**kwargs)
        self.config.validate()

        # 获取输出目录 | Get output directory
        output_dir = os.path.dirname(self.config.output_stat)

        # 初始化日志 | Initialize logging
        self.logger_manager = PopLDdecayLogger(output_dir)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, output_dir)

    def build_command(self) -> str:
        """构建PopLDdecay命令 | Build PopLDdecay command"""
        cmd_parts = [self.config.poplddecay_path]

        # 输入VCF文件 | Input VCF file
        cmd_parts.extend(["-InVCF", self.config.input_vcf])

        # 输出统计文件 | Output stat file
        cmd_parts.extend(["-OutStat", self.config.output_stat])

        # 亚群样本列表 | Subgroup sample list
        if self.config.sub_pop:
            cmd_parts.extend(["-SubPop", self.config.sub_pop])

        # 最大距离 | Max distance
        cmd_parts.extend(["-MaxDist", str(self.config.max_dist)])

        # 最小次等位基因频率 | Min minor allele frequency
        cmd_parts.extend(["-MAF", str(self.config.maf)])

        # 最大杂合率 | Max het ratio
        cmd_parts.extend(["-Het", str(self.config.het)])

        # 最大缺失率 | Max miss ratio
        cmd_parts.extend(["-Miss", str(self.config.miss)])

        # EHH分析 | EHH analysis
        if self.config.ehh:
            cmd_parts.extend(["-EHH", self.config.ehh])

        # 输出过滤SNP | Output filter SNP
        if self.config.out_filter_snp:
            cmd_parts.append("-OutFilterSNP")

        # 输出类型 | Output type
        cmd_parts.extend(["-OutType", str(self.config.out_type)])

        # 算法方法 | Algorithm method
        cmd_parts.extend(["-Method", str(self.config.method)])

        return " ".join(cmd_parts)

    def parse_sample_groups(self, sample_file: str) -> dict:
        """解析样本分组文件（两列格式：样本ID 分组名）| Parse sample group file (two columns: SampleID GroupName)"""
        groups = defaultdict(list)

        with open(sample_file, 'r') as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = line.split('\t') if '\t' in line else line.split()
                if len(parts) >= 2:
                    sample_id = parts[0]
                    group_name = parts[1]
                    groups[group_name].append(sample_id)
                elif len(parts) == 1:
                    # 只有一列，归入默认组 | Single column, put in default group
                    groups['default'].append(parts[0])

        return dict(groups)

    def run_group_analysis(self):
        """按分组运行LD衰减分析 | Run LD decay analysis by groups"""
        # 解析样本分组 | Parse sample groups
        self.logger.info("=" * 60)
        self.logger.info("PopLD: 连锁不平衡衰减分析（分组模式）")
        self.logger.info("PopLD: Linkage Disequilibrium Decay Analysis (Group Mode)")
        self.logger.info("=" * 60)

        self.logger.info(f"输入VCF | Input VCF: {self.config.input_vcf}")
        self.logger.info(f"样本分组文件 | Sample group file: {self.config.sub_pop}")

        groups = self.parse_sample_groups(self.config.sub_pop)

        if not groups:
            self.logger.error("无法解析样本分组 | Failed to parse sample groups")
            return False

        self.logger.info(f"检测到 {len(groups)} 个分组 | Found {len(groups)} groups:")
        for group_name, samples in groups.items():
            self.logger.info(f"  {group_name}: {len(samples)} 个样本 | samples")

        if len(groups) == 1 and 'default' in groups:
            self.logger.warning("样本文件只有一列，将作为单一群体处理 | Single column file, treating as single population")
            return self.run_single_analysis()

        # 创建临时目录存放分组样本列表 | Create temp directory for group sample lists
        output_dir = os.path.dirname(self.config.output_stat)
        temp_dir = os.path.join(output_dir, 'temp_groups')
        os.makedirs(temp_dir, exist_ok=True)

        # 为每个分组生成样本列表文件并运行PopLDdecay | Generate sample list for each group and run PopLDdecay
        group_results = []

        for group_name, samples in sorted(groups.items()):
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info(f"处理分组 | Processing group: {group_name} ({len(samples)} samples)")
            self.logger.info("=" * 60)

            # 生成分组样本列表文件 | Generate group sample list file
            group_sample_file = os.path.join(temp_dir, f"{group_name}_samples.txt")
            with open(group_sample_file, 'w') as f:
                for sample in samples:
                    f.write(f"{sample}\n")

            # 生成分组输出文件前缀 | Generate group output prefix
            group_output = f"{self.config.output_stat}.{group_name}"

            # 构建命令 | Build command
            cmd_parts = [
                self.config.poplddecay_path,
                "-InVCF", self.config.input_vcf,
                "-OutStat", group_output,
                "-SubPop", group_sample_file,
                "-MaxDist", str(self.config.max_dist),
                "-MAF", str(self.config.maf),
                "-Het", str(self.config.het),
                "-Miss", str(self.config.miss),
                "-OutType", str(self.config.out_type),
                "-Method", str(self.config.method)
            ]
            cmd = " ".join(cmd_parts)

            # 运行分析 | Run analysis
            start_time = time.time()
            success = self.cmd_runner.run(
                cmd,
                description=f"{group_name} LD衰减分析 | {group_name} LD decay analysis"
            )
            elapsed_time = time.time() - start_time

            if success:
                expected_output = f"{group_output}.stat.gz"
                if os.path.exists(expected_output):
                    file_size = os.path.getsize(expected_output)
                    self.logger.info(f"输出文件 | Output file: {expected_output}")
                    self.logger.info(f"文件大小 | File size: {file_size:,} bytes")
                    self.logger.info(f"运行时间 | Runtime: {elapsed_time:.2f} seconds")
                    group_results.append((group_name, expected_output))
                else:
                    self.logger.warning(f"预期输出文件未找到 | Expected output file not found: {expected_output}")
            else:
                self.logger.error(f"分组 {group_name} 分析失败 | Group {group_name} analysis failed")

        # 清理临时文件 | Clean up temporary files
        import shutil
        if os.path.exists(temp_dir):
            shutil.rmtree(temp_dir)

        # 生成分组合并列表文件用于绘图 | Generate merged list file for plotting
        if group_results:
            plot_list_file = f"{self.config.output_stat}.groups.list"
            with open(plot_list_file, 'w') as f:
                for group_name, result_file in group_results:
                    f.write(f"{result_file}\t{group_name}\n")

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("所有分组分析完成 | All group analyses completed")
            self.logger.info(f"成功分析 {len(group_results)} 个分组 | Successfully analyzed {len(group_results)} groups")
            self.logger.info(f"绘图列表文件 | Plot list file: {plot_list_file}")

            # 合并所有分组结果到一个文件 | Merge all group results into one file
            self.logger.info("")
            self.logger.info("合并分组结果 | Merging group results...")
            merged_file = self.merge_group_results(group_results)

            if merged_file:
                self.logger.info(f"合并结果文件 | Merged result file: {merged_file}")

            self.logger.info("")
            self.logger.info("使用以下命令绘图 | Use the following command to plot:")
            self.logger.info(f"  perl <PopLDdecay_path>/bin/Plot_MutiPop.pl -inList {plot_list_file} -output {self.config.output_stat}")
            self.logger.info("=" * 60)

            return len(group_results) > 0

        return False

    def merge_group_results(self, group_results: list) -> str:
        """合并分组结果到一个文件，在最后一列添加分组信息 | Merge group results into one file, add group info in last column"""
        import gzip

        merged_file = f"{self.config.output_stat}.merged.stat.gz"

        try:
            with gzip.open(merged_file, 'wt') as outf:
                # 写入表头（添加Group列）| Write header (add Group column)
                header = "#Dist\tMean_r^2\tMean_D'\tSum_r^2\tSum_D'\tNumberPairs\tGroup\n"
                outf.write(header)

                # 读取并合并每个分组的数据 | Read and merge data from each group
                for group_name, result_file in group_results:
                    try:
                        with gzip.open(result_file, 'rt') as inf:
                            # 跳过原表头 | Skip original header
                            header_line = inf.readline()
                            if not header_line:
                                continue

                            # 读取数据行并添加分组信息 | Read data lines and add group info
                            for line in inf:
                                line = line.strip()
                                if line:
                                    outf.write(f"{line}\t{group_name}\n")

                        self.logger.debug(f"  已合并 | Merged: {group_name} ({os.path.basename(result_file)})")
                    except Exception as e:
                        self.logger.warning(f"  跳过文件 | Skip file {result_file}: {e}")

            # 统计合并的记录数 | Count merged records
            with gzip.open(merged_file, 'rt') as f:
                line_count = sum(1 for _ in f) - 1  # 减去表头行 | Subtract header line

            self.logger.info(f"  合并总记录数 | Total merged records: {line_count:,}")
            return merged_file

        except Exception as e:
            self.logger.error(f"合并结果失败 | Failed to merge results: {e}")
            return None

    def run_single_analysis(self):
        """运行单个群体的LD衰减分析 | Run single population LD decay analysis"""
        self.logger.info("=" * 60)
        self.logger.info("PopLD: 连锁不平衡衰减分析")
        self.logger.info("PopLD: Linkage Disequilibrium Decay Analysis")
        self.logger.info("=" * 60)

        # 输出配置信息 | Output configuration
        self.logger.info(f"输入VCF | Input VCF: {self.config.input_vcf}")
        self.logger.info(f"输出文件 | Output file: {self.config.output_stat}")
        self.logger.info(f"软件路径 | Software path: {self.config.poplddecay_path}")

        if self.config.sub_pop:
            self.logger.info(f"亚群样本列表 | Subgroup sample list: {self.config.sub_pop}")

        self.logger.info(f"最大距离 | Max distance: {self.config.max_dist} kb")
        self.logger.info(f"最小MAF | Min MAF: {self.config.maf}")
        self.logger.info(f"最大杂合率 | Max het ratio: {self.config.het}")
        self.logger.info(f"最大缺失率 | Max miss ratio: {self.config.miss}")
        self.logger.info(f"输出类型 | Output type: {self.config.out_type}")
        self.logger.info(f"算法方法 | Algorithm method: {self.config.method}")
        self.logger.info("=" * 60)

        # 构建并执行命令 | Build and execute command
        cmd = self.build_command()
        self.logger.info("")
        self.logger.info("开始LD衰减分析 | Starting LD decay analysis")
        self.logger.info("=" * 60)

        start_time = time.time()

        success = self.cmd_runner.run(
            cmd,
            description="PopLD LD衰减分析 | PopLD LD decay analysis"
        )

        elapsed_time = time.time() - start_time

        if success:
            self.logger.info("=" * 60)
            self.logger.info("分析完成 | Analysis completed")
            self.logger.info(f"运行时间 | Runtime: {elapsed_time:.2f} seconds")

            # 检查输出文件 | Check output files
            expected_output = f"{self.config.output_stat}.stat.gz"
            if os.path.exists(expected_output):
                file_size = os.path.getsize(expected_output)
                self.logger.info(f"输出文件 | Output file: {expected_output}")
                self.logger.info(f"文件大小 | File size: {file_size:,} bytes")

                # 提供绘图命令 | Provide plotting command
                self.logger.info("")
                self.logger.info("使用以下命令绘图 | Use the following command to plot:")
                self.logger.info(f"  perl <PopLDdecay_path>/bin/Plot_OnePop.pl -inFile {expected_output} -output {self.config.output_stat}")
            else:
                self.logger.warning(f"预期输出文件未找到 | Expected output file not found: {expected_output}")

            self.logger.info("=" * 60)
            return True
        else:
            self.logger.error("=" * 60)
            self.logger.error("分析失败 | Analysis failed")
            self.logger.error("=" * 60)
            return False

    def run_analysis(self):
        """运行LD衰减分析（自动检测分组）| Run LD decay analysis (auto-detect groups)"""
        # 检查是否为两列格式的样本文件 | Check if it's a two-column sample file
        if self.config.sub_pop:
            try:
                with open(self.config.sub_pop, 'r') as f:
                    first_line = f.readline().strip()
                    if first_line and not first_line.startswith('#'):
                        parts = first_line.split('\t') if '\t' in first_line else first_line.split()
                        if len(parts) >= 2:
                            # 两列格式，使用分组分析 | Two-column format, use group analysis
                            return self.run_group_analysis()
            except Exception as e:
                self.logger.warning(f"无法检测样本文件格式，使用标准模式 | Cannot detect sample file format, using standard mode: {e}")

        # 标准单群体分析 | Standard single population analysis
        return self.run_single_analysis()


def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='PopLDdecay: 连锁不平衡(LD)衰减分析工具 | PopLDdecay: Linkage Disequilibrium Decay Analysis Tool',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='''
示例 | Examples:

  # 基本用法
  %(prog)s -i snp.vcf.gz -o LD_result

  # 使用亚群样本
  %(prog)s -i snp.vcf.gz -o LD_result -s subgroup.txt

  # 自定义参数
  %(prog)s -i snp.vcf.gz -o LD_result -d 500 -m 0.01 -t 2

  # 指定软件路径
  %(prog)s -i snp.vcf.gz -o LD_result -p /path/to/PopLDdecay
        '''
    )

    # 必需参数 | Required parameters
    required = parser.add_argument_group('必需参数 | Required arguments')
    required.add_argument('-i', '--input-vcf', required=True,
                         help='输入VCF文件路径 | Input VCF file path')
    required.add_argument('-o', '--output-stat', required=True,
                         help='输出统计文件前缀 | Output statistics file prefix')

    # 软件路径 | Software path
    software = parser.add_argument_group('软件配置 | Software configuration')
    software.add_argument('-p', '--poplddecay-path',
                         default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/poplddecay_v.3.43/bin/PopLDdecay',
                         help='PopLDdecay软件路径 | PopLDdecay software path')

    # 处理参数 | Processing parameters
    processing = parser.add_argument_group('处理参数 | Processing parameters')
    processing.add_argument('-s', '--sub-pop',
                           help='亚群样本列表文件 | Subgroup sample list file')
    processing.add_argument('-d', '--max-dist', type=int, default=300,
                           help='SNP间最大距离(kb) | Max distance between SNPs in kb (default: 300)')
    processing.add_argument('-m', '--maf', type=float, default=0.005,
                           help='最小次等位基因频率 | Min minor allele frequency (default: 0.005)')
    processing.add_argument('--het', type=float, default=0.88,
                           help='最大杂合率 | Max ratio of het allele (default: 0.88)')
    processing.add_argument('--miss', type=float, default=0.25,
                           help='最大缺失率 | Max ratio of miss allele (default: 0.25)')
    processing.add_argument('--ehh',
                           help='EHH起始位点 | EHH start site (NA for disabled)')
    processing.add_argument('--out-type', type=int, choices=[1, 2, 3, 4, 5, 6, 7, 8], default=1,
                           help='输出类型 | Output type (default: 1)')
    processing.add_argument('--method', type=int, choices=[1, 2], default=1,
                           help='算法方法(1=快速, 2=可能需要更多内存) | Algorithm method, 1=fast, 2=may use more memory (default: 1)')

    # 其他选项 | Other options
    other = parser.add_argument_group('其他选项 | Other options')
    other.add_argument('--out-filter-snp', action='store_true',
                      help='输出最终用于计算的SNP | Output final SNP used for calculation')
    other.add_argument('--log-file',
                      help='日志文件路径 | Log file path')
    other.add_argument('-v', '--verbose', action='count', default=0,
                      help='详细输出模式 | Verbose mode (-v: INFO, -vv: DEBUG)')

    args = parser.parse_args()

    try:
        # 创建运行器 | Create runner
        runner = PopLDRunner(
            input_vcf=args.input_vcf,
            output_stat=args.output_stat,
            poplddecay_path=args.poplddecay_path,
            sub_pop=args.sub_pop,
            max_dist=args.max_dist,
            maf=args.maf,
            het=args.het,
            miss=args.miss,
            ehh=args.ehh,
            out_filter_snp=args.out_filter_snp,
            out_type=args.out_type,
            method=args.method
        )

        # 运行分析 | Run analysis
        success = runner.run_analysis()

        if not success:
            sys.exit(1)

    except KeyboardInterrupt:
        runner.logger.warning("用户中断 | User interrupted")
        sys.exit(130)
    except Exception as e:
        if hasattr(runner, 'logger'):
            runner.logger.error(f"程序执行出错 | Program execution error: {str(e)}", exc_info=True)
        else:
            print(f"ERROR: {str(e)}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
