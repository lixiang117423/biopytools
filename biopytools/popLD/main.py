"""
PopLDdecay主程序模块 | PopLDdecay Main Module
用于连锁不平衡(LD)衰减分析 | Linkage Disequilibrium Decay Analysis
"""

import os
import sys
import argparse
import time
from .config import PopLDdecayConfig
from .utils import PopLDdecayLogger, CommandRunner


class PopLDdecayRunner:
    """PopLDdecay运行器类 | PopLDdecay Runner Class"""

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

    def run_analysis(self):
        """运行LD衰减分析 | Run LD decay analysis"""
        self.logger.info("=" * 60)
        self.logger.info("PopLDdecay: 连锁不平衡衰减分析")
        self.logger.info("PopLDdecay: Linkage Disequilibrium Decay Analysis")
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
            description="PopLDdecay LD衰减分析 | PopLDdecay LD decay analysis"
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
            else:
                self.logger.warning(f"预期输出文件未找到 | Expected output file not found: {expected_output}")

            self.logger.info("=" * 60)
            return True
        else:
            self.logger.error("=" * 60)
            self.logger.error("分析失败 | Analysis failed")
            self.logger.error("=" * 60)
            return False


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
        runner = PopLDdecayRunner(
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
