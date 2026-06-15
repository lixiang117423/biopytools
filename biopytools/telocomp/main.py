"""
TeloComp主程序模块|TeloComp Main Module
"""

import os
import sys
from pathlib import Path
from .config import TeloCompConfig
from .utils import TeloCompLogger, CondaCommandRunner, check_genome_index, create_genome_index


class TeloCompAnalyzer:
    """TeloComp端粒鉴定分析主类|TeloComp Telomere Identification Analyzer Main Class"""

    def __init__(self, **kwargs):
        """初始化|Initialize"""
        # 初始化配置|Initialize configuration
        self.config = TeloCompConfig(**kwargs)
        self.config.validate()

        # 初始化日志|Initialize logging
        log_file = os.path.join(self.config.output_dir, 'telocomp.log')
        self.logger_manager = TeloCompLogger(log_file)
        self.logger = self.logger_manager.get_logger()

        # 初始化命令执行器|Initialize command runner
        self.cmd_runner = CondaCommandRunner(self.logger, self.config.conda_env)

        # 设置环境变量|Set environment variables
        os.environ['TELOCOMP_BIN'] = self.config.telocomp_bin
        os.environ['GENOMESYN_BIN'] = self.config.genomesyn_bin

        # 输出文件路径|Output file paths
        self.filter1_ont_bam = os.path.join(self.config.output_dir, 'filter1_ont.bam')
        self.filter1_hifi_bam = os.path.join(self.config.output_dir, 'filter1_hifi.bam')
        self.filter2_output = os.path.join(self.config.output_dir, 'filter2_output')

        self.logger.info("TeloComp分析器初始化完成|TeloComp analyzer initialized")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")

    def run(self):
        """运行完整的端粒鉴定流程|Run complete telomere identification pipeline"""
        self.logger.info("="*60)
        self.logger.info("开始TeloComp端粒鉴定流程|Starting TeloComp telomere identification pipeline")
        self.logger.info("="*60)

        try:
            # 步骤1: 检查基因组索引|Step 1: Check genome index
            self._check_and_create_index()

            # 步骤2: Filter_1 - 过滤含端粒的reads|Step 2: Filter_1 - Filter telomere-containing reads
            if not self.config.skip_filter:
                self._run_filter_1()

                # 步骤3: Filter_2 - 检测和提取reads|Step 3: Filter_2 - Detect and extract reads
                self._run_filter_2()
            else:
                self.logger.info("跳过Filter步骤|Skipping Filter steps")

            # 步骤4: 可视化|Step 4: Visualization
            if self.config.run_visualization:
                self.logger.info("分析端粒位置|Analyzing telomere positions")
                self._analyze_telomere_positions()
                self._prepare_visualization_summary()
            else:
                self.logger.info("跳过可视化步骤|Skipping visualization step")

            self.logger.info("="*60)
            self.logger.info("TeloComp端粒鉴定流程完成|TeloComp pipeline completed")
            self.logger.info("="*60)

            return True

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {str(e)}")
            return False

    def _check_and_create_index(self):
        """检查并创建基因组索引|Check and create genome index"""
        self.logger.info("步骤1: 检查基因组索引|Step 1: Checking genome index")

        fai_file = f"{self.config.genome}.fai"

        if not os.path.exists(fai_file):
            self.logger.warning(f"基因组索引不存在|Genome index not found: {fai_file}")
            self.logger.info(f"自动创建基因组索引|Auto-creating genome index: {self.config.genome}")

            if not create_genome_index(self.config.genome, self.logger, self.config.conda_env):
                raise RuntimeError(f"无法创建基因组索引|Failed to create genome index: {self.config.genome}")

            self.logger.info("基因组索引创建成功|Genome index created successfully")
        else:
            self.logger.info("基因组索引已存在|Genome index exists")

    def _run_filter_1(self):
        """运行Filter_1步骤|Run Filter_1 step"""
        self.logger.info("步骤2: 运行Filter_1 - 过滤含端粒的reads|Step 2: Running Filter_1")

        # 构建参数|Build arguments
        args = [
            '--genome', self.config.genome,
            '--fai', f'{self.config.genome}.fai',
            '--threads', str(self.config.threads)
        ]

        # 添加ONT数据|Add ONT data
        if self.config.ont:
            args.extend(['--ont', self.config.ont])
            args.extend(['--Ob', self.filter1_ont_bam])

        # 添加HiFi数据|Add HiFi data
        if self.config.hifi:
            args.extend(['--hifi', self.config.hifi])
            args.extend(['--Hb', self.filter1_hifi_bam])

        # 添加motifs|Add motifs
        if self.config.motifs:
            args.extend(['--motifs'] + self.config.motifs)

        # 添加其他参数|Add other parameters
        args.extend([
            '--max_break', str(self.config.max_break),
            '--min_clip', str(self.config.min_clip)
        ])

        # 运行命令|Run command
        result = self.cmd_runner.run_telocomp_command('telocomp_Filter_1', args)

        if result.returncode != 0:
            raise RuntimeError("Filter_1执行失败|Filter_1 execution failed")

        self.logger.info("Filter_1执行完成|Filter_1 completed")

    def _run_filter_2(self):
        """运行Filter_2步骤|Run Filter_2 step"""
        # 检查是否同时有ONT和HiFi数据|Check if both ONT and HiFi data are available
        has_ont = self.config.ont and os.path.exists(self.filter1_ont_bam)
        has_hifi = self.config.hifi and os.path.exists(self.filter1_hifi_bam)

        if not (has_ont and has_hifi):
            self.logger.warning("="*60)
            self.logger.warning("Filter_2需要同时提供ONT和HiFi数据|Filter_2 requires both ONT and HiFi data")
            self.logger.warning("当前只有单一数据类型，跳过Filter_2步骤|Current has only single data type, skipping Filter_2")
            self.logger.warning("提示|Tip: Filter_1已完成端粒reads过滤，可直接使用filter1输出文件")
            self.logger.warning("="*60)

            # 创建空的filter2_output目录以保持结构一致性|Create empty filter2_output for consistency
            os.makedirs(self.filter2_output, exist_ok=True)
            return

        self.logger.info("步骤3: 运行Filter_2 - 检测和提取reads|Step 3: Running Filter_2")

        # 同时有ONT和HiFi数据，正常执行Filter_2|Both data types available, run Filter_2 normally
        ont_bam = self.filter1_ont_bam
        hifi_bam = self.filter1_hifi_bam

        # 构建参数|Build arguments
        args = [
            '-o', self.filter2_output,
            '-c', str(self.config.coverage),
            '-p', str(self.config.parallels),
            '--min_ratio', str(self.config.min_ratio),
            '--ont_bam', ont_bam,
            '--hifi_bam', hifi_bam
        ]

        # 运行命令|Run command
        result = self.cmd_runner.run_telocomp_command('telocomp_Filter_2', args)

        if result.returncode != 0:
            raise RuntimeError("Filter_2执行失败|Filter_2 execution failed")

        self.logger.info("Filter_2执行完成|Filter_2 completed")
        self._summarize_filter2_results()

    def _summarize_filter2_results(self):
        """汇总Filter_2结果|Summarize Filter_2 results"""
        self.logger.info("汇总Filter_2结果|Summarizing Filter_2 results")

        # 检查输出目录|Check output directories
        trim_l = os.path.join(self.filter2_output, 'trim_L')
        trim_r = os.path.join(self.filter2_output, 'trim_R')

        if os.path.exists(trim_l):
            l_files = len([f for f in os.listdir(trim_l) if f.endswith('.fa')])
            self.logger.info(f"左侧端粒reads|Left telomere reads: {l_files}个文件|files")

        if os.path.exists(trim_r):
            r_files = len([f for f in os.listdir(trim_r) if f.endswith('.fa')])
            self.logger.info(f"右侧端粒reads|Right telomere reads: {r_files}个文件|files")

    def _analyze_telomere_positions(self):
        """分析端粒在基因组上的位置|Analyze telomere positions on genome"""
        self.logger.info("开始分析端粒位置|Starting telomere position analysis")

        # 读取基因组信息|Read genome information
        try:
            # 直接读取FAI文件|Read FAI file directly
            chrom_lengths = {}
            with open(f"{self.config.genome}.fai", 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) >= 2:
                        chrom_lengths[fields[0]] = int(fields[1])
            self.logger.info(f"读取了{len(chrom_lengths)}条染色体信息|Read {len(chrom_lengths)} chromosomes")
        except Exception as e:
            self.logger.error(f"无法读取基因组索引|Failed to read genome index: {str(e)}")
            return

        # 统计端粒位置|Count telomere positions
        telomere_positions = {}

        # 分析ONT BAM|Analyze ONT BAM
        if self.config.ont and os.path.exists(self.filter1_ont_bam):
            # 创建BAM索引（如果不存在）|Create BAM index if not exists
            self._create_bam_index_if_needed(self.filter1_ont_bam)
            self._process_bam_for_positions(self.filter1_ont_bam, chrom_lengths, telomere_positions, "ONT")

        # 分析HiFi BAM|Analyze HiFi BAM
        if self.config.hifi and os.path.exists(self.filter1_hifi_bam):
            # 创建BAM索引（如果不存在）|Create BAM index if not exists
            self._create_bam_index_if_needed(self.filter1_hifi_bam)
            self._process_bam_for_positions(self.filter1_hifi_bam, chrom_lengths, telomere_positions, "HiFi")

        # 生成位置报告|Generate position report
        self._generate_position_report(telomere_positions, chrom_lengths)

    def _create_bam_index_if_needed(self, bam_file):
        """如果BAM文件没有索引，创建索引|Create BAM index if not exists"""
        try:
            import pysam
            import os

            bam_index = f"{bam_file}.bai"
            if not os.path.exists(bam_index):
                self.logger.info(f"创建BAM索引|Creating BAM index: {bam_file}")
                pysam.index(bam_file)
                self.logger.info(f"BAM索引创建完成|BAM index created")
        except Exception as e:
            self.logger.warning(f"创建BAM索引失败|Failed to create BAM index: {str(e)}")

    def _process_bam_for_positions(self, bam_file, chrom_lengths, telomere_positions, data_type):
        """处理BAM文件以统计端粒位置|Process BAM file to count telomere positions"""
        try:
            import pysam

            bam = pysam.AlignmentFile(bam_file, 'rb')
            total_reads = 0

            for read in bam.fetch():
                if read.is_unmapped:
                    continue

                total_reads += 1
                chrom = bam.get_reference_name(read.reference_id)
                pos = read.reference_start + 1  # 1-based position
                chrom_len = chrom_lengths.get(chrom, 0)

                # 使用cigartuples解析soft-clip|Use cigartuples to parse soft-clip
                if not read.cigartuples:
                    continue

                # 检查soft-clip|Check soft-clip
                left_clip = 0
                right_clip = 0

                # cigartuples: (operation, length), 4=S operation
                if read.cigartuples[0][0] == 4:  # 左侧soft-clip|Left soft-clip
                    left_clip = read.cigartuples[0][1]

                if len(read.cigartuples) > 1 and read.cigartuples[-1][0] == 4:  # 右侧soft-clip|Right soft-clip
                    right_clip = read.cigartuples[-1][1]

                # 判断是左端还是右端端粒|Determine if left or right telomere
                is_left = False
                is_right = False

                # 左端端粒：位置靠近起始位置且有左侧soft-clip|Left telomere: near start with left soft-clip
                if left_clip > 0 and pos <= 1000:
                    is_left = True

                # 右端端粒：位置靠近末端且有右侧soft-clip|Right telomere: near end with right soft-clip
                if right_clip > 0:
                    read_end = read.reference_end
                    if (chrom_len - read_end) <= 1000:
                        is_right = True

                # 记录端粒位置|Record telomere position
                if is_left or is_right:
                    if chrom not in telomere_positions:
                        telomere_positions[chrom] = {
                            'left': {'ONT': 0, 'HiFi': 0},
                            'right': {'ONT': 0, 'HiFi': 0}
                        }

                    if is_left:
                        telomere_positions[chrom]['left'][data_type] += 1
                    if is_right:
                        telomere_positions[chrom]['right'][data_type] += 1

            bam.close()
            self.logger.info(f"{data_type}: 处理了{total_reads}条reads|Processed {total_reads} reads")

        except Exception as e:
            self.logger.error(f"处理BAM文件失败|Failed to process BAM file {bam_file}: {str(e)}")
    def _generate_position_report(self, telomere_positions, chrom_lengths):
        """生成端粒位置报告|Generate telomere position report"""
        report_file = os.path.join(self.config.output_dir, 'telomere_positions.txt')

        try:
            with open(report_file, 'w') as f:
                f.write("="*80 + "\n")
                f.write("端粒位置检测报告|Telomere Position Detection Report\n")
                f.write("="*80 + "\n\n")

                if not telomere_positions:
                    f.write("未检测到端粒reads|No telomere reads detected\n")
                    f.write("请检查：\n")
                    f.write("1. 端粒motif参数是否正确（-m参数）\n")
                    f.write("2. 输入数据质量\n")
                    f.write("3. 基因组组装质量\n")
                else:
                    # 按染色体排序|Sort by chromosome
                    sorted_chroms = sorted(telomere_positions.keys(),
                                         key=lambda x: [int(c) if c.isdigit() else float('inf') for c in re.split(r'(\d+)', x)])

                    f.write(f"共在{len(sorted_chroms)}条染色体上检测到端粒\n")
                    f.write(f"Telomeres detected on {len(sorted_chroms)} chromosomes\n\n")
                    f.write("-"*80 + "\n\n")

                    for chrom in sorted_chroms:
                        chrom_len = chrom_lengths.get(chrom, 0)
                        positions = telomere_positions[chrom]

                        f.write(f"染色体|Chromosome: {chrom} (长度|Length: {chrom_len:,} bp)\n")
                        f.write("-"*80 + "\n")

                        # 左端端粒|Left telomere
                        left_ont = positions['left']['ONT']
                        left_hifi = positions['left']['HiFi']
                        left_total = left_ont + left_hifi

                        if left_total > 0:
                            f.write(f"  左端端粒|Left telomere: 1-{min(5000, chrom_len)} bp区域\n")
                            f.write(f"    ONT reads: {left_ont}\n")
                            f.write(f"    HiFi reads: {left_hifi}\n")
                            f.write(f"    总计|Total: {left_total}\n")
                        else:
                            f.write(f"  左端端粒|Left telomere: 未检测到|Not detected\n")

                        f.write("\n")

                        # 右端端粒|Right telomere
                        right_ont = positions['right']['ONT']
                        right_hifi = positions['right']['HiFi']
                        right_total = right_ont + right_hifi

                        if right_total > 0:
                            end_pos = max(chrom_len - 5000, 1)
                            f.write(f"  右端端粒|Right telomere: {end_pos}-{chrom_len} bp区域\n")
                            f.write(f"    ONT reads: {right_ont}\n")
                            f.write(f"    HiFi reads: {right_hifi}\n")
                            f.write(f"    总计|Total: {right_total}\n")
                        else:
                            f.write(f"  右端端粒|Right telomere: 未检测到|Not detected\n")

                        f.write("\n" + "="*80 + "\n\n")

                    # 统计摘要|Statistics summary
                    f.write("\n" + "="*80 + "\n")
                    f.write("统计摘要|Statistics Summary\n")
                    f.write("="*80 + "\n")

                    total_left = sum(p['left']['ONT'] + p['left']['HiFi'] for p in telomere_positions.values())
                    total_right = sum(p['right']['ONT'] + p['right']['HiFi'] for p in telomere_positions.values())
                    total_telomeres = sum(1 for p in telomere_positions.values()
                                          if (p['left']['ONT'] + p['left']['HiFi']) > 0)
                    total_telomeres += sum(1 for p in telomere_positions.values()
                                           if (p['right']['ONT'] + p['right']['HiFi']) > 0)

                    f.write(f"检测到左端端粒的染色体|Chromosomes with left telomere: {sum(1 for p in telomere_positions.values() if (p['left']['ONT'] + p['left']['HiFi']) > 0)}\n")
                    f.write(f"检测到右端端粒的染色体|Chromosomes with right telomere: {sum(1 for p in telomere_positions.values() if (p['right']['ONT'] + p['right']['HiFi']) > 0)}\n")
                    f.write(f"左端端粒reads总数|Total left telomere reads: {total_left}\n")
                    f.write(f"右端端粒reads总数|Total right telomere reads: {total_right}\n")
                    f.write(f"总端粒reads|Total telomere reads: {total_left + total_right}\n")

            self.logger.info(f"端粒位置报告已保存|Telomere position report saved to: {report_file}")

        except Exception as e:
            self.logger.error(f"生成端粒位置报告失败|Failed to generate position report: {str(e)}")

    def _prepare_visualization_summary(self):
        """准备可视化汇总|Prepare visualization summary"""
        self.logger.info("准备端粒鉴定结果汇总|Preparing telomere identification summary")

        summary_file = os.path.join(self.config.output_dir, 'telomere_summary.txt')

        try:
            with open(summary_file, 'w') as f:
                f.write("TeloComp端粒鉴定结果汇总|TeloComp Telomere Identification Summary\n")
                f.write("="*80 + "\n\n")
                f.write(f"基因组文件|Genome: {self.config.genome}\n")
                f.write(f"输出目录|Output: {self.config.output_dir}\n\n")

                # Filter结果|Filter results
                if not self.config.skip_filter:
                    f.write("过滤结果|Filter Results:\n")
                    f.write("-" * 40 + "\n")

                    # Filter_1 BAM文件统计|Filter_1 BAM file statistics
                    if self.config.ont and os.path.exists(self.filter1_ont_bam):
                        ont_count = self._count_bam_reads(self.filter1_ont_bam)
                        f.write(f"ONT端粒reads|ONT telomere reads: {ont_count}\n")
                        f.write(f"  文件|File: {self.filter1_ont_bam}\n")

                    if self.config.hifi and os.path.exists(self.filter1_hifi_bam):
                        hifi_count = self._count_bam_reads(self.filter1_hifi_bam)
                        f.write(f"HiFi端粒reads|HiFi telomere reads: {hifi_count}\n")
                        f.write(f"  文件|File: {self.filter1_hifi_bam}\n")

                    # Filter_2结果|Filter_2 results
                    trim_l = os.path.join(self.filter2_output, 'trim_L')
                    trim_r = os.path.join(self.filter2_output, 'trim_R')

                    if os.path.exists(trim_l) and os.listdir(trim_l):
                        l_files = [f for f in os.listdir(trim_l) if f.endswith('.fa')]
                        f.write(f"\n左侧端粒序列|Left telomere sequences: {len(l_files)}\n")
                        f.write(f"  目录|Directory: {trim_l}\n")

                    if os.path.exists(trim_r) and os.listdir(trim_r):
                        r_files = [f for f in os.listdir(trim_r) if f.endswith('.fa')]
                        f.write(f"\n右侧端粒序列|Right telomere sequences: {len(r_files)}\n")
                        f.write(f"  目录|Directory: {trim_r}\n")

                    # 可视化说明|Visualization note
                    f.write("\n" + "="*80 + "\n")
                    f.write("关于可视化|About Visualization:\n")
                    f.write("-"*80 + "\n")
                    f.write("完整的端粒可视化图表（telomere_plots）需要在Complement步骤生成。\n")
                    f.write("Complete telomere visualization plots are generated in the Complement step.\n\n")
                    f.write("Complement步骤需要额外的数据：\n")
                    f.write("The Complement step requires additional data:\n")
                    f.write("  - WGS reads (二代测序数据)\n")
                    f.write("  - NextPolish tool\n")
                    f.write("  - Assembled contigs from Assembly step\n\n")
                    f.write("当前模块仅实现了Filter步骤（端粒reads识别和过滤）。\n")
                    f.write("This module only implements Filter steps (telomere read identification and filtering).\n")
                    f.write("如需完整的可视化，请使用TeloComp的完整流程。\n")
                    f.write("For complete visualization, please use the full TeloComp pipeline.\n")
                    f.write("="*80 + "\n")

            self.logger.info(f"结果汇总已保存|Summary saved to: {summary_file}")
            self._print_summary_to_console()

        except Exception as e:
            self.logger.warning(f"无法创建结果汇总文件|Failed to create summary file: {str(e)}")

    def _count_bam_reads(self, bam_file):
        """统计BAM文件中的reads数量|Count reads in BAM file"""
        try:
            import pysam
            import os

            # 检查并创建BAM索引|Check and create BAM index
            bam_index = f"{bam_file}.bai"
            if not os.path.exists(bam_index):
                try:
                    pysam.index(bam_file)
                except:
                    pass  # 如果创建索引失败，继续尝试

            # 尝试使用count()（需要索引）|Try using count() (requires index)
            try:
                bam = pysam.AlignmentFile(bam_file, 'rb')
                count = bam.count()
                bam.close()
                return count
            except:
                # 如果count()失败，手动计数|If count() fails, count manually
                bam = pysam.AlignmentFile(bam_file, 'rb')
                count = 0
                for read in bam:
                    if not read.is_unmapped:
                        count += 1
                bam.close()
                return count

        except Exception as e:
            self.logger.warning(f"无法统计BAM文件|Failed to count BAM file: {str(e)}")
            return 0

    def _print_summary_to_console(self):
        """打印汇总信息到控制台|Print summary to console"""
        self.logger.info("")
        self.logger.info("="*60)
        self.logger.info("端粒鉴定结果摘要|Telomere Identification Summary")
        self.logger.info("="*60)

        if self.config.ont and os.path.exists(self.filter1_ont_bam):
            ont_count = self._count_bam_reads(self.filter1_ont_bam)
            self.logger.info(f"ONT端粒reads|ONT telomere reads: {ont_count}")

        if self.config.hifi and os.path.exists(self.filter1_hifi_bam):
            hifi_count = self._count_bam_reads(self.filter1_hifi_bam)
            self.logger.info(f"HiFi端粒reads|HiFi telomere reads: {hifi_count}")

        self.logger.info("="*60)


def main():
    """主函数|Main function"""
    import argparse

    parser = argparse.ArgumentParser(
        description='TeloComp端粒鉴定和可视化工具|TeloComp Telomere Identification and Visualization Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )

    # 必需参数|Required parameters
    parser.add_argument('-g', '--genome', required=True,
                       help='基因组FASTA文件|Genome FASTA file')
    parser.add_argument('-o', '--output-dir', required=True,
                       help='输出目录|Output directory')

    # 可选参数|Optional parameters
    parser.add_argument('--ont',
                       help='ONT数据文件|ONT data file')
    parser.add_argument('--hifi',
                       help='HiFi数据文件|HiFi data file')
    parser.add_argument('-t', '--threads', type=int, default=12,
                       help='线程数|Number of threads')
    parser.add_argument('-c', '--coverage', type=int, default=100,
                       help='覆盖度参数(0-100)|Coverage parameter (0-100)')
    parser.add_argument('-m', '--motif', default='CCCTAAA',
                       help='端粒重复序列|Telomeric repeat sequence')
    parser.add_argument('-M', '--motif-num', type=int, default=7,
                       help='端粒重复序列碱基数|Number of bases in telomere motif')

    # 流程控制|Pipeline control
    parser.add_argument('--skip-filter', action='store_true',
                       help='跳过Filter步骤|Skip filter steps')
    parser.add_argument('--no-visualization', action='store_true',
                       help='跳过可视化|Skip visualization')

    args = parser.parse_args()

    # 创建分析器并运行|Create analyzer and run
    analyzer = TeloCompAnalyzer(
        genome=args.genome,
        ont=args.ont,
        hifi=args.hifi,
        output_dir=args.output_dir,
        threads=args.threads,
        coverage=args.coverage,
        motif=args.motif,
        motif_num=args.motif_num,
        skip_filter=args.skip_filter,
        run_visualization=not args.no_visualization
    )

    success = analyzer.run()

    if not success:
        sys.exit(1)


if __name__ == "__main__":
    main()
