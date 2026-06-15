"""
OcBSA - F2/RILs群体BSA分析|OcBSA - F2/RILs Population BSA Analysis

支持SNP-index和ED(Euclidean Distance)两种分析方法。
"""

import gzip
import multiprocessing
import numpy as np


class F2bsaCalculator:
    """F2/RILs群体BSA分析计算器|F2/RILs Population BSA Analysis Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def _get_geno(self, geno_str):
        """提取基因型AD值|Extract genotype AD values"""
        num1 = geno_str.split(':')[1].split(',')[0]
        num2 = geno_str.split(':')[1].split(',')[1]
        return num1 + '|' + num2

    def read_vcf_data(self, infile, p1, p2, b1, b2):
        """读取VCF文件(F2群体)|Read VCF file for F2 population"""
        out_dict = {}
        open_fn = gzip.open if infile.endswith('.gz') else open
        with open_fn(infile, 'rt') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                line_list = line.strip().split()
                if ',' in line_list[3] or ',' in line_list[4]:
                    continue
                # F2群体过滤|F2 population filter
                if (line_list[p1].split(':')[0] == line_list[p2].split(':')[0]
                        or line_list[p2].startswith('./.')
                        or line_list[p1].startswith('./.')):
                    continue
                elif (line_list[b2].startswith('./.') or line_list[b1].startswith('./.')):
                    continue
                elif (line_list[p2].startswith('0/1') or line_list[p2].startswith('0|1')
                      or line_list[p1].startswith('0/1') or line_list[p1].startswith('0|1')):
                    continue

                p1_geno = self._get_geno(line_list[p1])
                if line_list[p2].startswith('./.'):
                    p2_geno = '0|0'
                else:
                    p2_geno = self._get_geno(line_list[p2])
                b1_geno = self._get_geno(line_list[b1])
                b2_geno = self._get_geno(line_list[b2])

                if line_list[0] not in out_dict:
                    out_dict[line_list[0]] = []
                out_dict[line_list[0]].append(
                    [int(line_list[1]), line_list[3], line_list[4],
                     p1_geno, p2_geno, b1_geno, b2_geno]
                )
        return out_dict

    def read_table_file(self, infile):
        """读取table格式文件|Read table format file"""
        out_dict = {}
        with open(infile, 'r') as table_file:
            for line in table_file:
                line_list = line.strip().split('\t')
                if line_list[0] not in out_dict:
                    out_dict[line_list[0]] = []
                ad_list = [int(line_list[1])] + line_list[2:]
                out_dict[line_list[0]].append(ad_list)
        return out_dict

    @staticmethod
    def pre_binom_test(rep):
        """预计算二项分布阈值(F2)|Pre-compute binomial thresholds for F2"""
        out_dict = {}
        for n in range(3, 4000, 1):
            num = np.percentile(np.random.binomial(n, 0.5, size=rep), 99.99)
            out_dict[n] = num
        return out_dict

    @staticmethod
    def _get_dep(geno):
        """提取覆盖度和基因型|Extract depth and genotype"""
        dep = int(geno.split('|')[0]) + int(geno.split('|')[1])
        geno_list = [int(geno.split('|')[0]), int(geno.split('|')[1])]
        return dep, geno_list

    def filter_data(self, key, chr_list, dep_lim, dep_high, dep_lim_pool,
                    dep_high_pool, binom_dict, snp_index, ED):
        """过滤SNP并计算SNP-index或ED值|Filter SNPs and calculate SNP-index or ED values"""
        out_snp_list = []

        for vcf_list in chr_list:
            p1_dep, p1_geno = self._get_dep(vcf_list[3])
            p2_dep, p2_geno = self._get_dep(vcf_list[4])
            b1_dep, b1_geno = self._get_dep(vcf_list[5])
            b2_dep, b2_geno = self._get_dep(vcf_list[6])

            # 过滤覆盖度|Filter by depth
            if (p1_dep < dep_lim or p1_dep > dep_high
                    or p2_dep < dep_lim or p2_dep > dep_high
                    or b1_dep < dep_lim_pool or b1_dep > dep_high_pool
                    or b2_dep < dep_lim_pool or b2_dep > dep_high_pool):
                continue

            if p1_geno[1] > p1_geno[0] and p2_geno[1] > p2_geno[0]:
                continue
            if p1_geno[1] < p1_geno[0] and p2_geno[1] < p2_geno[0]:
                continue
            if p1_geno[1] != 0 and p1_geno[0] != 0:
                continue
            if p2_geno[1] != 0 and p2_geno[0] != 0:
                continue

            if p2_geno[1] > p2_geno[0]:
                p1_num, p2_num = 0, 1
            else:
                p1_num, p2_num = 1, 0

            if len(vcf_list[1]) > 1 or len(vcf_list[2]) > 1:
                flag = 'indel'
            else:
                flag = 'snp'

            if snp_index:
                b1_index = b1_geno[p1_num] / b1_dep
                b2_index = b2_geno[p1_num] / b2_dep
                total_index = b1_geno[p1_num] / b1_dep - b2_geno[p1_num] / b2_dep
                num1 = binom_dict[int(b1_dep)]
                num2 = binom_dict[int(b2_dep)]
                virtual_index = abs(num1 / b1_dep - (b2_dep - num2) / b2_dep)
            elif ED:
                b1_index = (b1_geno[p1_num] / b1_dep - b2_geno[p1_num] / b2_dep) ** 2
                b2_index = (b1_geno[p2_num] / b1_dep - b2_geno[p2_num] / b2_dep) ** 2
                total_index = b1_index + b2_index
                num1 = binom_dict[int(b1_dep)]
                num2 = binom_dict[int(b2_dep)]
                index1 = (num1 / b1_dep - (b2_dep - num2) / b2_dep) ** 2
                index2 = ((b1_dep - num1) / b1_dep - num2 / b2_dep) ** 2
                virtual_index = index1 + index2

            out_snp_list.append(vcf_list + [virtual_index, b1_index, b2_index, total_index, flag])

        return [key, out_snp_list]

    @staticmethod
    def cal_dis_func(key, win_size, loop_value):
        """滑窗平滑计算(F2)|Sliding window smoothing for F2"""
        step = win_size / 10
        win_list = [0, win_size]
        total_list = []
        start_pos = 0

        while win_list[1] < loop_value[-1][0]:
            this_list = [[], [], [], []]
            for loop_list in loop_value[start_pos:]:
                i_list = [int(loop_list[0]), float(loop_list[-5]), float(loop_list[-4]),
                          float(loop_list[-3]), float(loop_list[-2])]
                if win_list[0] <= i_list[0] < win_list[1]:
                    this_list[0].append(i_list[1])
                    this_list[1].append(i_list[2])
                    this_list[2].append(i_list[3])
                    this_list[3].append(i_list[4])
                elif i_list[0] >= win_list[1]:
                    try:
                        total_list.append([
                            (win_list[0] + win_list[1]) / 2,
                            sum(this_list[0]) / len(this_list[0]),
                            sum(this_list[1]) / len(this_list[0]),
                            sum(this_list[2]) / len(this_list[0]),
                            sum(this_list[3]) / len(this_list[0]),
                            len(this_list[0]),
                        ])
                    except ZeroDivisionError:
                        total_list.append([(win_list[0] + win_list[1]) / 2, 0, 0, 0, 0, 0])
                    win_list[0] += step
                    win_list[1] += step
                    break
                elif i_list[0] < win_list[0]:
                    start_pos += 1

        return [key, total_list]

    def _detect_samples(self, infile):
        """从VCF header检测样本名和偏移量|Detect sample names and offset from VCF header"""
        open_fn = gzip.open if infile.endswith('.gz') else open
        sample_offset = 0
        sample_names = []
        with open_fn(infile, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    header_cols = line.strip().split('\t')
                    sample_offset = header_cols.index('FORMAT') + 1
                    sample_names = header_cols[sample_offset:]
                    break
        return sample_offset, sample_names

    def run(self):
        """运行F2群体BSA分析流程|Run F2 population BSA analysis pipeline"""
        cfg = self.config
        snp_index = cfg.method == "snpindex"
        ED = cfg.method == "ED"

        # 检测VCF样本列偏移|Detect VCF sample column offset
        sample_offset, sample_names = self._detect_samples(cfg.input_vcf)
        if not sample_names:
            self.logger.error("未检测到样本列|No sample columns detected in VCF header")
            return False
        total_samples = len(sample_names)
        for col_name, col_val in [('p1', cfg.parent1), ('p2', cfg.parent2),
                                   ('b1', cfg.pool1), ('b2', cfg.pool2)]:
            if col_val < 1 or col_val > total_samples:
                self.logger.error(f"{col_name}={col_val} 超出范围，共{total_samples}个样本|{col_name}={col_val} out of range, {total_samples} samples total")
                return False
        self.logger.info(f"检测到|Detected {total_samples} 个样本: {', '.join(sample_names)}")
        self.logger.info(f"列号|Column mapping: p1={cfg.parent1}({sample_names[cfg.parent1-1]}), "
                         f"p2={cfg.parent2}({sample_names[cfg.parent2-1]}), "
                         f"b1={cfg.pool1}({sample_names[cfg.pool1-1]}), "
                         f"b2={cfg.pool2}({sample_names[cfg.pool2-1]})")

        p1 = sample_offset + cfg.parent1 - 1
        p2 = sample_offset + cfg.parent2 - 1
        b1 = sample_offset + cfg.pool1 - 1
        b2 = sample_offset + cfg.pool2 - 1

        # 预计算二项分布阈值|Pre-compute binomial thresholds
        self.logger.info("预计算二项分布阈值|Pre-computing binomial thresholds")
        binom_dict = self.pre_binom_test(10000)

        # 读取VCF数据|Read VCF data
        self.logger.info(f"读取VCF文件|Reading VCF file: {cfg.input_vcf}")
        chr_geno_dict = self.read_vcf_data(cfg.input_vcf, p1, p2, b1, b2)
        self.logger.info(f"读取到|Read {len(chr_geno_dict)} 条染色体数据")

        # 过滤并计算|Filter and calculate
        method_name = "SNP-index" if snp_index else "ED"
        self.logger.info(f"开始计算{method_name}值|Starting {method_name} calculation")
        chr_geno_dict_filter = {}
        pool = multiprocessing.Pool(len(chr_geno_dict))
        results = []
        for key, value in chr_geno_dict.items():
            self.logger.info(f"处理染色体|Processing chromosome: {key}")
            results.append(pool.apply_async(
                self.filter_data,
                (key, value, cfg.parent_min_dep, cfg.parent_max_dep,
                 cfg.pool_min_dep, cfg.pool_max_dep, binom_dict, snp_index, ED)
            ))
        pool.close()
        pool.join()

        for result in results:
            chr_list = result.get()
            chr_geno_dict_filter[chr_list[0]] = chr_list[1]

        # 输出中间文件|Output intermediate file
        if snp_index:
            mid_file = str(cfg.output_path / "f2bsa.snpindex")
        else:
            mid_file = str(cfg.output_path / "f2bsa.ED")

        with open(mid_file, 'w') as out_file:
            for key, value in chr_geno_dict_filter.items():
                for i_list in value:
                    out_file.write(key + '\t')
                    for i in i_list:
                        out_file.write(str(i) + '\t')
                    out_file.write('\n')
        self.logger.info(f"{method_name}结果已保存|{method_name} results saved to: {mid_file}")

        # 滑窗平滑|Sliding window smoothing
        self.logger.info(f"开始滑窗平滑|Starting sliding window smoothing (window={cfg.window_size})")
        pool = multiprocessing.Pool(len(chr_geno_dict_filter))
        result = []
        for key, value in chr_geno_dict_filter.items():
            self.logger.info(f"平滑染色体|Smoothing chromosome: {key}")
            result.append(pool.apply_async(self.cal_dis_func, (key, cfg.window_size, value)))
        pool.close()
        pool.join()

        smoothed_file = str(cfg.output_path / "f2bsa.smoothed")
        with open(smoothed_file, 'w') as file3:
            for i_list in result:
                i_list = i_list.get()
                for i in i_list[1][1:]:
                    file3.write(i_list[0] + '\t' + str(int(i[0])) + '\t')
                    for x in i[1:]:
                        file3.write(str(x) + '\t')
                    file3.write('\n')
        self.logger.info(f"滑窗结果已保存|Smoothed results saved to: {smoothed_file}")

        # 输出R绘图用汇总数据|Output summary data for R plotting
        import re
        summary_file = str(cfg.output_path / "f2bsa.summary.tsv")
        chr_sorted = sorted(chr_geno_dict_filter.keys(), key=lambda x: int(re.findall(r'\d+', x)[0]))
        chr_len_dict = {}
        for ch in chr_sorted:
            chr_len_dict[ch] = chr_geno_dict_filter[ch][-1][0]

        with open(summary_file, 'w') as sf:
            sf.write('Chr\tPos\tChrLen\tVirtualIndex\tPool1Index\tPool2Index\tDeltaIndex\tMarkerNum\n')
            for i_list in result:
                i_list = i_list.get()
                chr_name = i_list[0]
                for row in i_list[1][1:]:
                    sf.write(f'{chr_name}\t{int(row[0])}\t{chr_len_dict[chr_name]}'
                             f'\t{row[1]}\t{row[2]}\t{row[3]}\t{row[4]}\t{row[5]}\n')
        self.logger.info(f"R绘图汇总数据已保存|Summary data for R plotting saved to: {summary_file}")

        # 自动绘图|Auto plotting
        self.logger.info("开始绘图|Starting figure plotting")
        from .config import BsaFigConfig
        from .fig import BsaFigPlotter

        fig_type = 'snpindex' if snp_index else 'ed'
        fig_file = str(cfg.output_path / f"f2bsa.{method_name}.png")
        fig_config = BsaFigConfig(
            input_file=smoothed_file,
            output_file=fig_file,
            plot_type=fig_type,
        )
        fig_config.validate()
        plotter = BsaFigPlotter(fig_config, self.logger)
        plotter.run()

        self.logger.info(f"F2群体{method_name}分析完成|F2 population {method_name} analysis completed")
        return True
