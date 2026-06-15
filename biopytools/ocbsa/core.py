"""
OcBSA - F1群体DHHP分析核心算法|OcBSA - F1 Population DHHP Analysis Core Algorithm

基于DHHP算法(Dominant Heterozygous Parents-based Haplotype Analysis)，
专门处理F1杂交群体中的显性基因定位。

参考|Reference: OcBSA: an NGS-based Bulk Segregant Analysis Tool for
Outcross Populations, Molecular Plant, 2024.
"""

import gzip
import multiprocessing
import numpy as np


class OcbsaCalculator:
    """OcBSA F1群体DHHP分析计算器|OcBSA F1 Population DHHP Analysis Calculator"""

    def __init__(self, config, logger):
        """初始化计算器|Initialize calculator

        Args:
            config: OcbsaConfig 配置对象|OcbsaConfig configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def _get_geno(self, geno_str):
        """提取基因型AD值|Extract genotype AD values"""
        num1 = geno_str.split(':')[1].split(',')[0]
        num2 = geno_str.split(':')[1].split(',')[1]
        return num1 + '|' + num2

    def read_vcf_data(self, infile, p1, p2, b1, b2):
        """读取VCF文件并提取基因型信息|Read VCF file and extract genotype info

        Args:
            infile: VCF文件路径|VCF file path
            p1: 显性亲本列号(0-based，已含VCF固定列偏移)|Dominant parent column (0-based, with VCF fixed column offset)
            p2: 隐性亲本列号(0-based)|Recessive parent column (0-based)
            b1: 显性混池列号(0-based)|Dominant pool column (0-based)
            b2: 隐性混池列号(0-based)|Recessive pool column (0-based)

        Returns:
            按染色体分组的基因型字典|Genotype dict grouped by chromosome
        """
        out_dict = {}
        open_fn = gzip.open if infile.endswith('.gz') else open
        with open_fn(infile, 'rt') as vcf_file:
            for line in vcf_file:
                if line.startswith('#'):
                    continue
                line_list = line.strip().split()
                if ',' in line_list[3] or ',' in line_list[4]:
                    continue
                # F1群体过滤|F1 population filter
                if (line_list[p1].split(':')[0] == line_list[p2].split(':')[0]
                        or line_list[p1].startswith('./.') or line_list[p2].startswith('0/1')
                        or line_list[p2].startswith('0|1') or line_list[b1].startswith('./.')
                        or line_list[b2].startswith('./.')):
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
    def pre_binom_test(rep, p):
        """预计算二项分布阈值|Pre-compute binomial distribution thresholds"""
        out_dict = {}
        for n in range(3, 4000, 1):
            num = np.percentile(np.random.binomial(n, 0.5, rep), p)
            out_dict[n] = num
        return out_dict

    @staticmethod
    def _get_dep(geno):
        """提取覆盖度和基因型|Extract depth and genotype"""
        dep = int(geno.split('|')[0]) + int(geno.split('|')[1])
        geno_list = [int(geno.split('|')[0]), int(geno.split('|')[1])]
        return dep, geno_list

    @staticmethod
    def _binom_index(b1_dep, b2_dep, binom_dict):
        """计算binom index|Calculate binom index"""
        num1 = binom_dict[int(b1_dep)]
        num2 = binom_dict[int(b2_dep)]
        index1 = (num1 - (b1_dep - num1)) / b1_dep
        index2 = (num2 - (b2_dep - num2)) / b2_dep
        return index1 + index2

    @staticmethod
    def _binom_index2(b1_dep, b2_dep, binom_dict):
        """计算binom index2|Calculate binom index2"""
        num1 = binom_dict[int(b1_dep)]
        num2 = binom_dict[int(b2_dep)]
        index1 = num1 / b1_dep - (b2_dep - num2) / b2_dep
        index2 = (b1_dep - num1) / b1_dep - num2 / b2_dep
        return index1 ** 4 + index2 ** 4

    def filter_data(self, key, chr_list, dep_lim, dep_high, dep_lim_pool,
                    dep_high_pool, binom_dict):
        """过滤SNP位点并计算DHHP值|Filter SNP sites and calculate DHHP values

        Args:
            key: 染色体名|Chromosome name
            chr_list: 该染色体的SNP列表|SNP list for this chromosome
            dep_lim: 亲本最低覆盖度|Parent min depth
            dep_high: 亲本最高覆盖度|Parent max depth
            dep_lim_pool: 混池最低覆盖度|Pool min depth
            dep_high_pool: 混池最高覆盖度|Pool max depth
            binom_dict: 二项分布阈值字典|Binomial threshold dict

        Returns:
            [染色体名, 过滤后的SNP列表]|[chromosome, filtered SNP list]
        """
        out_snp_list = []

        for vcf_list in chr_list:
            p1_dep, p1_geno = self._get_dep(vcf_list[3])
            p2_dep, p2_geno = self._get_dep(vcf_list[4])
            b1_dep, b1_geno = self._get_dep(vcf_list[5])
            b2_dep, b2_geno = self._get_dep(vcf_list[6])

            # 过滤覆盖度|Filter by depth
            if (p1_dep < dep_lim or p1_dep > dep_high or p2_dep > dep_high
                    or b1_dep < dep_lim_pool or b1_dep > dep_high_pool
                    or b2_dep < dep_lim_pool or b2_dep > dep_high_pool):
                continue

            # 过滤P1基因型|Filter by P1 genotype
            if p1_geno[0] > binom_dict[p1_dep] or p1_geno[1] > binom_dict[p1_dep]:
                continue

            # P2基因型分类|P2 genotype classification
            if p2_dep < 2:
                flag = 'PAV'
                b1_index = (b1_geno[0] - b1_geno[1]) / b1_dep
                b2_index = (b2_geno[0] - b2_geno[1]) / b2_dep
                total_index = abs(b1_index - b2_index)
            elif (p2_geno[0] == 0 or p2_geno[1] == 0
                  or p2_geno[1] / p2_geno[0] < 0.1
                  or p2_geno[1] / p2_geno[0] > 10):
                if p2_geno[1] > p2_geno[0]:
                    p1_num, p2_num = 0, 1
                else:
                    p1_num, p2_num = 1, 0

                if len(vcf_list[1]) > 1 or len(vcf_list[2]) > 1:
                    flag = 'indel'
                else:
                    flag = 'snp'

                b1_index = ((b1_geno[p2_num] - (b1_dep / 2)) / (b1_dep / 2)
                            - (b2_geno[p2_num] - (b2_dep / 2)) / (b2_dep / 2))
                b2_index = (b1_geno[p1_num] / (b1_dep / 2)
                            - b2_geno[p1_num] / (b2_dep / 2))

                if abs(b1_index) > 1.3 or abs(b2_index) > 1.3:
                    continue

                total_index = abs(b1_index ** 4 + b2_index ** 4)
            else:
                continue

            virtual_index = self._binom_index2(b1_dep / 2, b2_dep / 2, binom_dict)
            out_snp_list.append(vcf_list + [virtual_index, total_index, flag])

        return [key, out_snp_list]

    @staticmethod
    def cal_dis_func(key, win_size, loop_value):
        """滑窗平滑计算|Sliding window smoothing calculation"""
        step = win_size / 10
        win_list = [0, win_size]
        total_list = []
        start_pos = 0

        while win_list[1] < loop_value[-1][0]:
            this_list = [[], []]
            for loop_list in loop_value[start_pos:]:
                i_list = [int(loop_list[0]), float(loop_list[-3]), float(loop_list[-2])]
                if win_list[0] <= i_list[0] < win_list[1]:
                    this_list[0].append(i_list[1])
                    this_list[1].append(i_list[2])
                elif i_list[0] >= win_list[1]:
                    try:
                        total_list.append([
                            (win_list[0] + win_list[1]) / 2,
                            sum(this_list[0]) / len(this_list[0]),
                            sum(this_list[1]) / len(this_list[0]),
                            len(this_list[0])
                        ])
                    except ZeroDivisionError:
                        total_list.append([(win_list[0] + win_list[1]) / 2, 0, 0, 0])
                    win_list[0] += step
                    win_list[1] += step
                    break
                elif i_list[0] < win_list[0]:
                    start_pos += 1

        return [key, total_list]

    def run(self):
        """运行F1群体OcBSA分析流程|Run F1 population OcBSA analysis pipeline"""
        cfg = self.config

        # 检测VCF样本列偏移|Detect VCF sample column offset
        open_fn = gzip.open if cfg.input_vcf.endswith('.gz') else open
        sample_offset, sample_names = 0, []
        with open_fn(cfg.input_vcf, 'rt') as f:
            for line in f:
                if line.startswith('#CHROM'):
                    header_cols = line.strip().split('\t')
                    sample_offset = header_cols.index('FORMAT') + 1
                    sample_names = header_cols[sample_offset:]
                    break
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
        binom_dict = self.pre_binom_test(1000, cfg.pvalue)

        # 读取VCF数据|Read VCF data
        self.logger.info(f"读取VCF文件|Reading VCF file: {cfg.input_vcf}")
        chr_geno_dict = self.read_vcf_data(cfg.input_vcf, p1, p2, b1, b2)
        self.logger.info(f"读取到|Read {len(chr_geno_dict)} 条染色体数据")

        # 过滤数据并计算DHHP值|Filter data and calculate DHHP values
        self.logger.info("开始计算DHHP值|Starting DHHP value calculation")
        chr_geno_dict_filter = {}
        pool = multiprocessing.Pool(len(chr_geno_dict))
        results = []
        for key, value in chr_geno_dict.items():
            self.logger.info(f"处理染色体|Processing chromosome: {key}")
            results.append(pool.apply_async(
                self.filter_data,
                (key, value, cfg.parent_min_dep, cfg.parent_max_dep,
                 cfg.pool_min_dep, cfg.pool_max_dep, binom_dict)
            ))
        pool.close()
        pool.join()

        for result in results:
            chr_list = result.get()
            if len(chr_list[1]) < 1000:
                self.logger.info(f"跳过标记数不足的染色体|Skipping chromosome with insufficient markers: {chr_list[0]} ({len(chr_list[1])})")
                continue
            chr_geno_dict_filter[chr_list[0]] = chr_list[1]

        # 输出OcValue文件|Output OcValue file
        ocvalue_file = str(cfg.output_path / "ocbsa.OcValue")
        with open(ocvalue_file, 'w') as out_file:
            out_file.write('#chr.\tpos\tref\talt\tP1\tP2\tPool1\tPool2\tthresholds\tOcValue\tmarker_type\n')
            for key, value in chr_geno_dict_filter.items():
                for i_list in value:
                    out_file.write(key + '\t')
                    for i in i_list:
                        out_file.write(str(i) + '\t')
                    out_file.write('\n')
        self.logger.info(f"OcValue结果已保存|OcValue results saved to: {ocvalue_file}")

        # 滑窗平滑|Sliding window smoothing
        self.logger.info(f"开始滑窗平滑|Starting sliding window smoothing (window={cfg.window_size})")
        pool = multiprocessing.Pool(len(chr_geno_dict_filter))
        result = []
        for key, value in chr_geno_dict_filter.items():
            self.logger.info(f"平滑染色体|Smoothing chromosome: {key}")
            result.append(pool.apply_async(self.cal_dis_func, (key, cfg.window_size, value)))
        pool.close()
        pool.join()

        smoothed_file = str(cfg.output_path / "ocbsa.smoothed")
        with open(smoothed_file, 'w') as file3:
            file3.write('#chr.\tpos\tthresholds\tsmoothed_OcValue\tmarker_number\n')
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
        summary_file = str(cfg.output_path / "ocbsa.summary.tsv")
        chr_sorted = sorted(chr_geno_dict_filter.keys(), key=lambda x: int(re.findall(r'\d+', x)[0]))
        chr_len_dict = {}
        for ch in chr_sorted:
            chr_len_dict[ch] = chr_geno_dict_filter[ch][-1][0]

        with open(summary_file, 'w') as sf:
            sf.write('Chr\tPos\tChrLen\tOcValue\tThreshold\tMarkerNum\n')
            for i_list in result:
                i_list = i_list.get()
                chr_name = i_list[0]
                for row in i_list[1][1:]:
                    sf.write(f'{chr_name}\t{int(row[0])}\t{chr_len_dict[chr_name]}'
                             f'\t{row[2]}\t{row[1]}\t{row[3]}\n')
        self.logger.info(f"R绘图汇总数据已保存|Summary data for R plotting saved to: {summary_file}")

        # 自动绘图|Auto plotting
        self.logger.info("开始绘图|Starting figure plotting")
        from .config import BsaFigConfig
        from .fig import BsaFigPlotter

        fig_file = str(cfg.output_path / "ocbsa.OcValue.png")
        fig_config = BsaFigConfig(
            input_file=smoothed_file,
            output_file=fig_file,
            plot_type='ocvalue',
        )
        fig_config.validate()
        plotter = BsaFigPlotter(fig_config, self.logger)
        plotter.run()

        self.logger.info("F1群体OcBSA分析完成|F1 population OcBSA analysis completed")
        return True
