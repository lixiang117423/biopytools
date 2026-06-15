"""
OcBSA - BSA结果可视化绘图|OcBSA - BSA Result Visualization
"""

import re


class BsaFigPlotter:
    """BSA结果绘图器|BSA Result Plotter"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def plot_ocvalue(self):
        """绘制OcValue图|Plot OcValue figure"""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np

        cfg = self.config

        with open(cfg.input_file, 'r') as file1:
            if cfg.position:
                fig = plt.figure(figsize=(10, 5))
            else:
                fig = plt.figure(figsize=(30, 5))
            matplotlib.rcParams['pdf.fonttype'] = 42
            gs = matplotlib.gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0:1, 0:1], facecolor='white')

            snp_num_list = []
            chr_dict = {}

            if cfg.position:
                chr_name = cfg.position.split(',')[0]
                filter_start = int(cfg.position.split(',')[1])
                filter_end = int(cfg.position.split(',')[2])

            for line in file1:
                if line.startswith('#'):
                    continue
                line_list = line.strip().split('\t')
                snp_num_list.append(float(line_list[4]))

                if cfg.position:
                    if line_list[0] != chr_name:
                        continue
                    elif int(line_list[1]) < filter_start:
                        continue
                    elif int(line_list[1]) > filter_end:
                        continue

                if line_list[0] in chr_dict:
                    chr_dict[line_list[0]].append(
                        [float(line_list[2]), float(line_list[3]),
                         float(line_list[4]), int(line_list[1])])
                else:
                    chr_dict[line_list[0]] = [
                        [float(line_list[2]), float(line_list[3]),
                         float(line_list[4]), int(line_list[1])]]

            chr_list = sorted(chr_dict.keys(), key=lambda x: int(re.findall(r'\d+', x)[0]))
            total_list_x = []
            total_list_y = []
            total_list_y1 = []
            start = 0
            chr_pos_list = []
            chr_name_list = []
            chr_pos2_list = []
            color_list = []
            big_num = np.percentile(snp_num_list, 90)

            for chrom in chr_list:
                chr_name_list.append(chrom)
                chr_pos2_list.append(start + chr_dict[chrom][-1][-1] / 2)
                chr_pos_list.append([chrom, start + chr_dict[chrom][-1][-1] / 2, start])
                for i in chr_dict[chrom]:
                    total_list_x.append(i[3] + start)
                    total_list_y.append(i[0])
                    total_list_y1.append(i[1])
                    if i[2] > big_num:
                        color_list.append(big_num)
                    else:
                        color_list.append(i[2])
                start += chr_dict[chrom][-1][-1]

            for chrom in chr_pos_list:
                plt.axvline(chrom[2], ls="--", c="black", linewidth='0.3')

            axx = ax1.scatter(total_list_x, total_list_y1, s=8, c=color_list, cmap=cfg.color)

            if cfg.position:
                plt.xlim([total_list_x[0], start])
                step = round(int((filter_end - filter_start) / 20), -6)
                if step == 0:
                    step = round(int((filter_end - filter_start) / 10), -6)
                step = int(step)
                x_ticks = [i for i in range(filter_start, filter_end, step)]
                plt.xticks(x_ticks)
            else:
                sorted_y_list = sorted(total_list_y1)
                top95 = sorted_y_list[int(len(sorted_y_list) * 0.95)]
                top99 = sorted_y_list[int(len(sorted_y_list) * 0.99)]
                top999 = sorted_y_list[int(len(sorted_y_list) * 0.999)]
                ax1.axhline(y=top95, color='r', linestyle='--')
                ax1.axhline(y=top99, color='g', linestyle='--')
                ax1.axhline(y=top999, color='b', linestyle='--')
                plt.xlim([0, start])
                plt.xticks(chr_pos2_list, chr_name_list, fontsize=15)

            plt.yticks(fontsize=15)
            position = fig.add_axes([0.55, 0.01, 0.3, 0.02])
            fig.colorbar(axx, ax=ax1, cax=position, orientation='horizontal')
            plt.savefig(cfg.output_file, bbox_inches="tight", pad_inches=0.6, dpi=600)
            self.logger.info(f"OcValue图已保存|OcValue figure saved to: {cfg.output_file}")
            plt.close()

    def plot_ed(self):
        """绘制ED/SNP-index图|Plot ED/SNP-index figure"""
        import matplotlib
        matplotlib.use('Agg')
        import matplotlib.pyplot as plt
        import numpy as np

        cfg = self.config

        with open(cfg.input_file, 'r') as file1:
            fig = plt.figure(figsize=(30, 5))
            matplotlib.rcParams['pdf.fonttype'] = 42
            gs = matplotlib.gridspec.GridSpec(1, 1)
            ax1 = fig.add_subplot(gs[0:1, 0:1], facecolor='white')

            snp_num_list = []
            chr_dict = {}

            if cfg.position:
                chr_name = cfg.position.split(',')[0]
                filter_start = int(cfg.position.split(',')[1])
                filter_end = int(cfg.position.split(',')[2])

            for line in file1:
                if line.startswith('#'):
                    continue
                line_list = line.strip().split('\t')
                snp_num_list.append(float(line_list[-1]))

                if cfg.position:
                    if line_list[0] != chr_name:
                        continue
                    elif int(line_list[1]) < filter_start:
                        continue
                    elif int(line_list[1]) > filter_end:
                        continue

                if line_list[0] in chr_dict:
                    chr_dict[line_list[0]].append(
                        [float(line_list[2]), float(line_list[3]),
                         float(line_list[4]), float(line_list[5]),
                         int(line_list[1]), int(line_list[-1])])
                else:
                    chr_dict[line_list[0]] = [
                        [float(line_list[2]), float(line_list[3]),
                         float(line_list[4]), float(line_list[5]),
                         int(line_list[1]), int(line_list[-1])]]

            chr_list = sorted(chr_dict.keys(), key=lambda x: int(re.findall(r'\d+', x)[0]))
            total_list_x = []
            total_list_y = []
            total_list_y1 = []
            total_list_y2 = []
            total_list_y3 = []
            start = 0
            chr_pos_list = []
            chr_name_list = []
            chr_pos2_list = []
            color_list = []
            big_num = np.percentile(snp_num_list, 90)

            for chrom in chr_list:
                chr_name_list.append(chrom)
                chr_pos2_list.append(start + chr_dict[chrom][-1][-2] / 2)
                chr_pos_list.append([chrom, start + chr_dict[chrom][-1][-2] / 2, start])
                for i in chr_dict[chrom]:
                    total_list_x.append(i[-2] + start)
                    total_list_y.append(i[0])
                    total_list_y1.append(i[1])
                    total_list_y2.append(i[2])
                    total_list_y3.append(i[3])
                    if i[-1] > big_num:
                        color_list.append(big_num)
                    else:
                        color_list.append(i[-1])
                start += chr_dict[chrom][-1][-2]

            color_list[-1] = 0

            for chrom in chr_pos_list:
                plt.axvline(chrom[2], ls="--", c="black", linewidth='0.3')

            sorted_y_list = sorted(total_list_y1)
            top95 = sorted_y_list[int(len(sorted_y_list) * 0.95)]
            top99 = sorted_y_list[int(len(sorted_y_list) * 0.99)]
            top999 = sorted_y_list[int(len(sorted_y_list) * 0.999)]

            snpindex = cfg.plot_type == "snpindex"
            axx = ax1.scatter(total_list_x, total_list_y3, s=2, c=color_list, cmap=cfg.color)

            if snpindex:
                ax1.plot(total_list_x, total_list_y1, color='g', lw=0.3, label='Pool1 snp index')
                ax1.plot(total_list_x, total_list_y2, color='b', lw=0.3, label='Pool2 snp index')

            if cfg.position:
                plt.xlim([total_list_x[0], start])
                step = round(int((filter_end - filter_start) / 20), -6)
                if step == 0:
                    step = round(int((filter_end - filter_start) / 10), -6)
                x_ticks = [i for i in range(filter_start, filter_end, step)]
                plt.xticks(x_ticks)
            else:
                ax1.plot(total_list_x, total_list_y, color='r', lw=0.5, label='thresholds')
                plt.xlim([0, start])
                plt.xticks(chr_pos2_list, chr_name_list, fontsize=15)

            ax1.legend()
            ax1.legend(loc='upper right')
            plt.yticks(fontsize=15)
            position = fig.add_axes([0.55, 0.01, 0.3, 0.02])
            fig.colorbar(axx, ax=ax1, cax=position, orientation='horizontal')
            plt.savefig(cfg.output_file, bbox_inches="tight", pad_inches=0.6, dpi=600)
            self.logger.info(f"图已保存|Figure saved to: {cfg.output_file}")
            plt.close()

    def run(self):
        """运行绘图|Run plotting"""
        self.config.validate()
        self.logger.info(f"输入文件|Input file: {self.config.input_file}")
        self.logger.info(f"图表类型|Plot type: {self.config.plot_type}")

        if self.config.plot_type == "ocvalue":
            self.plot_ocvalue()
        else:
            self.plot_ed()

        self.logger.info("绘图完成|Plotting completed")
        return True
