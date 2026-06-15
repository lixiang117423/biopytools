import json
import math
import os
from multiprocessing import Pool, cpu_count
from functools import partial

import numpy as np
import pandas as pd
import statsmodels.api as sm
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from keras.models import load_model
from get_path import resource_path


def evaluate_frac(data, position):
    """对数据进行 5 k-fold CV 得到MSE最小的frac"""
    import rpy2.robjects as R
    print("get window ratio")
    data = np.array(data)
    position = np.array(position) / 1e6

    r_script = """
loess_as <- 
function(x, y, degree=1, criterion=c("aicc", "gcv"), 
	family = c("gaussian", "symmetric"), user.span=NULL, plot=FALSE, ...)
{
	criterion <- match.arg(criterion)
	family <- match.arg(family)
	x <- as.matrix(x)

	data.bind <- data.frame(x=x, y=y)
	if (ncol(x) == 1) {
		names(data.bind) <- c("x", "y")
	} else { names(data.bind) <- c("x1", "x2", "y") }

	opt.span <- function(model, criterion=c("aicc", "gcv"), span.range=c(.05, .95)){	
		as.crit <- function (x) {
			span <- x$pars$span
			traceL <- x$trace.hat
			sigma2 <- sum(x$residuals^2 ) / (x$n-1)
			aicc <- log(sigma2) + 1 + 2* (2*(traceL+1)) / (x$n-traceL-2)
			gcv <- x$n*sigma2 / (x$n-traceL)^2
			result <- list(span=span, aicc=aicc, gcv=gcv)
			return(result)
			}
		criterion <- match.arg(criterion)
		fn <- function(span) {
			mod <- update(model, span=span)
			as.crit(mod)[[criterion]]
		}
		result <- optimize(fn, span.range)
		return(list(span=result$minimum, criterion=result$objective))
		}

	if (ncol(x)==1) {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x, degree=degree, family = family, data=data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
		} else {
			span1 <- user.span
		}		
		fit <- loess(y ~ x, degree=degree, span=span1, family = family, data=data.bind, ...)
	} else {
		if (is.null(user.span)) {
			fit0 <- loess(y ~ x1 + x2, degree=degree,family = family, data.bind, ...)
			span1 <- opt.span(fit0, criterion=criterion)$span
		} else {
			span1 <- user.span
		}		
	}
	all_var = ls()
	rm(list = all_var[which(all_var!="span1")])
    gc()
	return(span1)
}
    """

    R.r(r_script)

    best_fraction = list(R.r['loess_as'](R.FloatVector(data), R.FloatVector(position)))[0]

    # rm_r_script = """
    # rm(list = ls())
    # gc()
    # """
    # R.r(rm_r_script)

    print("auto window ratio:", best_fraction)
    return round(best_fraction, 3)


def tri_cube_kernel_regression(data, window_size, step):
    # 计算距离权重的匿名函数
    weight_cal = lambda x: (1 - x ** 3) ** 3

    try:
        assert window_size % 2 == 1
        # 窗口大小为奇数
        center = int((window_size + 1) / 2)
        left = int(center - ((window_size - 1) / 2))
        right = int(center + ((window_size - 1) / 2))
    except:
        # 窗口大小为偶数
        center = int(window_size / 2)
        left = int(center - (window_size / 2) + 1)
        right = int(center + window_size / 2)
    left_offset = center - left
    right_offset = right - center
    # 数据扩充
    padding_num = 0  # data[0]
    head_padding = np.array([padding_num] * left_offset)
    padding_num = 0  # data[-1]
    tail_padding = np.array([padding_num] * (right_offset + step - (len(data) - 1) % step))
    padded_data = np.concatenate((head_padding, data, tail_padding))
    # 构造距离向量
    distance_vector = np.arange(len(padded_data))
    result = []
    for c_index in range(center - 1, len(padded_data) - right_offset - 1, step):
        window_data = padded_data[c_index - left_offset: c_index + right_offset + 1]
        window_distance = np.abs(
            distance_vector[c_index - left_offset: c_index + right_offset + 1] - distance_vector[c_index])
        distance_weight = weight_cal(window_distance / np.max(window_distance))
        k_weight = distance_weight / np.sum(distance_weight)
        result.append(np.sum(np.dot(k_weight, window_data)))
        QApplication.processEvents()
    return np.nan_to_num(result)


def smooth_function(data, function_name, window_size, step):
    try:
        assert function_name == "Tri-kernel-smooth"
        return tri_cube_kernel_regression(data, int(window_size * len(data)), step)
    except:
        pass
    try:
        assert function_name == "LOWESS"
        loess = sm.nonparametric.lowess
        z = loess(data, [i for i in range(len(data))], frac=window_size, it=0, delta=0)
        return z[:, 1]
    except:
        pass
    try:
        assert function_name == "Moving Average"
        return np.convolve(data, np.ones((window_size,)) / window_size, mode="same")
    except:
        pass


def oneD_peaks_finder(smooth_value, position, height):
    """1-D peaks finder"""
    # 获得所有大于阈值的点
    sets = np.where(smooth_value >= height)[0]
    if np.max(smooth_value) < height:
        return ["-"], ["-"], ["-"], ["-"]
    break_sets = []
    temp_set = []
    start = sets[0]
    for i in sets:
        try:
            assert i - start <= 1
            temp_set.append(i)
        except:
            break_sets.append(temp_set)
            temp_set = []
        start = i
        if i == sets[-1]:
            break_sets.append(temp_set)
    peaks, lefts, rights, values = [], [], [], []
    # 将连续的点划分为一个类
    for i in break_sets:
        maxi = np.argmax(smooth_value[i])
        left = position[i[0]]
        right = position[i[-1]]

        peaks.append(position[i[maxi]])
        lefts.append(left)
        rights.append(right)
        values.append(round(smooth_value[i[maxi]], 5))

    return lefts, peaks, rights, values


def peaks_finder(smooth_data, position, height, chrome_set, save_path):
    chromes, lefts, peaks, rights, values = [], [], [], [], []
    for chr_smooth_data, pos, chrome in zip(smooth_data, position, chrome_set):
        left, peak, right, value = oneD_peaks_finder(np.array(chr_smooth_data), pos, height)
        lefts = np.concatenate((lefts, left), axis=0)
        peaks = np.concatenate((peaks, peak), axis=0)
        rights = np.concatenate((rights, right), axis=0)
        values = np.concatenate((values, value), axis=0)
        chromes = np.concatenate((chromes, [chrome] * len(left)), axis=0)
    QTLs = [i + 1 for i in range(len(chromes))]
    # 按峰值降序排序
    index = np.argsort(values)[::-1]
    chromes = np.array(chromes)[index]
    lefts = np.array(lefts)[index]
    peaks = np.array(peaks)[index]
    rights = np.array(rights)[index]
    values = np.array(values)[index]
    QTLs = list(map(lambda x: round(x), QTLs))
    # chromes = list(map(lambda x: str(int(int(x))), chromes))
    dic = {"QTL": QTLs,
           "Chr": chromes,
           "Left": lefts,
           "Peak": peaks,
           "Right": rights,
           "Value": values}
    df = pd.DataFrame(dic)
    df.to_csv(save_path, index=False)
    return dic


def confidence_interval_plot(dict, data, position, chrome_set, save_path, func_name):
    import json

    this_dic = {}
    for chrome, peak in zip(dict["Chr"], dict["Peak"]):
        if peak == "-":
            break
        peak = float(peak) * 1e6
        index = list(chrome_set).index(chrome)
        chr_data = data[index]
        chr_pos = position[index] * 1e6
        # ed4_data = 4 * np.power(np.array(chr_data).T[-1] - np.array(chr_data).T[0], 4)
        ed4_data = np.array(chr_data)
        left_ptr = right_ptr = peak
        left_pos = []
        right_pos = []
        left_ed4 = []
        right_ed4 = []
        while 1:  # left_ptr > peak - 1e6:
            left_ptr -= 1e4
            if left_ptr < chr_pos[0]:
                break
            valid_index = list(set(np.where(left_ptr < chr_pos)[0]) & set(np.where(chr_pos < left_ptr + 1e4)[0]))
            if len(valid_index) == 0:
                continue
            valid_ed4 = ed4_data[valid_index]
            valid_pos = chr_pos[valid_index]
            left_ed4.append(np.mean(valid_ed4))
            left_pos.append(np.mean(valid_pos))
        while 1:  # right_ptr < peak + 1e6:
            right_ptr += 1e4
            if right_ptr > chr_pos[-1]:
                break
            valid_index = list(set(np.where(right_ptr - 1e4 < chr_pos)[0]) & set(np.where(chr_pos < right_ptr)[0]))
            if len(valid_index) == 0:
                continue
            valid_ed4 = ed4_data[valid_index]
            valid_pos = chr_pos[valid_index]
            right_ed4.append(np.mean(valid_ed4))
            right_pos.append(np.mean(valid_pos))

        # rank = np.argsort(left_pos)
        # left_pos = np.array(left_pos)[rank]
        # left_ed4 = np.array(left_ed4)[rank]
        if chrome in this_dic:      # 同一条染色体有多个位点的情况
            while 1:
                chrome = chrome + "*"
                if chrome not in this_dic:
                    break
            this_dic.update({
                chrome: {
                    "data": left_ed4[::-1] + right_ed4,
                    "pos": left_pos[::-1] + right_pos,
                    "peak": peak
                }
            })
        else:
            this_dic.update({
                chrome: {
                    "data": left_ed4[::-1] + right_ed4,
                    "pos": left_pos[::-1] + right_pos,
                    "peak": peak
                }
            })
    info_json = json.dumps(this_dic)
    with open(os.path.join(save_path, func_name + "_confidence_interval_data.json"), "w") as f:
        f.write(info_json)
    # np.save(os.path.join(save_path, "confidence_interval_data.npy"), this_dic)


def cal_threshold(data):
    from decimal import Decimal
    # data = data[np.where(data > 0.05)[0]]
    std = np.std(data)
    med = np.median(data)
    threshold = med + 3 * std

    if np.all(data < round(threshold, 4)):    # 不满足
        threshold = np.percentile(data, 90)
    return Decimal(threshold).quantize(Decimal("0.0001"), rounding="ROUND_HALF_UP")


def generator(x_data, batch_size=32, window=64, step=64):  # x_data,y_data are array types
    data = np.copy(x_data)

    x, y = [], []
    start, end = 0, window
    for i in range(batch_size):
        x.append(data[start: end])
        start += step
        end = start + window
    return np.array(x).astype("float64")


def cal_function(x, y):
    def fitting_function(x, k, b):  # 线性拟合
        return k * x + b

    from scipy import optimize
    k, b = optimize.curve_fit(fitting_function, x, y)[0]
    return abs(k)


def G_statistic(data):
    Gs = []
    for nn in data:
        G = 0
        n1, n2 = nn[0], 1 - nn[0]
        n3, n4 = nn[-1], 1 - nn[-1]
        tabel = np.array([[n1, n2],  # np.sum(tabel, axis=0) = [n1+n3, n2+n4]
                          [n3, n4]])  # np.sum(tabel, axis=1) = [n1+n2, n3+n4]
        for i in range(2):
            for j in range(2):
                ni = tabel[i][j]
                n_hat = np.sum(tabel, axis=i)[j] * np.sum(tabel, axis=j)[i] / np.sum(tabel)
                G += ni * np.log(ni / n_hat)
        Gs.append(2 * G)
    Gs = np.nan_to_num(Gs)
    return Gs


def LOD_statistic(ref_data, mut_data):
    ref_data = np.array(ref_data)
    mut_data = np.array(mut_data)
    lods = []
    for index in range(ref_data.shape[0]):
        n_AL, n_aL = ref_data[index][0], mut_data[index][0]
        n_AH, n_aH = ref_data[index][1], mut_data[index][1]
        n_L = n_aL + n_AL
        n_H = n_aH + n_AH
        p_L, p_H = n_AL / (n_AL + n_aL), n_AH / (n_aH + n_AH)
        # 化简后的计算公式
        lod = n_AL * np.log10(p_L) + n_aL * np.log10(1 - p_L) + \
              n_AH * np.log10(p_H) + n_aH * np.log10(1 - p_H) - \
              (n_L + n_H) * np.log10(1 / 2)
        lods.append(lod)
    lods = np.nan_to_num(lods)
    return lods


def ridit(ref_data, mut_data):
    def cal_ridit(a):
        a[a == 0] = 0.0001
        num = len(a)
        sum_total = 0
        sum_ref = 0
        sum_snp = 0
        sum = {}
        cumsum = {0: 0}
        r = {}
        cumfr = 0
        cumfrr = 0
        r_ref = 0
        r_snp = 0
        p = 0
        for i in range(num):
            sum_total += a[i]
            if i % 2 == 0:
                sum_ref += a[i]
            else:
                sum_snp += a[i]
        i = 0
        while i < num - 1:
            m = i / 2
            sum[m] = a[i] + a[i + 1]
            if m + 1 < num / 2:
                cumsum[m + 1] = cumsum[m] + sum[m]
            r[m] = (cumsum[m] + sum[m] / 2) / sum_total
            cumfr += sum[m] * r[m]
            cumfrr += sum[m] * r[m] * r[m]
            r_ref += (a[i] * r[m] / sum_ref)
            r_snp += (a[i + 1] * r[m] / sum_snp)
            i += 2
        var = (cumfrr - cumfr * cumfr / sum_total) / (sum_total - 1)
        try:
            n = abs(r_ref - r_snp) / math.sqrt(var * sum_total / (sum_ref * sum_snp))
        except:
            n = 101
        if n > 100 or n < -100:
            p = 0
        elif n <= 100 and n >= 1.9:
            for i in range(18, 0, -1):
                p = i / (abs(n) + p)
            p = math.exp(-.5 * abs(n) * abs(n)) / math.sqrt(2 * math.pi) / (abs(n) + p)
        elif n >= -100 and n <= -1.9:
            for i in range(18, 0, -1):
                p = i / (abs(n) + p)
            p = math.exp(-.5 * abs(n) * abs(n)) / math.sqrt(2 * math.pi) / (abs(n) + p)
            p = 1 - p
        elif n < 1.9 and n >= 0:
            p = (1 + abs(n) * (.049867347 + abs(n) * (.0211410061 + abs(n) * (.0032776263 + abs(n) * (
                    .0000380036 + abs(n) * (.0000488906 + abs(n) * .000005383)))))) ** -16 / 2
        elif n > -1.9 and n < 0:
            p = (1 + abs(n) * (.049867347 + abs(n) * (.0211410061 + abs(n) * (.0032776263 + abs(n) * (
                    .0000380036 + abs(n) * (.0000488906 + abs(n) * .000005383)))))) ** -16 / 2
            p = 1 - p
        p = float(p)
        if p == 0:
            p += 0.01
        result = -(math.log(abs(float(p))))
        return result

    read_data = np.zeros(shape=(ref_data.shape[0] * 2, ref_data.shape[1]))
    for i in range(ref_data.shape[0]):
        read_data[2 * i] = ref_data[i]
        read_data[2 * i + 1] = mut_data[i]
    read_data = read_data.T
    res = []
    for i in read_data:
        res.append(cal_ridit(i))
    return res


def get_data(func_name, data=None, ref_data=None, mut_data=None, num_pools=None):
    print("get data")
    if func_name == "DL":
        model = load_model(os.path.join(resource_path("Models"),
                                        "row_finetune" + str(num_pools) + "pool.h5"))
        a = len(data) % 64
        pad_n = [0] * num_pools
        pad_l = np.array([pad_n for _ in range(64 - a)])
        chrome_data = np.concatenate((data, pad_l), axis=0)
        batch_size = len(chrome_data) // 64
        x = generator(x_data=chrome_data, batch_size=batch_size)
        pre = model.predict(x)
        dl_data = pre.reshape(-1)[:a - 64]
        return dl_data
    elif func_name == "K":
        x = [i * (1 / (num_pools - 1)) for i in range(num_pools)]
        k_data = [cal_function(x, i) for i in data]
        return np.array(k_data)
    elif func_name == "ED4":
        ed4_data = 4 * np.power(np.array(data).T[-1] - np.array(data).T[0], 4)
        return ed4_data
    elif func_name == "SNP":
        snp_data = np.abs(np.array(data).T[-1] - np.array(data).T[0])
        return snp_data
    elif func_name == "SmoothG":
        g_data = G_statistic(data)
        return g_data
    elif func_name == "SmoothLOD":
        lod_data = LOD_statistic(ref_data, mut_data)
        return lod_data
    elif func_name == "Ridit":
        ridit_data = ridit(np.array(ref_data).T, np.array(mut_data).T)
        return ridit_data


def mean_data(data, pos):
    """每10kb计算平均"""
    mean_data = []
    mean_pos = []
    ptr = pos[0] + 1e4
    mean_chr_data = []
    mean_chr_pos = []
    while 1:
        chr_data = np.array(data)
        chr_pos = np.array(pos)
        valid_index = list(set(np.where(ptr - 1e4 < chr_pos)[0]) & set(np.where(chr_pos < ptr)[0]))
        if len(valid_index) != 0:
            mean_chr_data.append(np.mean(chr_data[valid_index]))
            mean_chr_pos.append(np.mean(chr_pos[valid_index]))

        ptr = ptr + 1e4
        if ptr > pos[-1]:
            break
    print(np.array(mean_chr_data).sum())
    return mean_chr_data, mean_chr_pos


class Statistic(object):

    def __init__(self,
                 func_name,
                 data_path,
                 ref_data_path,
                 mut_data_path,
                 pos_path,
                 chrome_set,
                 read_number,
                 smooth_func,
                 smooth_window_size,
                 threshold,
                 save_path):
        super(Statistic, self).__init__()
        self.func_name = func_name
        self.data_path = data_path
        self.ref_data_path = ref_data_path
        self.mut_data_path = mut_data_path
        self.pos_path = pos_path
        self.chrome_set = chrome_set
        self.read_number = read_number
        self.smooth_func = smooth_func
        self.smooth_window_size = smooth_window_size
        self.threshold = threshold
        self.default_threshold = threshold
        self.save_path = save_path

        self.data = np.load(self.data_path, allow_pickle=True)
        self.ref_data = np.load(self.ref_data_path, allow_pickle=True)
        self.mut_data = np.load(self.mut_data_path, allow_pickle=True)
        self.position = np.load(self.pos_path, allow_pickle=True)

        self.num_pools = np.array(self.data[0]).shape[1]

        self.auto_win = True if self.smooth_window_size == 0 else False

    def run(self):
        fig = plt.figure(figsize=(20, 5), dpi=300)
        colors = ["#b8c6d5", "#9fb2c6"]
        all_data_for_percentile = []
        all_data_for_plot = []
        smooth_data_for_plot = []
        smooth_for_threshold = []
        # 计算数据
        Mpos = []
        for chrome_index in range(self.data.shape[0]):
            print("chrome:",chrome_index+1)
            if len(self.data[chrome_index]) == 0:
                continue
            chr_data = get_data(self.func_name, self.data[chrome_index], self.ref_data[chrome_index],
                                self.mut_data[chrome_index], self.num_pools)
            #############   用10kb平均
            # chr_data, mean_pos = mean_data(chr_data, self.position[chrome_index])
            # print(len(chr_data), len(mean_pos))
            # Mpos.append(mean_pos)
            #############   用10kb平均
            if self.auto_win:
                self.smooth_window_size = evaluate_frac(chr_data, self.position[chrome_index]) # mean_pos
            windows_callback = int(len(chr_data) * self.smooth_window_size)
            y = smooth_function(chr_data, self.smooth_func, self.smooth_window_size, 1)
            all_data_for_percentile = np.concatenate((all_data_for_percentile, chr_data), axis=0)
            all_data_for_plot.append(chr_data)
            smooth_data_for_plot.append(y)
            smooth_for_threshold = np.concatenate((smooth_for_threshold, y), axis=0)
        # 画图
        pos = [np.array(i) / 1e6 for i in self.position]
        self.position = pos
        sort_pos = []
        sort_points = []
        for i in range(len(self.position)):
            index = np.argsort(self.position[i])
            sort_pos.append(self.position[i][index])
            sort_points.append(np.array(all_data_for_plot[i])[index])
        self.position = sort_pos
        all_data_for_plot = sort_points
        #############   用10kb平均
        # Mpos = [np.array(i) / 1e6 for i in Mpos]
        # self.position = Mpos
        #############   用10kb平均
        for index in range(self.data.shape[0]):
            if len(self.data[index]) == 0:
                continue
            ax = fig.add_subplot(1, self.data.shape[0], index + 1)
            # 散点图
            ax.scatter(self.position[index][:len(all_data_for_plot[index])], all_data_for_plot[index],
                       linewidths=0.1, color=colors[index % 2])
            plt.plot(self.position[index][:len(smooth_data_for_plot[index])], smooth_data_for_plot[index],
                     color="#ffa333", linewidth=2.0)
            # 阈值线
            if self.threshold == 0.0:
                self.threshold = cal_threshold(smooth_for_threshold)
            plt.hlines(self.threshold, xmin=0, xmax=max(self.position[index][:len(all_data_for_plot[index])]),
                       colors="#33ccff", linestyles="--", linewidth=2.0)

            plt.ylim(0, max(all_data_for_percentile) * 1.01)
            # 标上染色体号、隐藏横坐标
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
            # ax.set_xlabel(str(index + 1))
            ax.set_title(self.chrome_set[index])
            # ax.set_xticks([])

            frame = plt.gca()
            if index == 0:
                ax.spines["right"].set_visible(False)
            elif index == self.data.shape[0] - 1:
                ax.spines["left"].set_visible(False)
                frame.axes.get_yaxis().set_visible(False)
                ax.set_xlabel("10e6", loc="right")
            else:
                ax.spines["right"].set_visible(False)
                ax.spines["left"].set_visible(False)
                frame.axes.get_yaxis().set_visible(False)

        if self.default_threshold == 0:
            self.threshold = format(self.threshold, '.4f')
        else:
            self.threshold = format(self.threshold, '.2f')

        if self.auto_win:
            self.smooth_window_size = "auto"
        # 保存png用于展示
        plt.savefig(os.path.join(self.save_path,
                                 "{}-{}-{}-{}-{}.png".format(self.read_number, self.func_name, self.smooth_func,
                                                             self.smooth_window_size, self.threshold)),
                    bbox_inches='tight', dpi=fig.dpi, pad_inches=0.0)
        # 保存pdf
        plt.savefig(os.path.join(self.save_path,
                                 "{}-{}-{}-{}-{}.pdf".format(self.read_number, self.func_name, self.smooth_func,
                                                             self.smooth_window_size, self.threshold)),
                    bbox_inches='tight', dpi=fig.dpi, pad_inches=0.0)

        np.save(os.path.join(self.save_path, "all_data_for_percentile_{}.npy".format(self.func_name)),
                all_data_for_percentile)
        np.save(os.path.join(self.save_path, "all_data_for_plot_{}.npy".format(self.func_name)), all_data_for_plot)
        np.save(os.path.join(self.save_path, "smooth_data_for_plot_{}.npy".format(self.func_name)),
                smooth_data_for_plot)
        # 保存为txt
        f = open(os.path.join(self.save_path, "{} values.txt".format(self.func_name)), "w")
        for flag in range(len(self.chrome_set)):
            for txt_pos, txt_value in zip(self.position[flag], all_data_for_plot[flag]):
                txt_line = "{}\t{}\t{}\n".format(self.chrome_set[flag], int(txt_pos * 1e6), txt_value)
                f.writelines(txt_line)
        f.close()

        peaks_save_path = os.path.join(self.save_path,
                                       "{}-{}-{}-{}-{}.csv".format(self.read_number, self.func_name, self.smooth_func,
                                                                   self.smooth_window_size, self.threshold))
        dic = peaks_finder(smooth_data_for_plot, self.position, float(self.threshold), chrome_set=self.chrome_set,
                           save_path=peaks_save_path)
        confidence_interval_plot(dic, smooth_data_for_plot, self.position, self.chrome_set, self.save_path,
                                 self.func_name)

    def run_parallel(self, num_threads=None):
        """
        多线程版本的run方法
        """
        if num_threads is None:
            num_threads = cpu_count()
            num_threads = min(num_threads, self.data.shape[0], 16)

        print("Using {} threads for processing {} chromosomes".format(num_threads, self.data.shape[0]))

        # 准备参数
        args_list = [(i, self.data, self.ref_data, self.mut_data, self.position,
                      self.func_name, self.num_pools, self.auto_win, self.smooth_window_size, self.smooth_func)
                     for i in range(self.data.shape[0])]

        # 使用多进程池处理
        with Pool(processes=num_threads) as pool:
            results = pool.starmap(_process_single_chromosome, args_list)

        # 收集结果
        all_data_for_percentile = []
        all_data_for_plot = []
        smooth_data_for_plot = []
        smooth_for_threshold = []
        position_list = []
        valid_chromosomes = []

        for result in results:
            if result is not None:
                chrome_index, chr_data, smooth_y, pos = result
                all_data_for_percentile = np.concatenate((all_data_for_percentile, chr_data), axis=0)
                all_data_for_plot.append(chr_data)
                smooth_data_for_plot.append(smooth_y)
                smooth_for_threshold = np.concatenate((smooth_for_threshold, smooth_y), axis=0)
                position_list.append(pos)
                valid_chromosomes.append(chrome_index)

        # 更新chrome_set为有效染色体
        self.chrome_set = [self.chrome_set[i] for i in valid_chromosomes]

        # 画图
        fig = plt.figure(figsize=(20, 5), dpi=300)
        colors = ["#b8c6d5", "#9fb2c6"]

        pos = [np.array(i) / 1e6 for i in position_list]
        self.position = pos
        sort_pos = []
        sort_points = []
        for i in range(len(self.position)):
            index = np.argsort(self.position[i])
            sort_pos.append(self.position[i][index])
            sort_points.append(np.array(all_data_for_plot[i])[index])
        self.position = sort_pos
        all_data_for_plot = sort_points

        for idx in range(len(valid_chromosomes)):
            ax = fig.add_subplot(1, len(valid_chromosomes), idx + 1)
            ax.scatter(self.position[idx][:len(all_data_for_plot[idx])], all_data_for_plot[idx],
                       linewidths=0.1, color=colors[idx % 2])
            plt.plot(self.position[idx][:len(smooth_data_for_plot[idx])], smooth_data_for_plot[idx],
                     color="#ffa333", linewidth=2.0)
            if self.threshold == 0.0:
                self.threshold = cal_threshold(smooth_for_threshold)
            plt.hlines(self.threshold, xmin=0, xmax=max(self.position[idx][:len(all_data_for_plot[idx])]),
                       colors="#33ccff", linestyles="--", linewidth=2.0)

            plt.ylim(0, max(all_data_for_percentile) * 1.01)
            plt.subplots_adjust(left=None, bottom=None, right=None, top=None, wspace=0, hspace=0)
            ax.set_title(self.chrome_set[idx])

            frame = plt.gca()
            if idx == 0:
                ax.spines["right"].set_visible(False)
            elif idx == len(valid_chromosomes) - 1:
                ax.spines["left"].set_visible(False)
                frame.axes.get_yaxis().set_visible(False)
                ax.set_xlabel("10e6", loc="right")
            else:
                ax.spines["right"].set_visible(False)
                ax.spines["left"].set_visible(False)
                frame.axes.get_yaxis().set_visible(False)

        if self.default_threshold == 0:
            self.threshold = format(self.threshold, '.4f')
        else:
            self.threshold = format(self.threshold, '.2f')

        if self.auto_win:
            self.smooth_window_size = "auto"

        plt.savefig(os.path.join(self.save_path,
                                 "{}-{}-{}-{}-{}.png".format(self.read_number, self.func_name, self.smooth_func,
                                                             self.smooth_window_size, self.threshold)),
                    bbox_inches='tight', dpi=fig.dpi, pad_inches=0.0)
        plt.savefig(os.path.join(self.save_path,
                                 "{}-{}-{}-{}-{}.pdf".format(self.read_number, self.func_name, self.smooth_func,
                                                             self.smooth_window_size, self.threshold)),
                    bbox_inches='tight', dpi=fig.dpi, pad_inches=0.0)

        np.save(os.path.join(self.save_path, "all_data_for_percentile_{}.npy".format(self.func_name)),
                all_data_for_percentile)
        np.save(os.path.join(self.save_path, "all_data_for_plot_{}.npy".format(self.func_name)), all_data_for_plot)
        np.save(os.path.join(self.save_path, "smooth_data_for_plot_{}.npy".format(self.func_name)),
                smooth_data_for_plot)

        f = open(os.path.join(self.save_path, "{} values.txt".format(self.func_name)), "w")
        for flag in range(len(self.chrome_set)):
            for txt_pos, txt_value in zip(self.position[flag], all_data_for_plot[flag]):
                txt_line = "{}\t{}\t{}\n".format(self.chrome_set[flag], int(txt_pos * 1e6), txt_value)
                f.writelines(txt_line)
        f.close()

        peaks_save_path = os.path.join(self.save_path,
                                       "{}-{}-{}-{}-{}.csv".format(self.read_number, self.func_name, self.smooth_func,
                                                                   self.smooth_window_size, self.threshold))
        dic = peaks_finder(smooth_data_for_plot, self.position, float(self.threshold), chrome_set=self.chrome_set,
                           save_path=peaks_save_path)
        confidence_interval_plot(dic, smooth_data_for_plot, self.position, self.chrome_set, self.save_path,
                                 self.func_name)


# ========== 多线程支持函数 ==========

def _process_single_chromosome(chrome_index, data, ref_data, mut_data, position,
                               func_name, num_pools, auto_win, smooth_window_size, smooth_func):
    """处理单条染色体的辅助函数（用于多线程）"""
    if len(data[chrome_index]) == 0:
        return None

    print("Processing chromosome {}".format(chrome_index + 1))

    # 计算数据
    chr_data = get_data(func_name, data[chrome_index], ref_data[chrome_index],
                        mut_data[chrome_index], num_pools)

    # 自动窗口大小
    if auto_win:
        smooth_window_size = evaluate_frac(chr_data, position[chrome_index])

    # 平滑处理
    y = smooth_function(chr_data, smooth_func, smooth_window_size, 1)

    return (chrome_index, chr_data, y, position[chrome_index])
