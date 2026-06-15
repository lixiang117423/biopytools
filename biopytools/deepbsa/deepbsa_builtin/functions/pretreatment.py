import os

import numpy as np
import pandas as pd
from tqdm import tqdm

class Pretreatment(object):
    def __init__(self, step_one, step_two, step_three, file_path, file_name, save_path):
        super(Pretreatment, self).__init__()
        self.read_number = step_one
        self.step_two = step_two
        self.step_three = step_three
        self.file_path = file_path
        self.file_name = file_name
        self.save_path = save_path

    def first_check(self, A_pool, a_pool):
        """
        :param low_pool: 所有池A数据
        :param high_pool: 所有池a数据
        :return: 返回是否通过一次检验，通过返回True；否则反之
        """
        from scipy.stats import chisquare
        A_pool[A_pool == 0] = A_pool[A_pool == 0] + 0.0001
        a_pool[a_pool == 0] = a_pool[a_pool == 0] + 0.0001
        # 筛选条件1：每个池有一个read数的和（A+a）小于read number的都不要a
        try:
            assert np.any((A_pool + a_pool) < self.read_number)
            return False
        except:
            pass
        # 筛选条件2：每个池的卡方检验都大于0.05的要，或者第一个池、最后一个池都小于0.05并且分布在0.5两端
        p_values = []
        A_ratios = []
        for A, a in zip(A_pool, a_pool):
            p_value = chisquare(f_obs=[A, a], f_exp=[(A + a) / 2, (A + a) / 2])[1]
            p_values.append(p_value)
            A_ratios.append(A / (A + a))
        p_values = np.array(p_values)
        A_ratios = np.array(A_ratios)

        if np.all(p_values >= 0.05):
            return True
        elif (p_values[0] < 0.05 and p_values[-1] < 0.05) and (
                (A_ratios[0] > 0.5 and A_ratios[-1] < 0.5) or (A_ratios[0] < 0.5 and A_ratios[-1] > 0.5)):
            return True

    def second_check(self, A, a, freq, pos, is_process):
        """
        :param A: A的read数
        :param a: a的read数
        :param freq: A的频率
        :param pos:  位置
        :return: 返回经过筛选的4组数据的结果
        此筛选条件为：一个数据点每个池的A的频率与其前后的点的对应的每个池的频率差异小于0.1
        """
        result_A, result_a, result_freq, result_pos = [], [], [], []
        for index in range(len(A)):
            signal = False  # 设置是否通过检验的信号
            if index == 0:
                if np.all(abs(freq[index] - freq[index + 1]) < 0.1):
                    signal = True
            elif index == len(A) - 1:
                if np.all(abs(freq[index] - freq[index - 1]) < 0.1):
                    signal = True
            else:
                if np.all(abs(freq[index] - freq[index + 1]) < 0.1) and np.all(
                        abs(freq[index] - freq[index - 1]) < 0.1):
                    signal = True
            if is_process:
                signal = True
            if signal:
                result_A.append(A[index])
                result_a.append(a[index])
                result_freq.append(freq[index])
                result_pos.append(pos[index])
        return result_A, result_a, result_freq, result_pos

    def run(self, return_path):
        from collections import Counter

        ref_data_path = os.path.join(self.save_path, self.file_name + "_{}".format(self.read_number) + "_ref.npy")
        mut_data_path = os.path.join(self.save_path, self.file_name + "_{}".format(self.read_number) + "_mut.npy")
        freq_data_path = os.path.join(self.save_path, self.file_name + "_{}".format(self.read_number) + "_freq.npy")
        pos_data_path = os.path.join(self.save_path, self.file_name + "_{}".format(self.read_number) + "_pos.npy")

        data = pd.read_csv(self.file_path, header=None, sep=None, engine='python')
        chrome_set = list(Counter(list(data[0])).keys())
        max_chrome = len(chrome_set)
        A_data, a_data, freq_data, position = [], [], [], []  # 保存最后结果
        temp_A, temp_a, temp_freq, temp_pos = [], [], [], []  # 保存一次筛选后结果
        if return_path:
            return ref_data_path, mut_data_path, freq_data_path, pos_data_path, chrome_set
        for chrome in tqdm(range(int(max_chrome)), desc="pretreatment"):
            row_Ap, row_ap = [], []
            row_freq, row_pos = [], []
            chrome_data = data.loc[data[0] == chrome_set[chrome]]
            try:
                for _, row in chrome_data.iterrows():
                    Ap = np.array(row[4::2])
                    ap = np.array(row[5::2])
                    # 一次筛选
                    if self.first_check(Ap, ap):
                        row_Ap.append(Ap)
                        row_ap.append(ap)
                        row_freq.append(Ap / (ap + Ap))
                        row_pos.append(row[1])
                temp_A.append(row_Ap)
                temp_a.append(row_ap)
                temp_freq.append(row_freq)
                temp_pos.append(row_pos)
            except:
                print(chrome_data)
        # 二次筛选
        for A, a, freq, pos, chrome in tqdm(zip(temp_A, temp_a, temp_freq, temp_pos, [i for i in range(len(temp_a))]), desc="pretreatment"):
            chrome_A, chrome_a, chrome_freq, chrome_pos = self.second_check(A, a, freq, pos, self.step_three)
            A_data.append(chrome_A)
            a_data.append(chrome_a)
            freq_data.append(chrome_freq)
            position.append(chrome_pos)

        np.save(ref_data_path, A_data)
        np.save(mut_data_path, a_data)
        np.save(freq_data_path, freq_data)
        np.save(pos_data_path, position)

        return ref_data_path, mut_data_path, freq_data_path, pos_data_path, chrome_set
