import argparse
import os
import sys

from functions.vcf_handle import VCF2Excel
from functions.pretreatment import Pretreatment
from functions.none_pretreatment import NonePretreatment
from functions.Statistic_Methods import Statistic

def main(args):
    if args.p == 1:
        args.p = True
    else:
        args.p = False
    
    if args.p2 == 1:
        args.p2 = True
    else:
        args.p2 = False
    
    if args.p3 == 1:
        args.p3 = True
    else:
        args.p3 = False
    print("-" * 100)
    print("file path:{}\nmethod:{}\nis pretreatment:{}\nread number:{}\nChi-square test:{}\nContinuity test:{}\nsmooth method:{}\nsmooth window size:{}\nthreshold:{}".format(
        args.i, args.m, args.p, args.p1, args.p2, args.p3, args.s, args.w, args.t))
    print("-" * 100)
    # 新建文件夹
    root_path = os.getcwd()
    pretreatment_dir = os.path.join(root_path, "Pretreated_Files")
    nopretreatment_dir = os.path.join(root_path, "NoPretreatment")
    excel_path = os.path.join(root_path, "Excel_Files")
    for path in [pretreatment_dir, nopretreatment_dir, excel_path]:
        if not os.path.exists(path):
            os.mkdir(path)
    # 获得文件名、类型
    temp = args.i.split("/")[-1]
    file_name = temp.split(".")[0]
    file_type = temp.split(".")[1]

    if "vcf" in file_type:
        vcf2excel = VCF2Excel(args.i, file_name, excel_path)
        file_path = vcf2excel.run()
    else:
        file_path = args.i

    if args.p == True:
        return_path = os.path.exists(os.path.join(pretreatment_dir, file_name + "_{}".format(args.p1) + "_freq.npy"))
        pretreat = Pretreatment(args.p1, args.p2, args.p3, file_path, file_name, pretreatment_dir)
        ref_data_path, mut_data_path, freq_data_path, pos_data_path, chrome_set = pretreat.run(return_path)
        if not return_path:
            print("pretreatment & files do not exist")
        else:
            print("pretreatment & files exist")
    else:
        return_path = os.path.exists(os.path.join(nopretreatment_dir, file_name + "_ref.npy"))
        nopretreat = NonePretreatment(file_path, file_name, nopretreatment_dir)
        ref_data_path, mut_data_path, freq_data_path, pos_data_path, chrome_set = nopretreat.run(return_path)
        if not return_path:
            print("nonepretreatment & files do not exist")
        else:
            print("nonepretreatment & files exist")
    print(chrome_set)
    rsp = os.path.join(os.getcwd(), "Results")
    if not os.path.exists(rsp):
        os.mkdir(rsp)
    save_path = os.path.join(rsp, file_name)
    if not os.path.exists(save_path):
        os.mkdir(save_path)
    dl = Statistic(args.m, freq_data_path, ref_data_path, mut_data_path,
                    pos_data_path, chrome_set, args.p1, args.s,
                        args.w, args.t, save_path)

    # 根据参数选择单线程或多线程
    if hasattr(args, 'threads') and args.threads > 1:
        print("Running in multi-thread mode with {} threads".format(args.threads))
        dl.run_parallel(num_threads=args.threads)
    else:
        print("Running in single-thread mode")
        dl.run()

# example:
# python main.py --i /media/xaun/CXX/DeepBSA-terminal/bin/Excel_Files/nc-planthigh-pop1.csv --p False
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    # 载入数据及预处理
    parser.add_argument("--i", default=None, required=True, type=str, help="The input file path(vcf/csv).")
    parser.add_argument("--m", default="DL", required=False, type=str, help="The algorithm(DL/K/ED4/SNP/SmoothG/SmoothLOD/Ridit) used. Default is DL.")
    parser.add_argument("--p", default=1, type=int, help="Whether to pretreatment data(1[True] or 0[False]). Default is True.")
    parser.add_argument("--p1", default=0, type=int, help="Pretreatment step 1: Number of read thread, the SNP whose number lower than it will be filtered. Default is 0.")
    parser.add_argument("--p2", default=1, type=int, help="Pretreatment step 2: Chi-square test(1[True] or 0[False]). Default is 1[True].")
    parser.add_argument("--p3", default=1, type=int, help="Pretreatment step 3: Continuity test(1[True] or 0[False]). Default is 1[True].")
    # 方法选择等
    parser.add_argument("--s", default="LOWESS", type=str,
                    help="The function to smooth the result(Tri-kernel-smooth\LOWESS\Moving Average), Defalut is LOWESS")
    parser.add_argument("--w", default=0, type=int, help="Windows size of LOESS. The number is range from 0-1. 0 presents the best size for minimum AICc. Default is 0(auto).")
    parser.add_argument("--t", default=0, type=float, help="The threshold to find peaks(float). Default is 0(auto)")
    parser.add_argument("--threads", default=1, type=int, help="Number of threads for parallel processing. Default is 1 (single-thread). Use >1 for multi-thread mode.")

    args = parser.parse_args()

    main(args)
