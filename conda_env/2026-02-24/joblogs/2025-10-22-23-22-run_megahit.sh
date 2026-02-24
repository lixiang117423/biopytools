#!/bin/bash

# --- MEGAHIT Metagenome Assembly Script ---

# ** 脚本说明 **
#
# 该脚本使用MEGAHIT软件对双端测序数据进行宏基因组组装。
#
# ** 参数 **
#
# - `-1`: 双端测序的read 1文件路径.
# - `-2`: 双端测序的read 2文件路径.
# - `-o`: 输出结果的文件夹路径.
# - `-t`: 使用的CPU线程数 (可以根据服务器性能自行调整).
#
# --------------------------------------------------

# --- 请根据您的实际情况修改以下变量 ---

# 定义输入文件所在的文件夹
INPUT_DIR="/share/org/YZWL/yzwl_lixg/tmp/megahit/data"

# 定义输出结果的文件夹
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/tmp/megahit/megahit"

# 定义输入的Read 1和Read 2文件名
READ1="all_1.fq"
READ2="all_2.fq"

# 定义使用的CPU线程数 (建议使用服务器的核心数)
THREADS=$(nproc) # 自动获取当前系统的CPU核心数

# --- 脚本主体 ---

echo "--- 开始运行 MEGAHIT 宏基因组组装 ---"
echo "输入文件夹: ${INPUT_DIR}"
echo "输出文件夹: ${OUTPUT_DIR}"
echo "Read 1 文件: ${READ1}"
echo "Read 2 文件: ${READ2}"
echo "使用线程数: ${THREADS}"
echo "-----------------------------------------"

# 检查输出文件夹是否存在，如果不存在则创建
if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "输出文件夹不存在，正在创建: ${OUTPUT_DIR}"
    mkdir -p "${OUTPUT_DIR}"
fi

# 运行 MEGAHIT 命令
megahit -1 "${INPUT_DIR}/${READ1}" \
        -2 "${INPUT_DIR}/${READ2}" \
        -o "${OUTPUT_DIR}" \
        -t "${THREADS}"

# 检查MEGAHIT是否成功运行
if [ $? -eq 0 ]; then
    echo "--- MEGAHIT 组装成功完成 ---"
    echo "最终的contigs文件位于: ${OUTPUT_DIR}/final.contigs.fa"
    echo "---------------------------------"
else
    echo "--- MEGAHIT 运行失败 ---"
    echo "请检查错误信息。"
    echo "---------------------------"
fi
