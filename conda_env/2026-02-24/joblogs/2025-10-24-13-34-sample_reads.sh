#!/bin/bash

# --- 设置 ---
INPUT_DIR="/share/org/YZWL/yzwl_lixg/tmp/megahit/data"
OUTPUT_DIR="${INPUT_DIR}/subsampled" # 将抽样后的文件放入一个新文件夹

READ1_IN="${INPUT_DIR}/all_1.fq"
READ2_IN="${INPUT_DIR}/all_2.fq"

READ1_OUT="${OUTPUT_DIR}/all_1.sub50.fq"
READ2_OUT="${OUTPUT_DIR}/all_2.sub50.fq"

# 抽样比例 (0.5 = 50%)
FRACTION=0.4

# 随机种子 (一个固定的数字，确保两次抽样结果可复现且配对)
SEED=100

# --- 脚本主体 ---
echo "--- 开始对数据进行随机抽样 ---"
echo "抽样比例: ${FRACTION}"

# 创建输出文件夹
mkdir -p "${OUTPUT_DIR}"

# 对 R1 文件进行抽样
echo "正在抽样文件: ${READ1_IN}"
seqtk sample -s${SEED} ${READ1_IN} ${FRACTION} > ${READ1_OUT}

# 使用相同的种子对 R2 文件进行抽样
echo "正在抽样文件: ${READ2_IN}"
seqtk sample -s${SEED} ${READ2_IN} ${FRACTION} > ${READ2_OUT}

echo "--- 抽样完成 ---"
echo "输出文件:"
echo "${READ1_OUT}"
echo "${READ2_OUT}"
