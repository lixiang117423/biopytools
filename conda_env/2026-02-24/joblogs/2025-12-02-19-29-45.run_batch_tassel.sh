#!/bin/bash

# ================= 配置区域 =================
# 工作目录
WORK_DIR="/share/org/YZWL/yzwl_lixg/project/94.rice_gas/11.tassel_gwas"

# 绝对路径的输入文件
VCF_FILE="${WORK_DIR}/variation.filtered.snp.maf5.vcf.gz"
PHENO_SOURCE="${WORK_DIR}/17.微生物和甲烷菌用于GWAS的数据.txt"

# TASSEL 运行脚本路径 (使用 $HOME 代表 ~)
RUN_SCRIPT="$HOME/software/scripts/44.run_tassel.sh"

# 错误日志文件 (记录失败的表型名称)
FAIL_LOG="${WORK_DIR}/failed_phenotypes.txt"

# ================= 脚本逻辑 =================

# 1. 进入工作目录
if ! cd "${WORK_DIR}"; then
    echo "严重错误: 无法进入工作目录 ${WORK_DIR}，脚本终止。"
    exit 1
fi

# 初始化/清空错误日志
: > "${FAIL_LOG}"

# 2. 获取表型文件的总列数 (基于第一行表头)
TOTAL_COLS=$(head -n 1 "${PHENO_SOURCE}" | awk '{print NF}')

echo "检测到总列数: ${TOTAL_COLS}"
echo "开始遍历处理第 2 列到第 ${TOTAL_COLS} 列..."
echo "运行失败的表型将被记录在: ${FAIL_LOG}"

# 3. 循环遍历从第2列到最后一列
for (( i=2; i<=TOTAL_COLS; i++ )); do
    
    # 获取当前列的表头名称（作为表型ID和文件夹名）
    PHENO_NAME=$(head -n 1 "${PHENO_SOURCE}" | awk -v col=$i '{print $col}')
    
    # 去除名称中可能存在的特殊字符，例如 / 替换为 _
    DIR_NAME=${PHENO_NAME//\//_}
    
    echo "------------------------------------------------------"
    echo "正在处理表型: ${PHENO_NAME} (列索引: $i)"
    
    # 创建子文件夹
    if [ ! -d "${DIR_NAME}" ]; then
        mkdir -p "${DIR_NAME}"
    fi
    
    # 尝试进入子文件夹 (修改点：失败则跳过本次循环，不退出脚本)
    if ! cd "${DIR_NAME}"; then
        echo "警告: 无法进入目录 ${DIR_NAME}，跳过该表型。"
        echo "${PHENO_NAME} (Directory Error)" >> "${FAIL_LOG}"
        continue
    fi
    
    # 定义提取出来的单个表型文件名
    SINGLE_PHENO_FILE="${PHENO_NAME}.pheno.txt"
    
    # 提取第1列和第i列
    # 使用 awk 判断是否成功执行
    awk -v col=$i 'BEGIN{OFS="\t"} {print $1, $col}' "${PHENO_SOURCE}" > "${SINGLE_PHENO_FILE}"
    
    if [ $? -ne 0 ]; then
        echo "错误: 提取表型数据失败，跳过。"
        echo "${PHENO_NAME} (Extraction Error)" >> "${FAIL_LOG}"
        cd .. # 返回上一级
        continue
    fi
    
    # 运行 TASSEL 脚本 (修改点：捕获运行结果)
    echo "执行 TASSEL 脚本..."
    
    # 这里的 bash 命令被放在 if 条件中
    if bash "${RUN_SCRIPT}" -i "${VCF_FILE}" -p "${SINGLE_PHENO_FILE}" --keep-temp; then
        echo ">>> 表型 ${PHENO_NAME} 运行成功。"
    else
        echo ">>> 警告: 表型 ${PHENO_NAME} 运行失败！已记录到日志。"
        # 将失败的表型名称追加写入日志文件
        echo "${PHENO_NAME}" >> "${FAIL_LOG}"
        # 注意：这里没有 exit，脚本会继续向下运行
    fi
    
    # 返回上一级目录，准备处理下一个
    cd ..
    
done

echo "======================================================"
echo "所有表型处理流程结束。"
if [ -s "${FAIL_LOG}" ]; then
    echo "以下表型运行失败，请检查 ${FAIL_LOG} :"
    cat "${FAIL_LOG}"
else
    echo "完美！所有表型均成功运行。"
fi
echo "======================================================"