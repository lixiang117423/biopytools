#!/bin/bash

# --- 用户可配置变量 ---

# VCF 输入文件的路径
VCF_FILE="renamed_chromosomes.vcf.gz"

# PLINK 输出文件的前缀
PLINK_PREFIX="renamed_chromosomes"

# ADMIXTURE 运行的线程数
THREADS=8

# --- 脚本主体 ---

# 步骤 1: 使用 PLINK 将 VCF 文件转换为 .bed 格式
# --allow-extra-chr 选项允许非标准的染色体名称
# --set-missing-var-ids @:# 选项为没有ID的变体设置一个唯一的ID
echo "----------------------------------------------------"
echo "步骤 1: 正在将 VCF 文件转换为 PLINK .bed 格式..."
echo "----------------------------------------------------"
plink --vcf ${VCF_FILE} --allow-extra-chr --set-missing-var-ids @:# --make-bed --out ${PLINK_PREFIX}

# 检查 PLINK 是否成功生成了 .bed 文件
if [ ! -f "${PLINK_PREFIX}.bed" ]; then
    echo "错误: PLINK 未能生成 .bed 文件。请检查您的 VCF 文件或 PLINK 安装。"
    exit 1
fi
echo "PLINK .bed 文件已成功生成。"
echo ""

# 步骤 2: 使用 for 循环运行 ADMIXTURE，K 值从 2 到 6
# --cv 标志用于计算交叉验证误差，这有助于确定最佳的 K 值。
# 输出被重定向到日志文件，以便后续查看交叉验证误差。
echo "----------------------------------------------------"
echo "步骤 2: 正在为 K=2 到 6 运行 ADMIXTURE..."
echo "----------------------------------------------------"
for K in {2..6}
do
  echo "  正在运行 K=${K}..."
  admixture --cv ${PLINK_PREFIX}.bed ${K} -j${THREADS} | tee log.K=${K}.out
done
echo "ADMIXTURE 分析完成。"
echo ""

# 步骤 3: 提取交叉验证误差
# 这将帮助您评估哪个 K 值最适合您的数据。
# 交叉验证误差最低的 K 值通常被认为是最佳的。
echo "----------------------------------------------------"
echo "步骤 3: 正在从日志文件中提取交叉验证 (CV) 误差..."
echo "----------------------------------------------------"
grep -h "CV error" log.K=*.out

echo "----------------------------------------------------"
echo "脚本执行完毕。"
echo "请检查生成的 .Q 和 .P 文件以及 log.K=*.out 文件中的交叉验证误差。"
echo "----------------------------------------------------"
