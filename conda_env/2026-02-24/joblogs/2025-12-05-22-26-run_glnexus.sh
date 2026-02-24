#!/bin/bash

# ==============================================================================
# 脚本名称: run_glnexus_fixed.sh
# 功能: 自动修复 gVCF 中的负数 PL 值，并运行 GLnexus 联合变异检测
# ==============================================================================

# --- 错误处理设置 ---
# set -e: 遇到错误立即退出 (注意：在 if 中使用 grep 不会触发退出)
# set -o pipefail: 管道中任意命令失败则视为整体失败
set -e
set -o pipefail

# --- [用户配置区域] 请根据实际情况修改 ---
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/01.data/genome/genome.fa"
GVCF_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint/vcf"
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/21.野生大豆群体/03.glnexus_joint"
GLNEXUS_CONFIG="gatk"

# 线程设置
FIX_THREADS=8      # 修复单个 VCF 时的压缩线程数 (不宜过大，IO是瓶颈)
INDEX_THREADS=16   # 建立索引时的线程数
GLNEXUS_THREADS=32 # GLnexus 主程序线程数 (需要大量内存)
CONVERT_THREADS=64 # 最终转换格式的线程数

# ==============================================================================

# 0. 环境检查
if ! command -v bcftools &> /dev/null; then
    echo "Error: 未找到 bcftools，请先加载环境 (例如: conda activate your_env)"
    exit 1
fi
if ! command -v glnexus_cli &> /dev/null; then
    echo "Error: 未找到 glnexus_cli，请先加载环境"
    exit 1
fi

# echo "========================================"
# echo " 步骤 1: 扫描并修复异常 gVCF 文件"
# echo "========================================"
# cd "$GVCF_DIR" || exit

# # 遍历所有 .g.vcf.gz 文件
# for vcf_file in *.g.vcf.gz; do
#     # 检查文件是否存在
#     [ -e "$vcf_file" ] || continue

#     # 使用 zgrep 快速预检 (只读模式，不消耗过多资源)
#     # 逻辑：查找 "PL" 后跟任意字符直到出现 "-数字"
#     # 使用 if 语句包裹 zgrep，即使没找到也不会触发 set -e 退出
#     if zgrep -m 1 "PL.*-[0-9]" "$vcf_file" > /dev/null; then
#         echo "[发现异常] $vcf_file 含有负数 PL 值，正在修复..."
        
#         # 定义临时文件名
#         temp_vcf="${vcf_file}.fixing.tmp.gz"

#         # 使用 bcftools 剔除 PL < 0 的记录
#         # 注意：这里使用 overwrite 逻辑
#         bcftools view -e 'FORMAT/PL < 0' "$vcf_file" \
#             --threads "$FIX_THREADS" \
#             -Oz -o "$temp_vcf"

#         # 检查生成是否成功且非空
#         if [ -s "$temp_vcf" ]; then
#             # 覆盖原文件 (保持文件名不变)
#             mv "$temp_vcf" "$vcf_file"
            
#             # 删除旧索引 (强制后续步骤重新建立索引)
#             rm -f "${vcf_file}.tbi"
            
#             echo "   -> [修复完成] $vcf_file (旧索引已删除)"
#         else
#             echo "   -> [错误] bcftools 处理失败，保留原文件。请检查磁盘空间。"
#             rm -f "$temp_vcf"
#             exit 1
#         fi
#     else
#         # 仅显示进度，不刷屏
#         # echo "   [正常] $vcf_file" 
#         : # 空命令，相当于 pass
#     fi
# done

# echo "========================================"
# echo " 步骤 2: 检查并重建 VCF 索引"
# echo "========================================"
# # 这一步会为刚刚修复过的文件（没有.tbi）以及原本就缺索引的文件建立索引
# for file in *.g.vcf.gz; do
#     if [ ! -f "${file}.tbi" ]; then
#         echo "正在建立索引: $file"
#         tabix -@ "$INDEX_THREADS" -p vcf "$file"
#     fi
# done

echo "========================================"
echo " 步骤 3: 准备参考基因组 BED"
echo "========================================"
mkdir -p "${OUTPUT_DIR}"

BED_FILE="${OUTPUT_DIR}/genome.bed"
BCF_OUT="${OUTPUT_DIR}/merged_output.bcf"
VCF_OUT="${OUTPUT_DIR}/merged_output.vcf.gz"

# 检查参考基因组索引
if [ ! -f "${REF_GENOME}.fai" ]; then
    echo "正在为参考基因组创建 .fai 索引..."
    samtools faidx "${REF_GENOME}"
fi

# 生成 BED 文件
awk -v OFS='\t' '{print $1, 0, $2}' "${REF_GENOME}.fai" > "${BED_FILE}"
echo "BED 文件路径: ${BED_FILE}"


echo "========================================"
echo " 步骤 4: 清理旧数据库并运行 GLnexus"
echo "========================================"
DB_DIR="${OUTPUT_DIR}/GLnexus.DB"

# 强制清理旧数据库，防止 'duplicate' 报错
if [ -d "${DB_DIR}" ]; then
    echo "警告: 检测到旧数据库目录，正在删除: ${DB_DIR}"
    rm -rf "${DB_DIR}"
fi
# 防止当前目录下有残留
if [ -d "GLnexus.DB" ]; then
    rm -rf "GLnexus.DB"
fi

echo "启动 glnexus_cli (Threads: ${GLNEXUS_THREADS})..."

# 运行 GLnexus
glnexus_cli \
  --config "${GLNEXUS_CONFIG}" \
  --bed "${BED_FILE}" \
  --dir "${DB_DIR}" \
  --threads "${GLNEXUS_THREADS}" \
  "${GVCF_DIR}"/*.g.vcf.gz \
  > "${BCF_OUT}"

echo "GLnexus 运行完成，BCF 输出至: ${BCF_OUT}"


echo "========================================"
echo " 步骤 5: 转换为 VCF 并建立索引"
echo "========================================"

if [ -s "${BCF_OUT}" ]; then
    bcftools view --threads "${CONVERT_THREADS}" "${BCF_OUT}" | bgzip -c -@ "${CONVERT_THREADS}" > "${VCF_OUT}"
    bcftools index -t "${VCF_OUT}"
    
    echo "#################################################"
    echo " 所有任务圆满完成！"
    echo " 最终文件: ${VCF_OUT}"
    echo "#################################################"
else
    echo "Error: BCF 输出文件为空，GLnexus 可能运行失败，请检查上面的日志。"
    exit 1
fi