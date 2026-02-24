#!/bin/bash

# 🧬 完整的BSA分析流程脚本
# 包含：质量控制 → 单样品变异检测 → 分组合并BSA分析
# 作者: [Your Name]
# 日期: $(date)

# set -e  # 遇到错误立即退出
# set -u  # 使用未定义变量时报错

# 加载GTX环境
echo "🔧 加载GTX环境..."
source ~/.bashrc
module load gtx/2.2.1

# ============================================================================
# 🔧 参数配置区域 - 请在此处修改你的项目参数
# ============================================================================

# 📁 基本路径配置
RAW_DIR="~/project/16.荠菜/01.data/raw"                    # 原始fastq文件路径
CLEAN_DIR="~/project/16.荠菜/01.data/clean"               # fastp输出路径  
GENOME_FILE="~/project/16.荠菜/01.data/genome/genome.fa"   # 参考基因组路径
EACH_OUTPUT_DIR="~/project/16.荠菜/02.each"                # 单样品变异检测输出路径
BSA_OUTPUT_DIR="~/project/16.荠菜/03.bsa"                  # BSA分析输出路径

# 📋 分组信息文件
GROUP_INFO_FILE="~/project/16.荠菜/01.data/sample_groups.txt"      # 分组信息文件路径

# 🏷️ 文件命名格式
RAW_PATTERN="_1.fq.gz"                                     # 原始文件命名格式 (默认{sample}_1.fq.gz)

# ⚡ 线程配置
FASTP_THREADS=16                                          # fastp线程数
GTX_THREADS=80                                            # gtx线程数
GTX_JOINT_THREADS=80                                      # gtx联合调用线程数

# 🎯 过滤参数

# 小规模项目 (< 50个样本)
MAF=0.01              # 较宽松，保留更多低频变异
MAX_MISSING=0.9       # 允许20%缺失
HWE_PVALUE=0.0001     # 较宽松，BSA混池本身可能偏离HWE
MIN_DEPTH=5           # 较宽松的最小深度
MAX_DEPTH=500         # 根据测序深度调整
MIN_QUAL=20           # 标准设置

# 中等规模项目 (50-100个样本)
# MAF=0.05              # 适中设置
# MAX_MISSING=0.9       # 允许10%缺失  
# HWE_PVALUE=0.001      # 适中设置
# MIN_DEPTH=5           # 标准最小深度
# MAX_DEPTH=500         # 适中上限
# MIN_QUAL=20           # 标准设置

# 大规模项目 (> 100个样本)
# MAF=0.05              # 可以更严格
# MAX_MISSING=0.95      # 更严格，只允许5%缺失
# HWE_PVALUE=0.01       # 可以更严格
# MIN_DEPTH=8           # 更高的最小深度要求
# MAX_DEPTH=300         # 避免过高深度区域
# MIN_QUAL=30           # 更高质量要求

# ============================================================================
# 🚀 流程开始 - 请不要修改以下代码
# ============================================================================

# 展开波浪号路径
RAW_DIR=$(eval echo ${RAW_DIR})
CLEAN_DIR=$(eval echo ${CLEAN_DIR})
GENOME_FILE=$(eval echo ${GENOME_FILE})
EACH_OUTPUT_DIR=$(eval echo ${EACH_OUTPUT_DIR})
BSA_OUTPUT_DIR=$(eval echo ${BSA_OUTPUT_DIR})
GROUP_INFO_FILE=$(eval echo ${GROUP_INFO_FILE})

# 创建日志目录
LOG_DIR="$(dirname ${RAW_DIR})/logs"
mkdir -p ${LOG_DIR}

echo "🧬=============================================="
echo "🧬    单样品变异检测+BSA变异检测流程开始"
echo "🧬=============================================="
echo "🕒 开始时间: $(date)"
echo "📂 原始数据: ${RAW_DIR}"
echo "📂 清理数据: ${CLEAN_DIR}"
echo "📂 单样品结果: ${EACH_OUTPUT_DIR}"
echo "📂 BSA结果: ${BSA_OUTPUT_DIR}"
echo "📋 分组文件: ${GROUP_INFO_FILE}"

# 记录总开始时间
TOTAL_START=$(date +%s)

# # ============================================================================
# # 🧹 阶段1: 质量控制 (fastp)
# # ============================================================================

# echo ""
# echo "🧹=============================================="
# echo "🧹    阶段1: 质量控制 (fastp)"
# echo "🧹=============================================="
# echo "🕒 开始时间: $(date)"

# # 创建必要目录
# mkdir -p ${CLEAN_DIR}

# # 检查输入文件
# echo "🔍 检查输入文件..."
# if [ ! -d "${RAW_DIR}" ]; then
#     echo "❌ 错误: 原始数据目录不存在: ${RAW_DIR}"
#     exit 1
# fi

# # 根据模式查找文件
# RAW_PATTERN_1="${RAW_PATTERN}"
# RAW_PATTERN_2="${RAW_PATTERN/_1/_2}"  # 将_1替换为_2

# RAW_FILES=$(find ${RAW_DIR} -name "*${RAW_PATTERN_1}" | wc -l)
# echo "📊 找到原始数据文件: ${RAW_FILES} 个 (${RAW_PATTERN_1}格式)"

# if [ ${RAW_FILES} -eq 0 ]; then
#     echo "❌ 错误: 在 ${RAW_DIR} 中未找到 ${RAW_PATTERN_1} 格式的文件"
#     exit 1
# fi

# # 记录fastp开始时间
# FASTP_START=$(date +%s)

# # 运行fastp
# echo "🏃 运行fastp质量控制..."
# echo "⚡ 使用线程数: ${FASTP_THREADS}"

# biopytools fastp \
#     -i ${RAW_DIR} \
#     -o ${CLEAN_DIR} \
#     --read1-suffix ${RAW_PATTERN_1} \
#     --read2-suffix ${RAW_PATTERN_2} \
#     -t ${FASTP_THREADS}

# # 检查fastp结果
# if [ $? -ne 0 ]; then
#     echo "❌ 错误: fastp质量控制失败"
#     exit 1
# fi

# # 计算fastp耗时
# FASTP_END=$(date +%s)
# FASTP_TIME=$((FASTP_END - FASTP_START))

# echo "✅ fastp质量控制完成"
# echo "⏱️ 耗时: ${FASTP_TIME} 秒 ($(date -d@${FASTP_TIME} -u +%H:%M:%S))"

# # 统计clean数据
# CLEAN_FILES=$(find ${CLEAN_DIR} -name "*.clean.fq.gz" -o -name "*.clean.fq" | wc -l)
# echo "📊 生成清理数据文件: ${CLEAN_FILES} 个"

# # ============================================================================
# # 🔍 阶段2: 单样品变异检测 (GTX)
# # ============================================================================

# echo ""
# echo "🔍=============================================="
# echo "🔍    阶段2: 单样品变异检测 (GTX)"
# echo "🔍=============================================="
# echo "🕒 开始时间: $(date)"

# # 创建输出目录
# mkdir -p ${EACH_OUTPUT_DIR}

# # 检查参考基因组
# echo "🔍 检查参考基因组..."
# if [ ! -f "${GENOME_FILE}" ]; then
#     echo "❌ 错误: 参考基因组文件不存在: ${GENOME_FILE}"
#     exit 1
# fi

# echo "📊 参考基因组: $(basename ${GENOME_FILE})"
# echo "📊 参考基因组大小: $(ls -lh ${GENOME_FILE} | awk '{print $5}')"

# # # 验证GTX是否成功加载
# # if ! command -v biopytools &> /dev/null; then
# #     echo "❌ 错误: 无法找到biopytools命令，GTX模块加载失败"
# #     exit 1
# # fi

# # echo "✅ GTX环境加载成功"

# # 记录GTX开始时间
# GTX_EACH_START=$(date +%s)

# # 运行GTX单样品变异检测
# echo "🏃 运行GTX单样品变异检测..."
# echo "⚡ 使用线程数: ${GTX_THREADS}"
# echo "⚡ 联合线程数: ${GTX_JOINT_THREADS}"

# biopytools gtx \
#     -i ${CLEAN_DIR} \
#     -o ${EACH_OUTPUT_DIR} \
#     -r ${GENOME_FILE} \
#     -t ${GTX_THREADS} \
#     --joint-threads ${GTX_JOINT_THREADS} \
#     --enable-joint

# # 检查GTX结果
# if [ $? -ne 0 ]; then
#     echo "❌ 错误: GTX单样品变异检测失败"
#     exit 1
# fi

# # 计算GTX耗时
# GTX_EACH_END=$(date +%s)
# GTX_EACH_TIME=$((GTX_EACH_END - GTX_EACH_START))

# echo "✅ GTX单样品变异检测完成"
# echo "⏱️ 耗时: ${GTX_EACH_TIME} 秒 ($(date -d@${GTX_EACH_TIME} -u +%H:%M:%S))"

# # 验证输出文件
# EACH_VCF="${EACH_OUTPUT_DIR}/joint/merged_gtx.vcf.gz"
# if [ -f "${EACH_VCF}" ]; then
#     echo "✅ 单样品联合VCF文件生成成功"
    
#     if command -v bcftools &> /dev/null; then
#         EACH_VARIANTS=$(bcftools view -H ${EACH_VCF} | wc -l)
#         echo "📊 单样品变异总数: ${EACH_VARIANTS}"
#     fi
# else
#     echo "❌ 错误: 未找到预期的VCF输出文件: ${EACH_VCF}"
#     exit 1
# fi

# # ============================================================================
# # 🎯 阶段2.1: 单样品变异过滤
# # ============================================================================

# echo ""
# echo "🎯=============================================="
# echo "🎯    阶段2.1: 单样品变异过滤"
# echo "🎯=============================================="
# echo "🕒 开始时间: $(date)"

# EACH_FILTER_START=$(date +%s)

# # 提取单样品SNP
# echo "📤 提取单样品SNP变异..."
# gatk --java-options "-Xmx128g" SelectVariants \
#     -V ${EACH_VCF} \
#     -select-type SNP \
#     -O ${EACH_OUTPUT_DIR}/snp.raw.vcf.gz \
#     2>&1 | tee ${LOG_DIR}/select_snp_each.log

# EACH_RAW_SNPS=$(bcftools view -H ${EACH_OUTPUT_DIR}/snp.raw.vcf.gz | wc -l)
# echo "📊 单样品原始SNP数量: ${EACH_RAW_SNPS}"

# # 单样品SNP过滤
# echo "🎯 单样品SNP质量过滤..."
# vcftools \
#     --gzvcf ${EACH_OUTPUT_DIR}/snp.raw.vcf.gz \
#     --maf ${MAF} \
#     --max-missing ${MAX_MISSING} \
#     --hwe ${HWE_PVALUE} \
#     --min-meanDP ${MIN_DEPTH} \
#     --max-meanDP ${MAX_DEPTH} \
#     --minQ ${MIN_QUAL} \
#     --min-alleles 2 \
#     --max-alleles 2 \
#     --remove-indels \
#     --recode --recode-INFO-all \
#     --out ${EACH_OUTPUT_DIR}/snp.filtered \
#     2>&1 | tee ${LOG_DIR}/filter_snp_each.log

# # 压缩和索引单样品SNP文件
# echo "📦 压缩和索引过滤后的单样品SNP文件..."
# bgzip -f ${EACH_OUTPUT_DIR}/snp.filtered.recode.vcf
# tabix -f -p vcf ${EACH_OUTPUT_DIR}/snp.filtered.recode.vcf.gz

# EACH_FILTERED_SNPS=$(bcftools view -H ${EACH_OUTPUT_DIR}/snp.filtered.recode.vcf.gz | wc -l)
# echo "📊 单样品过滤后SNP数量: ${EACH_FILTERED_SNPS}"

# # 提取单样品INDEL（仅用于统计）
# echo "📤 提取单样品INDEL变异 (仅用于统计参考)..."
# gatk --java-options "-Xmx128g" SelectVariants \
#     -V ${EACH_VCF} \
#     -select-type INDEL \
#     -O ${EACH_OUTPUT_DIR}/indel.raw.vcf.gz \
#     2>&1 | tee ${LOG_DIR}/select_indel_each.log

# EACH_RAW_INDELS=$(bcftools view -H ${EACH_OUTPUT_DIR}/indel.raw.vcf.gz | wc -l)
# echo "📊 单样品原始INDEL数量: ${EACH_RAW_INDELS} (仅供参考)"

# # 创建单样品最终过滤文件 (仅SNP)
# echo "🎯 创建单样品最终过滤文件 (仅包含SNP)..."
# cp ${EACH_OUTPUT_DIR}/snp.filtered.recode.vcf.gz ${EACH_OUTPUT_DIR}/Each_final_SNP.vcf.gz
# cp ${EACH_OUTPUT_DIR}/snp.filtered.recode.vcf.gz.tbi ${EACH_OUTPUT_DIR}/Each_final_SNP.vcf.gz.tbi

# echo "📊 单样品最终过滤SNP数量: ${EACH_FILTERED_SNPS}"

# EACH_FILTER_END=$(date +%s)
# EACH_FILTER_TIME=$((EACH_FILTER_END - EACH_FILTER_START))

# echo "✅ 单样品变异过滤完成"
# echo "⏱️ 过滤耗时: ${EACH_FILTER_TIME} 秒 ($(date -d@${EACH_FILTER_TIME} -u +%H:%M:%S))"

# ============================================================================
# 🧬 阶段3: BSA分组合并分析
# ============================================================================

echo ""
echo "🧬=============================================="
echo "🧬    阶段3: BSA分组合并分析"
echo "🧬=============================================="
echo "🕒 开始时间: $(date)"

# 创建BSA输出目录
mkdir -p ${BSA_OUTPUT_DIR}

# 检查分组信息文件
echo "🔍 检查分组信息文件..."
if [ ! -f "${GROUP_INFO_FILE}" ]; then
    echo "❌ 错误: 分组信息文件不存在: ${GROUP_INFO_FILE}"
    exit 1
fi

# 读取分组信息
echo "📋 读取分组信息: $(basename ${GROUP_INFO_FILE})"

# 检查分组文件格式
if ! head -1 "${GROUP_INFO_FILE}" | grep -q "sample"; then
    echo "❌ 错误: 分组文件格式不正确，第一行应包含'sample'和'group'标题"
    exit 1
fi

# 获取所有分组
GROUPS=($(tail -n +2 "${GROUP_INFO_FILE}" | cut -f2 | sort | uniq))
echo "📊 发现分组: ${GROUPS[@]}"

# 统计每组样本数量
for group in "${GROUPS[@]}"; do
    group_count=$(awk -v grp="$group" '$2==grp {count++} END {print count+0}' "${GROUP_INFO_FILE}")
    echo "📊 ${group}组样本数量: ${group_count}"
done

# 记录BSA开始时间
BSA_START=$(date +%s)

# 3.1 按组合并fastq文件
echo ""
echo "📁 步骤3.1: 按组合并fastq文件"

for group in "${GROUPS[@]}"; do
    echo "🔄 处理${group}组..."
    
    # 删除已存在的合并文件
    [ -f "${BSA_OUTPUT_DIR}/${group}_1.clean.fq" ] && rm -f ${BSA_OUTPUT_DIR}/${group}_1.clean.fq
    [ -f "${BSA_OUTPUT_DIR}/${group}_2.clean.fq" ] && rm -f ${BSA_OUTPUT_DIR}/${group}_2.clean.fq
    
    # 获取当前组的样本列表
    group_samples=($(awk -v grp="$group" '$2==grp {print $1}' "${GROUP_INFO_FILE}"))
    
    echo "  📋 ${group}组样本: ${group_samples[@]}"
    
    # 合并每个样本的fastq文件
    for sample in "${group_samples[@]}"; do
        r1_file="${CLEAN_DIR}/${sample}_1.clean.fq.gz"
        r2_file="${CLEAN_DIR}/${sample}_2.clean.fq.gz"
        
        if [ -f "${r1_file}" ] && [ -f "${r2_file}" ]; then
            echo "    ➕ 添加样本 ${sample}"
            zcat ${r1_file} >> ${BSA_OUTPUT_DIR}/${group}_1.clean.fq
            zcat ${r2_file} >> ${BSA_OUTPUT_DIR}/${group}_2.clean.fq
        else
            echo "    ⚠️ 警告: 样本 ${sample} 的clean文件不存在，跳过..."
        fi
    done
    
    # 统计合并结果
    if [ -f "${BSA_OUTPUT_DIR}/${group}_1.clean.fq" ]; then
        reads_count=$(($(wc -l < ${BSA_OUTPUT_DIR}/${group}_1.clean.fq) / 4))
        echo "  ✅ ${group}组合并完成，reads数量: ${reads_count}"
    else
        echo "  ❌ 错误: ${group}组合并失败"
        exit 1
    fi
done

# 3.2 对pool文件进行变异检测
echo ""
echo "🔍 步骤3.2: 对pool文件进行GTX变异检测"

GTX_BSA_START=$(date +%s)

echo "🏃 运行BSA GTX变异检测..."
biopytools gtx \
    -i ${BSA_OUTPUT_DIR} \
    -o ${BSA_OUTPUT_DIR} \
    -r ${GENOME_FILE} \
    -t ${GTX_THREADS} \
    --joint-threads ${GTX_JOINT_THREADS} \
    --enable-joint \
    2>&1 | tee ${LOG_DIR}/gtx_bsa.log

# 检查GTX结果
if [ $? -ne 0 ]; then
    echo "❌ 错误: BSA GTX变异检测失败"
    exit 1
fi

GTX_BSA_END=$(date +%s)
GTX_BSA_TIME=$((GTX_BSA_END - GTX_BSA_START))

echo "✅ BSA GTX变异检测完成"
echo "⏱️ 耗时: ${GTX_BSA_TIME} 秒 ($(date -d@${GTX_BSA_TIME} -u +%H:%M:%S))"

# 验证BSA VCF文件
BSA_VCF="${BSA_OUTPUT_DIR}/joint/merged_gtx.vcf.gz"
if [ -f "${BSA_VCF}" ]; then
    echo "✅ BSA联合VCF文件生成成功"
    
    if command -v bcftools &> /dev/null; then
        BSA_VARIANTS=$(bcftools view -H ${BSA_VCF} | wc -l)
        echo "📊 BSA变异总数: ${BSA_VARIANTS}"
    fi
else
    echo "❌ 错误: 未找到BSA VCF输出文件: ${BSA_VCF}"
    exit 1
fi

# 3.3 过滤BSA VCF文件
echo ""
echo "🎯 步骤3.3: 过滤BSA VCF文件"

FILTER_START=$(date +%s)

# 提取SNP
echo "📤 提取BSA SNP变异..."
gatk --java-options "-Xmx128g" SelectVariants \
    -V ${BSA_VCF} \
    -select-type SNP \
    -O ${BSA_OUTPUT_DIR}/snp.raw.vcf.gz \
    2>&1 | tee ${LOG_DIR}/select_snp_bsa.log

RAW_SNPS=$(bcftools view -H ${BSA_OUTPUT_DIR}/snp.raw.vcf.gz | wc -l)
echo "📊 BSA原始SNP数量: ${RAW_SNPS}"

# SNP过滤
echo "🎯 BSA SNP质量过滤..."
vcftools \
    --gzvcf ${BSA_OUTPUT_DIR}/snp.raw.vcf.gz \
    --maf ${MAF} \
    --max-missing ${MAX_MISSING} \
    --hwe ${HWE_PVALUE} \
    --min-meanDP ${MIN_DEPTH} \
    --max-meanDP ${MAX_DEPTH} \
    --minQ ${MIN_QUAL} \
    --min-alleles 2 \
    --max-alleles 2 \
    --remove-indels \
    --recode --recode-INFO-all \
    --out ${BSA_OUTPUT_DIR}/snp.filtered \
    2>&1 | tee ${LOG_DIR}/filter_snp_bsa.log

# 压缩和索引
echo "📦 压缩和索引过滤后的BSA SNP文件..."
bgzip -f ${BSA_OUTPUT_DIR}/snp.filtered.recode.vcf
tabix -f -p vcf ${BSA_OUTPUT_DIR}/snp.filtered.recode.vcf.gz

FILTERED_SNPS=$(bcftools view -H ${BSA_OUTPUT_DIR}/snp.filtered.recode.vcf.gz | wc -l)
echo "📊 BSA过滤后SNP数量: ${FILTERED_SNPS}"

# 可选：提取INDEL（仅用于统计，不用于BSA分析）
echo "📤 提取BSA INDEL变异 (仅用于统计参考)..."
gatk --java-options "-Xmx128g" SelectVariants \
    -V ${BSA_VCF} \
    -select-type INDEL \
    -O ${BSA_OUTPUT_DIR}/indel.raw.vcf.gz \
    2>&1 | tee ${LOG_DIR}/select_indel_bsa.log

RAW_INDELS=$(bcftools view -H ${BSA_OUTPUT_DIR}/indel.raw.vcf.gz | wc -l)
echo "📊 BSA原始INDEL数量: ${RAW_INDELS} (仅供参考)"

# 创建最终BSA分析文件 (仅SNP)
echo "🎯 创建最终BSA分析文件 (仅包含SNP)..."
cp ${BSA_OUTPUT_DIR}/snp.filtered.recode.vcf.gz ${BSA_OUTPUT_DIR}/BSA_final_SNP.vcf.gz
cp ${BSA_OUTPUT_DIR}/snp.filtered.recode.vcf.gz.tbi ${BSA_OUTPUT_DIR}/BSA_final_SNP.vcf.gz.tbi

FINAL_SNP_COUNT=${FILTERED_SNPS}
echo "📊 最终BSA文件SNP数量: ${FINAL_SNP_COUNT}"

FILTER_END=$(date +%s)
FILTER_TIME=$((FILTER_END - FILTER_START))

echo "✅ BSA VCF过滤完成"
echo "⏱️ 过滤耗时: ${FILTER_TIME} 秒 ($(date -d@${FILTER_TIME} -u +%H:%M:%S))"

# 计算BSA总耗时
BSA_END=$(date +%s)
BSA_TIME=$((BSA_END - BSA_START))

# ============================================================================
# 📊 最终统计和报告
# ============================================================================

echo ""
echo "📊=============================================="
echo "📊    流程完成统计报告"
echo "📊=============================================="

# 计算总耗时
TOTAL_END=$(date +%s)
TOTAL_TIME=$((TOTAL_END - TOTAL_START))

# 生成详细报告
REPORT_FILE="${LOG_DIR}/complete_pipeline_report.txt"

cat > ${REPORT_FILE} << EOF
🧬============================================
    完整BSA分析流程报告
============================================
🕒 分析时间: $(date)
📂 项目目录结构:
   - 原始数据: ${RAW_DIR}
   - 清理数据: ${CLEAN_DIR}  
   - 单样品结果: ${EACH_OUTPUT_DIR}
   - BSA结果: ${BSA_OUTPUT_DIR}
   - 分组文件: ${GROUP_INFO_FILE}

⚙️ 参数配置:
   - 文件命名格式: ${RAW_PATTERN}
   - fastp线程数: ${FASTP_THREADS}
   - GTX线程数: ${GTX_THREADS}
   - GTX联合线程数: ${GTX_JOINT_THREADS}

🎯 过滤参数:
   - 最小等位基因频率: ${MAF}
   - 最大缺失率: ${MAX_MISSING}
   - Hardy-Weinberg平衡p值: ${HWE_PVALUE}
   - 深度范围: ${MIN_DEPTH}-${MAX_DEPTH}
   - 最小质量分数: ${MIN_QUAL}

📊 样本分组信息:
$(for group in "${GROUPS[@]}"; do
    group_count=$(awk -v grp="$group" '$2==grp {count++} END {print count+0}' "${GROUP_INFO_FILE}")
    echo "   - ${group}组: ${group_count}个样本"
done)

📈 变异统计:
   - 单样品原始变异总数: ${EACH_VARIANTS:-"N/A"}
   - 单样品原始SNP数量: ${EACH_RAW_SNPS:-"N/A"}
   - 单样品过滤后SNP数量: ${EACH_FILTERED_SNPS:-"N/A"} (保留率: $(echo "scale=1; ${EACH_FILTERED_SNPS:-0}*100/${EACH_RAW_SNPS:-1}" | bc)%)
   - 单样品原始INDEL数量: ${EACH_RAW_INDELS:-"N/A"} (仅供参考)
   - BSA原始变异数: ${BSA_VARIANTS:-"N/A"}
   - BSA原始SNP数量: ${RAW_SNPS}
   - BSA过滤后SNP数量: ${FILTERED_SNPS} (保留率: $(echo "scale=1; ${FILTERED_SNPS}*100/${RAW_SNPS}" | bc)%)
   - BSA原始INDEL数量: ${RAW_INDELS:-"N/A"} (仅供参考)
   - 最终BSA分析SNP数: ${FINAL_SNP_COUNT} (仅SNP，用于BSA分析)

⏱️ 运行时间统计:
   - 质量控制(fastp): ${FASTP_TIME}秒 ($(date -d@${FASTP_TIME} -u +%H:%M:%S))
   - 单样品变异检测: ${GTX_EACH_TIME}秒 ($(date -d@${GTX_EACH_TIME} -u +%H:%M:%S))
   - 单样品变异过滤: ${EACH_FILTER_TIME}秒 ($(date -d@${EACH_FILTER_TIME} -u +%H:%M:%S))
   - BSA变异检测: ${GTX_BSA_TIME}秒 ($(date -d@${GTX_BSA_TIME} -u +%H:%M:%S))
   - BSA变异过滤: ${FILTER_TIME}秒 ($(date -d@${FILTER_TIME} -u +%H:%M:%S))
   - 总运行时间: ${TOTAL_TIME}秒 ($(date -d@${TOTAL_TIME} -u +%H:%M:%S))

🎯 重要输出文件:
   - 单样品最终过滤文件: ${EACH_OUTPUT_DIR}/Each_final_SNP.vcf.gz (仅SNP)
   - 单样品过滤后SNP文件: ${EACH_OUTPUT_DIR}/snp.filtered.recode.vcf.gz
   - 单样品联合原始VCF: ${EACH_OUTPUT_DIR}/joint/merged_gtx.vcf.gz
   - 最终BSA分析文件: ${BSA_OUTPUT_DIR}/BSA_final_SNP.vcf.gz (仅SNP)
   - BSA过滤后SNP文件: ${BSA_OUTPUT_DIR}/snp.filtered.recode.vcf.gz
   - 详细日志目录: ${LOG_DIR}/

📋 后续分析建议:
   1. 单样品分析：使用 Each_final_SNP.vcf.gz 进行个体遗传变异分析
   2. BSA分析：使用 BSA_final_SNP.vcf.gz 进行BSA统计分析 (仅SNP，推荐)
   3. 可用软件：QTL-seq, BSA-seq等
   4. 进行SNP-index和Δ(SNP-index)计算
   5. 滑动窗口分析定位候选区间
   6. 注意：所有最终文件仅包含SNP，不含INDEL，符合标准分析流程

============================================
EOF

echo "📄 详细报告已生成: ${REPORT_FILE}"

# 显示最终总结
echo ""
echo "🎉=============================================="
echo "🎉    流程成功完成!"
echo "🎉=============================================="
echo "🕒 结束时间: $(date)"
echo "⏱️ 总耗时: ${TOTAL_TIME} 秒 ($(date -d@${TOTAL_TIME} -u +%H:%M:%S))"
echo ""
echo "✅ 主要输出文件:"
echo "   🧬 单样品最终过滤文件: ${EACH_OUTPUT_DIR}/Each_final_SNP.vcf.gz"
echo "   🧬 最终BSA文件 (仅SNP): ${BSA_OUTPUT_DIR}/BSA_final_SNP.vcf.gz"
echo "   📊 详细报告: ${REPORT_FILE}"
echo "   📁 所有结果: ${EACH_OUTPUT_DIR}/ 和 ${BSA_OUTPUT_DIR}/"
echo ""
echo "🔍 压缩pool文件节省空间..."
gzip ${BSA_OUTPUT_DIR}/*_1.clean.fq ${BSA_OUTPUT_DIR}/*_2.clean.fq 2>/dev/null || true

echo ""
echo "💡 下一步："
echo "   单样品分析：使用 Each_final_SNP.vcf.gz 进行个体分析"
echo "   BSA分析：使用 BSA_final_SNP.vcf.gz 进行BSA统计分析 (仅SNP，推荐)"
echo "   推荐软件：QTL-seq, BSA-seq等"
echo "   注意：所有最终文件仅包含过滤后的SNP，符合标准分析流程"
echo ""
echo "🎊 分析完成！享受你的结果吧！ 🎊"