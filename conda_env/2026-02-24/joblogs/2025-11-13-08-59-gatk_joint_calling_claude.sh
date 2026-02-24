#!/bin/bash

# 🧬 GATK联合分型脚本 - 大豆疫霉菌变异检测（染色体并行版）
# 作者: 李老师团队
# 日期: $(date +%Y-%m-%d)

# set -e  # 遇到错误立即退出
# set -u  # 使用未定义变量时报错
# set -o pipefail  # 管道命令中任何一个失败都会导致整个管道失败

# 记录开始时间
start_time=$(date +%s)

# ========================================
# 📁 路径设置（使用绝对路径）
# ========================================
PROJECT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌"
GVCF_DIR="${PROJECT_DIR}/02.mapping/vcf"
GENOME="${PROJECT_DIR}/01.data/genome/Phytophthora_sojae_JS2.fa"
OUTPUT_DIR="${PROJECT_DIR}/03.gatk_joint"
LOG_DIR="${OUTPUT_DIR}/logs"
CHR_DIR="${OUTPUT_DIR}/by_chromosome"
TMP_DIR="${OUTPUT_DIR}/tmp"

# 打印配置信息
echo "════════════════════════════════════════════════════════════"
echo "🧬 GATK联合分型 - 染色体并行版"
echo "════════════════════════════════════════════════════════════"
echo "📅 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "📁 项目目录: ${PROJECT_DIR}"
echo "📂 输出目录: ${OUTPUT_DIR}"
echo "🧬 参考基因组: ${GENOME}"
echo "📋 gVCF目录: ${GVCF_DIR}"
echo "════════════════════════════════════════════════════════════"
echo ""

# 创建输出目录
echo "📁 创建输出目录..."
mkdir -p "${OUTPUT_DIR}"
mkdir -p "${LOG_DIR}"
mkdir -p "${CHR_DIR}"
mkdir -p "${TMP_DIR}"

# 验证目录创建成功
if [ ! -d "${OUTPUT_DIR}" ]; then
    echo "❌ 错误: 无法创建输出目录 ${OUTPUT_DIR}"
    exit 1
fi

echo "✅ 输出目录已创建: ${OUTPUT_DIR}"
echo ""

# ========================================
# 🔧 GATK参数设置
# ========================================
GATK="gatk"  # 如果GATK不在PATH中,请修改为完整路径
MAX_THREADS=88    # 最大线程数
MAX_MEM="900G"    # 最大内存
THREADS_PER_CHR=4 # 每个染色体任务的线程数
MEM_PER_CHR="10G" # 每个染色体任务的内存

# 计算最大并发任务数
MAX_JOBS=$((MAX_THREADS / THREADS_PER_CHR))

echo "⚙️  资源配置:"
echo "   • 最大线程数: ${MAX_THREADS}"
echo "   • 最大内存: ${MAX_MEM}"
echo "   • 单染色体线程: ${THREADS_PER_CHR}"
echo "   • 单染色体内存: ${MEM_PER_CHR}"
echo "   • 最大并发任务: ${MAX_JOBS}"
echo ""

# ========================================
# 📝 生成样品列表文件
# ========================================
echo "📋 正在生成样品列表..."
SAMPLE_MAP="${OUTPUT_DIR}/sample_map.txt"

# 检查gVCF目录
if [ ! -d "${GVCF_DIR}" ]; then
    echo "❌ 错误: gVCF目录不存在: ${GVCF_DIR}"
    exit 1
fi

> "${SAMPLE_MAP}"  # 清空文件
sample_count=0

for gvcf in "${GVCF_DIR}"/*.g.vcf.gz; do
    if [ -f "$gvcf" ]; then
        sample_name=$(basename "${gvcf}" .g.vcf.gz)
        echo -e "${sample_name}\t${gvcf}" >> "${SAMPLE_MAP}"
        echo "  ✓ 添加样品: ${sample_name}"
        ((sample_count++))
    fi
done

# 检查是否有样品
if [ ${sample_count} -eq 0 ]; then
    echo "❌ 错误: 在 ${GVCF_DIR} 中未找到任何 *.g.vcf.gz 文件!"
    exit 1
fi

echo "✅ 共找到 ${sample_count} 个样品"
echo "📄 样品列表已保存到: ${SAMPLE_MAP}"
echo ""

# ========================================
# 🧬 提取染色体列表
# ========================================
echo "🧬 正在提取染色体信息..."
CHR_LIST="${OUTPUT_DIR}/chromosome_list.txt"

# 检查参考基因组文件
if [ ! -f "${GENOME}" ]; then
    echo "❌ 错误: 参考基因组文件不存在: ${GENOME}"
    exit 1
fi

# 从参考基因组fai文件获取染色体列表
if [ ! -f "${GENOME}.fai" ]; then
    echo "📊 正在生成参考基因组索引..."
    samtools faidx "${GENOME}"
    if [ $? -ne 0 ]; then
        echo "❌ 错误: 无法生成参考基因组索引"
        exit 1
    fi
fi

awk '{print $1}' "${GENOME}.fai" > "${CHR_LIST}"
chr_count=$(wc -l < "${CHR_LIST}")
echo "✅ 共找到 ${chr_count} 条染色体/scaffold"
echo "📄 染色体列表已保存到: ${CHR_LIST}"

# 显示染色体列表（最多显示10条）
echo ""
echo "📋 染色体列表（前10条）:"
head -10 "${CHR_LIST}" | while read chr; do
    echo "   • ${chr}"
done
if [ ${chr_count} -gt 10 ]; then
    echo "   ... 还有 $((chr_count - 10)) 条"
fi
echo ""

# ========================================
# 🗄️ 步骤1: GenomicsDBImport - 按染色体创建数据库
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "🗄️  步骤1: 按染色体创建GenomicsDB数据库"
echo "════════════════════════════════════════════════════════════"
echo "⚡ 最大并发任务数: ${MAX_JOBS}"
echo ""

# 创建任务计数器
job_count=0
failed_jobs=0

# 为每条染色体创建GenomicsDB
while read chr; do
    # 控制并发数
    while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
        sleep 2
    done
    
    echo "🔄 [$(date '+%H:%M:%S')] 启动任务: ${chr} (GenomicsDB)"
    
    (
        GENOMICS_DB="${CHR_DIR}/${chr}_genomicsdb"
        CHR_TMP="${TMP_DIR}/tmp_${chr}_genomicsdb"
        
        # 创建临时目录
        mkdir -p "${CHR_TMP}"
        
        # 删除旧数据库
        if [ -d "${GENOMICS_DB}" ]; then
            rm -rf "${GENOMICS_DB}"
        fi
        
        ${GATK} --java-options "-Xmx${MEM_PER_CHR} -Djava.io.tmpdir=${CHR_TMP}" GenomicsDBImport \
            --sample-name-map "${SAMPLE_MAP}" \
            --genomicsdb-workspace-path "${GENOMICS_DB}" \
            -L "${chr}" \
            --tmp-dir "${CHR_TMP}" \
            --reader-threads ${THREADS_PER_CHR} \
            --batch-size 50 \
            --genomicsdb-shared-posixfs-optimizations true \
            > "${LOG_DIR}/genomicsdb_${chr}.log" 2>&1
        
        exit_code=$?
        
        # 清理临时目录
        rm -rf "${CHR_TMP}"
        
        if [ ${exit_code} -eq 0 ]; then
            echo "✅ [$(date '+%H:%M:%S')] ${chr} - GenomicsDB创建成功"
        else
            echo "❌ [$(date '+%H:%M:%S')] ${chr} - GenomicsDB创建失败 (查看日志: ${LOG_DIR}/genomicsdb_${chr}.log)"
            exit 1
        fi
    ) &
    
    ((job_count++))
done < "${CHR_LIST}"

# 等待所有GenomicsDB任务完成
echo ""
echo "⏳ 等待所有GenomicsDB任务完成..."
wait
exit_code=$?

if [ ${exit_code} -ne 0 ]; then
    echo "❌ 有GenomicsDB任务失败,请检查日志文件"
    exit 1
fi

echo "✅ 所有染色体的GenomicsDB创建完成!"
echo ""

# ========================================
# 🧬 步骤2: GenotypeGVCFs - 按染色体联合分型
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "🧬 步骤2: 按染色体进行联合分型"
echo "════════════════════════════════════════════════════════════"
echo "⚡ 最大并发任务数: ${MAX_JOBS}"
echo ""

# 重置任务计数器
job_count=0

# 为每条染色体进行联合分型
while read chr; do
    # 控制并发数
    while [ $(jobs -r | wc -l) -ge ${MAX_JOBS} ]; do
        sleep 2
    done
    
    echo "🔄 [$(date '+%H:%M:%S')] 启动任务: ${chr} (GenotypeGVCFs)"
    
    (
        GENOMICS_DB="${CHR_DIR}/${chr}_genomicsdb"
        OUTPUT_VCF="${CHR_DIR}/${chr}.vcf.gz"
        CHR_TMP="${TMP_DIR}/tmp_${chr}_genotype"
        
        # 创建临时目录
        mkdir -p "${CHR_TMP}"
        
        ${GATK} --java-options "-Xmx${MEM_PER_CHR} -Djava.io.tmpdir=${CHR_TMP}" GenotypeGVCFs \
            -R "${GENOME}" \
            -V gendb://"${GENOMICS_DB}" \
            -O "${OUTPUT_VCF}" \
            --tmp-dir "${CHR_TMP}" \
            > "${LOG_DIR}/genotype_${chr}.log" 2>&1
        
        exit_code=$?
        
        # 清理临时目录
        rm -rf "${CHR_TMP}"
        
        if [ ${exit_code} -eq 0 ]; then
            echo "✅ [$(date '+%H:%M:%S')] ${chr} - 联合分型成功"
        else
            echo "❌ [$(date '+%H:%M:%S')] ${chr} - 联合分型失败 (查看日志: ${LOG_DIR}/genotype_${chr}.log)"
            exit 1
        fi
    ) &
    
    ((job_count++))
done < "${CHR_LIST}"

# 等待所有分型任务完成
echo ""
echo "⏳ 等待所有联合分型任务完成..."
wait
exit_code=$?

if [ ${exit_code} -ne 0 ]; then
    echo "❌ 有联合分型任务失败,请检查日志文件"
    exit 1
fi

echo "✅ 所有染色体的联合分型完成!"
echo ""

# ========================================
# 🔗 步骤3: 合并所有染色体的VCF
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "🔗 步骤3: 合并所有染色体的VCF文件"
echo "════════════════════════════════════════════════════════════"

# 生成染色体VCF列表
VCF_LIST="${OUTPUT_DIR}/vcf_list.txt"
> "${VCF_LIST}"

while read chr; do
    vcf_file="${CHR_DIR}/${chr}.vcf.gz"
    if [ -f "${vcf_file}" ]; then
        echo "${vcf_file}" >> "${VCF_LIST}"
    else
        echo "⚠️  警告: 染色体VCF文件不存在: ${vcf_file}"
    fi
done < "${CHR_LIST}"

# 检查VCF列表
vcf_file_count=$(wc -l < "${VCF_LIST}")
echo "📋 找到 ${vcf_file_count} 个染色体VCF文件"

if [ ${vcf_file_count} -eq 0 ]; then
    echo "❌ 错误: 没有找到任何染色体VCF文件"
    exit 1
fi

# 使用GATK GatherVcfs合并
FINAL_VCF="${OUTPUT_DIR}/joint_genotyped.vcf.gz"
GATHER_TMP="${TMP_DIR}/tmp_gather"
mkdir -p "${GATHER_TMP}"

echo ""
echo "🔗 正在合并VCF文件..."

${GATK} --java-options "-Xmx${MAX_MEM} -Djava.io.tmpdir=${GATHER_TMP}" GatherVcfs \
    -I "${VCF_LIST}" \
    -O "${FINAL_VCF}" \
    2>&1 | tee "${LOG_DIR}/gather_vcfs.log"

exit_code=$?

# 清理临时目录
rm -rf "${GATHER_TMP}"

if [ ${exit_code} -eq 0 ]; then
    echo "✅ VCF文件合并成功!"
    echo "📄 最终VCF文件: ${FINAL_VCF}"
else
    echo "❌ VCF文件合并失败,请检查日志: ${LOG_DIR}/gather_vcfs.log"
    exit 1
fi

echo ""

# ========================================
# 📊 步骤4: 生成统计信息
# ========================================
echo "════════════════════════════════════════════════════════════"
echo "📊 步骤4: 生成VCF统计信息"
echo "════════════════════════════════════════════════════════════"

# 统计最终合并文件
if [ -f "${FINAL_VCF}" ]; then
    echo "📊 正在统计最终VCF文件..."
    final_variants=$(zgrep -v "^#" "${FINAL_VCF}" | wc -l)
    final_snps=$(zgrep -v "^#" "${FINAL_VCF}" | awk '{if(length($4)==1 && length($5)==1) print}' | wc -l)
    final_indels=$(zgrep -v "^#" "${FINAL_VCF}" | awk '{if(length($4)!=length($5)) print}' | wc -l)
else
    echo "⚠️  警告: 最终VCF文件不存在"
    final_variants=0
    final_snps=0
    final_indels=0
fi

# 计算运行时间
end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(((runtime % 3600) / 60))
seconds=$((runtime % 60))

echo ""
echo "════════════════════════════════════════════════════════════"
echo "🎉 联合分型完成!"
echo "════════════════════════════════════════════════════════════"
echo "📅 完成时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "⏱️  总运行时间: ${hours}h ${minutes}m ${seconds}s"
echo ""
echo "📊 统计信息:"
echo "   • 样品数量: ${sample_count}"
echo "   • 染色体数量: ${chr_count}"
echo "   • 总变异数: ${final_variants}"
echo "   • SNPs数量: ${final_snps}"
echo "   • InDels数量: ${final_indels}"
echo ""
echo "📂 输出文件:"
echo "   • 最终VCF: ${FINAL_VCF}"
echo "   • 染色体VCF目录: ${CHR_DIR}/"
echo "   • 日志目录: ${LOG_DIR}/"
echo ""
echo "💡 后续操作建议:"
echo "   1. 过滤VCF: gatk VariantFiltration"
echo "   2. 质量控制: bcftools stats ${FINAL_VCF}"
echo "   3. 节省空间: rm -rf ${CHR_DIR}/ ${TMP_DIR}/"
echo "════════════════════════════════════════════════════════════"
echo ""