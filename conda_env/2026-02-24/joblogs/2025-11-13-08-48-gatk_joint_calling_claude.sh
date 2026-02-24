#!/bin/bash

# 🧬 GATK联合分型脚本 - 大豆疫霉菌变异检测（染色体并行版）
# 作者: 李老师团队
# 日期: $(date +%Y-%m-%d)

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错

# 记录开始时间
start_time=$(date +%s)

# ========================================
# 📁 路径设置
# ========================================
PROJECT_DIR="/share/org/YZWL/yzwl_lixg/project/19.大豆疫霉菌"
GVCF_DIR="${PROJECT_DIR}/02.mapping/vcf"
GENOME="${PROJECT_DIR}/01.data/genome/Phytophthora_sojae_JS2.fa"
OUTPUT_DIR="${PROJECT_DIR}/03.gatk_joint"
LOG_DIR="${OUTPUT_DIR}/logs"
CHR_DIR="${OUTPUT_DIR}/by_chromosome"

# 创建输出目录
mkdir -p ${OUTPUT_DIR}
mkdir -p ${LOG_DIR}
mkdir -p ${CHR_DIR}

# ========================================
# 🔧 GATK参数设置
# ========================================
GATK="gatk"  # 如果GATK不在PATH中,请修改为完整路径
MAX_THREADS=88    # 最大线程数
MAX_MEM="900G"    # 最大内存
THREADS_PER_CHR=4 # 每个染色体任务的线程数
MEM_PER_CHR="10G" # 每个染色体任务的内存

# ========================================
# 📝 生成样品列表文件
# ========================================
echo "📋 正在生成样品列表..."
SAMPLE_MAP="${OUTPUT_DIR}/sample_map.txt"

> ${SAMPLE_MAP}  # 清空文件
for gvcf in ${GVCF_DIR}/*.g.vcf.gz; do
    if [ -f "$gvcf" ]; then
        sample_name=$(basename ${gvcf} .g.vcf.gz)
        echo -e "${sample_name}\t${gvcf}" >> ${SAMPLE_MAP}
        echo "  ✓ 添加样品: ${sample_name}"
    fi
done

# 检查是否有样品
sample_count=$(wc -l < ${SAMPLE_MAP})
if [ ${sample_count} -eq 0 ]; then
    echo "❌ 错误: 未找到任何gVCF文件!"
    exit 1
fi

echo "✅ 共找到 ${sample_count} 个样品"

# ========================================
# 🧬 提取染色体列表
# ========================================
echo ""
echo "🧬 正在提取染色体信息..."
CHR_LIST="${OUTPUT_DIR}/chromosome_list.txt"

# 从参考基因组fai文件获取染色体列表
if [ ! -f "${GENOME}.fai" ]; then
    echo "📊 生成参考基因组索引..."
    samtools faidx ${GENOME}
fi

awk '{print $1}' ${GENOME}.fai > ${CHR_LIST}
chr_count=$(wc -l < ${CHR_LIST})
echo "✅ 共找到 ${chr_count} 条染色体/scaffold"

# 显示染色体列表（最多显示10条）
echo "📋 染色体列表（前10条）:"
head -10 ${CHR_LIST} | while read chr; do
    echo "   • ${chr}"
done
if [ ${chr_count} -gt 10 ]; then
    echo "   ... 还有 $((chr_count - 10)) 条"
fi

# ========================================
# 🗄️ 步骤1: GenomicsDBImport - 按染色体创建数据库
# ========================================
echo ""
echo "🗄️  步骤1: 按染色体创建GenomicsDB数据库..."
echo "⚡ 并行处理 - 最大并发任务数: $((MAX_THREADS / THREADS_PER_CHR))"

# 创建任务计数器
job_count=0
max_jobs=$((MAX_THREADS / THREADS_PER_CHR))

# 为每条染色体创建GenomicsDB
while read chr; do
    # 控制并发数
    while [ $(jobs -r | wc -l) -ge ${max_jobs} ]; do
        sleep 2
    done
    
    echo "  🔄 正在处理染色体: ${chr}"
    
    (
        GENOMICS_DB="${CHR_DIR}/${chr}_genomicsdb"
        
        # 删除旧数据库
        if [ -d "${GENOMICS_DB}" ]; then
            rm -rf ${GENOMICS_DB}
        fi
        
        ${GATK} --java-options "-Xmx${MEM_PER_CHR}" GenomicsDBImport \
            --sample-name-map ${SAMPLE_MAP} \
            --genomicsdb-workspace-path ${GENOMICS_DB} \
            -L ${chr} \
            --tmp-dir ${OUTPUT_DIR}/tmp_${chr} \
            --reader-threads ${THREADS_PER_CHR} \
            --batch-size 50 \
            > ${LOG_DIR}/genomicsdb_${chr}.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "    ✅ ${chr} - GenomicsDB创建成功"
        else
            echo "    ❌ ${chr} - GenomicsDB创建失败"
            exit 1
        fi
    ) &
    
    ((job_count++))
done < ${CHR_LIST}

# 等待所有GenomicsDB任务完成
echo "⏳ 等待所有GenomicsDB任务完成..."
wait
echo "✅ 所有染色体的GenomicsDB创建完成!"

# ========================================
# 🧬 步骤2: GenotypeGVCFs - 按染色体联合分型
# ========================================
echo ""
echo "🧬 步骤2: 按染色体进行联合分型..."
echo "⚡ 并行处理 - 最大并发任务数: $((MAX_THREADS / THREADS_PER_CHR))"

# 重置任务计数器
job_count=0

# 为每条染色体进行联合分型
while read chr; do
    # 控制并发数
    while [ $(jobs -r | wc -l) -ge ${max_jobs} ]; do
        sleep 2
    done
    
    echo "  🔄 正在分型染色体: ${chr}"
    
    (
        GENOMICS_DB="${CHR_DIR}/${chr}_genomicsdb"
        OUTPUT_VCF="${CHR_DIR}/${chr}.vcf.gz"
        
        ${GATK} --java-options "-Xmx${MEM_PER_CHR}" GenotypeGVCFs \
            -R ${GENOME} \
            -V gendb://${GENOMICS_DB} \
            -O ${OUTPUT_VCF} \
            --tmp-dir ${OUTPUT_DIR}/tmp_${chr} \
            > ${LOG_DIR}/genotype_${chr}.log 2>&1
        
        if [ $? -eq 0 ]; then
            echo "    ✅ ${chr} - 联合分型成功"
        else
            echo "    ❌ ${chr} - 联合分型失败"
            exit 1
        fi
    ) &
    
    ((job_count++))
done < ${CHR_LIST}

# 等待所有分型任务完成
echo "⏳ 等待所有联合分型任务完成..."
wait
echo "✅ 所有染色体的联合分型完成!"

# ========================================
# 🔗 步骤3: 合并所有染色体的VCF
# ========================================
echo ""
echo "🔗 步骤3: 合并所有染色体的VCF文件..."

# 生成染色体VCF列表
VCF_LIST="${OUTPUT_DIR}/vcf_list.txt"
> ${VCF_LIST}
while read chr; do
    echo "${CHR_DIR}/${chr}.vcf.gz" >> ${VCF_LIST}
done < ${CHR_LIST}

# 使用GATK GatherVcfs合并
FINAL_VCF="${OUTPUT_DIR}/joint_genotyped.vcf.gz"

${GATK} --java-options "-Xmx${MAX_MEM}" GatherVcfs \
    -I ${VCF_LIST} \
    -O ${FINAL_VCF} \
    2>&1 | tee ${LOG_DIR}/gather_vcfs.log

if [ $? -eq 0 ]; then
    echo "✅ VCF文件合并成功!"
else
    echo "❌ VCF文件合并失败,请检查日志: ${LOG_DIR}/gather_vcfs.log"
    exit 1
fi

# ========================================
# 📊 步骤4: 生成统计信息
# ========================================
echo ""
echo "📊 步骤4: 生成VCF统计信息..."

# 统计每条染色体的变异数量
echo "📋 各染色体变异统计:"
total_variants=0
total_snps=0
total_indels=0

while read chr; do
    chr_vcf="${CHR_DIR}/${chr}.vcf.gz"
    if [ -f "${chr_vcf}" ]; then
        variants=$(zgrep -v "^#" ${chr_vcf} | wc -l)
        snps=$(zgrep -v "^#" ${chr_vcf} | awk '{if(length($4)==1 && length($5)==1) print}' | wc -l)
        indels=$(zgrep -v "^#" ${chr_vcf} | awk '{if(length($4)!=length($5)) print}' | wc -l)
        
        total_variants=$((total_variants + variants))
        total_snps=$((total_snps + snps))
        total_indels=$((total_indels + indels))
        
        printf "   • %-20s: %10d 变异 (%d SNPs, %d InDels)\n" "${chr}" ${variants} ${snps} ${indels}
    fi
done < ${CHR_LIST}

# 统计最终合并文件
echo ""
echo "📊 最终VCF文件统计:"
final_variants=$(zgrep -v "^#" ${FINAL_VCF} | wc -l)
final_snps=$(zgrep -v "^#" ${FINAL_VCF} | awk '{if(length($4)==1 && length($5)==1) print}' | wc -l)
final_indels=$(zgrep -v "^#" ${FINAL_VCF} | awk '{if(length($4)!=length($5)) print}' | wc -l)

echo ""
echo "════════════════════════════════════════════════════════════"
echo "🎉 联合分型完成!"
echo "════════════════════════════════════════════════════════════"
echo "📁 输出文件: ${FINAL_VCF}"
echo "📊 统计信息:"
echo "   • 样品数量: ${sample_count}"
echo "   • 染色体数量: ${chr_count}"
echo "   • 总变异数: ${final_variants}"
echo "   • SNPs数量: ${final_snps}"
echo "   • InDels数量: ${final_indels}"
echo "════════════════════════════════════════════════════════════"
echo ""
echo "📂 中间文件目录:"
echo "   • 按染色体VCF: ${CHR_DIR}/"
echo "   • 日志文件: ${LOG_DIR}/"
echo ""
echo "💡 提示:"
echo "   - 使用了 ${MAX_THREADS} 个线程并行处理"
echo "   - 最大内存使用: ${MAX_MEM}"
echo "   - 可以使用bcftools或GATK进一步过滤VCF文件"
echo "   - 如需删除中间文件节省空间: rm -rf ${CHR_DIR}"
echo ""
echo "🚀 性能统计:"
end_time=$(date +%s)
runtime=$((end_time - start_time))
hours=$((runtime / 3600))
minutes=$(((runtime % 3600) / 60))
seconds=$((runtime % 60))
echo "   • 总运行时间: ${hours}h ${minutes}m ${seconds}s"
echo ""