#!/bin/bash

# 这是一个使用 NVIDIA Parabricks 和 biopytools 的混合变异检测流程脚本。
# 流程分为两个核心阶段：
# 1. (GPU加速) 使用Parabricks为每个样本单独生成GVCF，以获得最大速度。
# 2. (CPU) 使用 biopytools 封装的GATK工具对所有GVCF进行联合基因分型和过滤。
# 脚本已为四倍体生物进行了适配 (PLOIDY=4)。

# --- 安全设置 ---
set -e
set -o pipefail

# --- 用户配置区 ---

# 1. 输入文件路径
BAM_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads/bam"
REF_GENOME="/share/org/YZWL/yzwl_lixg/project/16.荠菜/01.data/genome/genome.fa"

# 2. 输出文件夹路径
OUTPUT_DIR="/share/org/YZWL/yzwl_lixg/project/16.荠菜/19.筛选唯一比对上的reads/vcf_parabricks" # 建议使用新目录

# 3. Parabricks 和系统资源配置
PARABRICKS_SIF="/share/apps/containers/parabricks.sif" # Parabricks容器的绝对路径
BASE_BIND_PATH="/share/org/YZWL/yzwl_lixg/project/16.荠菜" # 绑定一个包含所有路径的父目录
PB_RUNNER="apptainer exec --nv --bind ${BASE_BIND_PATH}:${BASE_BIND_PATH}"
THREADS=80 # 用于biopytools的线程数

# 4. 生物学参数
PLOIDY=2 # 【重要】为四倍体生物设置倍性

# --- 脚本主体 ---

# echo "================================================="
# echo "=== 混合流程启动 (Parabricks + biopytools) ==="
# echo "================================================="
# echo "输入 BAM 文件夹: $BAM_DIR"
# echo "输出主文件夹: $OUTPUT_DIR"
# echo "参考基因组: $REF_GENOME"
# echo "倍性 (Ploidy): $PLOIDY"
# echo "================================================="

# # --- 准备工作: 创建目录和索引 ---
# echo -e "\n--- 准备工作: 创建目录结构和准备参考基因组 ---"
# GVCF_DIR="$OUTPUT_DIR/gvcfs"
# mkdir -p $GVCF_DIR

# # 参考基因组索引（如果已存在则跳过）
# if [ ! -f "${REF_GENOME}.fai" ]; then samtools faidx $REF_GENOME; fi
# REF_DICT=$(echo $REF_GENOME | sed 's/\.fasta$/.dict/' | sed 's/\.fa$/.dict/')
# if [ ! -f "$REF_DICT" ]; then 
#     echo "创建参考基因组字典文件..."
#     gatk CreateSequenceDictionary -R $REF_GENOME -O $REF_DICT
# fi
# echo "准备工作完成。"

# # ===================================================================
# # === 阶段一: (GPU加速) 为每个样本生成GVCF (pbrun haplotypecaller) ===
# # ===================================================================
# echo -e "\n--- 阶段一: 开始为每个样本生成GVCF文件 (Ploidy=$PLOIDY) ---"

# for input_bam in "$BAM_DIR"/*.bam
# do
#     SAMPLE_NAME=$(basename "$input_bam" .bam)
#     echo -e "\n>>> 正在处理样本: $SAMPLE_NAME <<<"

#     gvcf="$GVCF_DIR/${SAMPLE_NAME}.g.vcf.gz"
    
#     if [ -f "$gvcf" ]; then
#         echo "样本 $SAMPLE_NAME 的GVCF文件已存在，跳过 haplotypecaller。"
#         continue
#     fi

#     if [ ! -f "${input_bam}.bai" ]; then
#         echo "为 $input_bam 创建索引..."
#         samtools index "$input_bam"
#     fi

#     # 运行 Parabricks haplotypecaller 生成GVCF
#     # 【【【核心命令: pbrun haplotypecaller --ploidy N --gvcf】】】
#     $PB_RUNNER $PARABRICKS_SIF pbrun haplotypecaller \
#         --ref $REF_GENOME \
#         --in-bam "$input_bam" \
#         --out-variants "$gvcf" \
#         --gvcf \
#         --ploidy $PLOIDY

#     echo "样本 $SAMPLE_NAME 的GVCF生成完毕。"
# done
# echo -e "\n--- 阶段一完成：所有样本的GVCF文件已生成于 $GVCF_DIR ---"


# ===================================================================
# === 阶段二: (CPU) 联合基因分型与过滤 (biopytools) ===
# ===================================================================
echo -e "\n--- 阶段二: 开始进行联合基因分型与过滤 ---"

# 步骤 2.1: 使用 biopytools gatk-joint 进行联合基因分型
echo -e "\n--- 步骤 2.1: 使用 biopytools gatk-joint 进行联合基因分型 ---"
# NOTE: 请确认 biopytools gatk-joint 是否会自动处理或需要额外参数来指定非二倍体生物。
# 该工具可能会默认使用二倍体。如果结果不符合预期，请查阅其帮助文档。
biopytools gatk-joint \
    -i $GVCF_DIR \
    -o $OUTPUT_DIR \
    -r $REF_GENOME

# 步骤 2.2: 使用 biopytools filter-snp-indel 进行变异过滤
echo -e "\n--- 步骤 2.2: 使用 biopytools filter-snp-indel 进行硬过滤 ---"
# 假设 gatk-joint 的输出文件名为 combined.g.vcf 或 combined.vcf.gz，请根据实际情况修改
RAW_VCF_NAME="combined.g.vcf" # 根据您的示例命名
RAW_VCF_PATH="$OUTPUT_DIR/$RAW_VCF_NAME" 

if [ ! -f "$RAW_VCF_PATH" ]; then
    echo "错误: 未找到预期的原始VCF文件: $RAW_VCF_PATH"
    echo "请检查 biopytools gatk-joint 的输出文件名是否为 '$RAW_VCF_NAME'"
    exit 1
fi

biopytools filter-snp-indel \
    -i $RAW_VCF_PATH \
    -o $OUTPUT_DIR \
    -t $THREADS

# 最终文件的命名可能需要根据 biopytools 的实际输出来确认
FINAL_VCF_PATH="$OUTPUT_DIR/$(basename $RAW_VCF_PATH .vcf).filtered.vcf"
echo "变异过滤完成。"
echo "最终分析就绪的VCF文件可能位于: $FINAL_VCF_PATH (请根据biopytools的实际输出确认)"
echo "================================================="
echo "=== 所有混合分析流程成功结束! ==="
echo "================================================="