#!/bin/bash

echo "========================================================"
echo "诊断开始：检查主机环境"
echo "主机上的 CUDA_VISIBLE_DEVICES: $CUDA_VISIBLE_DEVICES"
echo "主机上的 nvidia-smi 输出:"
nvidia-smi
echo "========================================================"

echo ""

echo "========================================================"
echo "诊断核心：在 Apptainer 容器内检查环境"
APPTAINER_IMAGE="/share/apps/containers/parabricks.sif" # 替换成您的容器路径

apptainer exec --nv "$APPTAINER_IMAGE" bash -c "
  echo '--- 进入容器内部 ---';
  echo '容器内的 CUDA_VISIBLE_DEVICES: \$CUDA_VISIBLE_DEVICES';
  echo '容器内的 nvidia-smi 输出:';
  nvidia-smi;
  echo '--- 退出容器内部 ---';
"
echo "========================================================"

echo ""

echo "准备运行 Python 主脚本..."

# 构建索引
bwa index /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/01.data/GCF_000149755.1_P.sojae_V3.0_genomic_modified.fna

# 主脚本
# python3 \
#     /share/org/YZWL/yzwl_lixg/software/scripts/run_parabricks.py \
#     --container /share/apps/containers/parabricks.sif \
#     --input-dir /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/01.data/clean \
#     --output-dir /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/03.mapping \
#     --ref-genome /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/01.data/GCF_000149755.1_P.sojae_V3.0_genomic_modified.fna \
#     --workflow fq2bam_and_call \
#     --num-threads 64 \
#     --max-parallel 1 \
#     --fastq-pattern "*_1.clean.fq.gz"

biopytools parabricks \
  -i /share/org/YZWL/yzwl_lixg/database/yimeijun/resequence/clean \
  -o /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/03.mapping \
  -r /share/org/YZWL/yzwl_lixg/project/17.大豆疫霉基因变异_李磊/01.data/GCF_000149755.1_P.sojae_V3.0_genomic_modified.fna \
  -t 64