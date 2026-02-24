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
bwa index /share/org/YZWL/yzwl_lixg/tmp/gaoyong/genome.fa

# 主脚本
python3 \
    /share/org/YZWL/yzwl_lixg/software/scripts/run_parabricks.py \
    --container /share/apps/containers/parabricks.sif \
    --input-dir /share/org/YZWL/yzwl_lixg/tmp/gaoyong/data \
    --output-dir /share/org/YZWL/yzwl_lixg/tmp/gaoyong/bam \
    --ref-genome /share/org/YZWL/yzwl_lixg/tmp/gaoyong/genome.fa \
    --workflow fq2bam_and_call \
    --num-threads 16 \
    --max-parallel 1 \
    --fastq-pattern "*_1.clean.fq.gz"