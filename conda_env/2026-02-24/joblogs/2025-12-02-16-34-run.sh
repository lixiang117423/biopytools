#!/bin/bash

# ================= 配置路径 =================
# 1. 设置 Singularity 执行程序路径
SINGULARITY_BIN="/share/org/YZWL/yzwl_lixg/miniforge3/envs/singularity_v.3.8.7/bin/singularity"

# 2. 设置 SIF 镜像文件路径
SIF_IMAGE="/share/org/YZWL/yzwl_lixg/software/singularity/gapit_v1.5.sif"

# 3. 设置 R 脚本名称
R_SCRIPT="run_gapit.R"

# ================= 运行设置 =================

# 获取当前 Shell 所在的目录作为工作目录
WORK_DIR=$(pwd)

echo "Singularity Path: $SINGULARITY_BIN"
echo "Image Path: $SIF_IMAGE"
echo "Working Directory: $WORK_DIR"

# 检查文件是否存在
if [ ! -f "$R_SCRIPT" ]; then
    echo "错误: 找不到 R 脚本: $R_SCRIPT"
    exit 1
fi

# ================= 执行命令 =================
# -B 参数用于挂载目录，确保容器能访问宿主机的文件
# 这里挂载 /share 目录，保证你的软件路径和数据路径都能被容器读取
# exec 命令在容器内执行 Rscript

$SINGULARITY_BIN exec \
    -B /share \
    -B "$WORK_DIR" \
    "$SIF_IMAGE" \
    Rscript "$R_SCRIPT"

echo "任务运行结束。"
