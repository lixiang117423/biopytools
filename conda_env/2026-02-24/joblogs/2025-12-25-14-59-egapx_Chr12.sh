#!/bin/bash
# EGAPx 运行脚本 - Chr12
# 自动生成时间: 2025-12-25 14:51:42

#cd "/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/33.EGAPx注释post_review后的基因组/each_chr/Chr12" || exit 1

# 激活conda环境
source ~/.bashrc
conda activate base

export JAVA_HOME=/share/org/YZWL/yzwl_lixg/miniforge3/envs/EGAPx_Chr12_v.0.4.0-alpha
export PATH=$JAVA_HOME/bin:$PATH
export PATH="/share/org/YZWL/yzwl_lixg/software:$PATH"

# 运行EGAPx_Chr12
python3 \
    ui/egapx.py \
    /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/33.EGAPx注释post_review后的基因组/each_chr/Chr12/Chr12.yaml \
    -e singularity \
    -w /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/33.EGAPx_Chr12注释post_review后的基因组/each_chr/Chr12/work \
    -o /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/33.EGAPx_Chr12注释post_review后的基因组/each_chr/Chr12/output \
    -lc /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/local_cache \
    -r EGAPx_Chr12

# 清理软链接
if [[ -f .egapx_symlinks ]]; then
    xargs -I {} rm -f {} < .egapx_symlinks 2>/dev/null || true
    rm -f .egapx_symlinks
    echo '[INFO] 已清理软链接'
fi

echo '[INFO] Chr12 任务完成'
