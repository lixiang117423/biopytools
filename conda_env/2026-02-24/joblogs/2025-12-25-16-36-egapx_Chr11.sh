#!/bin/bash
# 激活conda环境
source ~/.bashrc
conda activate base

export JAVA_HOME=/share/org/YZWL/yzwl_lixg/miniforge3/envs/EGAPx_v.0.4.0-alpha
export PATH=$JAVA_HOME/bin:$PATH
export PATH="/share/org/YZWL/yzwl_lixg/software:$PATH"

# 运行EGAPx
python3 \
    ui/egapx.py \
    /share/org/YZWL/yzwl_lixg/tmp/test_egapx/each_chr/Chr11/Chr11.yaml \
    -e singularity \
    -w /share/org/YZWL/yzwl_lixg/tmp/test_egapx/each_chr/Chr11/work \
    -o /share/org/YZWL/yzwl_lixg/tmp/test_egapx/each_chr/Chr11/output \
    -lc /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/local_cache \
    -r EGAPx_Chr11
