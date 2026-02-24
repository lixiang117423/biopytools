# #!/bin/bash

# 激活conda环境
source ~/.bashrc
conda activate base

export JAVA_HOME=/share/org/YZWL/yzwl_lixg/miniforge3/envs/EGAPx_v.0.4.0-alpha
export PATH=$JAVA_HOME/bin:$PATH
export PATH="/share/org/YZWL/yzwl_lixg/software:$PATH"

# 运行EGAPx
python3 \
    ui/egapx.py \
    /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/07.ncbi_egapx/YaHS/OV53_1.yaml \
    -e singularity \
    -w /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/07.ncbi_egapx/YaHS/work \
    -o /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/07.ncbi_egapx/YaHS/output \
    -lc /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/local_cache \
    -r zhugecai