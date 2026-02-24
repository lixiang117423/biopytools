
# 自动生成时间 | Auto generated time: 2025-12-22 15:43:39
# cd "./each_chr/Chr12" || exit 1

# 激活conda环境
source ~/.bashrc
conda activate base

export JAVA_HOME=/share/org/YZWL/yzwl_lixg/miniforge3/envs/EGAPx_v.0.4.0-alpha
export PATH=$JAVA_HOME/bin:$PATH
export PATH="/share/org/YZWL/yzwl_lixg/software:$PATH"

# 运行EGAPx
python3 \
    ui/egapx.py \
    ./Chr12.yaml \
    -e singularity \
    -w ./work \
    -o ./output \
    -lc /share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/local_cache \
    -r EGAPx_Chr12


# 清理软链接 | Cleanup symlinks
if [[ -f .egapx_symlinks ]]; then
    xargs -I {} rm -f {} < .egapx_symlinks 2>/dev/null || true
    rm -f .egapx_symlinks
    echo '[INFO] 已清理软链接 | Symlinks cleaned up'
fi

echo '[INFO] Chr12 任务完成 | Chr12 task completed'
