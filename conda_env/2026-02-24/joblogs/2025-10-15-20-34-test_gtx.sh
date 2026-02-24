#!/bin/bash

# 🧬 完整的BSA分析流程脚本
# 包含：质量控制 → 单样品变异检测 → 分组合并BSA分析
# 作者: [Your Name]
# 日期: $(date)

# set -e  # 遇到错误立即退出
# set -u  # 使用未定义变量时报错

# 加载GTX环境
echo "🔧 加载GTX环境..."
source ~/.bashrc
module load gtx/2.2.1

/share/software/GTX.CAT_2.2.1/bin/gtx index /share/org/YZWL/yzwl_lixg/project/16.荠菜/15.each_8_samples/genome.fa