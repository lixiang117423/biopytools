#!/bin/bash
#================================================
# 自动生成的作业脚本
# 原始命令文件: glnexus_each_chr_singularity.sh
# 原始行号: 12
# 生成时间: 2025-12-09 11:26:29
#================================================

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任何一个失败都返回失败

echo "🚀 作业开始: glnexus_12"
echo "⏰ 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "📍 执行命令: ~/software/scripts/50.运行Glnexus批量合并gVCF文件_singularity版本.sh -r ../../01.data/genome/genome.fa -i ../vcf -o ./ -c Chr12 --keep --no-validate"
echo "=========================================="

# 执行原始命令
~/software/scripts/50.运行Glnexus批量合并gVCF文件_singularity版本.sh -r ../../01.data/genome/genome.fa -i ../vcf -o ./ -c Chr12 --keep --no-validate

EXIT_CODE=$?

echo "=========================================="
echo "⏰ 结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✅ 作业成功: glnexus_12"
else
    echo "❌ 作业失败: glnexus_12 (退出码: $EXIT_CODE)"
fi

exit $EXIT_CODE
