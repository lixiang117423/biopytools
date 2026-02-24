#!/bin/bash
#================================================
# 自动生成的作业脚本
# 原始命令文件: 10.第二批bam文件转fastq.sh
# 原始行号: 11
# 生成时间: 2026-01-20 08:19:31
#================================================

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任何一个失败都返回失败

echo "🚀 作业开始: fastq2bam_11"
echo "⏰ 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "📍 执行命令: biopytools bam2fastq -i 01.data/raw/hifi/bam2/N3-6.hifi_reads.bam -o 01.data/raw/hifi/fastq -t 64"
echo "=========================================="

# 执行原始命令
biopytools bam2fastq -i 01.data/raw/hifi/bam2/N3-6.hifi_reads.bam -o 01.data/raw/hifi/fastq -t 64

EXIT_CODE=$?

echo "=========================================="
echo "⏰ 结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✅ 作业成功: fastq2bam_11"
else
    echo "❌ 作业失败: fastq2bam_11 (退出码: $EXIT_CODE)"
fi

exit $EXIT_CODE
