#!/bin/bash
#================================================
# 自动生成的作业脚本
# 原始命令文件: 12.fastq文件质控.sh
# 原始行号: 107
# 生成时间: 2026-01-20 11:59:55
#================================================

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任何一个失败都返回失败

echo "🚀 作业开始: fastp_107"
echo "⏰ 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "📍 执行命令: biopytools fastp -i 01.data/raw/hifi/fastq/N11-7.hifi_reads.fastq.gz -o 01.data/clean/hifi --read1-suffix .hifi_reads.fastq.gz --single-end"
echo "=========================================="

# 执行原始命令
biopytools fastp -i 01.data/raw/hifi/fastq/N11-7.hifi_reads.fastq.gz -o 01.data/clean/hifi --read1-suffix .hifi_reads.fastq.gz --single-end

EXIT_CODE=$?

echo "=========================================="
echo "⏰ 结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✅ 作业成功: fastp_107"
else
    echo "❌ 作业失败: fastp_107 (退出码: $EXIT_CODE)"
fi

exit $EXIT_CODE
