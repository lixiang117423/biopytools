#!/bin/bash
#================================================
# 自动生成的作业脚本
# 原始命令文件: all_jobs_submit.list.sh
# 原始行号: 3
# 生成时间: 2025-12-11 20:02:31
#================================================

set -e  # 遇到错误立即退出
set -u  # 使用未定义变量时报错
set -o pipefail  # 管道命令中任何一个失败都返回失败

echo "🚀 作业开始: EGAPx_3"
echo "⏰ 开始时间: $(date '+%Y-%m-%d %H:%M:%S')"
echo "📍 执行命令: /share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/26.测试批量EGAPx流程/each_chr/Chr03/egapx_Chr03.sh"
echo "=========================================="

# 执行原始命令
/share/org/YZWL/yzwl_lixg/project/06.longliuxing_BSA/68.三代数据组装和注释/26.测试批量EGAPx流程/each_chr/Chr03/egapx_Chr03.sh

EXIT_CODE=$?

echo "=========================================="
echo "⏰ 结束时间: $(date '+%Y-%m-%d %H:%M:%S')"
if [[ $EXIT_CODE -eq 0 ]]; then
    echo "✅ 作业成功: EGAPx_3"
else
    echo "❌ 作业失败: EGAPx_3 (退出码: $EXIT_CODE)"
fi

exit $EXIT_CODE
