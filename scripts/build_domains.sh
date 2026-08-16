#!/bin/bash
# =============================================================
# 域环境批量创建脚本|Domain env batch creation (在超算执行|run on HPC)
# 用法|Usage:
#   bash scripts/build_domains.sh              # 创建全部 14 个域环境
#   bash scripts/build_domains.sh align        # 只创建指定域
#   bash scripts/build_domains.sh align pop    # 创建多个域
# 说明|Notes:
#   - 已存在的环境自动跳过|Existing envs are skipped
#   - 不会删除/修改任何已有环境|Never touches existing envs
#   - mamba 路径可用环境变量覆盖|Override mamba path via env var:
#     MAMBA=/path/to/mamba bash scripts/build_domains.sh
# =============================================================
set -u

MAMBA="${MAMBA:-/share/org/YZWL/yzwl_lixg/miniforge3/bin/mamba}"
ENV_DIR="$(cd "$(dirname "$0")/../envs" && pwd)"
FILTER="$*"

echo "mamba: $MAMBA"
echo "envs目录|envs dir: $ENV_DIR"
echo ""

ok_n=0; fail_n=0; skip_n=0
for yml in "$ENV_DIR"/*.yml; do
  domain="$(basename "$yml" .yml)"
  if [ -n "$FILTER" ]; then
    echo " $FILTER " | grep -q " $domain " || continue
  fi
  if "$MAMBA" env list | awk '{print $1}' | grep -qx "$domain"; then
    echo "[SKIP] $domain 已存在|exists"
    skip_n=$((skip_n+1))
    continue
  fi
  echo "========================================"
  echo "[CREATE] $domain"
  if "$MAMBA" env create -f "$yml"; then
    echo "[OK] $domain"
    ok_n=$((ok_n+1))
  else
    echo "[FAIL] $domain 创建失败|creation failed"
    fail_n=$((fail_n+1))
  fi
  echo ""
done

echo "========================================"
echo "完成|Done: 成功|OK=$ok_n 失败|FAIL=$fail_n 跳过|SKIP=$skip_n"
