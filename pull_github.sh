#!/bin/bash
# =============================================================
# 超算端安全拉取 GitHub|Safe GitHub pull for the HPC side
#
# 用法|Usage:
#   bash pull_github.sh
#   bash pull_github.sh .gitignore docs/xxx.md   # 自定义要恢复的本地独有文件
#   bash pull_github.sh --force                  # 不在 main 分支也继续(慎用)
#
# 背景|Why:
#   超算端只写代码不 commit,本地未提交改动与 Mac 已推送提交内容一致时,
#   git pull --ff-only 会以"本地修改将被覆盖"拒绝。本脚本:
#     1) 展示本地改动与 diff
#     2) git stash push -u 暂存全部改动(含未跟踪文件,不含已忽略文件)
#     3) git pull --ff-only origin main
#     4) 从 stash 只捞回"本地独有"文件(默认 .gitignore):
#        - stash 内容与远端一致 -> 跳过(改动已在提交里)
#        - 与远端不同          -> 存为 <文件>.pre-pull-local 侧车文件,绝不覆盖
#        - 远端没有该文件      -> 原样恢复(本地新增文件)
#     5) 丢弃 stash,打印最终状态
#   pull 失败会自动 git stash pop 恢复现场,不丢任何东西。
#   本脚本只做 pull/stash/restore,不会在超算产生任何 commit。
# =============================================================
set -euo pipefail

FORCE=0
RESTORE_FILES=()
for a in "$@"; do
    case "$a" in
        --force) FORCE=1 ;;
        -h|--help) sed -n '2,21p' "$0"; exit 0 ;;
        *) RESTORE_FILES+=("$a") ;;
    esac
done
[ ${#RESTORE_FILES[@]} -eq 0 ] && RESTORE_FILES=(.gitignore)

# 仓库根目录:脚本位于仓库根,任何位置执行都能定位|Repo root via script location
REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$REPO_DIR"

# 分支检查:三角工作流只用 main|Branch check: triangle workflow uses main only
CURRENT_BRANCH="$(git symbolic-ref --short -q HEAD || true)"
if [ "$CURRENT_BRANCH" != "main" ]; then
    echo "[错误] 当前分支: [$CURRENT_BRANCH],请先 git checkout main(或 --force 继续)|Not on main"
    [ "$FORCE" -eq 1 ] || exit 1
fi

echo "============================================================"
echo "1) 本地未提交改动|Local changes"
echo "============================================================"
git status --short

echo ""
echo "============================================================"
echo "2) 本地独有改动的 diff(恢复清单: ${RESTORE_FILES[*]})|Diff of local-only changes"
echo "============================================================"
git diff -- "${RESTORE_FILES[@]}" || true

# 工作区干净时直接 pull,不需要 stash|Pull directly when the tree is clean
if [ -z "$(git status --porcelain)" ]; then
    echo ""
    echo "工作区干净,直接 pull|Working tree clean, pulling directly"
    git pull --ff-only origin main
    git status -sb
    exit 0
fi

echo ""
echo "============================================================"
echo "3) 暂存全部改动(含未跟踪)|Stash all changes (incl. untracked)"
echo "============================================================"
git stash push -u -m 'pre-pull'
STASH_REF="stash@{0}"
echo "已暂存|Stashed: $STASH_REF"

echo ""
echo "4) 快进拉取|Fast-forward pull"
if ! git pull --ff-only origin main; then
    echo "[失败] pull 失败,恢复现场|Pull failed, restoring stash"
    git stash pop
    exit 1
fi

echo ""
echo "============================================================"
echo "5) 从 stash 恢复本地独有文件|Restore local-only files from stash"
echo "============================================================"
for f in "${RESTORE_FILES[@]}"; do
    # 取 stash 中该文件的内容 blob(先查主提交树,再查未跟踪树)
    stash_blob=""
    if git cat-file -e "$STASH_REF:$f" 2>/dev/null; then
        stash_blob="$(git rev-parse "$STASH_REF:$f")"
    elif git cat-file -e "$STASH_REF^3:$f" 2>/dev/null; then
        stash_blob="$(git rev-parse "$STASH_REF^3:$f")"
    else
        echo "[跳过] $f: stash 中无此改动|not in stash"
        continue
    fi

    head_blob=""
    if git cat-file -e "HEAD:$f" 2>/dev/null; then
        head_blob="$(git rev-parse "HEAD:$f")"
    fi

    if [ -n "$head_blob" ] && [ "$stash_blob" = "$head_blob" ]; then
        echo "[跳过] $f: 本地改动与远端一致,无需恢复|identical to remote"
    elif [ -n "$head_blob" ]; then
        git cat-file blob "$stash_blob" > "$f.pre-pull-local"
        echo "[警告] $f: 本地与远端内容不同,已存为 $f.pre-pull-local,请手动合并"
        echo "       |local differs from remote, saved sidecar file, merge manually"
    else
        git cat-file blob "$stash_blob" > "$f"
        echo "[恢复] $f: 远端无此文件,恢复本地版本|remote lacks it, local restored"
    fi
done

echo ""
echo "============================================================"
echo "6) 以下 stash 内容将被丢弃(不含上一步已恢复的)|Stash contents to be dropped"
echo "============================================================"
git stash show --include-untracked --stat "$STASH_REF" | head -40 || true

git stash drop "$STASH_REF"

echo ""
echo "============================================================"
echo "7) 最终状态|Final status"
echo "============================================================"
git status -sb
git log --oneline -1
echo "完成|Done"
