#!/usr/bin/env python3
"""
CHANGELOG 合规检查器|CHANGELOG compliance checker

§10.1 版本控制规范配套检查——"历史只增不改,新改动开新版本段"。
|Companion check for §10.1 versioning convention — "history only grows;
new changes get a NEW version section; released sections are never edited".

校验(最近 RECENT_N 段 + 最近 RECENT_N 个 tag 严格校验,远古历史降级警告)|
|Checks (recent RECENT_N sections + recent RECENT_N tags enforced strictly;
ancient history downgraded to WARNING):
  1. 版本段按新→旧排序(顶部最新,严格递减)|sections newest-first (strictly decreasing)
  2. 最近段内无重复版本段|no duplicate version headers in the recent range
  3. 最新段版本 == pyproject.toml 的 version|top section version == pyproject version
  4. 近期已发布 tag 的版本段未被修改(与 git show <tag>:CHANGELOG.md 逐字比对)
     |recent released sections unchanged vs their git tags (byte-exact)
  5. 近期 tag 都有对应版本段|every recent tag has a CHANGELOG section

用法|Usage:
  python scripts/check_changelog.py          # 检查全部
  python scripts/check_changelog.py --quiet  # 只输出统计
"""

import re
import subprocess
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
CHANGELOG = REPO / "CHANGELOG.md"
PYPROJECT = REPO / "pyproject.toml"

# 最近多少个版本段 / tag 严格校验|how many most-recent sections/tags are strict
RECENT_N = 20

SECTION_RE = re.compile(r"^## \[(\d+\.\d+\.\d+)\](?: - (\d{4}-\d{2}-\d{2}))?")


def version_key(v: str):
    """'1.2.3' → (1,2,3)|version string to comparable tuple"""
    return tuple(int(x) for x in v.split("."))


def parse_sections(text: str):
    """解析 CHANGELOG → (ordered_headers, sections_dict)
    ordered_headers: 按文件顺序的版本号列表(含重复)
    sections: {version: 段文本} 取最后出现|last occurrence wins
    """
    headers = []
    secs = {}
    cur = None
    for line in text.splitlines():
        m = SECTION_RE.match(line)
        if m:
            v = m.group(1)
            headers.append(v)
            cur = v
            secs[cur] = [line]
        elif cur is not None:
            secs[cur].append(line)
    return headers, {k: chr(10).join(v) for k, v in secs.items()}


def git_show(path: str, gitref: str):
    """git show <gitref>:<path> → str|None on failure"""
    out = subprocess.run(
        ["git", "show", f"{gitref}:{path}"],
        capture_output=True, text=True, cwd=str(REPO))
    if out.returncode != 0:
        return None
    return out.stdout


def main():
    quiet = "--quiet" in sys.argv

    issues = []
    def error(msg):
        issues.append(("ERROR", msg))
    def warn(msg):
        issues.append(("WARN", msg))

    # ---- 解析|parse ----
    text = CHANGELOG.read_text(encoding="utf-8")
    headers, secs = parse_sections(text)
    if not headers:
        error("CHANGELOG 无任何版本段|no version sections found")
        return finish(issues, quiet)

    # 最近 N 段的版本号(文件顺序前 N 个)|recent range = first N in file order
    recent_set = set(headers[:RECENT_N])

    # ---- 1. 排序|ordering ----
    keys = [version_key(v) for v in headers]
    # 最近 N 段必须严格递减|recent must be strictly decreasing
    for i in range(len(headers) - 1):
        if i >= RECENT_N - 1:
            break
        if keys[i] <= keys[i + 1]:
            error(
                f"版本段顺序错误|section ordering broken near top: "
                f"{headers[i]} 应新→旧|must be newer than {headers[i + 1]}")
    # 全文件应非递增(含远古,只警告)|whole file non-increasing (warn for old)
    for i in range(RECENT_N - 1, len(headers) - 1):
        if keys[i] < keys[i + 1]:
            warn(
                f"版本段顺序异常|ordering anomaly (old): {headers[i]} 后跟 "
                f"{headers[i + 1]} (WARN)")

    # ---- 2. 重复|duplicates ----
    seen = set()
    dup_first = {}
    for idx, v in enumerate(headers):
        if v in seen:
            dup_first.setdefault(v, idx)
        seen.add(v)
    for v in sorted(dup_first, key=version_key, reverse=True):
        idx = headers.index(v)   # 首次出现位置|first occurrence position
        if idx < RECENT_N:
            error(f"重复版本段|duplicate section: {v}")
        else:
            warn(f"重复版本段|duplicate section: {v} (旧|old, WARN)")

    # ---- 3. 顶部版本 == pyproject|top == pyproject ----
    py_v = None
    m = re.search(r'^version\s*=\s*"([^"]+)"',
                  PYPROJECT.read_text(encoding="utf-8"), re.MULTILINE)
    if m:
        py_v = m.group(1)
    if py_v and py_v != headers[0]:
        error(
            f"顶部版本段 [{headers[0]}] ≠ pyproject version [{py_v}]|"
            f"top section != pyproject version")

    # ---- 4/5. 已发布 tag 的版本段未被修改|released sections unchanged ----
    tags = sorted(
        (t for t in subprocess.run(
            ["git", "tag", "-l", "v[0-9]*.[0-9]*.[0-9]*"],
            capture_output=True, text=True, cwd=str(REPO)).stdout.split()
         if re.fullmatch(r"v\d+\.\d+\.\d+", t)),
        key=lambda t: version_key(t[1:]))

    recent_tags = tags[-RECENT_N:] if len(tags) > RECENT_N else tags
    old_tags = tags[:-RECENT_N] if len(tags) > RECENT_N else []

    for tag in recent_tags:
        ver = tag[1:]
        if ver not in secs:
            error(f"tag {tag} 无对应版本段|no section for {tag}")
            continue
        tag_text = git_show("CHANGELOG.md", tag)
        if tag_text is None:
            continue
        tag_secs = parse_sections(tag_text)[1]
        if ver not in tag_secs:
            continue
        if tag_secs[ver].strip() != secs[ver].strip():
            error(
                f"已发布版本段被修改|released section modified: [{ver}] "
                f"与 tag {tag} 不一致|differs from tag")

    for tag in old_tags:
        ver = tag[1:]
        if ver not in secs:
            warn(f"tag {tag} 无对应版本段|no section for {tag} (旧|old, WARN)")
            continue
        tag_text = git_show("CHANGELOG.md", tag)
        if tag_text is None:
            continue
        tag_secs = parse_sections(tag_text)[1]
        if ver not in tag_secs:
            continue
        if tag_secs[ver].strip() != secs[ver].strip():
            warn(
                f"已发布版本段被修改|released section modified: [{ver}] "
                f"与 tag {tag} 不一致|differs from tag (旧|old, WARN)")

    return finish(issues, quiet, len(headers), len(tags), len(recent_tags))


def finish(issues, quiet, n_secs=0, n_tags=0, n_recent=0):
    errors = [m for s, m in issues if s == "ERROR"]
    warns = [m for s, m in issues if s == "WARN"]
    if not quiet:
        for m in issues:
            print(m)
    print(f"统计|Stats: 版本段|sections={n_secs}, "
          f"tags={n_tags}(近期|recent={n_recent}), "
          f"错误|errors={len(errors)}, 警告|warnings={len(warns)}")
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())
