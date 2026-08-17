#!/usr/bin/env python3
"""
§14 文档合规检查器|§14 doc compliance checker

检查每个 CLI 注册模块的 docs/<mod>.md 是否满足 §14 硬性要求:
  1. 存在自动生成的参数速查标记区块(BEGIN/END PARAMS:auto)
  2. 必备章节齐全: 功能概述/快速开始/输入/输出/依赖/常见问题
  3. 开头标题后应有"一句话人话"段落(非空文本行)
  4. 无 emoji(按常见 emoji 码位范围粗筛)

用法|Usage:
  python scripts/check_docs.py                 # 检查全部,输出违规清单
  python scripts/check_docs.py --only cim fst  # 只检查指定模块
  python scripts/check_docs.py --quiet         # 只输出统计
"""

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs"

REQUIRED_SECTIONS = [
    "## 功能概述",
    "## 快速开始",
    "## 输入",
    "## 输出",
    "## 依赖",
    "## 常见问题",
]

# 常见 emoji 码位粗筛(覆盖笑脸/手势/符号/装饰符,不含 CJK 扩展)
EMOJI_RE = re.compile(
    "[\U0001F000-\U0001FAFF\U00002600-\U000027BF\U0001F900-\U0001F9FF\U00002B00-\U00002BFF\uFE0F]"
)


def load_registry_mods():
    import importlib.util
    spec = importlib.util.spec_from_file_location(
        "g", str(REPO / "scripts" / "gen_docs_params.py"))
    g = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(g)
    return [m for m, _, _ in g.load_registry()]


def check_one(mod, text):
    issues = []
    if "<!-- BEGIN PARAMS:auto -->" not in text or "<!-- END PARAMS:auto -->" not in text:
        issues.append("缺参数速查标记块|missing PARAMS:auto block")
    for sec in REQUIRED_SECTIONS:
        if sec not in text:
            issues.append("缺章节|missing section: %s" % sec)
    lines = text.split(chr(10))
    intro_ok = False
    for ln in lines[:12]:
        s = ln.strip()
        if s and not s.startswith("#"):
            intro_ok = True
            break
    if not intro_ok:
        issues.append("缺开头一句话人话|missing one-line intro")
    if EMOJI_RE.search(text):
        issues.append("含emoji|contains emoji")
    return issues


def main():
    args = sys.argv[1:]
    quiet = "--quiet" in args
    only = set()
    if "--only" in args:
        i = args.index("--only")
        for a in args[i + 1:]:
            if a.startswith("--"):
                break
            only.add(a)
    mods = load_registry_mods()
    if only:
        mods = [m for m in mods if m in only]
    bad = 0
    for mod in mods:
        doc = DOCS / ("%s.md" % mod)
        if not doc.exists():
            if not quiet:
                print("[缺文件|missing] %s" % mod)
            bad += 1
            continue
        issues = check_one(mod, doc.read_text(encoding="utf-8"))
        if issues:
            bad += 1
            if not quiet:
                for it in issues:
                    print("[违规|violation] %s: %s" % (mod, it))
    print("统计|Stats: 检查 %d 个, 违规 %d 个" % (len(mods), bad))
    sys.exit(1 if bad else 0)


if __name__ == "__main__":
    main()
