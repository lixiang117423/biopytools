#!/usr/bin/env python3
"""
从 module_groups.md 生成 mkdocs.yml 的 nav 段|Generate mkdocs.yml nav from module_groups.md

分组表维护在 docs/module_groups.md(用户确认后的权威分组),
本脚本按组生成导航并写回 mkdocs.yml 的 nav 段,避免手改 200+ 行 YAML。

用法|Usage:
  python scripts/gen_docs_nav.py            # 重写 mkdocs.yml 的 nav 段
  python scripts/gen_docs_nav.py --dry-run  # 只打印生成的 nav,不写文件
"""

import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
GROUPS_MD = REPO / "docs" / "module_groups.md"
MK = REPO / "mkdocs.yml"

GROUP_RE = re.compile(r"^## (.+?)\((\d+)\)\|\d+ modules$")
ROW_RE = re.compile(r"^\| (\S+) \| (\S+) \|")


def parse_groups():
    groups = []
    current = None
    for line in GROUPS_MD.read_text(encoding="utf-8").splitlines():
        m = GROUP_RE.match(line)
        if m:
            current = {"title": m.group(1).split("|")[0].strip(),
                       "mods": []}
            groups.append(current)
            continue
        if current is None:
            continue
        r = ROW_RE.match(line)
        if r and r.group(1) not in ("模块", "------", "总覆盖|Total:"):
            current["mods"].append((r.group(1), r.group(2)))
    return groups


def nav_lines(groups):
    lines = ["nav:", "  - 首页: index.md", "  - 模块分组草案: module_groups.md", "  - 模块文档:"]
    for g in groups:
        lines.append("      - %s:" % g["title"])
        for mod, cmd in g["mods"]:
            lines.append("          - %s: %s.md" % (cmd, mod))
    return lines


def apply(nav_text):
    text = MK.read_text(encoding="utf-8")
    new_text = re.sub(
        r"# 导航:模块按功能域分组.*?\nnav:\n(?:.*?\n)*(?=\n#|\Z)",
        lambda m: nav_text + "\n",
        text, count=1, flags=re.DOTALL,
    )
    if new_text == text:
        # 没有旧注释段则直接替换 nav 块
        new_text = re.sub(
            r"^nav:\n(?:  - .*\n|      - .*\n|          - .*\n)*",
            nav_text + "\n", text, count=1, flags=re.MULTILINE,
        )
    return new_text


def main():
    dry = "--dry-run" in sys.argv
    groups = parse_groups()
    nav = "\n".join(nav_lines(groups))
    if dry:
        print(nav)
        return
    MK.write_text(apply(nav), encoding="utf-8")
    print("已写入 mkdocs.yml nav|nav written: %d 组, %d 模块" %
          (len(groups), sum(len(g[1]) for g in [(g["title"], g["mods"]) for g in groups])))


if __name__ == "__main__":
    main()
