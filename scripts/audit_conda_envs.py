#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
conda环境依赖审计脚本|Conda environment dependency audit script

只读扫描仓库|Read-only scan of the repository, 提取以下引用类型|extracts these reference types:
  1. env-path          : (miniforge3/)?envs/<NAME>/[bin|share]/<TOOL> 路径引用
  2. conda-run-literal : conda run -n <NAME> <TOOL> 字面量引用
  3. conda-activate    : conda activate <NAME>
  4. env-field         : xxx_env: str = "NAME" 配置字段默认值
  5. env-const         : conda_env = "NAME" 常量
  6. conda-run-var     : conda run -n {var} 动态引用(记录变量名)

输出|Output: docs/conda_env_audit.csv 明细 + stdout 汇总
用法|Usage: python scripts/audit_conda_envs.py [repo_root]
"""

import os
import re
import sys
import csv
from collections import Counter, defaultdict

REPO = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else ".")
OUT_CSV = os.path.join(REPO, "docs", "conda_env_audit.csv")

# 占位符环境名(非真实环境)|Placeholder env names (not real envs)
PLACEHOLDERS = {"xxx", "env", "en", "your_env", "tool_env",
                "custom_panman_env", "...", "my_env", "conda_env",
                "x", "bcftools_env", "none", "e"}

# 通用解释器词,不计入工具清单|Generic interpreter words, excluded from tool lists
GENERIC = {"bash", "sh", "which", "python", "python2", "python3",
           "perl", "r", "rscript", "java", "source", "conda"}

PATTERNS = [
    ("env-path",
     re.compile(r"(?:miniforge3/)?envs/([A-Za-z0-9][A-Za-z0-9_.+-]*)"
                r"(?:/(bin|share|libexec|lib)(?:/([A-Za-z0-9][A-Za-z0-9_.+-]*))?)?"),),
    ("conda-run-literal",
     re.compile(r"conda\s+run\s+-n\s+([A-Za-z0-9][A-Za-z0-9_.+-]*)"
                r"(?:\s+--no-capture-output)?(?:\s+([A-Za-z0-9][A-Za-z0-9_.+-]*))?"),),
    ("conda-activate",
     re.compile(r"conda\s+activate\s+([A-Za-z0-9][A-Za-z0-9_.+-]+)"),),
    ("env-field",
     re.compile(r"^\s*[A-Za-z_][A-Za-z0-9_]*_env[A-Za-z0-9_]*\s*:\s*str\s*=\s*['\"]([A-Za-z0-9][A-Za-z0-9_.+-]*)['\"]", re.M),),
    ("env-const",
     re.compile(r"^\s*(?:conda_env|CONDA_ENV)\s*=\s*['\"]([A-Za-z0-9][A-Za-z0-9_.+-]*)['\"]", re.M),),
    ("conda-run-var",
     re.compile(r"conda\s+run\s+-n\s+\{([A-Za-z_][A-Za-z0-9_]*)\}"),),
]

SUFFIXES = {".py", ".md", ".yml", ".yaml", ".sh", ".r", ".txt"}
SKIP_DIRS = {".git", "node_modules", "__pycache__", "build", "dist", ".venv", "venv"}
SKIP_FILES = {"conda_env_audit.csv", "conda_env_audit.md", "audit_conda_envs.py"}


def iter_files():
    for root, dirs, files in os.walk(REPO):
        dirs[:] = [d for d in dirs if d not in SKIP_DIRS and not d.startswith(".")]
        for fn in files:
            if fn in SKIP_FILES:
                continue
            ext = os.path.splitext(fn)[1].lower()
            if ext not in SUFFIXES:
                continue
            yield os.path.join(root, fn)


def classify(filepath):
    """scope|范围: code 或 docs"""
    return "docs" if filepath.endswith(".md") else "code"


def module_of(filepath, scope):
    rel = os.path.relpath(filepath, REPO)
    if scope == "code" and rel.startswith("biopytools"):
        parts = rel[len("biopytools/"):].split(os.sep)
        return parts[0] if parts else "?"
    return os.path.basename(rel).replace(".md", "")


def base_name(env):
    """去版本号+小写归一化|Strip version suffix and lowercase for grouping"""
    low = env.lower()
    m = re.match(r"^(.*?)(?:[_-]?v\.?\d.*)$", low)
    return m.group(1) if m else low


def main():
    refs = []  # (scope, file, module, line, style, env, tool, var, snippet)
    for fpath in iter_files():
        scope = classify(fpath)
        module = module_of(fpath, scope)
        try:
            with open(fpath, encoding="utf-8", errors="replace") as fh:
                lines = fh.readlines()
        except OSError:
            continue
        for lineno, line in enumerate(lines, 1):
            snippet = line.strip()[:160]
            for style, rx in PATTERNS:
                for m in rx.finditer(line):
                    env = m.group(1)
                    tool = None
                    var = None
                    if style == "env-path":
                        tool = m.group(3) if m.lastindex and m.lastindex >= 3 else None
                    elif style in ("conda-run-literal",):
                        tool = m.group(2) if m.lastindex and m.lastindex >= 2 else None
                        if tool and tool.lower() in GENERIC:
                            tool = None
                    elif style == "conda-run-var":
                        var = env
                        env = "<var>"
                    refs.append((scope, fpath, module, lineno, style, env, tool, var, snippet))

    # 去重(同文件同行同env同style)|Deduplicate
    seen = set()
    uniq = []
    for r in refs:
        key = (r[1], r[3], r[4], r[5])
        if key not in seen:
            seen.add(key)
            uniq.append(r)
    refs = uniq

    # 聚合|Aggregate per env
    envs = defaultdict(lambda: {"refs": 0, "files": set(), "modules": set(),
                                "styles": Counter(), "tools": set(),
                                "scope": Counter()})
    for scope, fpath, module, lineno, style, env, tool, var, snippet in refs:
        envs[env]["refs"] += 1
        envs[env]["files"].add(fpath)
        envs[env]["modules"].add(module)
        envs[env]["styles"][style] += 1
        envs[env]["scope"][scope] += 1
        if tool:
            envs[env]["tools"].add(tool)

    # 归一化分组|Canonical groups
    groups = defaultdict(list)
    for env in envs:
        groups[base_name(env)].append(env)
    canon = {}
    for base, variants in groups.items():
        best = sorted(variants,
                      key=lambda v: (envs[v]["scope"].get("code", 0), envs[v]["refs"], v),
                      reverse=True)[0]
        for v in variants:
            canon[v] = best

    # 写CSV|Write CSV
    os.makedirs(os.path.dirname(OUT_CSV), exist_ok=True)
    with open(OUT_CSV, "w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["canonical_env", "env_name", "is_alias", "scope", "module",
                    "file", "line", "style", "tool", "var", "snippet"])
        for scope, fpath, module, lineno, style, env, tool, var, snippet in sorted(refs):
            c = canon[env]
            w.writerow([c, env, "yes" if c != env else "", scope, module,
                        os.path.relpath(fpath, REPO), lineno, style, tool or "", var or "",
                        snippet])

    # 汇总|Summary
    placeholders = [e for e in envs if e.lower() in PLACEHOLDERS or e == "<var>"]
    real = [e for e in envs if e.lower() not in PLACEHOLDERS and e != "<var>"]
    canon_real = set(canon[e] for e in real)
    canon_code = set(canon[e] for e in real if envs[e]["scope"].get("code", 0) > 0)
    var_refs = [r for r in refs if r[4] == "conda-run-var"]

    style_total = Counter()
    for e in real:
        style_total.update(envs[e]["styles"])

    print("=== 统计总览|Summary ===")
    print(f"引用点总数|Total refs: {len(refs)}")
    print(f"环境名字符串|Raw env names: {len(envs)}")
    print(f"占位符|Placeholders: {len(placeholders)}")
    print(f"真实环境名|Real env names: {len(real)}")
    print(f"归一化后独立环境|Canonical envs: {len(canon_real)} (代码引用 {len(canon_code)})")
    print(f"仅文档提及|Docs-only canonical envs: {len(canon_real - canon_code)}")
    print(f"动态conda-run变量引用|Dynamic conda-run-var refs: {len(var_refs)}")
    print(f"涉及文件|Files: {len(set(r[1] for r in refs))}")
    print(f"涉及模块|Modules: {len(set(r[2] for r in refs))}")
    print("引用方式分布|Style breakdown:")
    for k, v in style_total.most_common():
        print(f"  {k}: {v}")

    print()
    print("=== 独立环境清单(按代码引用数排序)|Canonical envs ranked ===")
    ranked = sorted(canon_real, key=lambda c: (
        sum(envs[v]["scope"].get("code", 0) for v in groups[base_name(c)]),
        sum(envs[v]["refs"] for v in groups[base_name(c)])), reverse=True)
    for i, c in enumerate(ranked, 1):
        variants = groups[base_name(c)]
        aliases = [v for v in variants if v != c]
        code_n = sum(envs[v]["scope"].get("code", 0) for v in variants)
        docs_n = sum(envs[v]["scope"].get("docs", 0) for v in variants)
        files_n = len(set().union(*[envs[v]["files"] for v in variants]))
        mods_n = len(set().union(*[envs[v]["modules"] for v in variants]))
        tools = sorted(set().union(*[envs[v]["tools"] for v in variants]))[:8]
        styles = Counter()
        for v in variants:
            styles.update(envs[v]["styles"])
        style_s = ",".join(f"{k}:{n}" for k, n in styles.most_common(3))
        alias_s = (" |别名|alias:" + ",".join(aliases)) if aliases else ""
        tool_s = ",".join(tools) if tools else "-"
        print(f"{i:3d}. {c}  代码:{code_n} 文档:{docs_n} 文件:{files_n} 模块:{mods_n} "
              f"[{style_s}] 工具:{tool_s}{alias_s}")

    print()
    print("=== 同工具多环境(冲突信号)|Tool in multiple envs ===")
    tool_envs = defaultdict(set)
    for e in real:
        for t in envs[e]["tools"]:
            tool_envs[t].add(canon[e])
    for t, es in sorted(tool_envs.items(), key=lambda kv: -len(kv[1])):
        if len(es) >= 2:
            print(f"  {t}: {sorted(es)}")

    print()
    print("=== 版本钉死环境名|Version-pinned env names ===")
    pinned = [c for c in canon_real if re.search(r"v\.?\d", c)]
    print("  共|total:", len(pinned), "个")
    print("  ", ", ".join(sorted(pinned)))

    print()
    print("=== 占位符/待修复引用|Placeholder refs ===")
    for scope, fpath, module, lineno, style, env, tool, var, snippet in refs:
        if env.lower() in PLACEHOLDERS or env == "<var>":
            print(f"  {os.path.relpath(fpath, REPO)}:{lineno} [{style}] "
                  f"{('env=' + env) if env != '<var>' else ('var=' + var)} :: {snippet[:80]}")

    print()
    print("=== 模块依赖环境数(代码范围,前20)|Modules by env count ===")
    mod_envs = defaultdict(set)
    for e in real:
        for fpath in envs[e]["files"]:
            if classify(fpath) == "code":
                mod_envs[module_of(fpath, "code")].add(canon[e])
    for mod, es in sorted(mod_envs.items(), key=lambda kv: -len(kv[1]))[:20]:
        print(f"  {mod}: {len(es)} 个环境")

    print()
    print(f"CSV已写入|CSV written: {OUT_CSV}")


if __name__ == "__main__":
    main()
