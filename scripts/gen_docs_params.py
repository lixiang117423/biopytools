#!/usr/bin/env python3
"""
模块文档参数速查表自动生成器|Auto-generate parameter reference tables for module docs

依据 §14.4:参数表从 CLI 定义自动提取,禁止手写。
数据源|Sources(均为 AST 静态解析,不 import 模块,避免重依赖):
  1. biopytools/cli/commands/<mod>.py 的 @click.option 装饰器(主入口参数)
  2. biopytools/<mod>/main.py 的 argparse add_argument(模块直调补充参数)

行为|Behavior:
  - 在 docs/<module>.md 中维护一段带标记的区块:
      <!-- BEGIN PARAMS:auto --> ... <!-- END PARAMS:auto -->
  - 已有标记 → 原地替换;无标记 → 在「## 依赖/## 常见问题/## 输出」前插入,否则文末追加
  - 文档不存在 → 生成 §14 骨架 + 参数速查表
  - 人工手写的分组参数说明完全不受影响(生成内容只追加/替换自身标记区块)

用法|Usage:
  python scripts/gen_docs_params.py --dry-run              # 全量试运行,只报统计
  python scripts/gen_docs_params.py --apply                # 全量写入
  python scripts/gen_docs_params.py --module cim --apply   # 只处理指定模块(可多个)
"""

import ast
import re
import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
DOCS = REPO / "docs"
CLI_DIR = REPO / "biopytools" / "cli" / "commands"
REGISTRY_PY = REPO / "biopytools" / "cli" / "main.py"

BEGIN = "<!-- BEGIN PARAMS:auto -->"
END = "<!-- END PARAMS:auto -->"
BLOCK_RE = re.compile(
    re.escape(BEGIN) + r".*?" + re.escape(END), re.DOTALL
)

# 插入锚点(按优先级):命中第一个存在的标题,在其前插入
INSERT_ANCHORS = ["## 依赖", "## 常见问题", "## 输出", "## 结果解读"]


def literal(node):
    """安全求值常量表达式(含 f-string 常量、拼接)|Eval constant expressions"""
    try:
        return ast.literal_eval(node)
    except Exception:
        pass
    if isinstance(node, ast.JoinedStr):
        parts = []
        for v in node.values:
            if isinstance(v, ast.Constant):
                parts.append(str(v.value))
            elif isinstance(v, ast.FormattedValue):
                parts.append("…")
        return "".join(parts)
    if isinstance(node, ast.BinOp):
        l, r = literal(node.left), literal(node.right)
        if l is not None and r is not None:
            return "%s%s" % (l, r)
    if isinstance(node, ast.UnaryOp):
        v = literal(node.operand)
        if isinstance(v, (int, float)):
            return -v if isinstance(node.op, ast.USub) else v
    if isinstance(node, ast.Name) and node.id in ("True", "False"):
        return node.id == "True"
    return None


def type_of(node):
    """从 click type=/argparse type= 节点提取类型名|Extract type name"""
    if node is None:
        return ""
    if isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute):
        if node.func.attr == "Choice" and node.args:
            choices = literal(node.args[0])
            if isinstance(choices, (list, tuple)):
                return "/".join(str(c) for c in choices)
        return node.func.attr
    if isinstance(node, ast.Attribute):
        return node.attr
    if isinstance(node, ast.Name):
        return node.id
    return ""


def parse_click_option(call):
    """解析单个 @click.option 调用|Parse one @click.option call"""
    flags, dest = [], None
    kw = {}
    for a in call.args:
        if isinstance(a, ast.Constant) and isinstance(a.value, str):
            if a.value.startswith("-"):
                flags.append(a.value)
            else:
                dest = a.value
    for k in call.keywords:
        if k.arg in ("default", "help", "required", "type", "is_flag",
                     "show_default", "flag_value", "multiple", "count"):
            kw[k.arg] = k.value
    if not flags:
        return None
    return {"flags": flags, "dest": dest, "kw": kw}


def extract_click_options(path):
    """提取文件内全部 click option|Extract click options from file"""
    src = Path(path).read_text(encoding="utf-8")
    tree = ast.parse(src)
    out, seen = [], set()
    for node in ast.walk(tree):
        if not isinstance(node, ast.FunctionDef):
            continue
        for dec in node.decorator_list:
            if not (isinstance(dec, ast.Call) and isinstance(dec.func, ast.Attribute)
                    and dec.func.attr == "option"):
                continue
            opt = parse_click_option(dec)
            if opt is None:
                continue
            key = "|".join(opt["flags"])
            if key in seen:
                continue
            seen.add(key)
            out.append(opt)
    return out


def extract_argparse_options(path):
    """提取文件内 argparse add_argument(直调参数)|Extract argparse options"""
    p = Path(path)
    if not p.exists():
        return []
    tree = ast.parse(p.read_text(encoding="utf-8"))
    out, seen = [], set()
    for node in ast.walk(tree):
        if not (isinstance(node, ast.Call) and isinstance(node.func, ast.Attribute)
                and node.func.attr == "add_argument"):
            continue
        flags, kw = [], {}
        for a in node.args:
            v = literal(a)
            if isinstance(v, str):
                flags.append(v)
        for k in node.keywords:
            if k.arg in ("default", "help", "required", "type", "action", "choices"):
                kw[k.arg] = k.value
        key = "|".join(f for f in flags if f.startswith("-"))
        if key in seen or not key:
            continue
        seen.add(key)
        out.append({"flags": flags, "kw": kw})
    return out


def sanitize(text):
    """清理说明文本:表格内 | 会破坏解析,换成全角|Sanitize help text for table cells"""
    if not text:
        return ""
    text = text.replace("|", "｜").replace("\n", " ")
    return re.sub(r"\s+", " ", text).strip()


def row_of(flags, kw, is_argparse=False):
    """生成表格行|Build a table row"""
    params = "`" + ", ".join(f for f in flags if f.startswith("-")) + "`"
    required = bool(kw.get("required") and literal(kw["required"]))
    dflt = kw.get("default")
    default_val = literal(dflt) if dflt is not None else None
    is_flag = bool(literal(kw.get("is_flag"))) or kw.get("action") is not None
    if is_argparse and isinstance(kw.get("action"), ast.Constant):
        is_flag = kw["action"].value in ("store_true", "store_false", "version")
    if required:
        dflt_str = "必填"
    elif default_val is not None:
        if default_val is True:
            dflt_str = "`True`"
        elif default_val is False:
            dflt_str = "`False`"
        else:
            dflt_str = "`%s`" % default_val
    elif is_flag:
        dflt_str = "—"
    else:
        dflt_str = "—"
    type_str = type_of(kw.get("type"))
    if is_argparse and kw.get("choices") is not None:
        ch = literal(kw["choices"])
        if isinstance(ch, (list, tuple)):
            type_str = "/".join(str(c) for c in ch)
    if is_argparse and isinstance(kw.get("action"), ast.Constant):
        if type_str:
            type_str = type_str + ":" + kw["action"].value
        else:
            type_str = kw["action"].value
    help_txt = ""
    if kw.get("help") is not None:
        help_txt = sanitize(literal(kw["help"]) or "")
    return "| %s | %s | %s | %s |" % (params, dflt_str, type_str, help_txt)


def build_block(mod_name, cmd_name, description):
    """构建参数速查区块|Build the parameter reference block"""
    click_path = CLI_DIR / ("%s.py" % mod_name)
    click_opts = extract_click_options(click_path) if click_path.exists() else []
    arg_path = REPO / "biopytools" / mod_name / "main.py"
    arg_opts = extract_argparse_options(arg_path) if arg_path.exists() else []

    lines = [
        BEGIN, "",
        "## 参数速查 | Parameter reference", "",
        "> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改"
        "|Auto-generated from CLI definitions; do not edit by hand", "",
    ]
    if click_opts:
        lines += [
            "### 命令行参数 | CLI options", "",
            "| 参数 | 默认值 | 类型 | 说明 |",
            "|------|--------|------|------|",
        ]
        for o in click_opts:
            lines.append(row_of(o["flags"], o["kw"]))
        lines.append("")
    if arg_opts:
        lines += [
            "### 模块直调参数 | Direct invocation options", "",
            "| 参数 | 默认值 | 类型 | 说明 |",
            "|------|--------|------|------|",
        ]
        for o in arg_opts:
            lines.append(row_of(o["flags"], o["kw"], is_argparse=True))
        lines.append("")
    if not click_opts and not arg_opts:
        lines += ["_未找到 CLI 参数定义|No CLI definitions found_", ""]
    lines.append(END)
    return "\n".join(lines)


def insert_block(text, block):
    """有标记则替换,否则按锚点插入,再否则文末追加|Replace or insert block"""
    if BEGIN in text:
        return BLOCK_RE.sub(lambda m: block, text)
    for anchor in INSERT_ANCHORS:
        m = re.search("^" + re.escape(anchor) + r"\b.*$", text, re.MULTILINE)
        if m:
            return text[: m.start()] + block + "\n\n" + text[m.start():]
    return text.rstrip("\n") + "\n\n" + block + "\n"


def process(entries, apply, only):
    """主流程|Main loop"""
    created = updated = skipped = 0
    for mod_name, cmd_name, desc in entries:
        if only and mod_name not in only:
            continue
        doc = DOCS / ("%s.md" % mod_name)
        block = build_block(mod_name, cmd_name, desc)
        if doc.exists():
            text = doc.read_text(encoding="utf-8")
            new_text = insert_block(text, block)
            if new_text != text:
                updated += 1
                if apply:
                    doc.write_text(new_text, encoding="utf-8")
            else:
                skipped += 1
        else:
            created += 1
            if apply:
                skeleton = (
                    "# %s | %s\n\n" % (cmd_name, desc)
                    + "> 本文档骨架由 `scripts/gen_docs_params.py` 生成,待按 §14 模板补全"
                    + "|Skeleton generated by scripts/gen_docs_params.py; to be completed per §14\n\n"
                )
                doc.write_text(skeleton + block + "\n", encoding="utf-8")
    print("统计|Stats: 新建 %d, 更新 %d, 无变化 %d" % (created, updated, skipped))
    print("模式|Mode: %s" % ("已写入|APPLIED" if apply else "试运行|DRY-RUN"))


def load_registry():
    """解析 COMMAND_REGISTRY|Parse COMMAND_REGISTRY"""
    tree = ast.parse(REGISTRY_PY.read_text(encoding="utf-8"))
    for node in ast.walk(tree):
        if (isinstance(node, ast.Assign) and
                any(isinstance(t, ast.Name) and t.id == "COMMAND_REGISTRY"
                    for t in node.targets)):
            entries = []
            for elt in node.value.elts:
                if isinstance(elt, ast.Tuple) and len(elt.elts) >= 3:
                    entries.append(tuple(literal(e) for e in elt.elts[:3]))
            return entries
    return []


def main():
    args = sys.argv[1:]
    apply = "--apply" in args
    dry = "--dry-run" in args
    only = set()
    if "--module" in args:
        i = args.index("--module")
        for a in args[i + 1:]:
            if a.startswith("--"):
                break
            only.add(a)
    if not apply and not dry:
        print("请指定 --dry-run 或 --apply|Specify --dry-run or --apply")
        sys.exit(1)
    entries = load_registry()
    print("注册命令|Registry entries: %d" % len(entries))
    process(entries, apply, only)


if __name__ == "__main__":
    main()
