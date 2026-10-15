#!/usr/bin/env python3
"""发布前自检：找出「使用了 conda_runner 名称但既未导入也未定义」的模块。

背景
----
`cb73fe28`（全库 conda 实现去重）之后，多个模块把自写的 conda 命令构造换成了
`common.conda_runner` 的公共实现，但**漏了 import 行**：代码只有在真正跑到那一行时
才抛 `NameError`，静态语法检查（`py_compile`）看不出来。2026-10-15 实测即因此让
`biopytools rmvp` 自 1.72.4 起在 conda 环境模式下完全不可用（依赖检查第 1 步报
`name 'conda_env_run_prefix' is not defined`，再误报「rMVP包未安装」），一个多月未被发现。
共有 5 个文件中招（rmvp/utils.py、rmvp/analyzer.py、oomycete_anno/utils.py、
qiime2/utils.py、selective_sweep/utils.py），已在 1.74.1 修复。

用法
----
    python3 scripts/check_conda_imports.py [包根目录，默认 biopytools/]

返回码：0 = 无隐患；1 = 有隐患（会把文件与行号列出来）。
"""

from __future__ import annotations

import ast
import sys
from pathlib import Path

# common/conda_runner 对外提供的、供各模块直接调用的名称
TARGETS = {
    "conda_env_run_prefix",
    "conda_run_prefix",
    "build_conda_command",
    "build_conda_command_string",
    "get_conda_env",
    "preferred_env",
    "_resolve_env_prefix",
}


def analyse(path: Path) -> list[str]:
    """返回「用了但没导入/定义」的名称及行号。"""
    tree = ast.parse(path.read_text(encoding="utf-8"))
    imported: set[str] = set()
    defined: set[str] = set()
    used: dict[str, int] = {}

    empty_args = ast.arguments(posonlyargs=[], args=[], vararg=None,
                               kwonlyargs=[], kw_defaults=[], kwarg=None, defaults=[])
    for node in ast.walk(tree):
        if isinstance(node, (ast.Import, ast.ImportFrom)):
            for alias in node.names:
                imported.add((alias.asname or alias.name).split(".")[0])
        elif isinstance(node, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            defined.add(node.name)
            for arg in getattr(node, "args", empty_args).args:
                defined.add(arg.arg)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Store):
            defined.add(node.id)
        elif isinstance(node, ast.Name) and isinstance(node.ctx, ast.Load) and node.id in TARGETS:
            used[node.id] = node.lineno

    missing = {k: v for k, v in used.items() if k not in imported and k not in defined}
    return [f"{k} (行 {v})" for k, v in sorted(missing.items(), key=lambda kv: kv[1])]


def main() -> int:
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).parent.parent / "biopytools"
    bad = [(f, m) for f in sorted(root.rglob("*.py")) if "__pycache__" not in f.parts
           for m in [analyse(f)] if m]

    if not bad:
        print(f"conda 导入自检通过：{root} 下未发现「用了但没导入」的 conda_runner 名称")
        return 0

    print(f"conda 导入自检失败：{len(bad)} 个文件存在「用了但没导入」的名称")
    for f, miss in bad:
        print(f"  {f}：{'; '.join(miss)}")
    return 1


if __name__ == "__main__":
    sys.exit(main())
