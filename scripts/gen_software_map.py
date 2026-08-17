#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""生成 软件→conda环境 速查表|Generate software→env quick reference"""
import os, re, sys
sys.path.insert(0, ".")
from biopytools.common.env_map import TOOL_DOMAIN_MAP

DOMAINS = ["align","pop","asm","hic","annot","repeat","rna","protein","phylo","pan","viz","misc","r","busco"]
DOMAIN_CN = {
    "align": "比对与变异核心|Alignment & variant calling",
    "pop": "群体遗传|Population genetics",
    "asm": "基因组组装|Genome assembly",
    "hic": "Hi-C 三维|Hi-C & 3D genome",
    "annot": "注释与功能预测|Annotation",
    "repeat": "重复序列|Repeats",
    "rna": "转录组|RNA-seq",
    "protein": "蛋白分析|Protein",
    "phylo": "系统发育|Phylogenetics",
    "pan": "泛基因组|Pan-genome",
    "viz": "可视化|Visualization",
    "misc": "杂项工具|Misc utilities",
    "r": "R 生态|R ecosystem",
    "busco": "BUSCO 评估|BUSCO assessment",
}

NOISE = {"lib","_libgcc","_openmp","_x86_64","_r-mutex","alsa","brotli","bzip2",
    "c-ares","ca-certificates","cairo","certifi","charset-normalizer","click","colorama",
    "contourpy","cycler","cyrus-sasl","dbus","decorator","double-conversion","exceptiongroup",
    "expat","filelock","font-ttf","fontconfig","fonts-conda","freetype","fribidi","fsspec",
    "gdk-pixbuf","gettext","giflib","glib","graphite2","gst-plugins","gstreamer","gtk2","harfbuzz",
    "icu","idna","importlib-metadata","iniconfig","jpeg","kiwisolver","krb5","lcms2","lerc",
    "libdeflate","libffi","libgcc","libglib","libiconv","libidn2","libjpeg-turbo","libllvm",
    "libnghttp2","libpng","libssh2","libstdcxx","libtiff","libtool","libunistring","libuuid",
    "libwebp","libxcb","libxml2","libxslt","libzlib","lz4-c","lzo","markupsafe","munkres",
    "ncurses","nspr","nss","openjpeg","openssl","p11-kit","packaging","pango","pcre2","pillow",
    "pip","pixman","platformdirs","pluggy","poppler","psutil","pyparsing","pytest","python",
    "python-dateutil","python-tzdata","python_abi","pytz","pyyaml","readline","requests",
    "setuptools","six","snappy","sqlite","tk","tomli","tornado","tqdm","tzdata","urllib3",
    "wheel","xz","yaml","zipp","zlib","zstd","brotli-bin","bzip2","ca-certificates",
    "libcups","libedit","libev","libexpat","libnsl","libsqlite","libzlib","lz4-c"}

protect = [l.strip() for l in open("scripts/protect_list.txt") if l.strip()]
delete = [l.strip() for l in open("scripts/delete_list.txt") if l.strip()]
legacy = [e for e in protect if e not in DOMAINS and e not in ("base","biopytools")]

def parse_yaml(name):
    deps = set()
    fn = os.path.join("conda_envs_backup", name + ".yaml")
    if not os.path.exists(fn): return deps
    in_deps = False
    for line in open(fn, encoding="utf-8", errors="replace"):
        line = line.rstrip()
        if line == "dependencies:": in_deps = True; continue
        if in_deps and line.startswith("- "):
            item = line[2:].strip()
            if item.startswith("pip:") or item.startswith("pip =="): continue
            m = re.match(r"^([a-zA-Z0-9_.+-]+)=?", item)
            if m: deps.add(m.group(1).lower())
    return deps

def signature(name):
    tools = []
    for pkg in sorted(parse_yaml(name)):
        if pkg.startswith(NOISE): continue
        if pkg.startswith(("perl-","r-","bioconductor-","xorg-","gcc","gxx","binutils","gfortran","kernel-headers","ld_impl","libgfortran","libgomp","libcblas","libblas","liblapack","libopenblas","_openmp")): continue
        if re.match(r"^py[0-9]", pkg): continue
        tools.append(pkg)
        if len(tools) >= 12: break
    return tools

L = []
L.append("# Conda 环境软件速查表|Conda Environment Software Quick Reference")
L.append("")
L.append("> 用途|Purpose: 超算上的 AI/开发者找现成软件时, 按本表选择调用环境。")
L.append("> 更新|Updated: 2026-08-16")
L.append("")
L.append("## 使用规则|Usage Rules")
L.append("")
L.append("1. **优先**使用 14 个功能域环境(第一部分)|Prefer the 14 domain envs (Part 1)")
L.append("2. 域环境没有的软件, 查保留独立环境(第二部分)|Fall back to kept standalone envs (Part 2)")
L.append("3. **禁止**在 scripts/delete_list.txt(" + str(len(delete)) + " 个待退役环境)中找软件, 也不要让新模块依赖它们")
L.append("4. 新模块引入新软件时, 优先并入现有域环境(yaml 配方在 envs/*.yml), 不要新建环境")
L.append("5. 调用方式: conda run -n <env> <tool> --no-capture-output")
L.append("")
L.append("---")
L.append("")
L.append("## 第一部分: 14 个功能域环境|Part 1: The 14 Domain Envs")
for d in DOMAINS:
    tools = sorted({t for t, dm in TOOL_DOMAIN_MAP.items() if dm == d})
    L.append("")
    L.append("### " + d + " — " + DOMAIN_CN[d])
    L.append("")
    L.append("```bash")
    L.append("conda run -n " + d + " <工具> --no-capture-output")
    L.append("```")
    L.append("")
    L.append("| 工具|Tool | 工具|Tool | 工具|Tool |")
    L.append("|---|---|---|---|")
    for i in range(0, len(tools), 3):
        row = tools[i:i+3]
        L.append("| " + " | ".join("`" + t + "`" for t in row) + " |" + "  |" * (3 - len(row)))
L.append("")
L.append("---")
L.append("")
L.append("## 第二部分: 保留的独立环境(强依赖)|Part 2: Kept Standalone Envs")
L.append("")
L.append("| 环境|Env | 关键软件|Key software |")
L.append("|---|---|---|")
for e in sorted(legacy):
    sig = signature(e)
    L.append("| " + e + " | " + (", ".join("`" + t + "`" for t in sig) if sig else "(见 backup yaml)") + " |")
L.append("")
L.append("---")
L.append("")
L.append("## 第三部分: 项目运行环境|Part 3: Project Runtime Env")
L.append("")
L.append("| 环境|Env | 内容|Content |")
L.append("|---|---|")
L.append("| biopytools | 项目 CLI(pip 装 biopytools) + kmtricks/rocksdb/modelscope — 用户日常命令行入口 |")
L.append("| base | conda 基础环境(conda/mamba 本体) |")
L.append("")
L.append("---")
L.append("")
L.append("## 第四部分: 禁止使用|Part 4: Forbidden (待退役)")
L.append("")
L.append("scripts/delete_list.txt 中的 " + str(len(delete)) + " 个环境已被域环境取代或从未被模块引用, **不要**再调用它们, 也不要让新模块依赖它们; 如需其中的软件, 先并入对应域环境。")
L.append("")
L.append("> 完整新旧对应关系见 docs/env_migration_map.md; 域环境重建配方见 envs/*.yml, 保留环境重建配方见 envs/legacy/*.yaml。")

with open("docs/conda_env_software_map.md", "w", encoding="utf-8") as f:
    f.write("\n".join(L) + "\n")
print("OK 域:", len(DOMAINS), "legacy:", len(legacy), "退役:", len(delete))
