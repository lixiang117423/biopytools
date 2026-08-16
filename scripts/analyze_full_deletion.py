#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""全量删除清单计算|Compute full deletion list
标准|Rule: 模块(代码)未引用的环境一律列入删除清单
|Delete all envs not referenced by module code.
输出分组|Output groups: 已迁移51 / 版本分裂对 / 文档提及 / 纯孤儿 / 怪异名
"""
import os
import re

DOMAINS = {"align", "pop", "asm", "hic", "annot", "repeat", "rna",
           "protein", "phylo", "pan", "viz", "misc", "r", "busco"}
ALWAYS_KEEP = {"base"}

# ---- 1. 备份环境名 ----
backup = sorted(fn[:-5] for fn in os.listdir("conda_envs_backup")
                if fn.endswith(".yaml"))

# ---- 2. 代码引用(非 tests) ----
RE_PATH = re.compile(r"envs/([A-Za-z0-9_.+-]+)")
RE_FIELD = re.compile(r'(?:_env|conda_env|env_name)[^=]*=\s*["\']([A-Za-z0-9_.+-]+)["\']')
RE_DEFAULT = re.compile(r'default\s*=\s*["\']([A-Za-z0-9_.+-]+)["\']')
NOISE = {"INFO", "DEBUG", "none", "auto", "both", "all", "txt", "csv", "tsv",
         "fast", "slow", "eu", "no", "yes", "eukarya", "GMYN", "MYN", "YN",
         "NG", "MA", "GY", "LPB"}
referenced = set()
for root, dirs, files in os.walk("biopytools"):
    dirs[:] = [d for d in dirs if d not in ("__pycache__", "tests")]
    for fn in files:
        if not (fn.endswith(".py") or fn.endswith(".yml") or fn.endswith(".yaml")):
            continue
        txt = open(os.path.join(root, fn), encoding="utf-8", errors="replace").read()
        for m in RE_PATH.finditer(txt):
            referenced.add(m.group(1))
        for line in txt.splitlines():
            for m in RE_FIELD.finditer(line):
                referenced.add(m.group(1))
            for m in RE_DEFAULT.finditer(line):
                v = m.group(1)
                if re.match(r"^[A-Za-z0-9_.+-]+$", v) and v not in NOISE:
                    referenced.add(v)

# ---- 3. 文档提及(仅标记, 不改变删除结论) ----
docs_mentioned = set()
for root, dirs, files in os.walk("."):
    dirs[:] = [d for d in dirs if d not in ("__pycache__", "tests", "conda_envs_backup", ".git")]
    for fn in files:
        if not (fn.endswith(".md") or fn.endswith(".yml") or fn.endswith(".yaml")):
            continue
        p = os.path.join(root, fn)
        if p.startswith("./conda_envs_backup"):
            continue
        txt = open(p, encoding="utf-8", errors="replace").read()
        for m in RE_PATH.finditer(txt):
            docs_mentioned.add(m.group(1))
        for m in RE_FIELD.finditer(txt):
            docs_mentioned.add(m.group(1))

# ---- 4. 计算 ----
keeps = referenced | DOMAINS | ALWAYS_KEEP
DELETE = [e for e in backup if e not in keeps]

PREV51 = {
    "BUSCO_v.6.0.0", "BioinfTools", "GATK_v.4.6.2.0", "Genome_dedup", "K-mer",
    "RNA_Seq", "RSeQC_v.5.0.4", "agat_v.1.7.0", "bcftools_v.1.22", "canu_v.2.3",
    "eggnog-mapper_v.2.1.15", "eviann_v.2.0.5", "freebayes", "genomescope2_v.2.1.0",
    "genometools_v.1.6.5", "getorganelle_v.1.7.71", "gffcompare_v.0.12.10",
    "haphic", "hifiasm_v.0.25.0", "iqtree_v.3.0.1", "kakscalculator2_v.2.0.1",
    "kmindex_v.0.6.0", "kmtricks_v.1.5.1", "ltr_retriever_v.3.0.1", "mafft_v.7.525",
    "meme_v.5.5.9", "miniprot_v.0.18", "needle_v.1.0.3", "newick_utils_v.1.6",
    "orthofinder_v.3.1.5", "pairtools_v.1.1.3", "pan-blocks", "phobius_v.1.0.1",
    "pixy_v.2.0.0", "poplddecay_v.3.43", "primer3_v.2.6.1", "purge_dups_v.1.2.6",
    "pycirclize_v.1.10.1", "repeat_identiy", "repeatmodeler_v.2.0.7",
    "resistify_v.1.3.0", "samplot_v.1.3.0", "selective_sweep", "signalp6",
    "spades_v.4.3.0", "sv_calling", "tidk_v.0.2.65", "tmmhmm_v.2.0c",
    "trimal_v.1.5.0", "wgdi_v.0.75", "yahs_v.1.2.2",
}

# 版本分裂对(旧版 vs 被引用新版, 大小写不敏感归并判断)
def base_name(e):
    m = re.match(r"^(.*?)(?:[_-]?v\.?\d.*)$", e.lower())
    return m.group(1) if m else e.lower()
ref_bases = {base_name(e) for e in referenced}

groups = {"已迁移51": [], "版本分裂对": [], "文档提及": [], "纯孤儿": [], "怪异名": []}
for e in sorted(DELETE):
    if e == "Name":
        groups["怪异名"].append(e)
    elif e in PREV51:
        groups["已迁移51"].append(e)
    elif base_name(e) in ref_bases:
        groups["版本分裂对"].append(e)
    elif e in docs_mentioned:
        groups["文档提及"].append(e)
    else:
        groups["纯孤儿"].append(e)

print("=== 统计 ===")
print("备份环境总数: " + str(len(backup)))
print("代码引用(保留): " + str(len(referenced)))
print("域环境+base: " + str(len(DOMAINS) + 1))
print("待删除总数: " + str(len(DELETE)))
for k, v in groups.items():
    print("  " + k + ": " + str(len(v)))
print()
print("=== 代码引用但备份中不存在(死引用, 信息) ===")
dead = sorted(referenced - set(backup))
print(" ".join(dead))
print()
for k, v in groups.items():
    if k == "已迁移51":
        continue  # 不重复打印 51 个
    print("=== " + k + " ===")
    print(" ".join(v))
print()
print("=== 保留清单(供防呆名单使用) ===")
protect_full = sorted(DOMAINS | ALWAYS_KEEP | (set(backup) - set(DELETE)))
print(" ".join(protect_full))
print()
print("=== 删除清单总表(供脚本使用) ===")
print(" ".join(sorted(DELETE)))

# 落盘供删除脚本使用|Write lists for the deletion script
with open("scripts/delete_list.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(sorted(DELETE)) + "\n")
with open("scripts/protect_list.txt", "w", encoding="utf-8") as f:
    f.write("\n".join(protect_full) + "\n")
print()
print("已写 scripts/delete_list.txt (" + str(len(DELETE)) + ") 与 scripts/protect_list.txt (" + str(len(protect_full)) + ")")
