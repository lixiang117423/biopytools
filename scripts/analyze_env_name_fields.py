#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""扫描代码中的环境名引用(路径引用+名字段引用)|Scan env name refs in code
用于计算可安全删除的旧环境清单|Used to compute safe deletion list
"""
import os
import re

OLD_TO_DOMAIN = {
    "GATK_v.4.6.2.0": "align", "bcftools_v.1.22": "align", "sv_calling": "align",
    "freebayes": "align", "Genome_dedup": "align",
    "Population_genetics": "pop", "selective_sweep": "pop",
    "adamixture_v.1.0.2": "pop", "treemix_v.1.13": "pop",
    "pixy_v.2.0.0": "pop", "poplddecay_v.3.43": "pop",
    "canu_v.2.3": "asm", "hifiasm_v.0.25.0": "asm", "kmc_v.3.2.4": "asm",
    "K-mer": "asm", "merqury_v.1.3": "asm", "purge_dups_v.1.2.6": "asm",
    "genomescope_v.2.0.1": "asm", "genomescope2_v.2.1.0": "asm",
    "tidk_v.0.2.65": "asm", "getorganelle_v.1.7.71": "asm", "spades_v.4.3.0": "asm",
    "haphic": "hic", "pairtools_v.1.1.3": "hic", "yahs_v.1.2.2": "hic",
    "juicer_v.1.6": "hic",
    "Augustus_v.3.5.0": "annot", "agat_v.1.7.0": "annot",
    "gffcompare_v.0.12.10": "annot", "miniprot_v.0.18": "annot",
    "transdecoder_v.5.5.0": "annot", "eggnog-mapper_v.2.1.15": "annot",
    "orthofinder_v.3.1.5": "annot", "Blast_v.2.16.0": "annot",
    "genometools_v.1.6.5": "annot", "eviann_v.2.0.5": "annot",
    "repeatmodeler_v.2.0.7": "repeat", "repeat_identiy": "repeat",
    "ltr_retriever_v.3.0.1": "repeat",
    "RNA_Seq": "rna", "RSeQC_v.5.0.4": "rna",
    "signalp6": "protein", "resistify_v.1.3.0": "protein",
    "needle_v.1.0.3": "protein", "meme_v.5.5.9": "protein",
    "phobius_v.1.0.1": "protein", "tmmhmm_v.2.0c": "protein",
    "iqtree_v.3.0.1": "phylo", "mafft_v.7.525": "phylo", "trimal_v.1.5.0": "phylo",
    "newick_utils_v.1.6": "phylo", "wgdi_v.0.75": "phylo",
    "kakscalculator2_v.2.0.1": "phylo",
    "pggb_v.0.7.4": "pan", "kmtricks_v.1.5.1": "pan", "kmindex_v.0.6.0": "pan",
    "pan-blocks": "pan", "panman_v.0.1.4": "pan",
    "samplot_v.1.3.0": "viz", "pycirclize_v.1.10.1": "viz",
    "iseq_v.1.9.8": "misc", "primer3_v.2.6.1": "misc",
    "bbmap_v.39.81": "misc", "BioinfTools": "misc",
    "rMVP": "r", "WGCNA_v.1.73": "r",
    "BUSCO_v.6.0.0": "busco",
}

RE_PATH = re.compile(r"envs/([A-Za-z0-9_.+-]+)(?:/|$)")
RE_FIELD = re.compile(r'(?:_env|conda_env|env_name)[^=]*=\s*["\']([A-Za-z0-9_.+-]+)["\']')
RE_DEFAULT = re.compile(r'default\s*=\s*["\']([A-Za-z0-9_.+-]+)["\']')
NOISE = {"INFO", "DEBUG", "none", "auto", "both", "all", "txt", "csv", "tsv",
         "fast", "slow", "eu", "no", "yes", "eukarya", "GMYN", "MYN", "YN",
         "NG", "MA", "GY", "LPB"}

path_envs = {}
name_envs = {}
for root, dirs, files in os.walk("biopytools"):
    dirs[:] = [d for d in dirs if d not in ("__pycache__", "tests")]
    for fn in files:
        if not fn.endswith(".py"):
            continue
        p = os.path.join(root, fn)
        txt = open(p, encoding="utf-8", errors="replace").read()
        for m in RE_PATH.finditer(txt):
            e = m.group(1)
            path_envs.setdefault(e, set()).add(p)
        for line in txt.splitlines():
            for m in RE_FIELD.finditer(line):
                name_envs.setdefault(m.group(1), set()).add(p)
            for m in RE_DEFAULT.finditer(line):
                v = m.group(1)
                if re.match(r"^[A-Za-z0-9_.+-]+$", v) and v not in NOISE:
                    name_envs.setdefault(v, set()).add(p)

field_live = {e for e in name_envs if e in OLD_TO_DOMAIN}
path_live = {e for e in path_envs if e in OLD_TO_DOMAIN}
print("== 被名字段引用的 Tier1 旧环境(不能删) ==")
for e in sorted(field_live):
    print("  " + e)
print("== 名字段引用中疑似环境名(人工核对) ==")
for e in sorted(name_envs):
    if "v." in e or "_v" in e or e in OLD_TO_DOMAIN:
        print("  " + e)
print("== 仍有路径引用的 Tier1 旧环境(不能删) ==")
for e in sorted(path_live):
    print("  " + e + ": " + str(len(path_envs[e])) + " 文件")
print("== 可删除候选(Tier1 且无路径引用且无名字段引用) ==")
candidates = sorted(set(OLD_TO_DOMAIN) - path_live - field_live)
print(" ".join(candidates))
print("共 " + str(len(candidates)) + " 个")
