#!/bin/bash
# =============================================================
# 删除超算上已被域环境取代的旧 conda 环境|Delete obsolete legacy conda envs
#
# 删除前提|PREREQUISITES:
#   1. 新代码已同步到超算(github 拉取)且模块跑通验证
#   2. 确认 docs/env_migration_map.md 第二节「保留不动」清单之外
#
# 安全设计|Safety:
#   - 默认 DRY-RUN, 只列不删|dry-run by default
#   - 只删脚本内硬编码的 51 个白名单环境|only the 51 allowlisted envs
#   - 14 个域环境/base/其他一切环境永不触碰|never touches anything else
#
# 用法|Usage:
#   bash scripts/delete_obsolete_envs.sh              # dry-run 预览
#   bash scripts/delete_obsolete_envs.sh --delete     # 真正删除(需输入 yes)
# =============================================================
set -u

CONDA="${CONDA:-/share/org/YZWL/yzwl_lixg/miniforge3/bin/conda}"
ENVS_DIR="${ENVS_DIR:-/share/org/YZWL/yzwl_lixg/miniforge3/envs}"
DO_DELETE=0
if [ "${1:-}" = "--delete" ]; then DO_DELETE=1; fi

# 51 个可删除旧环境|51 deletable legacy envs
DELETE_LIST=(
  BUSCO_v.6.0.0 BioinfTools GATK_v.4.6.2.0 Genome_dedup K-mer RNA_Seq
  RSeQC_v.5.0.4 agat_v.1.7.0 bcftools_v.1.22 canu_v.2.3 eggnog-mapper_v.2.1.15
  eviann_v.2.0.5 freebayes genomescope2_v.2.1.0 genometools_v.1.6.5
  getorganelle_v.1.7.71 gffcompare_v.0.12.10 haphic hifiasm_v.0.25.0
  iqtree_v.3.0.1 kakscalculator2_v.2.0.1 kmindex_v.0.6.0 kmtricks_v.1.5.1
  ltr_retriever_v.3.0.1 mafft_v.7.525 meme_v.5.5.9 miniprot_v.0.18
  needle_v.1.0.3 newick_utils_v.1.6 orthofinder_v.3.1.5 pairtools_v.1.1.3
  pan-blocks phobius_v.1.0.1 pixy_v.2.0.0 poplddecay_v.3.43 primer3_v.2.6.1
  purge_dups_v.1.2.6 pycirclize_v.1.10.1 repeat_identiy repeatmodeler_v.2.0.7
  resistify_v.1.3.0 samplot_v.1.3.0 selective_sweep signalp6 spades_v.4.3.0
  sv_calling tidk_v.0.2.65 tmmhmm_v.2.0c trimal_v.1.5.0 wgdi_v.0.75
  yahs_v.1.2.2
)

# 保护名单|Protected (防呆)
PROTECT=(
  align pop asm hic annot repeat rna protein phylo pan viz misc r busco base
  qiime_v.2024.10.1 picrust_v.2.6.3 EDTA_v.2.2.2 EGAPx_v.0.4.0-alpha cphasing
  jcvi_v.1.5.7 telocomp rnaseq_val singularity_v.3.8.7 psvcp_v.1.0.1
  braker_v.3.0.8 fanc_v.0.9.23b deeptmhmm_v.1.0 DeepBSA HiC-Pro_v3.1.0
  SubPhaser pasa_v.2.5.3 gctb plothic_v.1.0.0 genomesyn2 alignoth a-liner
  kmeriaenv vcf2gwas_v.0.8.9 signalp_v.3.0b adamixture_v.1.0.2 vg_v.1.7.0
  kmc_v.3.2.4 merqury_v.1.3
)

echo "==== 前置检查|Preflight ===="
echo "conda: $CONDA"
echo "envs目录|envs dir: $ENVS_DIR"
[ -x "$CONDA" ] || { echo "错误|Error: conda 不存在: $CONDA"; exit 1; }
[ -d "$ENVS_DIR" ] || { echo "错误|Error: envs 目录不存在"; exit 1; }
if [ "$DO_DELETE" = "0" ]; then
  echo "模式|Mode: DRY-RUN (加 --delete 才真正删除)"
else
  echo "模式|Mode: DELETE !!"
fi

echo ""
echo "==== 白名单自检|Allowlist self-check ===="
BAD=0
for e in "${DELETE_LIST[@]}"; do
  for p in "${PROTECT[@]}"; do
    if [ "$e" = "$p" ]; then
      echo "错误|Conflict: $e 同时在删除和保护名单"; BAD=1
    fi
  done
done
[ "$BAD" = "1" ] && { echo "自检失败, 中止|Abort"; exit 1; }
echo "通过|OK: 删除清单与保护名单无冲突"

echo ""
echo "==== 扫描|Scanning ===="
TO_DELETE=()
for e in "${DELETE_LIST[@]}"; do
  d="$ENVS_DIR/$e"
  if [ -d "$d" ]; then
    size=$(du -sh "$d" 2>/dev/null | cut -f1)
    echo "  [存在|exists] $e ($size)"
    TO_DELETE+=("$e")
  else
    echo "  [不存在, 跳过|missing] $e"
  fi
done

echo ""
echo "==== 计划|Plan ===="
echo "  待删除|To delete: ${#TO_DELETE[@]} 个"

if [ "$DO_DELETE" = "0" ]; then
  echo "DRY-RUN 结束|Done. 未删除任何环境|Nothing deleted."
  echo "确认无误后执行|To execute: bash scripts/delete_obsolete_envs.sh --delete"
  exit 0
fi

read -r -p "警告|WARNING: 确认删除以上 ${#TO_DELETE[@]} 个环境? 输入 yes 继续: " ANS
if [ "$ANS" != "yes" ]; then echo "已取消|Cancelled"; exit 0; fi

echo ""
FAIL=0
for e in "${TO_DELETE[@]}"; do
  echo "  删除|Removing: $e"
  if "$CONDA" env remove -p "$ENVS_DIR/$e" -y >/dev/null 2>&1; then
    echo "    已删除|Removed: $e"
  else
    echo "    失败|Failed: $e"; FAIL=1
  fi
done
echo ""
echo "完成|Done. 失败|Failed: $FAIL"
