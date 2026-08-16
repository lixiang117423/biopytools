#!/bin/bash
# =============================================================
# conda 功能域合并求解预检|Conda domain merge solver preflight
# 只做 --dry-run 求解, 不创建任何环境|Solves only, creates nothing
#
# 用法|Usage:
#   bash scripts/dryrun_conda_domains.sh                       # 跑全部 14 个域
#   bash scripts/dryrun_conda_domains.sh ~/dryrun_logs         # 指定日志目录
#   bash scripts/dryrun_conda_domains.sh ~/dryrun_logs align pop   # 只跑指定域
#
# 镜像不可用时|When mirror is unreachable (CHANNELS 覆盖|override):
#   CHANNELS='-c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/conda-forge \
#             -c https://mirrors.tuna.tsinghua.edu.cn/anaconda/cloud/bioconda' \
#     bash scripts/dryrun_conda_domains.sh
#
# 结果判定|Result codes:
#   OK       求解成功(可合并)|Solvable
#   CONFLICT 依赖冲突(真实冲突,需调整分组)|Real dependency conflict
#   NOTFOUND 包名在channel中不存在(非冲突,包名需修正)|Package name not found
#   NETWORK  镜像连不上|Mirror unreachable
# =============================================================
set -u

OUT_DIR="${1:-$HOME/dryrun_logs}"
shift 2>/dev/null || true
DOMAIN_FILTER="$*"
mkdir -p "$OUT_DIR"

CHANNELS="${CHANNELS:--c conda-forge -c bioconda}"

# 域名|包列表 (空格分隔)|domain|packages
DOMAINS=(
  "align_jdk17|gatk4=4.6.2.0 openjdk=17 bcftools=1.22 samtools bwa freebayes bedtools"
  "align_jdk25|gatk4=4.6.2.0 bcftools=1.22 samtools bwa freebayes bedtools"
  "pop|plink vcftools admixture treemix pixy poplddecay"
  "asm|canu hifiasm kmc jellyfish merqury purge-dups genomescope2 tidk getorganelle spades"
  "hic|haphic pairtools yahs juicer"
  "annot|augustus pasa agat gffcompare miniprot transdecoder eggnog-mapper orthofinder blast genometools"
  "repeat|repeatmodeler repeatmasker ltr_retriever ltr_harvest_parallel ltr_finder_parallel"
  "rna|hisat2 stringtie fastp gffread rseqc"
  "protein|signalp6 resistify emboss meme phobius tmhmm"
  "phylo|iqtree mafft trimal newick_utils wgdi kakscalculator2 raxml-ng"
  "pan|pggb vg kmtricks kmindex panman mummer"
  "viz|samplot pycirclize"
  "misc|iseq primer3 bbmap seqkit"
  "r|r-qtl r-wgcna"
)

# 不预检的说明|Not preflighted (不在 conda / 分层留待后续):
#   gctb(自编译二进制)  mga(自研封装)  braker/fanc/deeptmhmm/DeepBSA(Tier2 legacy)
#   EDTA/cphasing/jcvi/qiime/picrust/EGAPx/BUSCO/telocomp/rnaseq_val/singularity(Tier3 例外)

SUMMARY="$OUT_DIR/00_summary.txt"
: > "$SUMMARY"
echo "时间|Time: $(date '+%F %T')" >> "$SUMMARY"
echo "域名|Domain | 结果|Result | 包数|Pkgs | 日志|Log" >> "$SUMMARY"
echo "----------------------------------------------------------------" >> "$SUMMARY"

run_domain() {
  local entry="$1"
  local domain="${entry%%|*}"
  local pkgs="${entry#*|}"
  local log="$OUT_DIR/dryrun_${domain}.log"
  local result="FAIL"
  local pkgnum="-"

  echo "" >&2
  echo "==== 求解|Solving: $domain ====" >&2
  echo "命令|Command: mamba create -n dryrun_${domain} --dry-run $CHANNELS $pkgs" >&2

  mamba create -n "dryrun_${domain}" --dry-run $CHANNELS $pkgs >"$log" 2>&1 || true

  if grep -q 'Dry run. Not executing the transaction' "$log"; then
    result="OK"
  elif grep -qi 'PackagesNotFoundError\|The following packages are not available\|not found in the channels' "$log"; then
    result="NOTFOUND"
  elif grep -qi 'Could not solve' "$log"; then
    result="CONFLICT"
  elif grep -q 'Download error' "$log"; then
    result="NETWORK"
  fi

  pkgnum=$(grep -m1 -oE '[0-9]+ packages' "$log" | grep -oE '[0-9]+' || true)
  [ -z "$pkgnum" ] && pkgnum="-"

  echo "$domain | $result | $pkgnum | $log" >> "$SUMMARY"
  echo ">> 结果|Result: $result (包数|pkgs: $pkgnum)" >&2

  if [ "$result" = "CONFLICT" ] || [ "$result" = "NOTFOUND" ]; then
    echo "----- 关键信息|Key lines -----" >&2
    grep -iE 'conflicts for|Could not solve|PackagesNotFoundError|not available|The following packages' "$log" | tail -8 >&2
    echo "----- 日志末尾|Log tail -----" >&2
    tail -6 "$log" >&2
  fi
}

ran=0
for entry in "${DOMAINS[@]}"; do
  d="${entry%%|*}"
  if [ -n "$DOMAIN_FILTER" ]; then
    echo " $DOMAIN_FILTER " | grep -q " $d " || continue
  fi
  run_domain "$entry"
  ran=$((ran+1))
done

echo "" >&2
echo "================ 汇总|Summary ================" >&2
cat "$SUMMARY" >&2
echo "" >&2
echo "完成|Done: $ran 个域。日志目录|Log dir: $OUT_DIR" >&2
