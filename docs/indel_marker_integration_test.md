# indel_marker 端到端集成测试|indel_marker End-to-End Integration Test

> 在计算节点运行|Run on a **compute node**, NOT the login node（CLAUDE.md §10.3）。
> 本文档涉及真实 bcftools / samtools / primer3 调用，登录节点禁止执行|This guide invokes real bcftools/samtools/primer3 binaries; login-node execution is forbidden.

---

## 0. 目的|Purpose

验证 `biopytools indel-marker` 在真实多样本 VCF + BAM + 参考基因组上的完整数据流：
提取 INDEL → 群体共显性判定 → 覆盖度质控与 deletion 骤降验证 → 侧翼序列提取 →
primer3 引物设计 → 候选主表/报告输出。重点确认人为注入的已知 deletion 能被回收为
`direction=resistant_specific` 且 `passes_deletion_drop=True` 的候选|Verify the full
dataflow on a real multi-sample VCF+BAM+reference: INDEL extraction → codominant
genotype calling → coverage QC + deletion-drop validation → flank extraction →
primer3 design → candidate table/report. The key assertion is that an injected
known deletion is recovered as `resistant_specific` with `passes_deletion_drop=True`.

---

## 1. 运行环境准备|Runtime Environment Setup

### 1.1 conda 环境|Conda environment

primer3-py 位于 conda 环境 `primer3_v.2.6.1`（python 3.13, primer3-py 2.3.0，
含 `primer3_core` 命令行）。该环境**缺失 pyyaml**，而 `biopytools/common/paths.py`
在读取 `~/.config/biopytools/config.yml` 时会 `import yaml`，因此必须补装：

```bash
# 进入计算节点后|after landing on a compute node
conda activate primer3_v.2.6.1
pip install pyyaml          # common/paths.py 依赖|required by common/paths.py
# 或在不方便装包时，用一个同时含 primer3-py + pyyaml + bcftools + samtools 的环境
```

> 替代方案|Alternative：若不想改 `primer3_v.2.6.1`，可新建一个环境：
> `conda create -n indel_marker_test -c bioconda -c conda-forge primer3-py pyyaml bcftools samtools bwa samtools python=3.13`

### 1.2 bcftools / samtools / bwa

需在 PATH 上可见（或通过环境变量 / 配置文件指定）|Must be on PATH or via env/config：

```bash
which bcftools samtools bwa          # 应均返回路径|all should resolve
bcftools --version | head -1
samtools --version | head -1
```

若不在默认 PATH，按优先级（高→低）配置|If not on default PATH, configure by priority (high→low)：

1. **环境变量|Env var**：`export BCFTOOLS_PATH=~/miniforge3/envs/.../bin/bcftools`
   `export SAMTOOLS_PATH=~/miniforge3/envs/.../bin/samtools`
2. **配置文件|Config file** `~/.config/biopytools/config.yml`：
   ```yaml
   tools:
     bcftools: ~/miniforge3/envs/<env>/bin/bcftools
     samtools: ~/miniforge3/envs/<env>/bin/samtools
   ```
3. 代码默认走 `shutil.which`，PATH 上有即可|Code default uses `shutil.which`.

### 1.3 影子模块陷阱|Shadow-module gotcha

仓库内存在 `biopytools/primer3/` 子包，当进程 cwd 恰好是**包目录**时，
会遮蔽 pip 安装的 `primer3-py`（`import primer3` 会命中本地子包）。务必从
**数据/工作目录**调用，不要 `cd` 进包目录|The repo ships a `biopytools/primer3/`
subpackage that shadows the pip `primer3-py` package when the process cwd is the
package dir. Always invoke from your **data/work directory**, never from inside
the package. `biopytools indel-marker` CLI 天然从调用 cwd 运行，正常使用即安全。

---

## 2. 准备小数据集|Prepare Small Dataset

目标：构造 4–6 个 Illumina 双端样本，在抗病组样本中人为注入一个已知 ~50bp deletion，
经 BWA→排序→索引→joint-call 得到多样本 VCF|Goal: build 4–6 Illumina PE samples,
inject a known ~50bp deletion into the resistant group only, then
BWA→sort→index→joint-call a multi-sample VCF.

### 2.1 取参考基因组片段|Reference fragment

```bash
mkdir -p ~/indel_test && cd ~/indel_test
# 取某条染色体约 1Mb 区段作为参考|take ~1Mb of a chromosome as reference
samtools faidx real_genome.fa chr01:1-1000000 > ref.fa
samtools faidx ref.fa                       # 建 faidx 索引|build faidx index
```

### 2.2 构造抗/感两套参考|Build R/S mutated references

选定一个注入位置，例如 `chr01:500000`，在抗病组参考里删除 50bp、可选插入 20bp|Pick
an injection site (e.g. chr01:500000); delete 50bp in the resistant reference and
optionally insert 20bp：

```bash
INJECT=500000; DEL_LEN=50
# 抗病组参考(删50bp)|resistant reference (50bp deletion)
awk -v inject=$INJECT -v dl=$DEL_LEN 'BEGIN{split("",a)} {
  if ($0 ~ /^>/) {print; next}
  s=$0; left=substr(s,1,inject-1); right=substr(s,inject+dl);
  print left right
}' ref.fa > ref.R_del.fa
# 感病组参考与 ref.fa 相同|susceptible reference == ref.fa
cp ref.fa ref.S.fa
samtools faidx ref.R_del.fa
```

### 2.3 模拟 reads 并人为注入突变|Simulate reads & inject the variant

用 `wgsim`（samtools 自带）或 `dwgsim` 为每组样本模拟 reads。抗病组从 `ref.R_del.fa`
模拟（自然携带 deletion），感病组从 `ref.S.fa` 模拟|Simulate reads per sample with
`wgsim`/`dwgsim`; R-group from `ref.R_del.fa` (carries the deletion), S-group from `ref.S.fa`：

```bash
# 每样本 ~30x，150bp PE，错误率 1%|~30x/sample, 150bp PE, 1% error
simulate_one () {
  local ref=$1 name=$2 seed=$3
  wgsim -N 200000 -1 150 -2 150 -e 0.01 -S $seed -r 0 $ref ${name}_R1.fq ${name}_R2.fq
}
simulate_one ref.R_del.fa R1 11
simulate_one ref.R_del.fa R2 22
simulate_one ref.S.fa    S1 33
simulate_one ref.S.fa    S2 44
# 可选第3个样本以提升组内一致性|optional 3rd sample to boost within-group consistency
# simulate_one ref.R_del.fa R3 55
# simulate_one ref.S.fa    S3 66
```

### 2.4 比对 / 排序 / 索引|Align / sort / index

用 `bwa`（本仓库 `biopytools bwa` 或裸 bwa 均可）比对到**原始 ref.fa**，确保 deletion
以比对缺口形式呈现|Align to the **original ref.fa** so the deletion shows up as a gap：

```bash
bwa index ref.fa
for s in R1 R2 S1 S2; do
  bwa mem -t 16 -R "@RG\tID:${s}\tSM:${s}" ref.fa ${s}_R1.fq ${s}_R2.fq \
    | samtools sort -@ 8 -o ${s}.bam -
  samtools index ${s}.bam
done
```

### 2.5 联合 call 出多样本 VCF|Joint-call multi-sample VCF

用 bcftools mpileup/call（或 GATK HaplotypeCaller → GenotypeGVCFs）联合 call|Joint-call
with bcftools mpileup/call (or GATK)：

```bash
bcftools mpileup -f ref.fa -R chr01:499000-501000 R1.bam R2.bam S1.bam S2.bam \
  | bcftools call -mv -Oz -o test.vcf.gz
tabix -p vcf test.vcf.gz
# 快速确认注入位点确实被 call 出 deletion|confirm the injected deletion was called
bcftools view -v indels test.vcf.gz | grep -A1 "chr01" | head
```

> 期望在 chr01:500000 附近看到长度 ~50 的 DEL 记录，且抗病组基因型为 `1/1`、感病组为 `0/0`|Expect
> a ~50bp DEL near chr01:500000 with R-group `1/1` and S-group `0/0`.

---

## 3. 准备 samplesheet|Prepare Samplesheet

格式：`sample_name <分隔> group <分隔> bam_path`，**表头可选**，分隔符默认 tab
（行内无 tab 时回退空格）；分组为任意两种不同标签，`R`/`S`、`resistant`/`susceptible`
等常见别名自动识别，其余按首次出现顺序映射为抗病(resistant)/感病(susceptible)|Format:
`sample_name<SEP>group<SEP>bam_path`, **header optional**, delimiter defaults to tab
(falls back to spaces when no tab); group = any two distinct labels — `R`/`S`,
`resistant`/`susceptible` are auto-recognized aliases, others map to R/S by first-seen order.

```bash
cat > ss.tsv <<'EOF'
sample_name	group	bam_path
R1	resistant	~/indel_test/R1.bam
R2	resistant	~/indel_test/R2.bam
S1	susceptible	~/indel_test/S1.bam
S2	susceptible	~/indel_test/S2.bam
EOF
```

> 注意|Notes：
> - 表头可选（首行第2列为 `group` 会自动识别跳过）；至少 3 列|header optional
>   (a first row whose 2nd column is `group` is auto-detected and skipped); ≥3 columns.
> - 分组须恰好两种不同标签|exactly 2 distinct group labels required.
> - 也支持无表头 + `R`/`S` 简写（见 `samplesheet.py` 别名表）|header-less +
>   `R`/`S` shorthand also supported (see alias table in `samplesheet.py`).
> - 每组样本数 ≥ `--min-samples-per-group`（默认 1）|sample count per group ≥
>   `--min-samples-per-group` (default 1).
> - bam 路径会被 `expand_path` 展开，支持 `~`|bam paths go through `expand_path`, `~` supported.

---

## 4. 运行|Run

从工作目录（`~/indel_test`）调用，避免影子模块陷阱|Invoke from your work directory
to avoid the shadow-module gotcha：

```bash
cd ~/indel_test
# 方式A：环境已 activate，biopytools 在 PATH|env already activated, biopytools on PATH
biopytools indel-marker -v test.vcf.gz -s ss.tsv -g ref.fa -o out/

# 方式B：用 conda run 包装（推荐，确保 primer3-py/pyyaml 可用）|wrap with conda run (recommended)
conda run -n primer3_v.2.6.1 --no-capture-output biopytools indel-marker \
  -v test.vcf.gz -s ss.tsv -g ref.fa -o out/
```

运行后在 `out/` 下生成标准目录结构|Standard output tree is created under `out/`：

```
out/
├── 00_pipeline_info/
├── 01_vcf_extract/      indels.gt_matrix.tsv
├── 02_genotype_call/
├── 03_coverage/         indels.coverage.tsv
├── 04_sequence/         indels.flank.fa
├── 05_primer/
├── 06_results/          indel_marker.candidates.tsv / .bed / summary.txt
└── 99_logs/             indel_marker.log
```

---

## 5. 预期结果|Expected Results

### 5.1 候选主表|Candidate master table

`out/06_results/indel_marker.candidates.tsv` 应包含注入的 ~50bp deletion 候选，
关键列断言|Must contain the injected ~50bp deletion with these key columns：

| 列|Column | 期望值|Expected |
|-------|--------|
| `candidate_id` | `chr01:500000-500050:DEL:R_spec`（格式 `{chrom}:{pos}-{end}:{DEL\|INS}:{R\|S}_spec`） |
| `indel_type` | `deletion` |
| `indel_size` | `50` |
| `direction` | `resistant_specific` |
| `present_group` | `resistant` |
| `R_hom_alt_rate` | ≥ `--min-group-consistency`（默认 0.9） |
| `S_hom_ref_rate` | ≥ `--min-group-consistency` |
| `passes_coverage_qc` | `True` |
| `passes_deletion_drop` | `True`（deletion 长度 ≥ `deletion_size_for_coverage_check=30` 才做骤降检查） |
| `deletion_depth_ratio` | < `--deletion-depth-ratio`（默认 0.3），即 present 组覆盖度显著低于 absent 组 |
| `primer_status` | `ok`（注入位点侧翼足够长，primer3 应成功设计） |

校验命令|Verify：

```bash
# 找到注入候选|locate the injected candidate
awk -F'\t' 'NR==1 || ($3==500000 && $7=="deletion" && $9=="resistant_specific")' \
  out/06_results/indel_marker.candidates.tsv

# 任何一个候选应 primer_status=ok|at least one candidate should have primer_status=ok
awk -F'\t' 'NR>1 && $24=="ok"' out/06_results/indel_marker.candidates.tsv
```

### 5.2 摘要|Summary

`out/06_results/indel_marker.summary.txt` 中 `成功设计引物|With primers` ≥ 1。

### 5.3 日志|Log

`out/99_logs/indel_marker.log` 应记录每步执行的完整命令（bcftools view、samtools depth
等，符合 CLAUDE.md §2.2.1 命令执行日志规范）|Each step's full command (bcftools view,
samtools depth, etc.) is logged at INFO per CLAUDE.md §2.2.1.

---

## 6. 反向验证|Sanity Checks

### 6.1 覆盖度骤降|Coverage drop

用 `samtools depth` 独立确认 R 组在 deletion 区覆盖度低于 S 组|Independently confirm
R-group coverage drops inside the deletion region vs the S-group：

```bash
# deletion 区间|deletion region chr01:500000-500050
for s in R1 R2 S1 S2; do
  mean=$(samtools depth -r chr01:500000-500050 ${s}.bam | awk '{s+=$3; n++} END{if(n) printf "%.1f", s/n; else print 0}')
  echo "$s mean_depth_in_del=$mean"
done
# 期望：R1/R2 ≈ 0（reads 跨越缺口映射不到参考该区段），S1/S2 ≈ 30x
# Expect R1/R2 ≈ 0 (reads span the gap, don't map to ref region), S1/S2 ≈ 30x
```

并与 `out/03_coverage/indels.coverage.tsv` 中记录的覆盖度矩阵核对一致|Cross-check
against the coverage matrix in `out/03_coverage/indels.coverage.tsv`.

### 6.2 引物 in silico PCR（可选）|Optional in silico PCR

用候选主表中的 `left_primer` / `right_primer` 对 `ref.R_del.fa` 与 `ref.S.fa` 做
`SequenceServer`/`isPcr`/`primer3_core --check_primers`，预期抗病组扩增产物比感病组
短约 50bp（deletion 长度）|With the designed primers, amplicon on `ref.R_del.fa` should
be ~50bp shorter than on `ref.S.fa` (the deletion length).

---

## 7. 可调参数参考|Tunable Parameters Reference

| 参数|Param | 默认|Default | 作用|Effect |
|------|--------|--------|--------|
| `--min-indel-size` | 10 | 最小 INDEL 长度（过滤过短）|min INDEL size |
| `--max-indel-size` | 100 | 最大 INDEL 长度（过滤过长，避开大型 SV）|max INDEL size |
| `--min-group-consistency` | 0.9 | 组内纯合一致比例（1.0=严格无杂合/缺失）|within-group homozygous consistency |
| `--min-samples-per-group` | 2 | 每组最少样本数|min samples per group |
| `--min-depth` | 10 | 覆盖度质控阈值|min depth QC threshold |
| `--deletion-depth-ratio` | 0.3 | deletion 骤降阈值（present/absent 平均覆盖度比）|deletion drop threshold |
| `--flank-length` | 300 | 引物设计侧翼长度|flank length for primer design |
| `--threads` | 12 | 线程数|threads |

调参提示|Tuning tips：
- 若注入的 deletion 未被回收，先放宽 `--min-group-consistency` 到 0.7 看是否为模拟噪声导致基因型不纯|If
  the injected deletion isn't recovered, relax `--min-group-consistency` to 0.7 first — simulation noise may make genotypes impure.
- 若 `passes_deletion_drop=False`，检查模拟覆盖度是否足够（wgsim `-N` 调高）或放宽 `--deletion-depth-ratio`|If
  `passes_deletion_drop=False`, raise coverage (`-N`) or relax `--deletion-depth-ratio`.
- insertion 类候选不做骤降检查（`deletion_depth_ratio=NA`，`passes_deletion_drop=NA`），属正常行为|Insertion
  candidates skip the drop check (NA) by design.

---

## 8. 故障排查|Troubleshooting

| 症状|Symptom | 原因|Cause | 处理|Fix |
|---------|--------|------|------|
| `ModuleNotFoundError: No module named 'yaml'` | primer3_v.2.6.1 缺 pyyaml | `conda run -n primer3_v.2.6.1 pip install pyyaml` |
| `ModuleNotFoundError: No module named 'primer3'` 但环境已装 | cwd 在包目录，影子模块 `biopytools/primer3/` 遮蔽 | 切到数据目录再运行（见 §1.3） |
| `命令|Command:` 日志里 `conda run` 缺 `-n <env>` | 传了命令名而非完整路径给 `build_conda_command`（本模块不应出现） | 确认未手动改 utils.py；用 CLI 入口调用 |
| 候选为空|No candidates | 注入位点未 call 出、或基因型不纯 | 放宽 `--min-group-consistency`；检查 `bcftools view test.vcf.gz` |
| `passes_deletion_drop=False` | 模拟覆盖度不足或 deletion 太短 | 调高 wgsim `-N`；确认 deletion 长度 ≥ 30 |
| 候选表 `primer_status` 全 fail | 侧翼序列不足或含 N | 调大 `--flank-length`；检查 `04_sequence/indels.flank.fa` |

---

## 9. 限制说明|Limitations

- 本文档不涉及 Git 操作（超算上禁止 commit，见 CLAUDE.md 顶部警告）。代码提交只在本地 Mac 进行|No
  Git operations on the supercomputer; commits happen on the local Mac only.
- 集成测试必须在计算节点运行，登录节点仅用于编辑/提交|Integration tests run on compute
  nodes only; login node is for editing/submitting.
- 本流程依赖外部 joint-call（bcftools/GATK）产出 VCF；上游 call 质量直接影响候选回收率|Upstream
  joint-call quality directly drives candidate recovery.
