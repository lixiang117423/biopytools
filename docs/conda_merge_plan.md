# Conda 环境合并方案|Conda Environment Consolidation Plan

> 基于|Based on: `docs/conda_env_audit.md`（代码引用清单）+ `conda_envs_backup/`（超算 214 个环境实测导出）
> 目标|Goal: 101 个代码引用环境 → **约 19~20 个受管环境**（减少 80%）

## 一、数据底座|Data Foundation

| 指标|Metric | 数值|Value |
|---|---|---|
| 超算实际环境|HPC actual envs | 214 |
| 代码引用环境（归一化）|Code-cited canonical envs | 112（code 102 / docs-only 10） |
| 名字匹配成功|Name-matched | 101 |
| 死引用|Dead refs | 21（真实问题 5 条，其余为占位符/测试/文档示例） |
| 孤儿环境（超算有、代码未引用）|Orphan envs | 113 |

## 二、死引用修复清单（第一优先级，独立于合并）|Dead-Ref Fixes

> 这些是「代码/文档写了，超算上没有这个环境」的实锤问题，无论是否合并都必须先修：

| # | 代码引用|Cited | 实际位置|Actual | 修复动作|Fix | 严重度|
|---|---|---|---|---|---|
| 1 | （已核实非bug）`vg_v.1.7.0` | 环境存在，内含 vg 1.70.0 | 不修；另有孤儿环境 `vg_v.1.67.0`（新版），列入版本分裂退役候选 | - |
| 2 | `kakscalculator` | `kakscalculator2_v.2.0.1` | kaks/config.py:63 → 改环境名 | 高 |
| 3 | `flye` | flye 在 `telocomp` 环境里 | common/config_template.yml:31 → 改默认值 | 中 |
| 4 | `raxml-ng_v.2.0.1` | raxml-ng 在 `orthofinder_v.3.1.5` | 仅 docs/superpowers 计划 + 测试引用，生产代码无引用 → 无需处理 | 低 |
| 5 | `sra-tools` | `sratoolkit_v.2.5.7` / iseq 环境 | docs/sra2fastq.md 是「环境创建教程」（自带 conda create 命令），自成体系 → 无需处理 | - |
| 6 | `minigraph` | 在 `swave_v.1.2` 环境 | docs/minigraph.md 同上为创建教程 → 无需处理 | - |
| 7 | `bwa_env`/`bwa_v0.7.17`/`busco_env` | bwa 在 9 个环境、busco 在 BUSCO_v.6.0.0 | 均为文档创建教程/示例 → 无需处理 | - |

> 占位符（`xxx`/`ENV`/`None`/`X`/`e`/`BCFTOOLS_ENV`/`<var>`/`your_env`/`tool_env`/`custom_panman_env`/`my_ltr_harvest`/`jcvi_custom`）均来自测试断言、文档示例与 nlr_annotator 注释，不在本次合并范围（nlr_annotator 的 `xxx` 见审计报告 4.4）。

## 三、合并架构：16 域 + 例外|Target Architecture

### Tier 1 —— 直接合并（成员版本兼容，无需重装）|Merge Directly

| 新域|Domain | 吸收成员|Members | 主要冲突点与处理|Conflicts |
|---|---|---|---|
| `align` | GATK_v.4.6.2.0(jdk17) + bcftools_v.1.22 + sv_calling + freebayes + Genome_dedup | jdk 统一 17（GATK 要求 17+）；sv_calling 的 py3.9 仅承载二进制工具，无需 python 版本统一 |
| `pop` | Population_genetics + selective_sweep + adamixture + gctb + poplddecay + treemix(R4.5.2) + pixy(py3.11) | pixy(228 pkgs) 与 treemix(R) 重装到统一 py3.13/R4.5 需 `--dry-run` 验证 |
| `asm` | canu(jdk25) + hifiasm + kmc + K-mer + merqury(R4.2.3,jdk24) + purge_dups + genomescope×2(R4.4.3) + tidk + mga(py3.10) + getorganelle + spades | **R 4.2.3 vs 4.4.3 统一**、**jdk 24 vs 25 统一**、mga 的 py3.10 python 生态重装验证 |
| `hic` | haphic + pairtools + yahs + juicer_v.1.6(jdk23, 含 python2 包) | jdk 统一 23；juicer 的 python2 独立包保留（不碰 python3） |
| `annot` | Augustus + pasa(R4.4.3) + agat(R4.4.3) + gffcompare + miniprot + transdecoder(R4.4.3) + eggnog-mapper(py3.11) + orthofinder(py3.14) + Blast + genometools + eviann | 都是新版本生态，统一 py3.13/R4.4.3 后求解验证；eggnog-mapper 的 py3.11 重装 |
| `repeat` | repeatmodeler + repeat_identiy + ltr_retriever + ltr_harvest_parallel + ltr_finder_parallel | repeatmodeler(py3.12.11) 升 py3.13 重装验证 |
| `rna` | RNA_Seq + RSeQC(R4.6.0) | R 4.6.0 与 annot 域 4.4.3 不同——RSeQC 的 R 依赖极轻，可保留 4.6 或降级验证 |
| `protein` | signalp6(py3.14) + resistify(py3.12.11) + needle + meme(py3.14) + phobius + tmmhmm + signalp_v.3.0b | resistify 与 signalp6 的 python 版本差 2 个小代际，重装验证；signalp_v.3.0b 为纯二进制，直接并 |
| `phylo` | iqtree + mafft + trimal + newick_utils + wgdi(py3.14) + **kakscalculator2_v.2.0.1**(孤儿收编) + **raxml-ng**(从 orthofinder 取) | 全轻量，死引用 2/4 在此域顺带解决 |
| `pan` | pan-blocks + **vg_v.1.67.0**(收编) + pggb(py3.12.13) + kmtricks + kmindex + panman(py3.11) + swave_v.1.2(含 minigraph) | pggb(176 pkgs) 重装验证；死引用 1/6 在此域顺带解决 |
| `viz` | samplot + alignoth + a-liner + genomesyn2 + plothic + pycirclize | 全轻量 python 生态，统一 py3.14 |
| `misc` | iseq + primer3 + BioinfTools + bbmap(jdk25) + kmeriaenv + **sratoolkit_v.2.5.7**(收编) | jdk25 统一；死引用 5 顺带解决 |
| `r` | Rqtl(R4.5.3) + rMVP(R4.5.3) + WGCNA(R4.4.3) | R 4.4.3 vs 4.5.3 统一（WGCNA 的 bioconductor 依赖需重装验证） |
| `busco` | BUSCO_v.6.0.0(jdk23, 267 pkgs) | 保留独立（内含 augustus+blast+hmmer+mafft 全套，并入 annot 会造成工具版本重复） |

### Tier 2 —— 重装升级后并入|Upgrade Then Merge

| 现状|Now | 升级到|Upgrade | 并入域|Into |
|---|---|---|---|
| braker_v.3.0.8（py3.7.12, jdk25） | py3.13 | annot |
| fanc_v.0.9.23b（py3.7.12） | py3.13 | hic |
| deeptmhmm_v.1.0（py3.8.20） | py3.13 | protein |
| DeepBSA（py3.7.12, R4.2.1） | R4.5 | r |
| HiC-Pro_v3.1.0 + SubPhaser（同为 py3.8 + R4.0.3 生态） | 合成一个 `hic-legacy` | （合并为一个，不再升级） |

### Tier 3 —— 例外保留|Keep As-Is

| 环境|Env | 原因|Reason |
|---|---|---|
| qiime_v.2024.10.1（762 pkgs, py3.10, R4.3.3, jdk22） | qiime2 生态锁死，谁碰谁炸 |
| picrust_v.2.6.3（235 pkgs, py3.12, R4.5.3） | 依赖 qiime2 版本对应关系 |
| EDTA_v.2.2.2（428 pkgs, py3.12, R4.3.3） | 巨无霸，单独观察，稳定后并入 repeat |
| EGAPx_v.0.4.0-alpha（96 pkgs） | 完整 pipeline 应用环境 |
| cphasing（394 pkgs, py3.12.0） | 巨无霸，单独观察，稳定后并入 asm |
| jcvi_v.1.5.7（354 pkgs, R4.5.2） | 巨无霸 R+python 生态，作 viz anchor 观察 |
| telocomp（105 pkgs, 含 flye + minimap2 + jdk25） | 杂货铺：flye 死引用根源，asm 域建成后退役观察 |
| rnaseq_val（291 pkgs, py3.10, jdk25） | 大环境，并入 rna 需单独验证 |
| singularity_v.3.8.7 / base | 基础设施 |

### 目标数字|Target Numbers

| 阶段|Stage | 受管环境数|Managed Envs |
|---|---|---|
| 现状|Current | 101（代码引用） |
| 修死引用 + Tier1 完成后|After Tier 1 | **14 域 + 9 例外 ≈ 23** |
| Tier2 重装完成后|After Tier 2 | **15 域 + 8 例外 ≈ 23**（braker/fanc/deeptmhmm/DeepBSA 并入，hic-legacy +1） |
| 观察期后（telocomp/cphasing/jcvi/EDTA 退役或并入）|Steady state | **≈ 19~20** |

> ⚠️ 诚实结论：压到 **15 个以内不现实**——py3.7/3.8 legacy 栈（7 个）、qiime2 锁死生态、EDTA/cphasing/jcvi 三个 350+ 包巨无霸是硬约束。
> 强行全并的代价是每次升级解依赖像拆炸弹，一崩全崩。**19~20 个受管环境**是数据支持的最优平衡点。


## 三-B、求解预检实测结果|Solver Preflight Results（2026-08-16 超算实测）

> 全部 15 个域用 `mamba create --dry-run` 实测，**无一真实依赖冲突**。
> 预检脚本：`scripts/dryrun_conda_domains.sh`（本次由 AI 通过 ssh 在 online 节点直跑，未写超算任何文件）

| 域|Domain | 结果|Result | 包数|Pkgs | 下载|Download | 备注|Note |
|---|---|---|---|---|---|
| align（jdk17 钉版） | OK | 104 | ~647MB | GATK 4.6.2.0 + openjdk 17 可解 |
| align（jdk25 默认版） | OK | 104 | ~596MB | 求解器默认选 jdk25，两者皆可 |
| pop | OK | 253 | - | py3.13.15 + R4.5.2 共存无冲突 |
| asm | OK | 276 | 908MB | 需包名 `purge_dups`（非 purge-dups）；merqury 拉大体积 |
| hic | OK | 161 | 439MB | haphic 需 `-c bioconda-legacy`、juicer 需 `-c hcc` |
| annot | OK | 336 | - | 求解器拖入 python 2.7.15 + R 4.0.3（pasa/augustus 依赖），建好后需功能验证 |
| repeat | OK | 221 | - | py3.12.13 |
| rna | OK | 182 | - | py3.13.15 + R4.6.1 |
| protein | OK | 311 | 823MB | signalp6/phobius/tmhmm 需 `-c predector` 渠道（授权软件包装包） |
| phylo | OK | 128 | - | kakscalculator2/raxml-ng 均在此域可解（死引用修复顺带验证通过） |
| pan | OK | 184 | - | py3.11.15 |
| viz | OK | 107 | - | py3.13.15 |
| misc | OK | 127 | - | openjdk25 |
| r | OK | 209 | - | R4.5.3 |

**三个「假冲突」的真相|The 3 "conflicts" debunked:**

1. **asm**：`purge-dups` 是错误包名，bioconda 实际叫 `purge_dups`（下划线）
2. **hic**：`haphic` 只存在于 **bioconda-legacy** 渠道，`juicer` 只存在于第三方 **hcc** 渠道
3. **protein**：`signalp6`/`phobius`/`tmhmm` 只存在于 **predector** 渠道（三者均为需 license 的授权软件）

**结论|Conclusion:** 按功能域合并**求解器层面完全可行**，你担心的「有一类冲突」在这个分组下没有出现。
真正的风险不在求解，而在**运行时兼容**（annot 域拖入 python2、GATK+jdk25 未官测）——
建好环境后每个域跑一次最小冒烟测试即可确认。

## 四、执行顺序|Execution Order

1. **修死引用（立即，5 条代码级）**——不影响环境、不依赖合并，先止血
2. **超算求解预检**——对 14 个域逐一跑（不建环境，只求解）：
   ```bash
   mamba create -n dryrun_align --dry-run -c conda-forge -c bioconda \
     gatk4=4.6.2.0 bcftools bedtools samtools freebayes ...
   ```
3. **建域环境 + 写环境定义文件**——`envs/<domain>.yml` 进仓库，可复现
4. **代码层改造**——`common/` 增加 工具→域环境 映射表；模块默认路径改为 `get_tool_path()` 查表；`get_conda_env()` 加旧环境名回退（过渡期不破）
5. **模块分批切换**——先切共享度最高的（Population_genetics 15 模块、GATK 16 模块），后切单模块
6. **退役观察**——版本分裂对（JCVI_v.1.5.6 vs v.1.5.7、LTR_retriever_v.3.0.1 vs v.3.0.4、Orthofinder_v.3.0.1b1 vs v.3.1.5、kakscalculator vs kakscalculator2、wgdi_v.0.74 vs v.0.75、vg_v.1.7.0 vs v.1.67.0、EDTA vs EDTA_v.2.2.2 vs edta_v.2.3.0）旧名字停用

## 五、孤儿环境处理原则|Orphan Env Policy

- **113 个孤儿环境一律不删**：conda 安装在共享路径 `/share/org/YZWL`，其他同事可能直接使用
- 其中约 11 个其实是代码引用的「版本分裂对」（见上节），随合并退役
- 其余约 102 个（aria2/checkm/quast/salmon/trinity 等）不属于 biopytools 管理范围，保持现状
- 未来新模块优先复用域环境，不再新建专属环境

## 六、分工|Division of Work

| 谁|Who | 做什么|What |
|---|---|---|
| 你（超算） | 跑 14 个域的 `mamba create --dry-run` 求解预检，把失败域报回来 |
| 我（代码） | 修 5 条死引用 → 建 common 映射表 → 模块分批切换 |

## 七、代码迁移进度|Code Migration Progress（2026-08-16）

| 批次|Batch | 内容|Content | 提交|Commit |
|---|---|---|---|
| 2A | common/env_map.py + get_domain_tool_path()（域解析+旧路径回退）+ 8 模块 30 处调用点切换 + 8 个单测 | 333eb8f |
| 2B | 迁移脚本 scripts/migrate_env_paths.py 按工具级映射批量替换 103 个文件 219 处默认路径（存在性核验，缺工具自动跳过） | 81eb430 |

**迁移后代码中的引用分布|Post-migration refs:**

- **约 208 处指向 14 个域环境**（目标状态）
- **约 40 处 Tier3 例外**（singularity/cphasing/psvcp/picrust/qiime/HiC-Pro/juicer/gctb 等）——有意保留
- **约 40 处 Tier2 legacy**（fanc/plothic/genomesyn2/alignoth/Rqtl/transdecoder/deeptmhmm/DeepBSA 等）——工具未入域或需升级验证
- **4 处 signalp_v.3.0b**——黑名单（SignalP3 与 SignalP6 同名不同软件）
- **3 处 xxx 占位符**——nlr_annotator 的 Java 注释欠账

**迁移中发现的遗留问题|Issues found during migration:**

1. adamixture_v.1.0.2/bin/adamixture 引用（4 处）：二进制名疑似笔误（bioconda 装的是 admixture），需超算验证旧环境 bin 名
2. 域环境待补工具（后续由用户安装后即可再迁移）：minimap2→align、merqury、samblaster、matlock、TransDecoder、RAiSD/xpclr、mstmap、ltr_harvest/finder_parallel、panman
3. 回归测试：231/231 通过（含新 8 个单测与受影响断言更新）

