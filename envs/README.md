# 域环境定义|Domain Environment Definitions

> 按功能域合并 conda 环境的定义文件（Phase 1 产物）。
> 求解预检全部通过（见 `docs/conda_merge_plan.md` 三-B 节）；Tier2/Tier3 例外环境不在此目录。


## 旧环境备份|Legacy Env Backups

历史 214 个旧环境的导出快照保存在本地 `conda_envs_backup/`（git 忽略，不入库），
仅用于回溯与故障排查。新旧环境对应关系见 `docs/env_migration_map.md`。

> 曾入库的 `conda_env/` 目录已于 2026-08-16 移除（与 `conda_envs_backup/` 完全重叠，backup 为其超集）。
## 创建方式|How to Build

```bash
# 单个域|Single domain
mamba env create -f envs/align.yml

# 批量(跳过已存在)|Batch (skips existing)
bash scripts/build_domains.sh
```

## 域清单|Domain List

| 域|Domain | 吸收的旧环境|Absorbed Old Envs | 特殊渠道|Special Channels | 备注|Note |
|---|---|---|---|---|
| align | GATK_v.4.6.2.0, bcftools_v.1.22, sv_calling, freebayes, Genome_dedup | - | jdk17 钉版 |
| pop | Population_genetics, selective_sweep, adamixture_v.1.0.2, treemix_v.1.13, pixy_v.2.0.0, poplddecay_v.3.43 | - | gctb 自编译, 不入 env |
| asm | canu_v.2.3, hifiasm_v.0.25.0, kmc_v.3.2.4, K-mer, merqury_v.1.3, purge_dups_v.1.2.6, genomescope_v.2.0.1, genomescope2_v.2.1.0, tidk_v.0.2.65, getorganelle_v.1.7.71, spades_v.4.3.0 | - | 包名 `purge_dups`(下划线) |
| hic | haphic, pairtools_v.1.1.3, yahs_v.1.2.2, juicer_v.1.6 | bioconda-legacy, hcc | haphic 在 legacy、juicer 在 hcc |
| annot | Augustus_v.3.5.0, pasa_v.2.5.3, agat_v.1.7.0, gffcompare_v.0.12.10, miniprot_v.0.18, transdecoder_v.5.5.0, eggnog-mapper_v.2.1.15, orthofinder_v.3.1.5, Blast_v.2.16.0, genometools_v.1.6.5 | - | 建后跑冒烟验证(python2 被拖入) |
| repeat | repeatmodeler_v.2.0.7, repeat_identiy, ltr_retriever_v.3.0.1, ltr_harvest_parallel_v.1.2, ltr_finder_parallel_v.1.3 | - | |
| rna | RNA_Seq, RSeQC_v.5.0.4 | - | |
| protein | signalp6, resistify_v.1.3.0, needle_v.1.0.3(emboss), meme_v.5.5.9, phobius_v.1.0.1, tmmhmm_v.2.0c | predector | 授权软件, 建后按旧环境方式补 license/模型 |
| phylo | iqtree_v.3.0.1, mafft_v.7.525, trimal_v.1.5.0, newick_utils_v.1.6, wgdi_v.0.75, kakscalculator2_v.2.0.1, raxml-ng(原在 orthofinder) | - | 顺带修复 2 条死引用 |
| pan | pggb_v.0.7.4, vg_v.1.67.0, kmtricks_v.1.5.1, kmindex_v.0.6.0, panman_v.0.1.4, pan-blocks, swave_v.1.2 | - | |
| viz | samplot_v.1.3.0, pycirclize_v.1.10.1 | - | jcvi(354pkgs 巨无霸)观察期不并 |
| misc | iseq_v.1.9.8, primer3_v.2.6.1, bbmap_v.39.81, BioinfTools(seqkit) | - | |
| r | Rqtl, rMVP, WGCNA_v.1.73 | - | DeepBSA(py3.7)在 Tier2 |
| busco | BUSCO_v.6.0.0 | - | 独立保留(内置全套注释工具) |

## 冒烟验证|Smoke Tests（建好后每个域跑一遍）

```bash
M=/share/org/YZWL/yzwl_lixg/miniforge3/bin/mamba
$M run -n align  gatk --version
$M run -n pop    plink --version
$M run -n asm    canu --version
$M run -n hic    haphic --version
$M run -n annot  augustus --version
$M run -n repeat RepeatMasker -v
$M run -n rna    hisat2 --version
$M run -n protein signalp6 --help
$M run -n phylo  iqtree --version
$M run -n pan    vg version
$M run -n viz    samplot --help
$M run -n misc   seqkit version
$M run -n r      Rscript -e 'library(qtl)'
$M run -n busco  busco --version
```

## 注意事项|Caveats

1. **旧环境保留**：代码切换到新域环境之前，旧环境（101 个）一律不动，双轨运行过渡。
2. **授权软件**：protein 域的 signalp6/phobius/tmhmm 建好后，需按现有 signalp6/phobius_v.1.0.1/tmmhmm_v.2.0c 环境的方式补装 license 文件。
3. **渠道依赖**：hic/protein 依赖 hcc/predector/bioconda-legacy 渠道名，若报 unknown channel，检查超算 `~/.condarc` 是否已定义（旧环境创建时已用过，应当已配好）。
4. **版本漂移**：yml 未钉死全部版本（除 gatk4/jdk/bcftools/busco），随渠道更新会漂移；全部锁定可在建好后 `mamba env export` 生成精确版本 yml 再入库。

## 冒烟终态|Final Smoke Status（2026-08-16）

**14 个域全部通过**。验证方式：`conda run -n <域> <工具> --version`（直接调 bin 二进制会因 PATH 缺 env/bin 误报，勿用）。

| 域 | 关键版本 | 状态 |
|---|---|---|
| align | gatk 4.6.2.0 / samtools 1.23.1 / bcftools 1.22 | ✅ |
| pop | plink 1.9 / pixy / vcftools | ✅ |
| asm | canu 2.3 / hifiasm 0.25 / jellyfish 2.2.10 / spades 4.3 | ✅ |
| hic | haphic 1.0.7 / juicer 1.6 / pairtools / yahs | ✅ |
| annot | augustus 3.5.0 / orthofinder 3.1.5 / gt 1.6.5 / py3.11 | ✅ |
| repeat | RepeatMasker 4.2.4 / LTR_retriever | ✅ |
| rna | hisat2 2.2.3 / stringtie 3.0.3 / fastp 1.3.6 | ✅ |
| protein | meme 5.5.9 / emboss / resistify / signalp6 / phobius / tmhmm | ✅（授权软件已注册） |
| phylo | iqtree 3.1.3 / mafft 7.526 / raxml-ng 2.0.2 / KaKs_Calculator | ✅ |
| pan | vg 1.63.1 / pggb / minigraph 0.21 / mummer | ✅ |
| viz | samplot 1.3.0 / pycirclize 1.10.1 | ✅ |
| misc | seqkit 2.13.0 / axel / iseq / primer3 / bbmap | ✅ |
| r | R 4.5.3 + qtl + WGCNA 加载成功 | ✅ |
| busco | BUSCO 6.0.0 | ✅ |

**建好后还差的「数据/license 类」准备**（与旧环境相同，非环境缺陷）：

1. protein 域：signalp6/phobius/tmhmm 补 license 与模型（按旧 signalp6/phobius_v.1.0.1/tmmhmm_v.2.0c 环境的方式）
2. annot 域：eggnog-mapper 数据库 `conda run -n annot download_eggnog_data.py`（按需）
3. asm 域：tidk 的 clade 数据库 `conda run -n asm tidk build`（按需）



## 授权软件注册手册|License Registration Guide（protein 域）

> 三个授权软件在超算上的 tarball 位置已核实（2026-08-16），注册命令可直接照抄。
> 机制：predector 渠道的包装包 + `*-register <tarball>` 把授权内容解进环境。

| 软件|Tool | tarball 路径|Tarball Path | 注册命令|Register Command |
|---|---|---|---|
| SignalP6 | `~/software/SignalP_v.6.0/signalp-6.0i.fast.tar.gz`（6.0i fast 模型，与旧环境一致） | `conda run -n protein signalp6-register ~/software/SignalP_v.6.0/signalp-6.0i.fast.tar.gz` |
| Phobius | `~/software/Phobius/phobius101_linux.tgz` | `conda run -n protein phobius-register ~/software/Phobius/phobius101_linux.tgz` |
| TMHMM | `~/software/tmhmm/tmhmm-2.0c.Linux.tar.gz` | `conda run -n protein tmhmm2-register ~/software/tmhmm/tmhmm-2.0c.Linux.tar.gz` |

注册后验证|Verify after register:

```bash
/share/org/YZWL/yzwl_lixg/miniforge3/bin/conda run -n protein signalp6 --help | head -3
/share/org/YZWL/yzwl_lixg/miniforge3/bin/conda run -n protein phobius.pl 2>&1 | head -2
/share/org/YZWL/yzwl_lixg/miniforge3/bin/conda run -n protein tmhmm -h 2>&1 | head -3
```

注意|Notes:
- `signalp6-register` 会提示「文件名不匹配 signalp-6.0h.fast.tar.gz」——正常，旧环境就是这么注册的（包版本 6.0h，模型用 6.0i fast），继续即可
- 旧环境的 license 文件在 `~/software/SignalP_v.6.0/signalp-6.0_license.txt` 和 `~/software/tmhmm/tmhmm-2.0c_license.txt`，未动

## 按需数据准备|On-Demand Data（旧环境同样没有，非回归）

```bash
# eggnog-mapper 数据库（约 20-40GB，做 eggnog 注释时才下）
/share/org/YZWL/yzwl_lixg/miniforge3/bin/conda run -n annot download_eggnog_data.py

# tidk clade 数据库（几 MB，找端粒时才建）
/share/org/YZWL/yzwl_lixg/miniforge3/bin/conda run -n asm tidk build
```

## 建环境踩坑实录|Build Pitfalls（2026-08-16 实测）

| # | 坑|Pitfall | 解决|Fix |
|---|---|---|---|
| 1 | annot 域被 `pasa`/老 `genometools` 包拖入 python2.7 生态（augustus 退到 3.2.2、orthofinder 退到 2.5.4、scipy 1.2.0） | pasa 留在旧 pasa_v.2.5.3 环境；genometools 包名改为 `genometools-genometools`（老包名是 py2.7 版） |
| 2 | asm 域 `jellyfish` 解到 conda-forge 的 1.2.1（python 字符串库，无 CLI） | 装 `jellyfish=2.2.10`（bioconda），命令需 bioconda 渠道在前：`mamba install -n asm -c bioconda -c conda-forge jellyfish=2.2.10` |
| 3 | pan 域缺 `minigraph`（yml 漏写） | 补装：`mamba install -n pan -c conda-forge -c bioconda minigraph` |
| 4 | misc 域 `iseq` 缺下载工具（axel/pigz/sra-tools/aspera-cli/wget） | 补装：`mamba install -n misc -c conda-forge -c bioconda axel pigz sra-tools aspera-cli wget` |
| 5 | busco/embossversion 直接调 bin 报错（PATH 缺 env/bin） | 假警报：经 `conda run -n <env>` 调用正常（仓库代码本来就是 conda run 包装） |
| 6 | `mamba run` 对 `--` 参数解析有 bug | 仓库代码用 `conda run`（验证正常），不受影响 |

