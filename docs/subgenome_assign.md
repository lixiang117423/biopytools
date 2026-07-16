# 亚基因组归属 | Subgenome Assignment

**异源多倍体基因组的亚基因组归属, 经亲本比对判定每条染色体的来源 | Assign subgenomes in an allopolyploid by aligning chromosomes against parental haplotypes**

## 功能概述 | Overview

subgenome_assign 模块将异源多倍体基因组的每条染色体归属到来源亲本。流程:

1. **检查依赖** — minimap2 / samtools
2. **合并亲本 hap FASTA** — 各亲本的多个单倍型合并
3. **minimap2 比对** — 目标基因组 vs 各亲本(`-x asm10` 等预设)
4. **归属判定** — 解析 PAF, 按比对一致性将每条染色体判给最优亲本, 低于 `--min-conf` 标记 `LOW_CONFIDENCE`
5. **拆分 FASTA** — 按归属结果把目标基因组拆成各亲本的 FASTA

## 快速开始 | Quick Start

```bash
# 二倍体亲本(每亲本2个单倍型)
biopytools subgenome-assign -i polyploid.fa \
  --parent Ca:Ca_hap1.fa,Ca_hap2.fa \
  --parent Ch:Ch_hap1.fa,Ch_hap2.fa \
  -o out_dir/

# 只做归属判定, 不拆分 FASTA
biopytools subgenome-assign -i target.fa --parent A:A.fa --parent B:B.fa -o out/ --no-split
```

## 参数说明 | Parameters

### 必需参数 | Required

| 参数 | 描述 |
|------|------|
| `-i, --target` | 目标多倍体基因组 FASTA |
| `--parent` | 亲本 `NAME:hap1.fa,hap2.fa`(可重复指定多个亲本) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output-dir` | `./subgenome_assign_output` | 输出目录 |
| `--preset` | `asm10` | minimap2 预设(asm5/asm10/asm20/asm25) |
| `-t, --threads` | `12` | 线程数 |
| `--minimap2-secondary` | `False` | 保留次要比对(默认 `--secondary=no`) |
| `--min-conf` | `0.65` | 置信度阈值, 低于此值标记 LOW_CONFIDENCE |
| `--no-split` | `False` | 不输出拆分的 FASTA |
| `--no-keep-unassigned` | `False` | 不输出未归属染色体的 FASTA |
| `--minimap2-path` | `~/miniforge3/envs/cphasing/bin/minimap2` | minimap2 路径 |
| `--samtools-path` | `~/.local/bin/samtools` | samtools 路径 |

(运行 `biopytools subgenome-assign -h` 查看完整参数列表)

## 输出 | Output

- 归属目录: `subgenome_assignment.tsv`(每条染色体 → 归属亲本 + 置信度)
- 拆分 FASTA 目录: 按亲本拆分的基因组 FASTA
- `00_pipeline_info/software_versions.yml`: 软件版本与运行参数
- `00_pipeline_info/pipeline_params.yaml`: 含归属结果汇总(各亲本染色体数)
- `99_logs/`: 运行日志

## 依赖 | Dependencies

- **minimap2**: 目标 vs 亲本比对
- **samtools**: 构建 `.fai` 索引

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
