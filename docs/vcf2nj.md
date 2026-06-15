# VCF构建NJ树 | VCF Neighbor-Joining Tree (vcf2nj)

**基于VCF2Dis计算距离矩阵并构建邻接(NJ)系统发育树, 支持外群重根化 | Distance-matrix and neighbor-joining tree from VCF**

## 功能概述 | Overview

vcf2nj 模块封装了 VCF2Dis 和 newick_utils, 用于从VCF变异文件计算样本间遗传距离矩阵, 并使用邻接(Neighbor-Joining, NJ)算法构建系统发育树。相比最大似然法, NJ算法速度快, 适合群体水平的大样本快速建树场景。模块支持指定外群进行重根化(rerooting), 也支持直接输入已有距离矩阵跳过VCF2Dis步骤。

## 快速开始 | Quick Start

```bash
# 从VCF构建NJ树
biopytools vcf2nj -i wild.snp.vcf -o output_dir -p wild_snp

# 指定外群重根化
biopytools vcf2nj -i wild.snp.vcf -o output_dir -p wild_snp --outgroup outgroup_sample

# 已有距离矩阵直接建树
biopytools vcf2nj -d distance.mat -o output_dir -p my_tree --skip-vcf2dis
```

## 参数说明 | Parameters

### 输入参数(二选一) | Input (one required)

| 参数 | 描述 |
|------|------|
| `-i, --input` | 输入VCF文件路径 |
| `-d, --distance-matrix` | 已有距离矩阵文件路径(配合`--skip-vcf2dis`) |

### 常用可选参数 | Common Options

| 参数 | 默认值 | 描述 |
|------|--------|------|
| `-o, --output` | `.` | 输出目录 |
| `-p, --prefix` | `vcf2nj` | 输出文件前缀 |
| `-t, --tree-output` | `None` | 系统发育树输出文件路径 |
| `--outgroup` | `None` | 外群样本标签(逗号分隔) |
| `--nw-reroot-path` | `~/miniforge3/envs/newick_utils_v.1.6/bin/nw_reroot` | nw_reroot程序路径 |
| `--vcf2dis-path` | `VCF2Dis` | VCF2Dis程序路径 |
| `-w, --working-dir` | `.` | 工作目录 |
| `--skip-vcf2dis` | `False` | 跳过VCF2Dis步骤(需配合`-d`) |

(运行 `biopytools vcf2nj -h` 查看完整参数列表)

## 输出 | Output

```
output_dir/
├── {prefix}.dis.mat        # VCF2Dis生成的距离矩阵
├── {prefix}.nwk            # NJ树(Newick格式)
└── {prefix}.rerooted.nwk   # 重根化后的树(若指定外群)
```

## 依赖 | Dependencies

- **VCF2Dis**: VCF距离矩阵计算 (https://github.com/BGI-shenzhen/VCF2Dis)
- **newick_utils (nw_reroot)**: 树重根化 (https://github.com/tjunier/newick_utils)
- **PHYLIP neighbor**: NJ建树(VCF2Dis内部调用)

## 引用 | Citation

- Liu H., et al. VCF2Dis: a simple and fast distance matrix calculation software for VCF format. (github.com/BGI-shenzhen/VCF2Dis)
- Saitou N., Nei M. (1987) The neighbor-joining method: a new method for reconstructing phylogenetic trees. Molecular Biology and Evolution. 4(4):406-425.

## 相关链接 | References

- [项目主页](https://github.com/lixiang117423/biopytools)
