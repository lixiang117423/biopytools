# SplitsTree6 建网建树 | splitstree6

给一份群体基因型数据(VCF)或者比对好的序列,自动构建进化网络/进化树——比常规建树多一层"网"的表达,能直观展示重组、杂交、水平转移等树状结构解释不了的演化信号,结果同时输出 Newick/Nexus/GML 等多种格式。

## 功能概述 | Overview { #overview }

- 封装官方 **SplitsTree6** 命令行(workflow-run),headless 运行(Xvfb 虚拟显示),无需图形界面
- **默认输入 VCF**:模块自动将样本基因型转成 p-distance 距离矩阵喂给 SplitsTree6;也支持 fasta/nexus/phylip/newick 等格式直传
- 内置 **NeighborNet**(网状网络)+ Hamming Distance 标准工作流;可自定义 `.stree6` 工作流
- 默认导出 **Newick + Nexus + GML**,可选 PlainText/Phylip/FastA/Clustal 等全格式
- `-n` 参数指定导出节点(Splits 网络/距离矩阵/Taxa 等),灵活取工作流中间产物
- `software_versions.yml` 记录版本与参数;日志三分离(.log/.out.log/.err.log)

## 快速开始 | Quick Start { #quick-start }

```bash
biopytools splitstree6 -i variants.vcf -o splitstree6_out/
```

读 VCF → 转 p-distance → NeighborNet 建网 → 输出 Newick/Nexus/GML 三种格式到 `splitstree6_out/01_workflow/`。

## 零基础概念速览 | Concepts in plain words { #concepts }

| 术语<br>Term | 通俗解释<br>In plain words |
|------|------|
| 进化网络(网络图) | 比树更泛化的演化图:树上每个点只有一个"爹",网络允许多个——正好用来表达重组、杂交、水平基因转移这类树装不下的事件 |
| NeighborNet | SplitsTree 家族招牌算法:先算 pairwise 距离,再拆成一组相容"splits"(分割),画出网状图 |
| local bootstrap | 工作流里输出的分支支持度;没有显式 bootstrap 时看 fit 值(R² 类指标,Nexus 头部 `fit=99.9` 表示拟合极好) |
| p-distance | 两两样本之间不同位点的比例,最朴素的遗传距离 |
| splits | 把所有样本分成两堆的一种切法;NeighborNet 的输出就是一堆"合理切法"的集合 |

## 输入 | Input { #input }

**首选:VCF 变异文件**(`.vcf` / `.vcf.gz`,二等位 SNP 效果最佳)

- 模块读取基因型列,按样本对计算 p-distance,生成 SplitsTree6 的 CSV 距离矩阵输入
- 缺失基因型(`.` 或 `./.`)在该位点比较中跳过,不中断运行
- 转换后的矩阵保留在 `01_workflow/*.distances.csv`,可与下游复用

**同样支持 SplitsTree6 原生格式**(直接透传 `-f`):FastA(多序列比对)、Nexus、Phylip、Newick、Stockholm、MSF、GML、CSV 距离矩阵

## 参数说明 | Parameters { #parameters }

**通俗理解|In plain words:** `-e/--export-formats` 决定你拿到什么格式的文件(默认 Newick+Nexus+GML 三件套);`-w/--workflow` 是进阶口——想换算法(如 Consensus Network、Hybridization Network)时,用 SplitsTree6 GUI 保存一个 .stree6 传入即可;其他情况一律不动。

<!-- BEGIN PARAMS:auto -->

## 参数速查 | Parameter reference

> 本表由 `scripts/gen_docs_params.py` 从 CLI 定义自动生成,勿手改|Auto-generated from CLI definitions; do not edit by hand

### 命令行参数 | CLI options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入数据 .vcf/.vcf.gz(自动转距离矩阵)或其他 SplitsTree6 格式｜Input: VCF (auto-converted) or any SplitsTree6-readable file |
| `-o, --output-dir` | `./splitstree6_output` |  | 输出目录(默认./splitstree6_output)｜Output directory |
| `-e, --export-formats` | — |  | 输出格式,逗号分隔(默认 Newick,Nexus,GML)｜Export formats comma-separated (default Newick,Nexus,GML)。可选｜valid: Newick,Nexus,GML,PlainText,Phylip,FastA,Clustal |
| `-w, --workflow` | — |  | 自定义 .stree6 工作流(默认内置 NeighborNet 模板)｜Custom .stree6 workflow (default built-in NeighborNet template) |
| `-n, --node-name` | `Splits` |  | 导出节点名(默认 Splits 网络节点)｜Node to export (default Splits) |
| `--input-format` | `` |  | 指定输入格式(默认自动识别)｜Input format (default auto-detect) |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Threads (default 12) |
| `--tools-dir` | — |  | splitstree6-tools 目录(jars 所在)｜splitstree6-tools dir |
| `--xvfb-path` | — |  | Xvfb 路径｜Xvfb path (JavaFX requires a display) |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

### 模块直调参数 | Direct invocation options

| 参数 | 默认值 | 类型 | 说明 |
|------|--------|------|------|
| `-i, --input` | 必填 |  | 输入数据:.vcf/.vcf.gz(默认自动转距离矩阵)或其他 SplitsTree6 支持格式(fasta/nexus/phylip/newick/…)｜Input: VCF (auto-converted) or any SplitsTree6-readable format |
| `-o, --output-dir` | `./splitstree6_output` |  | 输出目录(默认./splitstree6_output)｜Output directory |
| `-e, --export-formats` | — |  | 输出格式,逗号分隔(默认 Newick,Nexus,GML)｜Export formats comma-separated (default: Newick,Nexus,GML)。可选｜valid: … |
| `--input-format` | `` |  | 指定输入格式(默认自动识别)｜Input format (default auto-detect) |
| `-w, --workflow` | `` |  | 自定义 .stree6 工作流(默认内置 NeighborNet 模板)｜Custom .stree6 workflow (default: built-in NeighborNet template) |
| `-n, --node-name` | `Splits` |  | 导出节点名(默认 Splits 网络节点)｜Node to export (default: Splits) |
| `-t, --threads` | `12` | int | 线程数(默认12)｜Threads (default 12) |
| `--tools-dir` | — |  | splitstree6-tools 目录(jars 所在)｜splitstree6-tools dir |
| `--xvfb-path` | — |  | Xvfb 路径(JavaFX 需要虚拟显示)｜Xvfb path |
| `--log-level` | `INFO` |  | 日志级别(默认INFO)｜Log level (default INFO) |

<!-- END PARAMS:auto -->

## 输出 | Output { #output }

```
output/
├── 00_pipeline_info/
│   └── software_versions.yml     # 版本/参数/耗时
├── 01_workflow/
│   ├── <workflow 文件名>.stree6   # 本次使用的工作流副本(留档)
│   ├── *.distances.csv           # (VCF 输入时)p-distance 距离矩阵
│   ├── Splits.newick             # Newick 格式网络/树
│   ├── Splits.nexus              # Nexus 格式(含 splits block 最完整)
│   └── Splits.gml                # GML 图格式(Gephi/Cytoscape 可读)
└── 99_logs/
    ├── splitstree6.log
    └── run_<format>.log          # 每个导出格式的运行日志
```

## 结果解读 | Interpreting Results { #interpreting-results }

- **`Splits.nexus`**:信息最全的产物,`BEGIN SPLITS` 块含全部 splits 与权重;文件头部的 `fit=xx.x` 是 NeighborNet 对距离数据的拟合度,**>90 说明网络能很好表达数据,<80 要小心结论**
- **`Splits.newick`**:把网状结构投影成树的近似,兼容所有吃 Newick 的下游工具;注意:网里的"平行边"在 Newick 里会有损
- **平行边标记 `<x|y`**:Newick 中 `<1|tax6` 这种表示该节点有多重身份(splits 重叠),是网络的正常特征,不是解析错误
- **`*.distances.csv`**:转换出的距离矩阵,人工可查;p-distance ≈ 0 的两个样本过近时要怀疑测序混样
- **GML**:给 Cytoscape/Gephi 画图用;输出为 0 字节说明该网络不含可视化坐标(正常现象,用 Nexus 即可)

## 参数选择建议 | Parameter Guidance { #parameter-guidance }

| 场景<br>Scenario | 推荐参数<br>Recommended |
|------|------|
| VCF 直接建网 | 默认即可(`-i xxx.vcf`) |
| 多格式都要 | `-e Newick,Nexus,GML,PlainText` |
| 只要树不要网 | 自定义 workflow(去掉 Neighbor Net 节点)后 `-w` 传入 |
| 换算法(Hybridization 等) | GUI 里搭好后 File→Save 存 `.stree6`,`-w` 传入 |
| 大 VCF(>10k 位点) | 先抽样位点再跑;内存需求随位点数线性增长 |

**通俗理解|In plain words:** 平时一条命令完事;要特殊玩法(GUI 定制流程)才碰 `-w`;VCF 太大就先减位点。

## 依赖 | Dependencies { #dependencies }

| 依赖<br>Dependency | 说明<br>Note |
|------|------|
| splitstree6-tools 完整包(lib/jars) | 必需,默认 `~/software/splitstree6/splitstree6-tools-bin`(可用 `--tools-dir`/`SPLITSTREE6_TOOLS_DIR` 指定);已预装  |
| Xvfb | JavaFX 需虚拟显示;conda `xvfb_test` 环境已含(可用 `--xvfb-path`/`XVFB_PATH` 指定) 已确认可用 |
| java 17+(任意发行版) | 系统 PATH 可见即可(`~/.local/bin/java` 已确认) 已确认可用 |
| libcrypto.so.10 | Xvfb 运行需要;已放 `~/tmp/x11libs`(LD_LIBRARY_PATH 自动注入) |

## 常见问题 | FAQ { #faq }

- **能直接输出 VCF 吗?** 不能——SplitsTree6 只有读入没有 VCF 导出器(源码 ExportManager 无此 writer)。本模块"默认 VCF"指**输入**侧;输出走 Newick/Nexus 等标准格式
- **为什么跑的时候有 Xvfb?** SplitsTree6 主程序基于 JavaFX 启动,X 显示系统不可缺;Xvfb 提供"假屏幕"让它安心当命令行工具干活
- **报 `Unable to open DISPLAY`?** 旧版 bug 已规避:模块强制注入 DISPLAY=:99 并自动拉起 Xvfb;手动排查时先 `pgrep -f Xvfb` 看进程是否活着
- **GML 为什么是 0 字节?** `-n Splits` 导出的 network 无坐标属性,GML 仅落坐标类内容;请改用 Nexus
- **想换算法怎么办?** 用 SplitsTree6 桌面版搭好流程保存为 `.stree6`,`-w my_workflow.stree6` 一键替换内置模板
