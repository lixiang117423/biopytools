"""
RAxML系统发育分析命令 | RAxML Phylogenetic Analysis Command
"""

import click
from ...raxml.main import RAxMLPhylogeneticAnalyzer


@click.command(short_help='RAxML系统发育分析工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
# 必需参数 | Required arguments
@click.option('--sequence-file', '-s',
              required=True,
              type=click.Path(exists=True),
              help='输入序列文件 (PHYLIP格式) | Input sequence file (PHYLIP format)')
@click.option('--output-name', '-n',
              required=True,
              type=str,
              help='输出文件名称 | Output file name')

# 模型参数 | Model parameters
@click.option('--model', '-m',
              default='GTRGAMMA',
              type=str,
              help='替换模型 | Substitution model (GTRGAMMA, PROTGAMMAWAG, etc.) (default: GTRGAMMA)')
@click.option('--categories', '-c',
              default=25,
              type=int,
              help='速率异质性类别数 | Number of rate heterogeneity categories (default: 25)')
@click.option('--likelihood-epsilon', '-e',
              default=0.1,
              type=float,
              help='似然优化精度 | Likelihood optimization precision (default: 0.1)')

# 算法参数 | Algorithm parameters
@click.option('--algorithm', '-f',
              default='d',
              type=str,
              help='算法类型 | Algorithm type (d=rapid hill-climbing, a=rapid bootstrap, etc.) (default: d)')
@click.option('--parsimony-seed', '-p',
              type=int,
              help='简约法随机种子 | Parsimony random seed')
@click.option('--runs', '-N',
              default='1',
              type=str,
              help='运行次数或bootstrap次数 | Number of runs or bootstrap replicates (default: 1)')

# Bootstrap参数 | Bootstrap parameters
@click.option('--bootstrap-seed', '-b',
              type=int,
              help='Bootstrap随机种子 | Bootstrap random seed')
@click.option('--rapid-bootstrap-seed', '-x',
              type=int,
              help='快速bootstrap随机种子 | Rapid bootstrap random seed')
@click.option('--bootstrap-convergence', '-I',
              type=click.Choice(['autoFC', 'autoMR', 'autoMRE', 'autoMRE_IGN']),
              help='Bootstrap收敛标准 | Bootstrap convergence criterion')
@click.option('--bootstop-threshold', '-B',
              default=0.03,
              type=float,
              help='Bootstrap停止阈值 | Bootstrap stop threshold (default: 0.03)')
@click.option('--bootstop-perms',
              default=100,
              type=int,
              help='Bootstrap停止检验次数 | Bootstrap stop test permutations (default: 100)')
@click.option('--print-bootstrap-trees', '-k',
              is_flag=True,
              help='输出带分支长度的bootstrap树 | Print bootstrap trees with branch lengths')

# 树参数 | Tree parameters
@click.option('--starting-tree', '-t',
              type=click.Path(exists=True),
              help='起始树文件 | Starting tree file')
@click.option('--constraint-tree', '-g',
              type=click.Path(exists=True),
              help='约束树文件 | Constraint tree file')
@click.option('--outgroup', '-o',
              type=str,
              help='外群名称 (逗号分隔多个) | Outgroup name(s) (comma-separated)')

# 性能参数 | Performance parameters
@click.option('--threads', '-T',
              default=88,
              type=int,
              help='线程数 | Number of threads (default: 88)')
@click.option('--memory-saving', '-U',
              is_flag=True,
              help='启用内存节省模式 | Enable memory saving mode')
@click.option('--ml-search-convergence', '-D',
              is_flag=True,
              help='启用ML搜索收敛标准 | Enable ML search convergence criterion')

# 高级参数 | Advanced parameters
@click.option('--random-starting-tree', '-d',
              is_flag=True,
              help='使用随机起始树 | Use random starting tree')
@click.option('--disable-rate-heterogeneity', '-V',
              is_flag=True,
              help='禁用速率异质性模型 | Disable rate heterogeneity model')
@click.option('--gamma-median', '-u',
              is_flag=True,
              help='使用GAMMA模型中位数 | Use median for GAMMA model')
@click.option('--disable-pattern-compression', '-H',
              is_flag=True,
              help='禁用模式压缩 | Disable pattern compression')

# 输出和工具参数 | Output and tool parameters
@click.option('--output-dir', '-w',
              default='./raxml_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./raxml_output)')
@click.option('--raxml-path',
              default='raxmlHPC-PTHREADS',
              type=str,
              help='RAxML程序路径 | RAxML program path (default: raxmlHPC-PTHREADS)')

# 质量控制参数 | Quality control parameters
@click.option('--no-seq-check',
              is_flag=True,
              help='跳过序列检查 | Skip sequence checking')
@click.option('--silent',
              is_flag=True,
              help='静默模式 | Silent mode')
def raxml(sequence_file, output_name, model, categories, likelihood_epsilon,
                algorithm, parsimony_seed, runs, bootstrap_seed, rapid_bootstrap_seed,
                bootstrap_convergence, bootstop_threshold, bootstop_perms, print_bootstrap_trees,
                starting_tree, constraint_tree, outgroup, threads, memory_saving,
                ml_search_convergence, random_starting_tree, disable_rate_heterogeneity,
                gamma_median, disable_pattern_compression, output_dir, raxml_path,
                no_seq_check, silent):
    """
    RAxML系统发育分析脚本 (模块化版本)
    
    基于RAxML的高性能最大似然法系统发育树构建工具，支持多种进化模型、
    Bootstrap分析和并行计算，适用于分子进化和系统发育研究。
    
    功能特点 | Features:
    - 多种进化模型支持 (DNA/蛋白质)
    - 快速Bootstrap分析和收敛检测
    - 高效的并行计算优化
    - 灵活的约束树和外群设置
    - 自动序列格式验证
    - 详细的分析日志记录
    - 模块化可扩展架构
    
    分析流程 | Analysis Pipeline:
    1. 序列文件格式验证和统计
    2. 进化模型参数设置
    3. 最大似然法树搜索
    4. Bootstrap置信度分析
    5. 结果文件整理和验证
    6. 系统发育分析报告生成
    
    应用场景 | Use Cases:
    - 分子系统发育重建
    - 进化关系分析
    - 基因家族进化研究
    - 物种分化时间估算
    - 分子钟分析
    - 保守性和选择压力分析
    
    示例 | Examples:
    
    \b
    # 基本系统发育分析
    biopytools raxml-phylo -s alignment.phy -n my_tree
    
    \b
    # 完整的bootstrap分析
    biopytools raxml-phylo -s data.phy -n phylo_bootstrap \\
        -m GTRGAMMA -f a -x 12345 -N 1000
    
    \b
    # 蛋白质序列分析
    biopytools raxml-phylo -s proteins.phy -n protein_tree \\
        -m PROTGAMMAWAG -T 64
    
    \b
    # 使用约束树分析
    biopytools raxml-phylo -s alignment.phy -n constrained_tree \\
        -g constraint.tre -T 88
    
    \b
    # 快速bootstrap收敛分析
    biopytools raxml-phylo -s sequences.phy -n auto_bootstrap \\
        -f a -m GTRGAMMA -x 12345 -I autoMRE -T 64
    
    \b
    # 大数据集内存优化分析
    biopytools raxml-phylo -s large_alignment.phy -n large_tree \\
        -m GTRGAMMA -U -H -T 88 --silent
    
    \b
    # 指定外群的完整分析
    biopytools raxml-phylo -s data.phy -n outgroup_tree \\
        -o "Outgroup1,Outgroup2" -f a -m GTRGAMMA -x 12345 -N 100
    
    进化模型说明 | Evolutionary Models:
    
    DNA序列模型:
    - GTRGAMMA: GTR+GAMMA (最常用)
    - GTRCAT: GTR+CAT (快速近似)
    - GTRCATI: GTR+CAT+Inv (含不变位点)
    - HKY85GAMMA: HKY85+GAMMA
    
    蛋白质序列模型:
    - PROTGAMMAWAG: WAG+GAMMA (推荐)
    - PROTGAMMALG: LG+GAMMA
    - PROTGAMMAJTT: JTT+GAMMA
    - PROTCATWAG: WAG+CAT (快速)
    
    算法类型详解 | Algorithm Types:
    
    -f d: 新分析 (默认)
    - 快速爬山法搜索最佳树
    - 适用于初步分析和快速结果
    
    -f a: 快速Bootstrap分析
    - 结合快速bootstrap和最佳树搜索
    - 推荐用于发表级别的分析
    
    -f b: 仅Bootstrap分析
    - 从给定的起始树进行bootstrap
    - 用于已有最佳树的置信度评估
    
    -f x: 最大似然搜索
    - 从随机起始树开始
    - 用于复杂数据集的深度搜索
    
    Bootstrap分析选项 | Bootstrap Analysis Options:
    
    固定次数Bootstrap:
    - -N 100: 执行100次bootstrap
    - -N 1000: 执行1000次bootstrap (推荐发表)
    
    自动收敛Bootstrap:
    - -I autoFC: 频率停止标准
    - -I autoMR: 多数规则停止标准  
    - -I autoMRE: 扩展多数规则 (推荐)
    - -I autoMRE_IGN: 忽略分支长度的扩展多数规则
    
    性能优化建议 | Performance Optimization:
    
    线程设置:
    - 小数据集 (<1000序列): 8-16线程
    - 中等数据集 (1000-5000序列): 32-64线程
    - 大数据集 (>5000序列): 64-88线程
    
    内存优化:
    - --memory-saving: 适用于内存受限环境
    - --disable-pattern-compression: 禁用可节省内存但影响速度
    - 大对齐序列建议使用高内存机器
    
    算法选择策略:
    - 快速预览: -f d (几分钟到几小时)
    - 发表质量: -f a -N 1000 (几小时到几天)
    - 高精度: -f a -I autoMRE (自动收敛)
    
    输出文件说明 | Output Files:
    
    主要结果文件:
    - RAxML_bestTree.{name}: 最佳树 (无支持值)
    - RAxML_bipartitions.{name}: 带bootstrap支持值的最佳树
    - RAxML_bootstrap.{name}: 所有bootstrap树
    - RAxML_info.{name}: 详细运行信息
    
    中间和日志文件:
    - RAxML_log.{name}: 运行日志
    - RAxML_result.{name}: 最终似然值结果
    - RAxML_parsimonyTree.{name}: 简约法起始树
    
    高级应用技巧 | Advanced Usage Tips:
    
    约束树分析:
    - 使用-g参数指定拓扑约束
    - 约束树应为Newick格式
    - 约束树可以是部分拓扑结构
    
    外群设置:
    - 使用-o参数指定外群
    - 多个外群用逗号分隔
    - 确保外群名称与序列标签完全匹配
    
    大数据集处理:
    - 考虑使用RAxML-NG替代标准RAxML
    - 预先检查序列对齐质量
    - 移除gap过多的序列或位点
    
    模型选择指导:
    - DNA数据: 通常使用GTRGAMMA
    - 蛋白质数据: WAG模型适用于大多数情况
    - 特殊数据: 可使用ProtTest或ModelTest选择最佳模型
    """
    
    try:
        # 处理特殊参数逻辑
        rate_het_model = not disable_rate_heterogeneity
        pattern_compression = not disable_pattern_compression
        
        # 直接创建分析器，避免参数重构的性能损失
        analyzer = RAxMLPhylogeneticAnalyzer(
            sequence_file=sequence_file,
            output_name=output_name,
            model=model,
            output_dir=output_dir,
            algorithm=algorithm,
            parsimony_seed=parsimony_seed,
            bootstrap_seed=bootstrap_seed,
            rapid_bootstrap_seed=rapid_bootstrap_seed,
            runs=runs,
            starting_tree=starting_tree,
            constraint_tree=constraint_tree,
            outgroup=outgroup,
            categories=categories,
            rate_het_model=rate_het_model,
            gamma_median=gamma_median,
            threads=threads,
            likelihood_epsilon=likelihood_epsilon,
            ml_search_convergence=ml_search_convergence,
            random_starting_tree=random_starting_tree,
            memory_saving=memory_saving,
            pattern_compression=pattern_compression,
            bootstrap_convergence=bootstrap_convergence,
            bootstop_threshold=bootstop_threshold,
            bootstop_perms=bootstop_perms,
            print_bootstrap_trees=print_bootstrap_trees,
            raxml_path=raxml_path,
            no_seq_check=no_seq_check,
            silent=silent
        )
        
        # 直接运行分析，无额外包装
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        click.echo("\nRAxML系统发育分析被用户中断 | RAxML phylogenetic analysis interrupted by user", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"RAxML系统发育分析失败 | RAxML phylogenetic analysis failed: {e}", err=True)
        raise click.Abort()