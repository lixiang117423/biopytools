"""
VCF系统发育分析命令 | VCF Phylogenetic Analysis Command
"""

import click
import sys
from ...vcf_phylo.main import main as vcf_phylo_main


@click.command(short_help='VCF文件系统发育树构建和分析',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              type=click.Path(exists=True),
              help='输入VCF文件路径 | Input VCF file path')
@click.option('--distance-matrix', '-d',
              type=click.Path(exists=True),
              help='已有的距离矩阵文件路径（用于跳过VCF2Dis步骤）| Existing distance matrix file path (for skipping VCF2Dis step)')
@click.option('--output', '-o',
              default='phylo_analysis',
              type=str,
              help='输出文件前缀 | Output file prefix (default: phylo_analysis)')
@click.option('--tree-output', '-t',
              type=click.Path(),
              help='系统发育树输出文件路径（默认: OUTPUT.nwk）| Phylogenetic tree output file path (default: OUTPUT.nwk)')
@click.option('--vcf2dis-path',
              default='VCF2Dis',
              type=str,
              help='VCF2Dis程序路径 | VCF2Dis program path (default: VCF2Dis)')
@click.option('--working-dir', '-w',
              default='.',
              type=click.Path(),
              help='工作目录 | Working directory (default: current directory)')
@click.option('--skip-vcf2dis',
              is_flag=True,
              help='跳过VCF2Dis步骤，直接从距离矩阵构建树 | Skip VCF2Dis step, build tree directly from distance matrix')
def vcf_nj_tree(input, distance_matrix, output, tree_output, vcf2dis_path, working_dir, skip_vcf2dis):
    """
    VCF系统发育分析工具 v2.0
    
    基于VCF文件构建系统发育树，使用距离矩阵方法和邻接法(NJ)算法。
    支持从VCF文件计算遗传距离，构建高质量的系统发育关系图谱。
    
    功能特点 | Features:
    - VCF文件遗传距离计算
    - 邻接法(NJ)系统发育树构建
    - 灵活的输入输出选项
    - 结果验证和质量评估
    - 详细的分析报告生成
    - 支持已有距离矩阵输入
    
    分析流程 | Analysis Pipeline:
    1. VCF文件解析和预处理
    2. 遗传距离矩阵计算
    3. 邻接法系统发育树构建
    4. 结果验证和质量检查
    5. 分析报告和可视化输出
    
    示例 | Examples:
    
    \b
    # 基本系统发育分析
    biopytools vcf-phylo -i wild.snp.vcf -o wild_snp_phylo
    
    \b
    # 指定完整输出路径
    biopytools vcf-phylo --input data/wild.snp.vcf \\
        --output results/wild_snp_analysis \\
        --tree-output results/wild_snp.nwk
    
    \b
    # 从已有距离矩阵构建树
    biopytools vcf-phylo --distance-matrix existing_matrix.txt \\
        --tree-output phylo_tree.nwk --skip-vcf2dis
    
    \b
    # 自定义工作目录和工具路径
    biopytools vcf-phylo -i variants.vcf -o analysis \\
        --vcf2dis-path /opt/vcf2dis/VCF2Dis \\
        --working-dir /tmp/phylo_work/
    
    \b
    # 复杂群体系统发育分析
    biopytools vcf-phylo -i population_variants.vcf \\
        -o population_phylogeny \\
        -w ./phylo_analysis/ \\
        -t population_tree.newick
    
    输入文件要求 | Input Requirements:
    
    VCF文件格式:
    - 标准VCF 4.0+格式
    - 包含多个样本的基因型数据
    - 变异位点质量控制后的数据
    - 建议预先进行MAF和缺失率过滤
    
    距离矩阵文件格式（可选）:
    - 方阵格式，样本间遗传距离
    - 制表符或空格分隔
    - 第一行和第一列为样本名称
    
    输出文件 | Output Files:
    - {prefix}.dis: 遗传距离矩阵文件
    - {prefix}.nwk: Newick格式系统发育树
    - {prefix}_report.txt: 详细分析报告
    - {prefix}.log: 分析过程日志
    
    系统发育方法说明 | Phylogenetic Methods:
    
    距离计算:
    - 基于SNP差异的遗传距离
    - 考虑缺失数据的校正
    - 标准化距离矩阵
    
    树构建算法:
    - 邻接法(Neighbor-Joining)
    - 基于最小进化原理
    - 适用于大样本量分析
    
    应用场景 | Use Cases:
    - 物种进化关系研究
    - 群体结构和亲缘关系分析
    - 品种/品系分类和鉴定
    - 保护遗传学评估
    - 育种材料遗传背景分析
    
    参数使用指南 | Parameter Guidelines:
    
    --input vs --distance-matrix:
    - 使用--input从VCF开始完整分析
    - 使用--distance-matrix跳过距离计算步骤
    - 两者必须提供其中之一
    
    --skip-vcf2dis标志:
    - 与--distance-matrix配合使用
    - 适用于已有距离矩阵的情况
    - 可以节省计算时间
    
    工作目录设置:
    - 建议为大型分析设置专用目录
    - 确保有足够的磁盘空间
    - 临时文件会在此目录生成
    
    性能优化建议 | Performance Tips:
    - 大数据集建议预先进行质量控制
    - 使用SSD存储提升I/O性能
    - 合理设置工作目录避免路径问题
    - 监控内存使用情况
    
    质量控制建议 | Quality Control:
    - 输入VCF文件建议预过滤
    - 移除低质量和高缺失位点
    - 考虑MAF阈值设置
    - 检查样本标识一致性
    
    结果解读 | Result Interpretation:
    
    距离矩阵:
    - 数值越小表示遗传关系越近
    - 对角线为0（自身距离）
    - 矩阵对称性检验
    
    系统发育树:
    - 枝长代表进化距离
    - 拓扑结构显示亲缘关系
    - 支持Newick标准格式
    
    常见问题 | Troubleshooting:
    - 确保VCF2Dis工具正确安装
    - 检查输入文件路径和格式
    - 验证样本名称不含特殊字符
    - 监控内存和磁盘空间使用
    
    依赖软件 | Dependencies:
    - VCF2Dis: 遗传距离计算
    - Python科学计算库
    - 标准Unix工具
    
    注意事项 | Notes:
    - 系统发育分析结果需要生物学验证
    - 建议结合其他方法验证树拓扑
    - 大样本量可能需要较长计算时间
    - 结果可用于后续可视化和注释
    """
    
    # 验证参数逻辑
    if skip_vcf2dis:
        if not distance_matrix:
            raise click.ClickException("使用 --skip-vcf2dis 时必须指定 --distance-matrix | --distance-matrix must be specified when using --skip-vcf2dis")
    else:
        if not input:
            raise click.ClickException("必须指定输入VCF文件 | Input VCF file must be specified")
    
    # 构建参数列表传递给原始main函数
    args = ['vcf_nj_tree.py']
    
    # 输入文件参数
    if input:
        args.extend(['-i', input])
    
    if distance_matrix:
        args.extend(['-d', distance_matrix])
    
    # 输出文件参数
    if output != 'phylo_analysis':
        args.extend(['-o', output])
    
    if tree_output:
        args.extend(['-t', tree_output])
    
    # 工具参数
    if vcf2dis_path != 'VCF2Dis':
        args.extend(['--vcf2dis-path', vcf2dis_path])
    
    if working_dir != '.':
        args.extend(['-w', working_dir])
    
    # 行为控制参数
    if skip_vcf2dis:
        args.append('--skip-vcf2dis')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        vcf_phylo_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n系统发育分析被用户中断 | Phylogenetic analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"系统发育分析失败 | Phylogenetic analysis failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv