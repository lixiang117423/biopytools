"""
EDTA植物基因组TE注释工具命令 | EDTA Plant Genome TE Annotation Command
"""

import click
import sys
from pathlib import Path


@click.command(short_help='EDTA植物基因组转座元素注释工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genome', '-g',
              required=True,
              type=click.Path(exists=True),
              help='输入基因组FASTA文件路径 | Input genome FASTA file path')
@click.option('--output-dir', '-o',
              default='./edta_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./edta_output)')
@click.option('--species', '-s',
              default='others',
              type=click.Choice(['Rice', 'Maize', 'others']),
              help='物种类型 | Species for TIR candidate identification (default: others)')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='使用的线程数 | Number of threads to use (default: 88)')
@click.option('--step',
              default='all',
              type=click.Choice(['all', 'filter', 'final', 'anno']),
              help='指定运行的步骤 | Specify which steps to run (default: all)')
@click.option('--sensitive',
              default=0,
              type=click.Choice([0, 1]),
              help='使用RepeatModeler增强敏感性 | Use RepeatModeler to enhance sensitivity (default: 0)')
@click.option('--anno',
              default=1,
              type=click.Choice([0, 1]),
              help='执行全基因组TE注释 | Perform whole-genome TE annotation (default: 1)')
@click.option('--cds',
              type=click.Path(exists=True),
              help='提供CDS序列文件 | Provide CDS sequences file')
@click.option('--genome-list',
              type=click.Path(exists=True),
              help='包含多个基因组文件路径的列表文件 | File containing list of genome file paths')
@click.option('--batch-mode',
              is_flag=True,
              help='启用批量处理模式 | Enable batch processing mode')
@click.option('--resume',
              is_flag=True,
              help='从上次中断的地方继续分析 | Resume analysis from last interruption')
@click.option('--generate-plots',
              is_flag=True,
              default=True,
              help='生成可视化图表 | Generate visualization plots')
@click.option('--compare-results',
              is_flag=True,
              help='启用结果比较 | Enable results comparison')
@click.option('--check-dependencies',
              is_flag=True,
              help='检查依赖软件后退出 | Check dependencies and exit')
def edta(genome, output_dir, species, threads, step, sensitive, anno, cds,
         genome_list, batch_mode, resume, generate_plots, compare_results,
         check_dependencies):
    """
    EDTA植物基因组转座元素注释工具 (模块化版本)
    
    基于EDTA的植物基因组转座元素识别和注释工具，支持全自动化的TE分析流程，
    适用于植物基因组学研究中的重复序列分析和基因组注释。
    
    功能特点 | Features:
    - 全自动TE识别和分类
    - 多物种优化算法支持
    - 高效的并行处理能力
    - 断点续跑机制
    - 批量基因组处理
    - 丰富的可视化输出
    - 详细的注释报告
    
    分析流程 | Analysis Pipeline:
    1. 输入文件验证和预处理
    2. LTR逆转录转座子识别
    3. TIR转座子识别和分类
    4. Helitron转座子检测
    5. 重复序列库构建
    6. 全基因组TE注释
    7. 结果整理和可视化
    
    应用场景 | Use Cases:
    - 植物基因组重复序列分析
    - 转座元素进化研究
    - 基因组组装质量评估
    - 比较基因组学研究
    - 基因组注释流程
    - TE插入位点分析
    
    示例 | Examples:
    
    \b
    # 基本TE注释分析
    biopytools edta --g plant.fa -o edta_results
    
    \b
    # 水稻基因组专用分析
    biopytools edta -g rice_genome.fa -o rice_results \\
        --species Rice --threads 64
    
    \b
    # 增强敏感性分析
    biopytools edta --g plant.fa --cds proteins.fa \\
        --sensitive 1 --threads 88
    
    \b
    # 批量基因组处理
    biopytools edta --genome-list genomes.txt \\
        --batch-mode -o batch_results
    
    \b
    # 断点续跑分析
    biopytools edta -g plant.fa -o results \\
        --resume --generate-plots
    
    \b
    # 仅检查依赖软件
    biopytools edta -g dummy.fa --check-dependencies
    
    \b
    # 特定步骤分析
    biopytools edta --g plant.fa --step filter \\
        -o filtered_results --threads 32
    
    输入文件要求 | Input File Requirements:
    
    基因组文件:
    - FASTA格式 (.fa, .fasta, .fna)
    - 建议使用高质量基因组组装
    - 文件大小支持从几MB到几GB
    - 序列标识符应简洁明了
    
    CDS文件 (可选):
    - FASTA格式的编码序列
    - 用于提高TE识别准确性
    - 避免将基因误识别为TE
    - 建议使用高质量基因预测结果
    
    基因组列表文件:
    - 每行一个基因组文件路径
    - 支持绝对路径和相对路径
    - 以#开头的行为注释行
    - 用于批量处理模式
    
    物种参数选择 | Species Parameter Selection:
    
    --species Rice:
    - 针对水稻基因组优化
    - 使用水稻特异的TIR识别参数
    - 适用于禾本科植物
    
    --species Maize:
    - 针对玉米基因组优化
    - 处理大基因组和高TE含量
    - 适用于大基因组禾本科植物
    
    --species others:
    - 通用参数设置
    - 适用于大多数植物物种
    - 平衡敏感性和特异性
    
    分析步骤说明 | Analysis Step Description:
    
    --step all (推荐):
    - 完整的TE注释流程
    - 从头开始到最终注释
    - 包含所有质量控制步骤
    
    --step filter:
    - 仅执行TE序列过滤
    - 适用于已有初步结果的情况
    - 提高TE库质量
    
    --step final:
    - 执行最终整理步骤
    - 生成最终TE库文件
    - 进行序列分类和命名
    
    --step anno:
    - 仅执行全基因组注释
    - 需要预先存在TE库文件
    - 快速获得注释结果
    
    敏感性选项 | Sensitivity Options:
    
    --sensitive 0 (默认):
    - 标准敏感性模式
    - 平衡速度和准确性
    - 适用于大多数应用场景
    
    --sensitive 1:
    - 高敏感性模式
    - 使用RepeatModeler增强检测
    - 处理时间显著增加
    - 适用于新物种或复杂基因组
    
    性能优化建议 | Performance Optimization:
    
    线程配置:
    - 小基因组 (<500MB): 16-32线程
    - 中等基因组 (500MB-2GB): 32-64线程
    - 大基因组 (>2GB): 64-88线程
    
    内存需求:
    - 基础内存: 16GB+
    - 大基因组: 64GB+
    - 批量处理: 128GB+
    
    存储需求:
    - 临时空间: 输入文件的5-10倍
    - 最终结果: 输入文件的2-3倍
    - 建议使用SSD提高I/O性能
    
    输出文件说明 | Output Files Description:
    
    核心结果文件:
    - {genome}.EDTA.TElib.fa: 最终TE库文件
    - {genome}.EDTA.TEanno.gff3: TE注释GFF3文件
    - {genome}.EDTA.intact.fa: 完整TE序列
    - {genome}.EDTA.intact.gff3: 完整TE注释
    
    统计和报告文件:
    - {genome}.EDTA.TElib.nov.fa: 新发现的TE
    - EDTA.log: 详细分析日志
    - summary_report.html: 可视化总结报告
    - statistics.json: 详细统计数据
    
    中间文件:
    - LTR/: LTR逆转录转座子结果
    - TIR/: TIR转座子结果  
    - Helitron/: Helitron转座子结果
    - combine/: 组合结果文件
    
    断点续跑机制 | Resume Mechanism:
    
    自动检查点:
    - 每个主要步骤完成后保存状态
    - 断电或意外中断后可恢复
    - 自动识别已完成的步骤
    
    续跑策略:
    - 检查输出目录中的中间文件
    - 从最后一个完整步骤开始
    - 避免重复计算已完成的部分
    
    故障排除 | Troubleshooting:
    
    常见错误:
    1. 依赖软件未安装
       - 运行 --check-dependencies 检查
       - 按照EDTA官方文档安装依赖
    
    2. 内存不足
       - 减少线程数
       - 使用更大内存的机器
       - 分批处理大基因组
    
    3. 磁盘空间不足
       - 清理临时文件
       - 使用更大存储空间
       - 监控磁盘使用情况
    
    4. 分析中断
       - 使用 --resume 选项继续
       - 检查日志文件获取详细信息
       - 必要时重新开始分析
    
    最佳实践 | Best Practices:
    
    数据准备:
    - 使用高质量基因组组装
    - 确保基因组文件完整性
    - 提供高质量CDS序列(推荐)
    
    参数选择:
    - 根据物种选择合适的species参数
    - 新物种或复杂基因组考虑使用--sensitive 1
    - 合理分配计算资源
    
    结果验证:
    - 检查TE库的完整性
    - 验证注释结果的合理性
    - 与已知数据库进行比较
    
    高级应用技巧 | Advanced Usage Tips:
    
    自定义分析:
    - 根据需要选择特定分析步骤
    - 结合其他TE分析工具验证结果
    - 进行TE进化分析
    
    结果后处理:
    - 使用R或Python进行统计分析
    - 整合基因组注释信息
    - 进行比较基因组分析
    """
    
    # 懒加载：只在需要时导入重量级模块
    def lazy_import_edta():
        try:
            from ...edta import EDTAAnalyzer
            return EDTAAnalyzer
        except ImportError as e:
            click.echo(f"无法导入EDTA模块 | Failed to import EDTA module: {e}", err=True)
            click.echo("请确保已正确安装EDTA相关依赖 | Please ensure EDTA dependencies are properly installed", err=True)
            sys.exit(1)
    
    # 处理基因组列表文件
    genome_list_data = None
    if genome_list:
        try:
            with open(genome_list, 'r') as f:
                genome_list_data = [line.strip() for line in f if line.strip() and not line.startswith('#')]
            batch_mode = True
        except FileNotFoundError:
            click.echo(f"基因组列表文件不存在 | Genome list file not found: {genome_list}", err=True)
            sys.exit(1)
    
    try:
        # 延迟导入：只有在真正需要执行分析时才加载模块
        EDTAAnalyzer = lazy_import_edta()
        
        # 直接创建分析器，避免参数重构的性能损失
        analyzer = EDTAAnalyzer(
            genome=genome,
            genome_list=genome_list_data,
            output_dir=output_dir,
            species=species,
            threads=threads,
            step=step,
            sensitive=sensitive,
            anno=anno,
            cds=cds,
            batch_mode=batch_mode,
            resume=resume,
            generate_plots=generate_plots,
            compare_results=compare_results,
            check_dependencies=check_dependencies
        )
        
        # 如果只是检查依赖，执行检查后退出
        if check_dependencies:
            analyzer.check_dependencies()
            click.echo("所有依赖检查通过 | All dependencies check passed")
            return
        
        # 运行完整分析
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        click.echo("\nEDTA分析被用户中断 | EDTA analysis interrupted by user", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"EDTA分析失败 | EDTA analysis failed: {e}", err=True)
        raise click.Abort()


# 性能优化的懒加载函数示例
def get_edta_analyzer_class():
    """懒加载获取EDTAAnalyzer类 | Lazy load EDTAAnalyzer class"""
    if not hasattr(get_edta_analyzer_class, '_cached_class'):
        from ...edta import EDTAAnalyzer
        get_edta_analyzer_class._cached_class = EDTAAnalyzer
    return get_edta_analyzer_class._cached_class


# 额外的性能优化工具
class LazyModuleLoader:
    """懒加载模块加载器 | Lazy module loader"""
    
    def __init__(self):
        self._modules = {}
    
    def get_module(self, module_name):
        """获取模块，首次访问时才加载 | Get module, load on first access"""
        if module_name not in self._modules:
            if module_name == 'edta_analyzer':
                from ...edta_te import EDTAAnalyzer
                self._modules[module_name] = EDTAAnalyzer
            # 可以添加更多模块的懒加载逻辑
        
        return self._modules[module_name]


# 全局懒加载器实例
_lazy_loader = LazyModuleLoader()