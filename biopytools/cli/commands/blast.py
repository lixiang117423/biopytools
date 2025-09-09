# """
# BLAST比对分析命令 | BLAST Alignment Analysis Command
# """

# import click
# import sys


# @click.command(short_help='BLAST序列比对分析工具',
#                context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
# @click.option('--input', '-i',
#               type=click.Path(exists=True),
#               help='输入文件或目录路径（支持单个文件或包含多个文件的目录）| Input file or directory path')
# @click.option('--sample-map-file', '-s',
#               type=click.Path(exists=True),
#               help='样品映射文件，格式：文件路径<TAB>样品名称 | Sample mapping file, format: filepath<TAB>samplename')
# @click.option('--target-file', '-t',
#               required=True,
#               type=click.Path(exists=True),
#               help='目标基因序列文件 | Target gene sequence file')
# @click.option('--output-dir', '-o',
#               default='./blast_output',
#               type=click.Path(),
#               help='输出目录 | Output directory (default: ./blast_output)')
# @click.option('--blast-type',
#               default='blastn',
#               type=click.Choice(['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']),
#               help='BLAST程序类型 | BLAST program type (default: blastn)')
# @click.option('--evalue', '-e',
#               default=1e-5,
#               type=float,
#               help='E-value阈值 | E-value threshold (default: 1e-5)')
# @click.option('--max-target-seqs',
#               default=10,
#               type=int,
#               help='最大目标序列数 | Maximum target sequences (default: 10)')
# @click.option('--word-size',
#               default=11,
#               type=int,
#               help='词大小 (仅适用于blastn/tblastx) | Word size (for blastn/tblastx only) (default: 11)')
# @click.option('--threads', '-j',
#               default=88,
#               type=int,
#               help='线程数 | Number of threads (default: 88)')
# @click.option('--input-suffix',
#               default='*.fa',
#               type=str,
#               help='输入文件后缀模式（当输入为目录时使用）| Input file suffix pattern (default: *.fa)')
# @click.option('--target-db-type',
#               default='nucl',
#               type=click.Choice(['nucl', 'prot']),
#               help='目标数据库类型 | Target database type (default: nucl)')
# @click.option('--min-identity',
#               default=70.0,
#               type=float,
#               help='最小序列相似度 (%%) | Minimum sequence identity (%%) (default: 70.0)')
# @click.option('--min-coverage',
#               default=50.0,
#               type=float,
#               help='最小覆盖度 (%%) | Minimum coverage (%%) (default: 50.0)')
# @click.option('--high-quality-evalue',
#               default=1e-10,
#               type=float,
#               help='高质量比对E-value阈值 | High quality alignment E-value threshold (default: 1e-10)')
# @click.option('--auto-detect-samples',
#               is_flag=True,
#               default=True,
#               help='自动检测样品名称（仅当使用-i时有效）| Auto-detect sample names (only when using -i)')
# @click.option('--sample-name-pattern',
#               default=r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$',
#               type=str,
#               help='样品名称提取正则表达式 | Sample name extraction regex')
# @click.option('--makeblastdb-path',
#               default='makeblastdb',
#               type=str,
#               help='makeblastdb程序路径 | makeblastdb program path (default: makeblastdb)')
# @click.option('--blastn-path',
#               default='blastn',
#               type=str,
#               help='blastn程序路径 | blastn program path (default: blastn)')
# @click.option('--blastp-path',
#               default='blastp',
#               type=str,
#               help='blastp程序路径 | blastp program path (default: blastp)')
# @click.option('--blastx-path',
#               default='blastx',
#               type=str,
#               help='blastx程序路径 | blastx program path (default: blastx)')
# @click.option('--tblastn-path',
#               default='tblastn',
#               type=str,
#               help='tblastn程序路径 | tblastn program path (default: tblastn)')
# @click.option('--tblastx-path',
#               default='tblastx',
#               type=str,
#               help='tblastx程序路径 | tblastx program path (default: tblastx)')

"""
🧬 BLAST比对分析命令 | BLAST Alignment Analysis Command
"""

import click
import sys


@click.command(short_help='🔍 BLAST序列比对分析工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              type=click.Path(exists=True),
              help='📁 输入文件或目录路径（支持单个文件或包含多个文件的目录）| Input file or directory path')
@click.option('--sample-map-file', '-s',
              type=click.Path(exists=True),
              help='🗺️  样品映射文件，格式：文件路径<TAB>样品名称 | Sample mapping file, format: filepath<TAB>samplename')
@click.option('--target-file', '-t',
              required=True,
              type=click.Path(exists=True),
              help='🎯 目标基因序列文件 | Target gene sequence file')
@click.option('--output-dir', '-o',
              default='./blast_output',
              type=click.Path(),
              help='📂 输出目录 | Output directory (default: ./blast_output)')
@click.option('--blast-type',
              default='blastn',
              type=click.Choice(['blastn', 'blastp', 'blastx', 'tblastn', 'tblastx']),
              help='⚡ BLAST程序类型 | BLAST program type (default: blastn)')
@click.option('--evalue', '-e',
              default=1e-5,
              type=float,
              help='📊 E-value阈值 | E-value threshold (default: 1e-5)')
@click.option('--max-target-seqs',
              default=10,
              type=int,
              help='🔢 最大目标序列数 | Maximum target sequences (default: 10)')
@click.option('--word-size',
              default=11,
              type=int,
              help='📏 词大小 (仅适用于blastn/tblastx) | Word size (for blastn/tblastx only) (default: 11)')
@click.option('--threads', '-j',
              default=88,
              type=int,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--input-suffix',
              default='*.fa',
              type=str,
              help='📄 输入文件后缀模式（当输入为目录时使用）| Input file suffix pattern (default: *.fa)')
@click.option('--target-db-type',
              default='nucl',
              type=click.Choice(['nucl', 'prot']),
              help='💾 目标数据库类型 | Target database type (default: nucl)')
@click.option('--min-identity',
              default=70.0,
              type=float,
              help='🎯 最小序列相似度 (%%) | Minimum sequence identity (%%) (default: 70.0)')
@click.option('--min-coverage',
              default=50.0,
              type=float,
              help='📐 最小覆盖度 (%%) | Minimum coverage (%%) (default: 50.0)')
@click.option('--high-quality-evalue',
              default=1e-10,
              type=float,
              help='⭐ 高质量比对E-value阈值 | High quality alignment E-value threshold (default: 1e-10)')
@click.option('--auto-detect-samples',
              is_flag=True,
              default=True,
              help='🤖 自动检测样品名称（仅当使用-i时有效）| Auto-detect sample names (only when using -i)')
@click.option('--sample-name-pattern',
              default=r'([^/]+?)(?:\.fa|\.fasta|\.fna)?$',
              type=str,
              help='🔍 样品名称提取正则表达式 | Sample name extraction regex')
@click.option('--makeblastdb-path',
              default='makeblastdb',
              type=str,
              help='🛠️  makeblastdb程序路径 | makeblastdb program path (default: makeblastdb)')
@click.option('--blastn-path',
              default='blastn',
              type=str,
              help='🔬 blastn程序路径 | blastn program path (default: blastn)')
@click.option('--blastp-path',
              default='blastp',
              type=str,
              help='🧪 blastp程序路径 | blastp program path (default: blastp)')
@click.option('--blastx-path',
              default='blastx',
              type=str,
              help='🔄 blastx程序路径 | blastx program path (default: blastx)')
@click.option('--tblastn-path',
              default='tblastn',
              type=str,
              help='🔀 tblastn程序路径 | tblastn program path (default: tblastn)')
@click.option('--tblastx-path',
              default='tblastx',
              type=str,
              help='🔁 tblastx程序路径 | tblastx program path (default: tblastx)')

def blast(input, sample_map_file, target_file, output_dir, blast_type, evalue,
          max_target_seqs, word_size, threads, input_suffix, target_db_type,
          min_identity, min_coverage, high_quality_evalue, auto_detect_samples,
          sample_name_pattern, makeblastdb_path, blastn_path, blastp_path,
          blastx_path, tblastn_path, tblastx_path):
    """
    BLAST比对分析脚本 (模块化版本)
    
    基于BLAST+套件的序列比对分析工具，支持多种BLAST算法和批量序列处理，
    适用于基因功能注释、同源性分析和序列比较研究。
    
    功能特点 | Features:
    - 支持5种BLAST算法类型
    - 批量序列自动处理
    - 灵活的样品映射管理
    - 智能结果过滤和排序
    - 高质量比对结果筛选
    - 详细统计报告生成
    - 多线程并行加速
    
    分析流程 | Analysis Pipeline:
    1. 样品文件自动发现或映射解析
    2. 目标数据库构建和索引
    3. BLAST比对参数优化执行
    4. 比对结果解析和质量评估
    5. 多维度过滤和结果排序
    6. 高质量结果集筛选
    7. 统计报告和可视化生成
    
    应用场景 | Use Cases:
    - 基因功能注释
    - 同源基因识别
    - 序列相似性分析
    - 进化关系研究
    - 基因家族分析
    - 转录组功能注释
    - 蛋白质结构域预测
    
    示例 | Examples:
    
    \b
    # 目录批量核酸比对 (自动生成样品映射)
    biopytools blast -i sequences/ -t nlr_genes.fa -o blast_results
    
    \b
    # 单文件比对分析
    biopytools blast -i single_sequence.fa -t targets.fa -o results
    
    \b
    # 使用现有样品映射文件
    biopytools blast -s sample_map.txt -t nlr_genes.fa -o results
    
    \b
    # 蛋白质序列比对
    biopytools blast -i proteins/ -t targets.fa -o results \\
        --blast-type blastp --target-db-type prot --threads 32
    
    \b
    # 高严格性过滤分析
    biopytools blast -i sequences/ -t targets.fa -o results \\
        --min-identity 85 --min-coverage 80 --evalue 1e-10
    
    \b
    # 跨物种同源基因搜索
    biopytools blast -i query_genes.fa -t genome_proteins.fa \\
        -o homolog_search --blast-type blastx \\
        --evalue 1e-6 --max-target-seqs 20
    
    \b
    # 转录组功能注释
    biopytools blast -i transcripts/ -t functional_database.fa \\
        -o annotation_results --blast-type blastx \\
        --min-identity 60 --min-coverage 70
    
    BLAST算法选择指南 | BLAST Algorithm Selection Guide:
    
    blastn (核酸-核酸):
    - 查询序列: DNA/RNA
    - 数据库: DNA/RNA
    - 用途: 基因定位、转录本比对、SNP分析
    
    blastp (蛋白质-蛋白质):
    - 查询序列: 蛋白质
    - 数据库: 蛋白质
    - 用途: 蛋白质功能预测、结构域分析
    
    blastx (核酸-蛋白质):
    - 查询序列: DNA/RNA (6框翻译)
    - 数据库: 蛋白质
    - 用途: 编码基因功能注释、ORF预测
    
    tblastn (蛋白质-核酸):
    - 查询序列: 蛋白质
    - 数据库: DNA/RNA (6框翻译)
    - 用途: 基因发现、假基因识别
    
    tblastx (核酸-核酸翻译):
    - 查询序列: DNA/RNA (6框翻译)
    - 数据库: DNA/RNA (6框翻译)
    - 用途: 进化分析、假基因比较
    
    输入数据组织 | Input Data Organization:
    
    方式1: 目录输入 (推荐批量处理)
    ```
    sequences/
    ├── sample1.fa
    ├── sample2.fa
    └── sample3.fa
    ```
    
    方式2: 样品映射文件 (精确控制)
    ```
    /path/to/seq1.fa    Sample1
    /path/to/seq2.fa    Sample2_treated
    /path/to/seq3.fa    Sample3_control
    ```
    
    方式3: 单文件输入
    ```
    single_query.fa -> 直接比对
    ```
    
    参数优化策略 | Parameter Optimization Strategy:
    
    E-value设置:
    - 1e-3: 宽松，发现远程同源性
    - 1e-5: 标准，平衡敏感性和特异性
    - 1e-10: 严格，高置信度比对
    
    Word size调优:
    - blastn: 7(敏感) - 28(快速)
    - blastp: 2(敏感) - 6(快速)
    - 短序列建议使用小word size
    
    身份度和覆盖度:
    - 身份度 >90%: 近缘物种比较
    - 身份度 70-90%: 同科或同目比较
    - 身份度 50-70%: 远缘同源性检测
    - 覆盖度 >80%: 全长比对
    - 覆盖度 50-80%: 部分结构域比对
    
    性能调优建议 | Performance Tuning Recommendations:
    
    线程配置:
    - 小数据集 (<1GB): 8-16线程
    - 中等数据集 (1-10GB): 32-64线程
    - 大数据集 (>10GB): 64-88线程
    
    内存管理:
    - 小数据库 (<1GB): 4-8GB内存
    - 大数据库 (>10GB): 32GB+内存
    - 考虑数据库分割策略
    
    I/O优化:
    - 使用SSD存储提高速度
    - 临时文件放在本地高速存储
    - 避免网络存储瓶颈
    
    结果过滤和质控 | Result Filtering and Quality Control:
    
    自动过滤机制:
    - E-value阈值过滤
    - 序列身份度筛选
    - 覆盖度质量控制
    - 比对长度验证
    
    高质量结果筛选:
    - 更严格的E-value (默认1e-10)
    - 更高的身份度要求
    - 更好的覆盖度标准
    - 适用于关键基因鉴定
    
    结果排序策略:
    - 主要: 覆盖度降序
    - 次要: E-value升序
    - 第三: 身份度降序
    - 便于识别最佳比对
    
    输出文件说明 | Output Files Description:
    
    核心结果文件:
    - blast_results_sorted.tsv: 排序的比对结果
    - high_quality_results.tsv: 高质量比对筛选
    - blast_statistics.txt: 详细统计报告
    - sample_mapping.txt: 样品映射文件
    
    统计信息包含:
    - 总查询序列数量
    - 成功比对序列数量
    - 平均身份度分布
    - E-value分布统计
    - 覆盖度质量评估
    
    工作流程选择 | Workflow Selection:
    
    场景1: 快速基因注释
    ```bash
    # 使用默认参数，快速获得功能注释
    biopytools blast -i transcripts.fa -t nr_database.fa -o annotation
    ```
    
    场景2: 严格同源性分析
    ```bash
    # 高严格性参数，确保比对可靠性
    biopytools blast -i genes.fa -t targets.fa -o homologs \\
        --evalue 1e-15 --min-identity 90 --min-coverage 85
    ```
    
    场景3: 批量样品处理
    ```bash
    # 自动化处理多个样品文件
    biopytools blast -i samples_dir/ -t database.fa -o batch_results \\
        --threads 64 --input-suffix "*.fasta"
    ```
    
    常见问题解决 | Troubleshooting:
    
    输入文件问题:
    - 确保FASTA格式规范
    - 检查文件编码 (推荐UTF-8)
    - 验证序列标识符唯一性
    - 确认文件权限可读
    
    性能问题:
    - 监控内存使用情况
    - 调整线程数避免过载
    - 检查磁盘I/O性能
    - 考虑数据库大小限制
    
    结果异常:
    - 验证BLAST程序版本兼容性
    - 检查数据库构建是否成功
    - 确认参数设置合理性
    - 查看详细错误日志
    
    无结果或结果过少:
    - 放宽E-value阈值
    - 降低身份度要求
    - 减少覆盖度限制
    - 检查查询序列质量
    
    高级应用技巧 | Advanced Usage Tips:
    
    自定义过滤策略:
    - 根据研究目标调整阈值
    - 结合多个质量指标
    - 考虑生物学意义
    
    结果后处理:
    - 使用R/Python进行深度分析
    - 整合GO/KEGG功能注释
    - 进行网络分析和可视化
    
    大规模数据处理:
    - 分批处理策略
    - 并行任务调度
    - 结果合并和去重
    - 存储优化方案
    """
    
    # 参数验证
    if not input and not sample_map_file:
        click.echo("错误：必须指定输入路径(-i)或样品映射文件(-s)中的一个", err=True)
        sys.exit(1)
    
    # 懒加载：只在需要时导入重量级模块
    def lazy_import_blast():
        try:
            from ...blast_analysis import BLASTAnalyzer
            return BLASTAnalyzer
        except ImportError as e:
            click.echo(f"无法导入BLAST模块 | Failed to import BLAST module: {e}", err=True)
            click.echo("请确保已正确安装BLAST+相关依赖 | Please ensure BLAST+ dependencies are properly installed", err=True)
            sys.exit(1)
    
    try:
        # 延迟导入：只有在真正需要执行分析时才加载模块
        BLASTAnalyzer = lazy_import_blast()
        
        # 直接创建分析器，避免参数重构的性能损失
        analyzer = BLASTAnalyzer(
            input_path=input,
            target_file=target_file,
            output_dir=output_dir,
            sample_map_file=sample_map_file,
            blast_type=blast_type,
            evalue=evalue,
            max_target_seqs=max_target_seqs,
            word_size=word_size,
            threads=threads,
            input_suffix=input_suffix,
            target_db_type=target_db_type,
            min_identity=min_identity,
            min_coverage=min_coverage,
            high_quality_evalue=high_quality_evalue,
            auto_detect_samples=auto_detect_samples,
            sample_name_pattern=sample_name_pattern,
            makeblastdb_path=makeblastdb_path,
            blastn_path=blastn_path,
            blastp_path=blastp_path,
            blastx_path=blastx_path,
            tblastn_path=tblastn_path,
            tblastx_path=tblastx_path
        )
        
        # 直接运行分析，无额外包装
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        click.echo("\nBLAST分析被用户中断 | BLAST analysis interrupted by user", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"BLAST分析失败 | BLAST analysis failed: {e}", err=True)
        raise click.Abort()


# 性能优化的懒加载管理器
class BLASTLazyLoader:
    """BLAST懒加载管理器 | BLAST Lazy Loader Manager"""
    
    def __init__(self):
        self._analyzer_class = None
        self._loaded_modules = {}
    
    def get_analyzer_class(self):
        """获取分析器类，首次调用时加载 | Get analyzer class, load on first call"""
        if self._analyzer_class is None:
            from ...blast_analyzer import BLASTAnalyzer
            self._analyzer_class = BLASTAnalyzer
        return self._analyzer_class
    
    def get_module(self, module_name):
        """获取指定模块，支持细粒度懒加载 | Get specific module with fine-grained lazy loading"""
        if module_name not in self._loaded_modules:
            if module_name == 'config':
                from ...blast_analyzer.config import BLASTConfig
                self._loaded_modules[module_name] = BLASTConfig
            elif module_name == 'utils':
                from ...blast_analyzer import utils
                self._loaded_modules[module_name] = utils
            # 可以添加更多模块的细粒度加载
        
        return self._loaded_modules.get(module_name)
    
    def is_analyzer_loaded(self):
        """检查分析器是否已加载 | Check if analyzer is loaded"""
        return self._analyzer_class is not None
    
    def clear_cache(self):
        """清理缓存，释放内存 | Clear cache and free memory"""
        self._analyzer_class = None
        self._loaded_modules.clear()


# 全局懒加载器实例
_blast_loader = BLASTLazyLoader()