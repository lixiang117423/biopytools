"""
GenomeThreader基因预测分析命令 | GenomeThreader Gene Prediction Analysis Command
"""

import click
import sys


@click.command(short_help='GenomeThreader基因预测分析工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--genomic', '-g',
              required=True,
              type=click.Path(exists=True),
              help='输入基因组序列文件 (FASTA格式) | Input genomic sequence file (FASTA format)')
@click.option('--cdna', '-c',
              type=click.Path(exists=True),
              help='cDNA/转录本序列文件 | cDNA/transcript sequence file')
@click.option('--protein', '-p',
              type=click.Path(exists=True),
              help='蛋白质序列文件 | Protein sequence file')
@click.option('--est', '-e',
              type=click.Path(exists=True),
              help='EST序列文件 | EST sequence file')
@click.option('--output', '-o',
              default='./gth_output',
              type=click.Path(),
              help='输出目录 | Output directory (default: ./gth_output)')
@click.option('--species', '-s',
              type=click.Choice(['human', 'mouse', 'rat', 'chicken', 'drosophila', 'nematode', 
                                'fission_yeast', 'aspergillus', 'arabidopsis', 'maize', 'rice', 'medicago']),
              help='物种模型 | Species model')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='线程数 | Number of threads (default: 88)')
@click.option('--min-score',
              default=0.65,
              type=float,
              help='最小比对得分阈值 | Minimum alignment score threshold (default: 0.65)')
@click.option('--gc-min-coverage',
              default=50,
              type=int,
              help='全局链最小覆盖度百分比 | Global chain minimum coverage percentage (default: 50)')
@click.option('--paralogs',
              is_flag=True,
              help='计算旁系同源基因 | Compute paralogous genes')
@click.option('--intron-cutout',
              is_flag=True,
              help='启用内含子切除技术 | Enable intron cutout technique')
@click.option('--auto-intron-cutout',
              default=0,
              type=int,
              help='自动内含子切除矩阵大小(MB) | Auto intron cutout matrix size (MB) (default: 0)')
@click.option('--fastdp',
              is_flag=True,
              help='使用快速DP算法 | Use fast DP algorithm')
@click.option('--forward-only',
              is_flag=True,
              help='只分析正链 | Analyze forward strand only')
@click.option('--reverse-only',
              is_flag=True,
              help='只分析负链 | Analyze reverse strand only')
@click.option('--cdna-forward',
              is_flag=True,
              help='只比对cDNA正链 | Align cDNA forward strand only')
@click.option('--xml-output',
              is_flag=True,
              help='输出XML格式 | Output XML format')
@click.option('--intermediate',
              is_flag=True,
              help='输出中间结果用于后续处理 | Output intermediate results for further processing')
@click.option('--skip-alignment-out',
              is_flag=True,
              help='跳过比对输出 | Skip alignment output')
@click.option('--score-matrix',
              default='BLOSUM62',
              type=str,
              help='氨基酸替换评分矩阵 | Amino acid substitution scoring matrix (default: BLOSUM62)')
@click.option('--translation-table',
              default=1,
              type=int,
              help='密码子翻译表 | Codon translation table (default: 1)')
@click.option('--bssm-file',
              type=click.Path(exists=True),
              help='BSSM参数文件路径 | BSSM parameter file path')
@click.option('--from-pos',
              default=1,
              type=int,
              help='起始位置 | Starting position (default: 1)')
@click.option('--to-pos',
              type=int,
              help='结束位置 | Ending position')
@click.option('--width',
              default=1000000,
              type=int,
              help='序列处理宽度 | Sequence processing width (default: 1000000)')
@click.option('--first-alignments',
              default=0,
              type=int,
              help='每个基因组序列最大比对数 (0=无限制) | Max alignments per genomic sequence (0=unlimited) (default: 0)')
@click.option('--gth-path',
              default='gth',
              type=str,
              help='GenomeThreader程序路径 | GenomeThreader program path (default: gth)')
@click.option('--gthconsensus-path',
              default='gthconsensus',
              type=str,
              help='gthconsensus程序路径 | gthconsensus program path (default: gthconsensus)')
def genomethreader(genomic, cdna, protein, est, output, species, threads, min_score,
                   gc_min_coverage, paralogs, intron_cutout, auto_intron_cutout, fastdp,
                   forward_only, reverse_only, cdna_forward, xml_output, intermediate,
                   skip_alignment_out, score_matrix, translation_table, bssm_file,
                   from_pos, to_pos, width, first_alignments, gth_path, gthconsensus_path):
    """
    GenomeThreader基因预测分析工具 (模块化版本)
    
    基于GenomeThreader的基因结构预测工具，通过比对cDNA、蛋白质或EST序列
    到基因组序列来预测基因结构，适用于真核基因组注释和基因发现。
    
    功能特点 | Features:
    - 高精度基因结构预测
    - 多类型序列输入支持 (cDNA/蛋白质/EST)
    - 物种特异性参数优化
    - 并行计算加速处理
    - 旁系同源基因检测
    - 内含子边界精确识别
    - 多种输出格式支持
    
    预测流程 | Prediction Pipeline:
    1. 输入序列验证和统计
    2. 序列比对和剪接位点识别
    3. 外显子边界确定
    4. 基因结构模型构建
    5. 旁系同源基因分析 (可选)
    6. 结果质量评估
    7. 多格式输出生成
    
    应用场景 | Use Cases:
    - 新基因组注释
    - 基因结构验证
    - 选择性剪接分析
    - 同源基因预测
    - 基因家族分析
    - 比较基因组研究
    - 转录组辅助注释
    
    示例 | Examples:
    
    \b
    # cDNA序列基因预测
    biopytools genome-threader -g genome.fa -c cdna.fa -o gth_results
    
    \b
    # 蛋白质序列基因预测 (指定物种)
    biopytools genome-threader -g genome.fa -p proteins.fa \\
        -s arabidopsis -o results
    
    \b
    # 多类型序列联合预测
    biopytools genome-threader -g genome.fa -c cdna.fa -p proteins.fa \\
        -s human --paralogs -o comprehensive_results
    
    \b
    # 快速预测 (使用快速DP算法)
    biopytools genome-threader -g genome.fa -c cdna.fa \\
        --fastdp --intron-cutout -o fast_results
    
    \b
    # 高质量EST预测
    biopytools genome-threader -g genome.fa -e est.fa \\
        --intermediate --min-score 0.8 -o high_quality
    
    \b
    # 特定区域预测
    biopytools genome-threader -g genome.fa -c cdna.fa \\
        --from-pos 1000000 --to-pos 2000000 -o region_results
    
    \b
    # 单链分析
    biopytools genome-threader -g genome.fa -c cdna.fa \\
        --forward-only --xml-output -o forward_results
    
    输入文件要求 | Input File Requirements:
    
    基因组文件 (必需):
    - FASTA格式 (.fa, .fasta, .fna)
    - 高质量基因组组装
    - 染色体或scaffold级别组装
    - 序列标识符应简洁清晰
    
    序列文件 (至少需要一个):
    - cDNA文件: 全长转录本序列，最高预测准确性
    - 蛋白质文件: 同源蛋白质序列，适合远缘物种
    - EST文件: 表达序列标签，覆盖度可能不完整
    
    物种模型选择 | Species Model Selection:
    
    动物模型:
    - human: 人类基因组优化
    - mouse: 小鼠基因组特化
    - drosophila: 果蝇模型，短内含子
    - nematode: 线虫模型，紧凑基因组
    
    植物模型:
    - arabidopsis: 拟南芥模型，标准双子叶植物
    - maize: 玉米模型，大基因组单子叶
    - rice: 水稻模型，紧凑单子叶基因组
    - medicago: 苜蓿模型，豆科植物
    
    微生物模型:
    - fission_yeast: 裂殖酵母
    - aspergillus: 曲霉菌
    
    分析参数详解 | Analysis Parameters:
    
    基础参数:
    --min-score: 控制比对质量阈值
    - 0.5-0.6: 宽松，发现更多候选基因
    - 0.65: 默认值，平衡敏感性和特异性
    - 0.8+: 严格，高质量预测
    
    --gc-min-coverage: 全局链覆盖度
    - 30-40: 允许部分覆盖，适合EST数据
    - 50: 默认值，中等覆盖要求
    - 70+: 严格覆盖，适合高质量cDNA
    
    高级分析选项:
    --paralogs: 旁系同源基因检测
    - 识别基因家族成员
    - 分析基因复制事件
    - 适用于基因家族研究
    
    --intron-cutout: 内含子切除技术
    - 提高长内含子处理效率
    - 减少内存使用
    - 适合大基因组分析
    
    --fastdp: 快速动态规划
    - 显著提高计算速度
    - 略微降低精度
    - 适合大规模预测
    
    链方向控制:
    --forward-only / --reverse-only:
    - 单链分析，减少计算量
    - 用于方向性cDNA库
    - 提高特定方向的敏感性
    
    输出格式说明 | Output Format Description:
    
    标准输出文件:
    - genes.gff3: GFF3格式基因注释
    - transcripts.gff3: 转录本结构注释
    - analysis_summary.txt: 分析统计报告
    - analysis_summary.json: JSON格式统计
    
    XML输出 (--xml-output):
    - 详细的比对信息
    - 置信度分数
    - 剪接位点信息
    
    中间结果 (--intermediate):
    - 原始比对结果
    - 用于后续consensus分析
    - 支持多轮迭代优化
    
    性能优化策略 | Performance Optimization:
    
    线程配置:
    - 小基因组 (<100MB): 8-16线程
    - 中等基因组 (100MB-1GB): 32-64线程
    - 大基因组 (>1GB): 64-88线程
    
    内存管理:
    - 基础需求: 4-8GB
    - 大基因组: 16-32GB
    - 启用--intron-cutout可减少内存使用
    
    算法选择:
    - 高精度: 默认设置
    - 平衡模式: --fastdp + --intron-cutout
    - 高速模式: --fastdp + 较低--min-score
    
    质量控制建议 | Quality Control Recommendations:
    
    输入质量检查:
    - 验证FASTA格式完整性
    - 检查序列标识符唯一性
    - 确认物种模型匹配
    
    参数调优:
    - 根据序列类型调整阈值
    - 考虑基因组复杂度
    - 平衡敏感性与特异性
    
    结果验证:
    - 检查预测基因的合理性
    - 验证剪接位点准确性
    - 与已知注释比较
    
    常见问题解决 | Troubleshooting:
    
    序列格式问题:
    - 确保FASTA格式规范
    - 检查序列长度合理性
    - 验证字符编码正确
    
    内存不足:
    - 启用--intron-cutout
    - 减少线程数
    - 分区域处理大基因组
    
    预测质量差:
    - 调整--min-score阈值
    - 选择合适的物种模型
    - 检查输入序列质量
    
    运行速度慢:
    - 启用--fastdp
    - 使用--intron-cutout
    - 增加线程数
    
    高级应用技巧 | Advanced Usage Tips:
    
    多轮迭代预测:
    1. 初步预测获得候选基因
    2. 提取预测转录本作为新的cDNA输入
    3. 重新预测提高准确性
    
    结合其他工具:
    - 与AUGUSTUS等从头预测工具比较
    - 结合RNA-seq数据验证
    - 使用MAKER进行综合注释
    
    特殊应用:
    - 假基因识别: 降低阈值检测
    - 新基因发现: 使用多物种蛋白质
    - 剪接变体分析: 多cDNA库比对
    """
    
    # 参数验证
    if not any([cdna, protein, est]):
        click.echo("错误：必须至少提供一个序列文件 (--cdna, --protein, 或 --est)", err=True)
        sys.exit(1)
    
    if forward_only and reverse_only:
        click.echo("错误：不能同时指定 --forward-only 和 --reverse-only", err=True)
        sys.exit(1)
    
    # 懒加载：只在需要时导入重量级模块
    def lazy_import_gth():
        try:
            from ...genomethreader import GenomeThreaderAnalyzer
            return GenomeThreaderAnalyzer
        except ImportError as e:
            click.echo(f"无法导入GenomeThreader模块 | Failed to import GenomeThreader module: {e}", err=True)
            click.echo("请确保已正确安装GenomeThreader相关依赖 | Please ensure GenomeThreader dependencies are properly installed", err=True)
            sys.exit(1)
    
    try:
        # 延迟导入：只有在真正需要执行分析时才加载模块
        GenomeThreaderAnalyzer = lazy_import_gth()
        
        # 直接创建分析器，避免参数重构的性能损失
        analyzer = GenomeThreaderAnalyzer(
            genomic_file=genomic,
            output_dir=output,
            cdna_file=cdna,
            protein_file=protein,
            est_file=est,
            species=species,
            threads=threads,
            min_alignment_score=min_score,
            gc_min_coverage=gc_min_coverage,
            paralogs=paralogs,
            intron_cutout=intron_cutout,
            auto_intron_cutout=auto_intron_cutout,
            fastdp=fastdp,
            forward_only=forward_only,
            reverse_only=reverse_only,
            cdna_forward=cdna_forward,
            xml_output=xml_output,
            intermediate=intermediate,
            skip_alignment_out=skip_alignment_out,
            score_matrix=score_matrix,
            translation_table=translation_table,
            bssm_file=bssm_file,
            from_pos=from_pos,
            to_pos=to_pos,
            width=width,
            first_alignments=first_alignments,
            gth_path=gth_path,
            gthconsensus_path=gthconsensus_path
        )
        
        # 直接运行分析，无额外包装
        analyzer.run_analysis()
        
    except KeyboardInterrupt:
        click.echo("\nGenomeThreader分析被用户中断 | GenomeThreader analysis interrupted by user", err=True)
        raise click.Abort()
    except Exception as e:
        click.echo(f"GenomeThreader分析失败 | GenomeThreader analysis failed: {e}", err=True)
        raise click.Abort()


# 性能优化的懒加载管理器
class GTHLazyLoader:
    """GenomeThreader懒加载管理器 | GenomeThreader Lazy Loader Manager"""
    
    def __init__(self):
        self._analyzer_class = None
    
    def get_analyzer_class(self):
        """获取分析器类，首次调用时加载 | Get analyzer class, load on first call"""
        if self._analyzer_class is None:
            from ...genomethreader import GenomeThreaderAnalyzer
            self._analyzer_class = GenomeThreaderAnalyzer
        return self._analyzer_class
    
    def is_loaded(self):
        """检查模块是否已加载 | Check if module is loaded"""
        return self._analyzer_class is not None


# 全局懒加载器实例
_gth_loader = GTHLazyLoader()