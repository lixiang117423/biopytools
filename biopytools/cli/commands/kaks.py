"""
Ka/Ks计算分析命令 | Ka/Ks Calculation Analysis Command
"""

import click
import sys
from ...kakscalc.main import main as kaks_main


@click.command(short_help = '计算Ks,Ks及其比值',
               context_settings=dict(help_option_names=['-h', '--help'],max_content_width=120))
@click.option('--fasta1', '-1',
              required=True,
              type=click.Path(exists=True),
              help='🧬 第一个FASTA文件 (物种1 CDS序列) | First FASTA file (species 1 CDS sequences)')
@click.option('--fasta2', '-2',
              required=True,
              type=click.Path(exists=True),
              help='🧬 第二个FASTA文件 (物种2 CDS序列) | Second FASTA file (species 2 CDS sequences)')
@click.option('--pairs', '-p',
              required=True,
              type=click.Path(exists=True),
              help='🔗 序列配对文件 (TSV/CSV格式) | Sequence pair file (TSV/CSV format)')
@click.option('--output', '-o',
              required=True,
              type=click.Path(),
              help='📂 输出目录 | Output directory for results')
@click.option('--method', '-m',
              default='GMYN',
              type=click.Choice([
                  'GMYN', 'MYN', 'YN', 'NG', 'LWL', 'LPB', 'MLWL', 'MLPB',
                  'GY', 'MS', 'MA', 'GNG', 'GLWL', 'GLPB', 'GMLWL', 'GMLPB', 'GYN'
              ]),
              help='🧮 计算方法 | Calculation method (default: GMYN)')
@click.option('--kaks-path',
              default='KaKs_Calculator',
              type=str,
              help='🛠️ KaKs_Calculator可执行文件路径 | Path to KaKs_Calculator executable')
@click.option('--verbose', '-v',
              is_flag=True,
              help='🔍 启用详细日志记录 | Enable verbose logging')
@click.option('--temp-dir',
              type=click.Path(),
              help='🗂️ 自定义临时目录 | Custom temporary directory')
@click.option('--keep-temp',
              is_flag=True,
              help='🗑️ 保留临时文件 (用于调试) | Keep temporary files (for debugging)')
def kaks(fasta1, fasta2, pairs, output, method, kaks_path, verbose, temp_dir, keep_temp):
    """
    计算Ka/Ks
    
    计算同源基因对的Ka/Ks比值，用于进化选择压力分析。
    Ka代表非同义替换率，Ks代表同义替换率，Ka/Ks(ω)反映选择压力强度。
    
    输入要求 | Input Requirements:
    - FASTA文件需包含CDS序列，且长度为3的倍数
    - 配对文件格式：序列ID1<tab>序列ID2 或 序列ID1,序列ID2
    - 序列ID必须在对应的FASTA文件中存在
    
    输出文件 | Output Files:
    - kaks_detailed.csv: CSV格式详细结果
    - kaks_detailed.tsv: TSV格式详细结果  
    - kaks_summary.xlsx: Excel格式汇总文件
    - summary_stats.json: JSON格式统计数据
    - kaks_calculator.log: 分析日志文件
    
    示例 | Examples:
    
    \b
    # 基本用法
    biopytools kaks -1 species1.fasta -2 species2.fasta \\
        -p pairs.txt -o results/
    
    \b
    # 指定计算方法
    biopytools kaks --fasta1 human.fa --fasta2 mouse.fa \\
        --pairs orthologs.tsv --output analysis/ --method YN
    
    \b
    # 详细模式分析
    biopytools kaks -1 seqs1.fa -2 seqs2.fa -p pairs.csv \\
        -o out/ -m GMYN --verbose
    
    \b
    # 自定义KaKs_Calculator路径
    biopytools kaks --fasta1 cds1.fa --fasta2 cds2.fa \\
        --pairs pairs.txt --output results/ \\
        --kaks-path /usr/local/bin/KaKs_Calculator
    
    \b
    # 使用自定义临时目录并保留临时文件
    biopytools kaks -1 seq1.fa -2 seq2.fa -p pairs.txt \\
        -o output/ --temp-dir /tmp/kaks_work --keep-temp
    
    \b
    # 大批量分析（推荐gamma-MYN方法）
    biopytools kaks -1 genome1_cds.fa -2 genome2_cds.fa \\
        -p ortholog_pairs.tsv -o large_analysis/ \\
        --method GMYN --verbose
    
    计算方法说明 | Method Description:
    - GMYN: Gamma-MYN方法，推荐用于大多数分析
    - YN: Yang-Nielsen方法，经典方法
    - NG: Nei-Gojobori方法，简单快速
    - MYN: Modified YN方法，改进的YN方法
    - 其他方法适用于特定研究需求
    
    生物学解释 | Biological Interpretation:
    - ω < 1: 纯化选择 (负选择)
    - ω = 1: 中性进化
    - ω > 1: 正选择 (适应性进化)
    """
    
    # 构建参数列表传递给原始main函数
    args = ['kaks.py']
    
    # 必需参数
    args.extend(['-1', fasta1])
    args.extend(['-2', fasta2])
    args.extend(['-p', pairs])
    args.extend(['-o', output])
    
    # 可选参数（只在非默认值时添加）
    if method != 'GMYN':
        args.extend(['-m', method])
    
    if kaks_path != 'KaKs_Calculator':
        args.extend(['--kaks-path', kaks_path])
    
    if verbose:
        args.append('--verbose')
    
    if temp_dir:
        args.extend(['--temp-dir', temp_dir])
    
    if keep_temp:
        args.append('--keep-temp')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        kaks_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n❌ 分析被用户中断 | Analysis interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"❌ 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv