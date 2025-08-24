"""
FASTA ID分割命令 | FASTA ID Splitting Command
"""

import click
import sys
from ...split_fasta_id.main import main as split_fasta_id_main


@click.command(short_help='FASTA序列ID分割和清理工具',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入FASTA文件路径 | Input FASTA file path')
@click.option('--output', '-o',
              default='output.fasta',
              type=click.Path(),
              help='输出FASTA文件路径 | Output FASTA file path (default: output.fasta)')
@click.option('--position', '-p',
              default=0,
              type=int,
              help='提取位置 (0表示第一个元素) | Extract position (0 means first element) (default: 0)')
@click.option('--delimiter', '-d',
              default='auto',
              type=str,
              help='分隔符类型 | Delimiter type: "auto"(自动检测), "space"(空格), "tab"(制表符), "both"(空格和制表符), 或任意字符如","、"|"等 (default: auto)')
@click.option('--keep-original',
              is_flag=True,
              help='保留原始文件作为备份 | Keep original file as backup')
@click.option('--no-skip-empty',
              is_flag=True,
              help='不跳过空的序列名称行 | Do not skip empty sequence name lines')
@click.option('--preserve-comments',
              is_flag=True,
              help='保留序列名称行中的注释部分 | Preserve comments in sequence name lines')
def split_fasta_id(input, output, position, delimiter, keep_original, 
                   no_skip_empty, preserve_comments):
    """
    FASTA序列ID分割脚本 (模块化版本)
    
    用于清理和标准化FASTA文件中的序列ID，支持多种分隔符和提取位置，
    可以有效处理复杂的序列标识符，生成简洁规范的序列名称。
    
    功能特点 | Features:
    - 智能分隔符自动检测
    - 灵活的位置提取选项
    - 批量序列ID处理
    - 原始文件备份保护
    - 注释信息保留选项
    - 空行和异常处理
    
    应用场景 | Use Cases:
    - BLAST数据库序列ID清理
    - 基因组注释文件标准化
    - 系统发育分析序列准备
    - 多序列比对前处理
    - 序列数据格式转换
    
    示例 | Examples:
    
    \b
    # 基本用法 - 提取第一个元素
    biopytools split-fasta-id -i input.fasta -o output.fasta -p 0
    
    \b
    # 提取第二个元素（索引为1）
    biopytools split-fasta-id -i sequences.fa -o clean_ids.fa -p 1 -d tab
    
    \b
    # 自动检测分隔符并保留原文件
    biopytools split-fasta-id -i data.fasta -o result.fasta -p 0 \\
        -d auto --keep-original
    
    \b
    # 使用逗号分隔符提取第三个元素
    biopytools split-fasta-id -i file.fasta -o out.fasta -p 2 -d ","
    
    \b
    # 使用管道符分隔并保留注释
    biopytools split-fasta-id -i seqs.fa -o clean.fa -p 1 -d "|" \\
        --preserve-comments
    
    \b
    # 处理复杂ID格式不跳过空行
    biopytools split-fasta-id -i complex.fasta -o simplified.fasta \\
        -p 0 -d "both" --no-skip-empty --keep-original
    
    分隔符类型说明 | Delimiter Types:
    - "auto": 自动检测常见分隔符（空格、制表符、|、等）
    - "space": 空格字符
    - "tab": 制表符
    - "both": 空格和制表符
    - 任意字符: 如 ",", "|", ":", "_", 等
    
    输入格式示例 | Input Format Examples:
    
    原始复杂ID:
    >gi|123456|ref|NM_001234.1| Homo sapiens gene
    >sp|P12345|PROT_HUMAN Protein name OS=Homo sapiens
    >gene_001 transcript_variant_1 chromosome_1
    
    处理后简化ID (position=0):
    >gi
    >sp
    >gene_001
    
    处理后简化ID (position=1, delimiter="|"):
    >123456
    >P12345
    >gene_001
    
    位置索引说明 | Position Index:
    对于ID: "gene_001|transcript_1|chr1"
    - position=0: 提取 "gene_001"
    - position=1: 提取 "transcript_1"  
    - position=2: 提取 "chr1"
    - position=-1: 提取最后一个元素 "chr1"
    
    处理选项详解 | Processing Options:
    
    --keep-original:
    - 创建原文件的 .bak 备份
    - 防止意外数据丢失
    - 便于结果比较验证
    
    --preserve-comments:
    - 保留ID行中除分割部分外的其他信息
    - 适用于需要保持部分注释的场景
    
    --no-skip-empty:
    - 处理可能产生空ID的特殊情况
    - 用于调试和特殊数据处理
    
    常见使用场景 | Common Use Cases:
    
    1. NCBI序列清理:
       >gi|123456|ref|NM_001234.1| → >NM_001234.1
       使用: -p 3 -d "|"
    
    2. UniProt序列标准化:
       >sp|P12345|PROT_HUMAN → >P12345  
       使用: -p 1 -d "|"
    
    3. 转录组数据清理:
       >ENST00000123456.1|ENSG00000654321.2 → >ENST00000123456.1
       使用: -p 0 -d "|"
    
    4. 基因组注释清理:
       >gene_001 transcript_variant_1 → >gene_001
       使用: -p 0 -d "space"
    
    错误处理 | Error Handling:
    - 自动检测并处理格式异常
    - 保护原始数据完整性
    - 详细的处理日志输出
    - 异常情况优雅降级
    
    性能特点 | Performance Features:
    - 内存高效的流式处理
    - 适用于大型FASTA文件
    - 快速的模式识别算法
    - 最小化磁盘I/O操作
    """
    
    # 构建参数列表传递给原始main函数
    args = ['split_fasta_id.py']
    
    # 必需参数
    args.extend(['-i', input])
    
    # 可选参数（只在非默认值时添加）
    if output != 'output.fasta':
        args.extend(['-o', output])
    
    if position != 0:
        args.extend(['-p', str(position)])
    
    if delimiter != 'auto':
        args.extend(['-d', delimiter])
    
    # 处理选项
    if keep_original:
        args.append('--keep-original')
    
    if no_skip_empty:
        args.append('--no-skip-empty')
    
    if preserve_comments:
        args.append('--preserve-comments')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        split_fasta_id_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n分割流程被用户中断 | Splitting pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"分割失败 | Splitting failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv