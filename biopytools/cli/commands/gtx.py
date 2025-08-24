"""
GTX WGS分析命令 | GTX WGS Analysis Command
"""

import click
import sys
from ...gtx.main import main as gtx_main


@click.command(short_help = '运行GTX WGS流程',
               context_settings=dict(help_option_names=['-h', '--help'],max_content_width=120))
@click.option('--input-dir', '-i',
              required=True,
              type=click.Path(exists=True),
              help='📂 输入目录路径 (包含clean FASTQ文件) | Input directory path (containing clean FASTQ files)')
@click.option('--output-dir', '-o',
              required=True,
              type=click.Path(),
              help='📤 输出目录路径 | Output directory path')
@click.option('--reference', '-r',
              required=True,
              type=click.Path(exists=True),
              help='🧬 参考基因组文件路径 | Reference genome file path')
@click.option('--threads', '-t',
              default=88,
              type=int,
              help='🧵 线程数 | Number of threads (default: 88)')
@click.option('--gtx-path',
              default='/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx',
              type=click.Path(),
              help='💻 GTX程序路径 | GTX program path (default: /share/apps/gtx/GTX.CAT_2.2.1/bin/gtx)')
@click.option('--tmp-dir',
              type=click.Path(),
              help='🗑️ 临时目录路径 (默认使用输出目录下的tmp) | Temporary directory path (default: tmp under output directory)')
@click.option('--min-confidence',
              default=30,
              type=int,
              help='🎯 最小置信度阈值 | Minimum confidence threshold (default: 30)')
@click.option('--min-base-quality',
              default=20,
              type=int,
              help='✨ 最小碱基质量阈值 | Minimum base quality threshold (default: 20)')
@click.option('--ploidy',
              default=2,
              type=int,
              help='🧬 倍性 | Ploidy (default: 2)')
@click.option('--pcr-indel-model',
              default='CONSERVATIVE',
              type=str,
              help='🔬 PCR indel模型 | PCR indel model (default: CONSERVATIVE)')
@click.option('--read1-pattern',
              default='*_1.fq.gz',
              type=str,
              help='📄 R1文件匹配模式 | R1 file pattern (default: *_1.fq.gz)')
@click.option('--read2-pattern',
              default='*_2.fq.gz',
              type=str,
              help='📄 R2文件匹配模式 | R2 file pattern (default: *_2.fq.gz)')
def gtx(input_dir, output_dir, reference, threads, gtx_path, tmp_dir,
        min_confidence, min_base_quality, ploidy, pcr_indel_model,
        read1_pattern, read2_pattern):
    """
    🧬 GTX WGS批处理分析工具
    
    使用GTX软件对FASTQ文件进行全基因组测序分析，
    包括比对、变异检测和质量控制等步骤。
    
    示例 | Examples:
    
    \b
    # 🚀 基本分析
    biopytools gtx -i /path/to/clean/data -o /path/to/output -r /path/to/reference.fa
    
    \b
    # ⚡ 高线程数分析
    biopytools gtx -i ./clean_data -o ./gtx_results -r ./reference.fa -t 64
    
    \b
    # 💻 自定义GTX路径
    biopytools gtx -i /data/clean -o /results -r /genome/ref.fa --gtx-path /custom/gtx/path
    
    \b
    # 🔬 自定义质量控制参数
    biopytools gtx -i ./data -o ./results -r ./ref.fa \\
        --min-confidence 25 --min-base-quality 15 --ploidy 2
    
    \b
    # 📁 自定义文件匹配模式
    biopytools gtx -i ./data -o ./results -r ./ref.fa \\
        --read1-pattern "*_R1.fq.gz" --read2-pattern "*_R2.fq.gz"
    
    \b
    # 🗑️ 指定临时目录
    biopytools gtx -i ./data -o ./results -r ./ref.fa \\
        --tmp-dir /tmp/gtx_analysis --threads 32
    
    \b
    # 🔧 完整参数示例
    biopytools gtx -i /data/fastq -o /results/gtx \\
        -r /genome/reference.fa -t 128 \\
        --gtx-path /opt/gtx/bin/gtx \\
        --min-confidence 40 --min-base-quality 25 \\
        --pcr-indel-model AGGRESSIVE \\
        --read1-pattern "*_1.clean.fq.gz" \\
        --read2-pattern "*_2.clean.fq.gz"
    """
    
    # 构建参数列表传递给原始main函数
    args = ['gtx.py']
    
    # 必需参数
    args.extend(['-i', input_dir])
    args.extend(['-o', output_dir])
    args.extend(['-r', reference])
    
    # 可选参数（只在非默认值时添加）
    if threads != 88:
        args.extend(['-t', str(threads)])
    
    if gtx_path != '/share/apps/gtx/GTX.CAT_2.2.1/bin/gtx':
        args.extend(['--gtx-path', gtx_path])
    
    if tmp_dir:
        args.extend(['--tmp-dir', tmp_dir])
    
    # 质量控制参数
    if min_confidence != 30:
        args.extend(['--min-confidence', str(min_confidence)])
    
    if min_base_quality != 20:
        args.extend(['--min-base-quality', str(min_base_quality)])
    
    if ploidy != 2:
        args.extend(['--ploidy', str(ploidy)])
    
    if pcr_indel_model != 'CONSERVATIVE':
        args.extend(['--pcr-indel-model', pcr_indel_model])
    
    # 文件模式参数
    if read1_pattern != '*_1.fq.gz':
        args.extend(['--read1-pattern', read1_pattern])
    
    if read2_pattern != '*_2.fq.gz':
        args.extend(['--read2-pattern', read2_pattern])
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        gtx_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n⚠️ 用户中断操作 | User interrupted", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"💥 运行错误 | Runtime error: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv