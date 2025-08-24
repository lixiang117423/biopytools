"""
VCF筛选命令 | VCF Filtering Command
"""

import click
import sys
from ...vcf_filter.main import main as vcf_filter_main


@click.command(short_help='VCF文件高性能筛选和格式转换',
               context_settings=dict(help_option_names=['-h', '--help'], max_content_width=120))
@click.option('--input', '-i',
              required=True,
              type=click.Path(exists=True),
              help='输入VCF文件路径 | Input VCF file path')
@click.option('--output', '-o',
              type=click.Path(),
              help='输出VCF文件路径 | Output VCF file path')
@click.option('--chr', '-c',
              required=True,
              type=str,
              help='染色体名称 (支持逗号分隔的多个染色体) | Chromosome name(s) (comma-separated for multiple)')
@click.option('--start', '-s',
              type=int,
              help='起始位置 | Start position')
@click.option('--end', '-e',
              type=int,
              help='结束位置 | End position')
@click.option('--convert-format',
              is_flag=True,
              help='使用PLINK进行格式转换 | Use PLINK for format conversion')
@click.option('--plink-path',
              default='plink',
              type=str,
              help='PLINK可执行文件路径 | PLINK executable path (default: plink)')
@click.option('--allow-extra-chr',
              is_flag=True,
              default=True,
              help='允许额外染色体 | Allow extra chromosomes (default: True)')
@click.option('--maf',
              type=float,
              help='最小等位基因频率 | Minimum allele frequency')
@click.option('--max-missing',
              type=float,
              help='最大缺失率 | Maximum missing rate')
@click.option('--quality-threshold',
              type=float,
              help='质量阈值 | Quality threshold')
@click.option('--min-depth',
              type=int,
              help='最小深度 | Minimum depth')
@click.option('--max-depth',
              type=int,
              help='最大深度 | Maximum depth')
@click.option('--keep-samples',
              type=str,
              help='保留的样本名称 (逗号分隔) | Sample names to keep (comma-separated)')
@click.option('--remove-samples',
              type=str,
              help='移除的样本名称 (逗号分隔) | Sample names to remove (comma-separated)')
@click.option('--keep-ids',
              type=str,
              help='保留的变异位点ID (逗号分隔) | Variant IDs to keep (comma-separated)')
@click.option('--remove-ids',
              type=str,
              help='移除的变异位点ID (逗号分隔) | Variant IDs to remove (comma-separated)')
@click.option('--biallelic-only',
              is_flag=True,
              help='只保留双等位基因位点 | Keep only biallelic sites')
@click.option('--remove-indels',
              is_flag=True,
              help='移除插入缺失变异 | Remove indel variants')
@click.option('--skip-validation',
              is_flag=True,
              default=True,
              help='跳过输入验证以提高速度（默认开启）| Skip input validation for speed (default enabled)')
@click.option('--force-validation',
              is_flag=True,
              help='强制执行输入验证 | Force input validation')
@click.option('--verbose', '-v',
              is_flag=True,
              help='显示详细信息 | Show verbose information')
def vcf_filter(input, output, chr, start, end, convert_format, plink_path, allow_extra_chr,
               maf, max_missing, quality_threshold, min_depth, max_depth,
               keep_samples, remove_samples, keep_ids, remove_ids,
               biallelic_only, remove_indels, skip_validation, force_validation, verbose):
    """
    VCF文件筛选工具 (模块化高性能版本)
    
    高效筛选VCF文件，支持染色体、位置、样本、质量等多种筛选条件，
    可进行PLINK格式转换，针对大数据集优化性能。
    
    功能特点 | Features:
    - 高性能VCF文件处理
    - 多种筛选条件组合
    - PLINK格式转换支持
    - 质量控制和样本管理
    - 内存优化和快速处理
    - 灵活的输出格式选项
    
    筛选条件 | Filter Conditions:
    - 染色体和基因组位置筛选
    - 等位基因频率和缺失率过滤
    - 测序深度和质量阈值
    - 样本包含/排除管理
    - 双等位基因和SNP专用模式
    
    示例 | Examples:
    
    \b
    # 基本区域筛选
    biopytools vcf-filter -i input.vcf -c chr1 -s 1000 -e 2000
    
    \b
    # PLINK格式转换与质控
    biopytools vcf-filter -i input.vcf -c chr1 --convert-format --maf 0.05
    
    \b
    # 多染色体双等位基因筛选
    biopytools vcf-filter -i input.vcf -c "chr1,chr2,chr3" \\
        --biallelic-only --verbose
    
    \b
    # 样本筛选
    biopytools vcf-filter -i input.vcf -c chr1 \\
        --keep-samples "sample1,sample2,sample3"
    
    \b
    # 综合质控筛选
    biopytools vcf-filter -i variants.vcf -c chr1 \\
        --maf 0.01 --max-missing 0.1 --quality-threshold 30 \\
        --min-depth 10 --max-depth 200 --remove-indels
    
    \b
    # 高性能大数据集处理
    biopytools vcf-filter -i large_dataset.vcf \\
        -c "1,2,3,4,5" --skip-validation \\
        --biallelic-only --maf 0.05 -o filtered.vcf
    
    \b
    # 特定变异位点筛选
    biopytools vcf-filter -i input.vcf -c chr1 \\
        --keep-ids "rs123456,rs789012" \\
        --quality-threshold 40
    
    \b
    # 样本排除和格式转换
    biopytools vcf-filter -i cohort.vcf -c chr22 \\
        --remove-samples "outlier1,outlier2" \\
        --convert-format --plink-path /usr/local/bin/plink
    
    筛选参数详解 | Filter Parameters Explained:
    
    位置筛选:
    - --chr: 染色体名称，支持多个(逗号分隔)
    - --start/--end: 基因组坐标范围
    
    质量控制:
    - --maf: 最小等位基因频率，过滤稀有变异
    - --max-missing: 最大缺失率，移除数据质量差的位点
    - --quality-threshold: QUAL字段阈值
    - --min-depth/--max-depth: DP字段深度范围
    
    样本管理:
    - --keep-samples: 保留指定样本
    - --remove-samples: 排除指定样本
    
    变异筛选:
    - --keep-ids/--remove-ids: 按变异ID筛选
    - --biallelic-only: 只保留双等位基因位点
    - --remove-indels: 移除插入缺失，只保留SNP
    
    性能优化:
    - --skip-validation: 跳过验证提升速度(默认)
    - --force-validation: 强制验证确保数据完整性
    
    格式转换 | Format Conversion:
    使用 --convert-format 可以将VCF转换为PLINK格式：
    - 生成 .bed, .bim, .fam 文件
    - 适用于GWAS和群体遗传学分析
    - 需要系统安装PLINK软件
    
    输出文件 | Output Files:
    - 标准模式: 筛选后的VCF文件
    - PLINK模式: .bed/.bim/.fam文件集
    - 日志文件: 详细的筛选统计信息
    
    常见使用场景 | Common Use Cases:
    
    1. GWAS数据预处理:
       移除稀有变异和低质量位点
       biopytools vcf-filter --maf 0.05 --max-missing 0.1
    
    2. 群体遗传学分析:
       提取特定染色体的双等位基因位点
       biopytools vcf-filter -c chr22 --biallelic-only
    
    3. 区域关联分析:
       提取候选基因区域变异
       biopytools vcf-filter -c chr1 -s 1000000 -e 2000000
    
    4. 样本质控:
       排除异常样本后重新分析
       biopytools vcf-filter --remove-samples "outlier1,outlier2"
    
    性能建议 | Performance Tips:
    - 大文件使用 --skip-validation 提升速度
    - 按染色体分批处理超大数据集
    - 使用SSD存储提高I/O性能
    - 合理设置质控参数平衡数据量和质量
    
    错误处理 | Error Handling:
    - 自动检测和处理格式问题
    - 详细的筛选统计报告
    - 优雅的异常处理和错误提示
    - 数据完整性验证选项
    """
    
    # 构建参数列表传递给原始main函数
    args = ['vcf_filter.py']
    
    # 必需参数
    args.extend(['-i', input])
    args.extend(['-c', chr])
    
    # 可选参数（只在提供时添加）
    if output:
        args.extend(['-o', output])
    
    if start is not None:
        args.extend(['-s', str(start)])
    
    if end is not None:
        args.extend(['-e', str(end)])
    
    # 格式转换参数
    if convert_format:
        args.append('--convert-format')
    
    if plink_path != 'plink':
        args.extend(['--plink-path', plink_path])
    
    if not allow_extra_chr:  # 默认为True，只在False时添加
        pass  # allow_extra_chr默认为True，这里不需要处理
    
    # 质量控制参数
    if maf is not None:
        args.extend(['--maf', str(maf)])
    
    if max_missing is not None:
        args.extend(['--max-missing', str(max_missing)])
    
    if quality_threshold is not None:
        args.extend(['--quality-threshold', str(quality_threshold)])
    
    if min_depth is not None:
        args.extend(['--min-depth', str(min_depth)])
    
    if max_depth is not None:
        args.extend(['--max-depth', str(max_depth)])
    
    # 样本筛选参数
    if keep_samples:
        args.extend(['--keep-samples', keep_samples])
    
    if remove_samples:
        args.extend(['--remove-samples', remove_samples])
    
    # 变异位点筛选参数
    if keep_ids:
        args.extend(['--keep-ids', keep_ids])
    
    if remove_ids:
        args.extend(['--remove-ids', remove_ids])
    
    # 变异筛选参数
    if biallelic_only:
        args.append('--biallelic-only')
    
    if remove_indels:
        args.append('--remove-indels')
    
    # 性能优化参数
    if not skip_validation or force_validation:
        if force_validation:
            args.append('--force-validation')
    
    if verbose:
        args.append('--verbose')
    
    # 保存并恢复sys.argv
    original_argv = sys.argv
    sys.argv = args
    
    try:
        # 调用原始的main函数
        vcf_filter_main()
    except SystemExit as e:
        # 处理程序正常退出
        if e.code != 0:
            sys.exit(e.code)
    except KeyboardInterrupt:
        click.echo("\n筛选流程被用户中断 | Filtering pipeline interrupted by user", err=True)
        sys.exit(1)
    except Exception as e:
        click.echo(f"筛选失败 | Filtering failed: {e}", err=True)
        sys.exit(1)
    finally:
        sys.argv = original_argv