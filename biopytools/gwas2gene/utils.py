"""
GWAS2Gene工具函数模块|GWAS2Gene Utility Functions Module
"""

import sys
from collections import defaultdict


def parse_gff(gff_path):
    """
    解析GFF3文件|Parse GFF3 file

    Args:
        gff_path: GFF3文件路径|GFF3 file path

    Returns:
        tuple: (genes列表, gene_to_transcripts字典, chrom_bounds字典)
    """
    genes = []
    gene_to_transcripts = defaultdict(list)
    chrom_bounds = defaultdict(lambda: [float('inf'), 0])

    with open(gff_path, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            line = line.rstrip('\n')
            if not line:
                continue

            parts = line.split('\t')
            if len(parts) < 9:
                continue

            chrom, _, feature, start, end, _, strand, _, attrs = parts
            start, end = int(start), int(end)

            # 更新染色体范围|Update chromosome bounds
            if start < chrom_bounds[chrom][0]:
                chrom_bounds[chrom][0] = start
            if end > chrom_bounds[chrom][1]:
                chrom_bounds[chrom][1] = end

            attr_dict = parse_attributes(attrs)

            if feature == 'gene':
                gene_id = attr_dict.get('ID', '')
                gene_name = attr_dict.get('Name', gene_id)
                genes.append({
                    'chrom': chrom,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'gene_id': gene_id,
                    'gene_name': gene_name,
                })

            elif feature in ('mRNA', 'transcript'):
                transcript_id = attr_dict.get('ID', '')
                parent = attr_dict.get('Parent', '')
                if parent and transcript_id:
                    gene_to_transcripts[parent].append(transcript_id)

    # 转换chrom_bounds为普通dict|Convert chrom_bounds to regular dict
    chrom_bounds = {c: (v[0], v[1]) for c, v in chrom_bounds.items()}
    return genes, gene_to_transcripts, chrom_bounds


def parse_attributes(attrs_str):
    """
    解析GFF3属性字段|Parse GFF3 attributes field

    Args:
        attrs_str: 属性字符串|Attribute string

    Returns:
        dict: key->value 字典
    """
    result = {}
    for item in attrs_str.split(';'):
        item = item.strip()
        if '=' in item:
            k, v = item.split('=', 1)
            result[k.strip()] = v.strip()
    return result


def parse_gwas(gwas_path, pval_col, threshold, logger=None):
    """
    解析GWAS结果，筛选P值小于阈值的SNP|Parse GWAS results and filter significant SNPs

    Args:
        gwas_path: GWAS文件路径|GWAS file path
        pval_col: P值列名或列索引|P-value column name or index
        threshold: P值阈值|P-value threshold
        logger: 日志器|Logger

    Returns:
        list: 显著SNP列表|List of significant SNPs
    """
    significant = []

    with open(gwas_path, 'r') as f:
        # 读取并自动检测分隔符|Read and detect separator
        header_line = f.readline().rstrip('\n')
        sep = detect_separator(header_line)
        headers = header_line.split(sep)
        headers = [h.strip() for h in headers]

        # 确定P值列索引|Determine P-value column index
        pval_idx = resolve_column(headers, pval_col, 'P值')

        # 尝试推断SNP ID列、染色体列、位置列|Try to guess SNP ID, chromosome, position columns
        snp_idx = guess_column(headers, ['snp', 'snpid', 'marker', 'rs', 'id', 'name'])
        chrom_idx = guess_column(headers, ['chr', 'chrom', 'chromosome', 'scaffold'])
        pos_idx = guess_column(headers, ['pos', 'position', 'bp', 'ps', 'snp_pos'])

        if chrom_idx is None or pos_idx is None:
            sys.exit(
                '错误：无法自动识别染色体列或位置列，'
                '请确保列名包含 chr/chrom 和 pos/position/bp 等关键词。'
            )

        if logger:
            logger.info(
                f"识别到列|Columns detected: "
                f"SNP={headers[snp_idx] if snp_idx is not None else '(无，用行号)'}, "
                f"Chrom={headers[chrom_idx]}, Pos={headers[pos_idx]}, Pval={headers[pval_idx]}"
            )

        for i, line in enumerate(f, start=2):
            line = line.rstrip('\n')
            if not line:
                continue
            cols = line.split(sep)
            cols = [c.strip() for c in cols]

            try:
                pval = float(cols[pval_idx])
            except (ValueError, IndexError):
                continue

            if pval < threshold:
                snp_id = cols[snp_idx] if snp_idx is not None else f'SNP_line{i}'
                try:
                    chrom = cols[chrom_idx]
                    pos = int(float(cols[pos_idx]))
                except (ValueError, IndexError):
                    if logger:
                        logger.warning(f"第{i}行坐标解析失败，跳过|Line {i}: coordinate parsing failed, skipped")
                    continue

                significant.append({
                    'snp_id': snp_id,
                    'chrom': chrom,
                    'pos': pos,
                    'pval': pval,
                })

    return significant


def detect_separator(line):
    """
    检测分隔符|Detect separator

    Args:
        line: 文本行|Text line

    Returns:
        str: 分隔符|Separator
    """
    if '\t' in line:
        return '\t'
    elif ',' in line:
        return ','
    elif ' ' in line:
        # 检测空格分隔（优先检测多个连续空格）|Detect space separator (multiple spaces preferred)
        if '  ' in line:  # 两个或更多空格|Two or more spaces
            return '  '  # 返回双空格作为分隔符|Return double space as separator
        else:
            return ' '  # 单个空格|Single space
    else:
        return None


def resolve_column(headers, col_spec, label):
    """
    通过列名或列索引找到列下标|Find column index by name or index

    Args:
        headers: 列名列表|List of column names
        col_spec: 列名或列索引|Column name or index (1-based)
        label: 标签（用于错误信息）|Label for error message

    Returns:
        int: 列索引|Column index
    """
    # 尝试作为整数（1-based）|Try as integer (1-based)
    try:
        idx = int(col_spec) - 1
        if 0 <= idx < len(headers):
            return idx
        else:
            sys.exit(f'错误：{label}列索引 {col_spec} 超出范围（共{len(headers)}列）')
    except ValueError:
        pass

    # 尝试作为列名（大小写不敏感）|Try as column name (case-insensitive)
    lower_headers = [h.lower() for h in headers]
    col_lower = col_spec.lower()
    if col_lower in lower_headers:
        return lower_headers.index(col_lower)

    sys.exit(f'错误：找不到{label}列 "{col_spec}"，请检查列名或使用列序号（1-based）')


def guess_column(headers, keywords):
    """
    根据关键词猜测列索引|Guess column index by keywords

    Args:
        headers: 列名列表|List of column names
        keywords: 关键词列表|List of keywords

    Returns:
        int or None: 列索引|Column index
    """
    lower = [h.lower() for h in headers]
    for kw in keywords:
        for i, h in enumerate(lower):
            if kw in h:
                return i
    return None


def parse_function_file(func_path, logger=None):
    """
    解析功能注释文件|Parse function annotation file

    Args:
        func_path: 功能注释文件路径|Function annotation file path
        logger: 日志器|Logger

    Returns:
        dict: id -> function 字典
    """
    func_map = {}
    with open(func_path, 'r') as f:
        for line in f:
            line = line.rstrip('\n')
            if not line or line.startswith('#'):
                continue
            sep = detect_separator(line) or '\t'
            parts = line.split(sep, 1)
            if len(parts) < 2:
                continue
            gene_id = parts[0].strip()
            func = parts[1].strip()
            func_map[gene_id] = func
    return func_map


def calc_distance(snp_pos, gene_start, gene_end, strand):
    """
    计算SNP到基因的距离，考虑strand方向|Calculate SNP-to-gene distance considering strand

    Args:
        snp_pos: SNP位置|SNP position
        gene_start: 基因起始位置|Gene start position
        gene_end: 基因终止位置|Gene end position
        strand: 正负链|Strand

    Returns:
        int: 距离（0=基因内，负=上游，正=下游）|Distance (0=in gene, negative=upstream, positive=downstream)
    """
    if gene_start <= snp_pos <= gene_end:
        return 0

    if strand == '+' or strand == '.':
        if snp_pos < gene_start:
            return snp_pos - gene_start  # 负值，上游|Negative, upstream
        else:
            return snp_pos - gene_end  # 正值，下游|Positive, downstream
    else:  # strand == '-'
        if snp_pos > gene_end:
            return gene_end - snp_pos  # 负值，上游|Negative, upstream
        else:
            return gene_start - snp_pos  # 正值，下游|Positive, downstream


def detect_chromosome_format(chrom_list):
    """
    检测染色体名称格式|Detect chromosome name format
    
    Args:
        chrom_list: 染色体名称列表|List of chromosome names
    
    Returns:
        str: 格式类型|Format type ('chr_prefix', 'numeric', 'other')
    """
    if not chrom_list:
        return 'other'
    
    # 检查是否有Chr前缀（不区分大小写）|Check for Chr prefix (case-insensitive)
    has_chr_prefix = any(chrom.lower().startswith('chr') for chrom in chrom_list if chrom.lower() != 'chr')
    
    # 检查是否主要是数字|Check if mainly numeric
    is_numeric = any(chrom.isdigit() or (chrom.startswith('0') and chrom[1:].isdigit()) for chrom in chrom_list)
    
    if has_chr_prefix:
        return 'chr_prefix'
    elif is_numeric:
        return 'numeric'
    else:
        return 'other'


def normalize_chromosome(chrom, target_format):
    """
    标准化染色体名称|Normalize chromosome name
    
    Args:
        chrom: 染色体名称|Chromosome name
        target_format: 目标格式|Target format ('chr_prefix', 'numeric')
    
    Returns:
        str: 标准化后的染色体名称|Normalized chromosome name
    """
    chrom_lower = chrom.lower()
    
    if target_format == 'chr_prefix':
        # 添加Chr前缀|Add Chr prefix
        if chrom_lower.startswith('chr'):
            return chrom
        else:
            return f'Chr{chrom}'
    
    elif target_format == 'numeric':
        # 移除Chr前缀|Remove Chr prefix
        if chrom_lower.startswith('chr'):
            # 处理 Chr1, Chr01, ChrX 等|Handle Chr1, Chr01, ChrX, etc.
            suffix = chrom[3:] if chrom_lower.startswith('chr') else chrom[3:]
            return suffix
        else:
            return chrom
    
    return chrom
