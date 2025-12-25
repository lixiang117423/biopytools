"""
文件解析模块 | File Parsing Module
"""

class GFFParser:
    """GFF文件解析器 | GFF File Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse_gene_location(self, gff_file: str, gene_id: str):
        """
        从GFF文件中解析目标基因的位置信息
        Parse target gene location from GFF file
        """
        self.logger.info(f"🔍 正在解析GFF文件 | Parsing GFF file: {gff_file}")
        self.logger.info(f"🎯 查找基因ID | Searching for gene ID: {gene_id}")
        
        with open(gff_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                chrom = fields[0]
                feature_type = fields[2]
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                attributes = fields[8]
                
                # 查找目标基因 | Search target gene
                if feature_type == 'gene' and gene_id in attributes:
                    # 验证是否是完全匹配 | Verify exact match
                    attr_dict = {}
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key] = value
                    
                    if attr_dict.get('ID') == gene_id or attr_dict.get('Name') == gene_id.replace('gene-', ''):
                        self.logger.info(f"✅ 找到目标基因 | Found target gene:")
                        self.logger.info(f"   📍 染色体 | Chromosome: {chrom}")
                        self.logger.info(f"   📍 起始位置 | Start: {start}")
                        self.logger.info(f"   📍 终止位置 | End: {end}")
                        self.logger.info(f"   📍 链方向 | Strand: {strand}")
                        return chrom, start, end, strand
        
        self.logger.error(f"❌ 未找到基因 | Gene not found: {gene_id}")
        return None, None, None, None

class VariantParser:
    """变异文件解析器 | Variant File Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse_variant_line(self, line: str, is_exonic: bool = False):
        """
        解析变异文件的一行 | Parse a line from variant file
        智能检测格式：如果第一个字段以"line"开头，则按exonic格式解析
        """
        fields = line.strip().split('\t')
        
        try:
            # 智能检测：如果第一个字段是"lineXXX"格式，按exonic处理
            if fields[0].startswith('line'):
                # 强制使用exonic格式解析
                if len(fields) < 8:
                    return None
                return {
                    'line_number': fields[0],
                    'variant_type': fields[1],
                    'gene_info': fields[2],
                    'chrom': fields[3],
                    'position': int(fields[4]),
                    'end_pos': int(fields[5]),
                    'ref': fields[6],
                    'alt': fields[7],
                    'raw_line': line.strip()
                }
            
            # 否则按照指定的格式解析
            if is_exonic:
                # exonic_variant_function格式 | exonic_variant_function format
                if len(fields) < 8:
                    return None
                return {
                    'line_number': fields[0],
                    'variant_type': fields[1],
                    'gene_info': fields[2],
                    'chrom': fields[3],
                    'position': int(fields[4]),
                    'end_pos': int(fields[5]),
                    'ref': fields[6],
                    'alt': fields[7],
                    'raw_line': line.strip()
                }
            else:
                # variant_function格式 | variant_function format
                if len(fields) < 7:
                    return None
                return {
                    'function': fields[0],
                    'gene_info': fields[1],
                    'chrom': fields[2],
                    'position': int(fields[3]),
                    'end_pos': int(fields[4]),
                    'ref': fields[5],
                    'alt': fields[6],
                    'raw_line': line.strip()
                }
        except (ValueError, IndexError) as e:
            # 记录导致错误的行以便调试
            self.logger.debug(f"跳过无法解析的行 | Skipping unparseable line: {line.strip()[:100]}")
            return None
