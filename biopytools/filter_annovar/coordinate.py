"""
坐标转换模块 | Coordinate Conversion Module
"""

class CoordinateConverter:
    """基因组坐标转换器 | Genomic Coordinate Converter"""
    
    @staticmethod
    def calculate_gene_coordinate(genomic_pos: int, gene_start: int, gene_end: int, strand: str):
        """
        将基因组坐标转换为基因坐标（1-based）
        Convert genomic coordinate to gene coordinate (1-based)
        
        返回 | Returns: (gene_coord, location_type)
            gene_coord: 基因坐标或相对距离 | Gene coordinate or relative distance
            location_type: 'gene', 'upstream', 'downstream'
        
        符号约定 | Sign convention:
            上游(upstream): 负数 | Negative numbers
            下游(downstream): 正数 | Positive numbers
        
        链方向处理 | Strand handling:
            正链(+): 上游在左侧(< gene_start), 下游在右侧(> gene_end)
            负链(-): 上游在右侧(> gene_end), 下游在左侧(< gene_start)
        """
        if genomic_pos < gene_start:
            # 基因左侧 | Left of gene
            distance = gene_start - genomic_pos
            if strand == '+':
                # 正链: 左侧是上游 | Forward strand: left is upstream
                return f"-{distance}", "upstream"
            else:
                # 负链: 左侧是下游 | Reverse strand: left is downstream
                return f"+{distance}", "downstream"
        elif genomic_pos > gene_end:
            # 基因右侧 | Right of gene
            distance = genomic_pos - gene_end
            if strand == '+':
                # 正链: 右侧是下游 | Forward strand: right is downstream
                return f"+{distance}", "downstream"
            else:
                # 负链: 右侧是上游 | Reverse strand: right is upstream
                return f"-{distance}", "upstream"
        else:
            # 基因内部 | Within gene
            if strand == '+':
                # 正链: 从5'端(左侧)开始计数 | Forward: count from 5' end (left)
                gene_coord = genomic_pos - gene_start + 1
            else:
                # 负链: 从5'端(右侧)开始计数 | Reverse: count from 5' end (right)
                gene_coord = gene_end - genomic_pos + 1
            return str(gene_coord), "gene"
