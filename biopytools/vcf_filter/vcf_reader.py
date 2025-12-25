"""
VCF文件读取模块 | VCF File Reading Module
"""

import gzip
import os
from typing import List, Dict, Tuple, Iterator

class VCFReader:
    """VCF文件读取器 | VCF File Reader"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def read_vcf_header(self) -> Tuple[List[str], List[str]]:
        """读取VCF文件头信息 | Read VCF file header"""
        header_lines = []
        samples = []
        
        if self.config.verbose:
            self.logger.info(f"正在读取VCF文件头 | Reading VCF file header: {self.config.vcf_file}")
        
        opener = gzip.open if self.config.vcf_file.endswith('.gz') else open
        
        with opener(self.config.vcf_file, 'rt') as f:
            for line in f:
                if line.startswith('##'):
                    header_lines.append(line.strip())
                elif line.startswith('#CHROM'):
                    header_lines.append(line.strip())
                    # 提取样本名称 | Extract sample names
                    fields = line.strip().split('\t')
                    if len(fields) > 9:
                        samples = fields[9:]
                    break
        
        if self.config.verbose:
            self.logger.info(f"读取到 {len(samples)} 个样本 | Read {len(samples)} samples")
        
        return header_lines, samples
    
    def iterate_variants(self) -> Iterator[Dict]:
        """迭代读取变异位点 | Iterate through variants"""
        opener = gzip.open if self.config.vcf_file.endswith('.gz') else open
        
        with opener(self.config.vcf_file, 'rt') as f:
            # 跳过头部信息 | Skip header
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 8:
                    continue
                
                variant = {
                    'chrom': fields[0],
                    'pos': int(fields[1]),
                    'id': fields[2],
                    'ref': fields[3],
                    'alt': fields[4],
                    'qual': fields[5],
                    'filter': fields[6],
                    'info': fields[7],
                    'format': fields[8] if len(fields) > 8 else '',
                    'samples': fields[9:] if len(fields) > 9 else [],
                    'raw_line': line.strip()
                }
                
                yield variant
    
    def count_variants(self) -> int:
        """统计变异位点数量 | Count variant sites"""
        count = 0
        for _ in self.iterate_variants():
            count += 1
        return count
