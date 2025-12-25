"""
基因序列提取核心模块 🧬 | Gene Sequence Extraction Core Module
"""

import os
from typing import Dict, List, Tuple
from .utils import SequenceUtils, FileValidator

class GenomeLoader:
    """基因组加载器 📖 | Genome Loader"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def load_genome(self, genome_file: str) -> Dict[str, str]:
        """加载基因组序列 🧬 | Load genome sequences"""
        self.logger.info(f"📖 正在加载基因组文件 | Loading genome file: {genome_file}")
        
        if not FileValidator.validate_fasta_file(genome_file):
            raise ValueError(f"❌ 无效的FASTA格式 | Invalid FASTA format: {genome_file}")
        
        sequences = {}
        current_seq_id = None
        current_seq = []
        
        try:
            with open(genome_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue
                    
                    if line.startswith('>'):
                        # 保存前一个序列 | Save previous sequence
                        if current_seq_id and current_seq:
                            sequences[current_seq_id] = ''.join(current_seq)
                        
                        # 开始新序列 | Start new sequence
                        current_seq_id = line[1:].split()[0]  # 取第一个空格前的部分作为ID
                        current_seq = []
                    else:
                        # 序列行 | Sequence line
                        current_seq.append(line.upper())
                
                # 保存最后一个序列 | Save last sequence
                if current_seq_id and current_seq:
                    sequences[current_seq_id] = ''.join(current_seq)
        
        except Exception as e:
            raise Exception(f"❌ 读取基因组文件失败 | Failed to read genome file: {e}")
        
        self.logger.info(f"✅ 成功加载 {len(sequences)} 个染色体/contig | Successfully loaded {len(sequences)} chromosomes/contigs")
        return sequences

class GFFParser:
    """GFF文件解析器 📄 | GFF File Parser"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def parse_gff(self, gff_file: str, feature_type: str = 'gene') -> List[Dict]:
        """解析GFF文件 📄 | Parse GFF file"""
        self.logger.info(f"📖 正在解析GFF文件 | Parsing GFF file: {gff_file}")
        
        if not FileValidator.validate_gff_file(gff_file):
            raise ValueError(f"❌ 无效的GFF格式 | Invalid GFF format: {gff_file}")
        
        features = []
        
        try:
            with open(gff_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    
                    # 跳过注释行和空行 | Skip comment lines and empty lines
                    if not line or line.startswith('#'):
                        continue
                    
                    # 解析GFF行 | Parse GFF line
                    fields = line.split('\t')
                    if len(fields) != 9:
                        continue
                    
                    seqid, source, ftype, start, end, score, strand, phase, attributes = fields
                    
                    # 只提取指定类型的特征 | Only extract specified feature type
                    if ftype.lower() != feature_type.lower():
                        continue
                    
                    # 解析属性 | Parse attributes
                    attr_dict = {}
                    for attr in attributes.split(';'):
                        if '=' in attr:
                            key, value = attr.split('=', 1)
                            attr_dict[key] = value
                    
                    # 创建特征字典 | Create feature dictionary
                    feature = {
                        'seqid': seqid,
                        'source': source,
                        'type': ftype,
                        'start': int(start),
                        'end': int(end),
                        'score': score,
                        'strand': strand,
                        'phase': phase,
                        'attributes': attr_dict,
                        'line_num': line_num
                    }
                    
                    features.append(feature)
        
        except Exception as e:
            raise Exception(f"❌ 解析GFF文件失败 | Failed to parse GFF file: {e}")
        
        self.logger.info(f"✅ 成功解析 {len(features)} 个 {feature_type} 特征 | Successfully parsed {len(features)} {feature_type} features")
        return features

class GeneExtractor:
    """基因序列提取器 🧬 | Gene Sequence Extractor"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def extract_genes(self, genome_seqs: Dict[str, str], features: List[Dict], 
                     min_length: int = 0, verbose: bool = False) -> List[Tuple[str, str]]:
        """提取基因序列 🧬 | Extract gene sequences"""
        self.logger.info(f"🔬 正在提取基因序列 | Extracting gene sequences...")
        
        extracted_genes = []
        skipped_count = 0
        
        for i, feature in enumerate(features, 1):
            seqid = feature['seqid']
            start = feature['start']
            end = feature['end']
            strand = feature['strand']
            attributes = feature['attributes']
            
            # 检查染色体是否存在 | Check if chromosome exists
            if seqid not in genome_seqs:
                if verbose:
                    self.logger.warning(f"⚠️ 跳过基因 {i}: 染色体 {seqid} 不存在于基因组中 | Skip gene {i}: chromosome {seqid} not found in genome")
                skipped_count += 1
                continue
            
            # 检查坐标是否有效 | Check if coordinates are valid
            seq_length = len(genome_seqs[seqid])
            if start > seq_length or end > seq_length:
                if verbose:
                    self.logger.warning(f"⚠️ 跳过基因 {i}: 坐标超出染色体范围 | Skip gene {i}: coordinates exceed chromosome range ({start}-{end} > {seq_length})")
                skipped_count += 1
                continue
            
            # 提取序列 | Extract sequence
            gene_seq = SequenceUtils.extract_sequence_region(genome_seqs[seqid], start, end, strand)
            gene_length = len(gene_seq)
            
            # 长度过滤 | Length filtering
            if gene_length < min_length:
                if verbose:
                    self.logger.warning(f"⚠️ 跳过基因 {i}: 长度太短 | Skip gene {i}: too short ({gene_length} < {min_length})")
                skipped_count += 1
                continue
            
            # 构建基因ID | Build gene ID
            gene_id = attributes.get('ID', f"gene_{i}")
            gene_name = attributes.get('Name', gene_id)
            
            # 构建FASTA头部 | Build FASTA header
            header = f"{gene_id}"
            if gene_name != gene_id:
                header += f" {gene_name}"
            header += f" {seqid}:{start}-{end}({strand}) length={gene_length}"
            
            extracted_genes.append((header, gene_seq))
            
            if verbose and i % 100 == 0:
                self.logger.info(f"📊 已处理 {i}/{len(features)} 个基因 | Processed {i}/{len(features)} genes...")
        
        self.logger.info(f"✅ 成功提取 {len(extracted_genes)} 个基因序列 | Successfully extracted {len(extracted_genes)} gene sequences")
        if skipped_count > 0:
            self.logger.info(f"⚠️ 跳过 {skipped_count} 个基因 | Skipped {skipped_count} genes")
        
        return extracted_genes

class FASTAWriter:
    """FASTA文件写入器 💾 | FASTA File Writer"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def write_fasta(self, sequences: List[Tuple[str, str]], output_file: str, line_width: int = 60):
        """写入FASTA文件 💾 | Write FASTA file"""
        self.logger.info(f"💾 正在写入输出文件 | Writing output file: {output_file}")
        
        try:
            with open(output_file, 'w') as f:
                for header, seq in sequences:
                    f.write(f">{header}\n")
                    # 按指定宽度写入序列 | Write sequence with specified width
                    formatted_seq = SequenceUtils.format_fasta_sequence(seq, line_width)
                    f.write(f"{formatted_seq}\n")
        
        except Exception as e:
            raise Exception(f"❌ 写入输出文件失败 | Failed to write output file: {e}")
        
        self.logger.info(f"✅ 成功写入 {len(sequences)} 个基因序列到 {output_file} | Successfully wrote {len(sequences)} gene sequences to {output_file}")
