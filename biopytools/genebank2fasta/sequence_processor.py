"""
序列处理模块 🧬 | Sequence Processing Module
"""

from Bio import SeqIO
from Bio.Seq import Seq
import os
from collections import defaultdict

class SequenceExtractor:
    """序列提取器 🔬 | Sequence Extractor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def process_single_file(self, gb_file):
        """处理单个GenBank文件 📄 | Process single GenBank file"""
        filename = os.path.basename(gb_file)
        sample_name = os.path.splitext(filename)[0]
        
        sample_sequences = {
            'sample': sample_name,
            'file': filename,
            'genes_extracted': 0,
            'genes_list': [],
            'cds_sequences': [],
            'protein_sequences': []
        }
        
        try:
            # 解析GenBank文件 | Parse GenBank file
            for record in SeqIO.parse(gb_file, "genbank"):
                organism = record.description or sample_name
                
                # 提取CDS特征 | Extract CDS features
                for feature in record.features:
                    if feature.type == "CDS":
                        seq_info = self.extract_cds_sequence(feature, record, sample_name, organism)
                        if seq_info:
                            sample_sequences['cds_sequences'].append(seq_info)
                            sample_sequences['protein_sequences'].append(seq_info)
                            sample_sequences['genes_extracted'] += 1
                            sample_sequences['genes_list'].append(seq_info['gene'])
        
        except Exception as e:
            self.logger.error(f"处理文件 {gb_file} 时出错 ❌: {e}")
        
        return sample_sequences
    
    def extract_cds_sequence(self, feature, record, sample_name, organism):
        """提取CDS序列 🧪 | Extract CDS sequence"""
        # 获取基因信息 | Get gene information
        gene_name = feature.qualifiers.get('gene', ['unknown'])[0]
        product = feature.qualifiers.get('product', ['unknown product'])[0]
        
        # 跳过unknown基因 | Skip unknown genes
        if self.config.skip_unknown_genes and gene_name == 'unknown':
            return None
        
        try:
            # 提取CDS序列 | Extract CDS sequence
            cds_sequence = feature.location.extract(record.seq)
            
            # 检查序列完整性 | Check sequence integrity
            if len(cds_sequence) % 3 != 0:
                self.logger.warning(f"⚠️ {gene_name} CDS长度不是3的倍数 ({len(cds_sequence)} bp)")
            
            # 获取翻译表 | Get translation table
            transl_table = feature.qualifiers.get('transl_table', ['11'])[0]
            
            # 翻译成蛋白质序列 | Translate to protein sequence
            try:
                protein_sequence = cds_sequence.translate(table=int(transl_table), to_stop=True)
            except:
                protein_sequence = cds_sequence.translate(table=11, to_stop=True)
            
            # 质量检查 | Quality check
            if len(protein_sequence) < self.config.min_protein_length:
                self.logger.warning(f"⚠️ {gene_name} 蛋白质序列过短 ({len(protein_sequence)} aa)")
                return None
            
            # 存储序列信息 | Store sequence information
            seq_info = {
                'sample': sample_name,
                'organism': organism,
                'gene': gene_name,
                'product': product,
                'cds_seq': str(cds_sequence),
                'protein_seq': str(protein_sequence),
                'cds_length': len(cds_sequence),
                'protein_length': len(protein_sequence)
            }
            
            return seq_info
            
        except Exception as e:
            self.logger.error(f"提取 {gene_name} 序列时出错 ❌: {e}")
            return None
