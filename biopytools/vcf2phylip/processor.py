"""
VCF处理核心模块|VCF Processing Core Module
"""

import os
from pathlib import Path
from .utils import FileHandler

class VCFProcessor:
    """VCF处理器|VCF Processor"""
    
    def __init__(self, config, logger, parser, writer):
        self.config = config
        self.logger = logger
        self.parser = parser
        self.writer = writer
        self.temp_files = []
    
    def process_vcf(self, sample_names, num_samples):
        """处理VCF文件|Process VCF file"""
        self.logger.info(" 开始处理VCF基因型数据|Starting VCF genotype processing")
        
        # 准备临时文件|Prepare temporary files
        outfile_base = str(Path(self.config.output_path, self.writer.get_output_prefix()))
        temp_file = None
        temp_bin_file = None
        used_sites_file = None
        
        #  核苷酸矩阵临时文件|Nucleotide matrix temporary file
        if not self.config.phylip_disable or self.config.fasta or self.config.nexus:
            temp_file = outfile_base + ".tmp"
            self.temp_files.append(temp_file)
            temp_handle = open(temp_file, "w")
        
        #  二进制矩阵临时文件|Binary matrix temporary file
        if self.config.nexus_binary:
            temp_bin_file = outfile_base + ".bin.tmp"
            self.temp_files.append(temp_bin_file)
            temp_bin_handle = open(temp_bin_file, "w")
        
        #  使用位点记录文件|Used sites record file
        if self.config.write_used_sites:
            used_sites_file = outfile_base + ".used_sites.tsv"
            used_sites_handle = open(used_sites_file, "w")
            used_sites_handle.write("#CHROM\tPOS\tNUM_SAMPLES\n")
        
        # 处理统计变量|Processing statistics
        stats = {
            'snp_num': 0,
            'snp_accepted': 0,
            'snp_shallow': 0,
            'mnp_num': 0,
            'snp_biallelic': 0
        }
        
        #  处理VCF文件|Process VCF file
        with FileHandler.safe_open(self.config.input_file) as vcf:
            while True:
                # 分块读取大文件|Read large chunks of file
                vcf_chunk = vcf.readlines(50000)
                if not vcf_chunk:
                    break
                
                for line in vcf_chunk:
                    line = line.strip()
                    
                    if line and not line.startswith("#"):
                        record = line.split("\t")
                        stats['snp_num'] += 1
                        
                        # 每50万行打印进度|Print progress every 500k lines
                        if stats['snp_num'] % 500000 == 0:
                            self.logger.info(f" 已处理基因型|Genotypes processed: {stats['snp_num']}")
                        
                        # 检查记录完整性|Check record integrity
                        if self.parser.is_anomalous(record, num_samples):
                            self.logger.warning(f" 跳过格式错误行|Skipping malformed line: {line[:100]}...")
                            continue
                        
                        # 检查最小样本数|Check minimum sample count
                        num_samples_locus = self.parser.num_genotypes(record, num_samples)
                        if num_samples_locus < self.config.min_samples_locus:
                            stats['snp_shallow'] += 1
                            continue
                        
                        # 检查是否为SNP|Check if it's a SNP
                        if self.parser.is_snp(record):
                            #  处理核苷酸矩阵|Process nucleotide matrices
                            if temp_file:
                                site_tmp = self.parser.get_matrix_column(record, num_samples)
                                if site_tmp == "malformed":
                                    self.logger.warning(f" 跳过格式错误行|Skipping malformed line")
                                    continue
                                
                                stats['snp_accepted'] += 1
                                temp_handle.write(site_tmp + "\n")
                                
                                if used_sites_file:
                                    used_sites_handle.write(f"{record[0]}\t{record[1]}\t{num_samples_locus}\n")
                            
                            #  处理二进制NEXUS|Process binary NEXUS
                            if temp_bin_file and len(record[4]) == 1:
                                stats['snp_biallelic'] += 1
                                binsite_tmp = self.parser.get_matrix_column_bin(record, num_samples)
                                temp_bin_handle.write(binsite_tmp + "\n")
                        else:
                            stats['mnp_num'] += 1
        
        # 关闭临时文件句柄|Close temporary file handles
        if temp_file:
            temp_handle.close()
        if temp_bin_file:
            temp_bin_handle.close()
        if used_sites_file:
            used_sites_handle.close()
            self.logger.info(f" 使用位点已保存|Used sites saved: {used_sites_file}")
        
        # 输出统计信息|Output statistics
        self.logger.info(f" 处理统计|Processing Statistics:")
        self.logger.info(f"  总基因型数|Total genotypes: {stats['snp_num']}")
        self.logger.info(f"  缺失数据过多被排除|Excluded (missing data): {stats['snp_shallow']}")
        self.logger.info(f"  MNP被排除|Excluded (MNPs): {stats['mnp_num']}")
        self.logger.info(f"  通过筛选的SNP|SNPs passed filters: {stats['snp_accepted']}")
        if self.config.nexus_binary:
            self.logger.info(f"  二等位SNP(用于二进制NEXUS)|Biallelic SNPs (binary NEXUS): {stats['snp_biallelic']}")
        
        return stats, temp_file, temp_bin_file
    
    def cleanup_temp_files(self):
        """清理临时文件|Cleanup temporary files"""
        for temp_file in self.temp_files:
            if os.path.exists(temp_file):
                Path(temp_file).unlink()
                self.logger.debug(f" 临时文件已删除|Temporary file deleted: {temp_file}")
