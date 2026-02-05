"""
矩阵文件写入模块|Matrix File Writing Module
"""

from pathlib import Path

class MatrixWriter:
    """矩阵文件写入器|Matrix File Writer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.output_files = {}
    
    def initialize_output_files(self, sample_names, snp_accepted, snp_biallelic):
        """初始化输出文件|Initialize output files"""
        outfile_base = str(Path(self.config.output_path, self.get_output_prefix()))
        
        #  PHYLIP格式|PHYLIP format
        if not self.config.phylip_disable:
            self.output_files['phylip'] = open(outfile_base + ".phy", "w")
            self.output_files['phylip'].write(f"{len(sample_names)} {snp_accepted}\n")
            self.logger.info(" PHYLIP输出已初始化|PHYLIP output initialized")
        
        #  FASTA格式|FASTA format
        if self.config.fasta:
            self.output_files['fasta'] = open(outfile_base + ".fasta", "w")
            self.logger.info(" FASTA输出已初始化|FASTA output initialized")
        
        #  NEXUS格式|NEXUS format
        if self.config.nexus:
            self.output_files['nexus'] = open(outfile_base + ".nexus", "w")
            self.output_files['nexus'].write(
                f"#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={len(sample_names)} NCHAR={snp_accepted};\n"
                f"\tFORMAT DATATYPE=DNA MISSING=N GAP=- ;\nMATRIX\n"
            )
            self.logger.info(" NEXUS输出已初始化|NEXUS output initialized")
        
        #  二进制NEXUS格式|Binary NEXUS format
        if self.config.nexus_binary:
            self.output_files['nexus_binary'] = open(outfile_base + ".bin.nexus", "w")
            self.output_files['nexus_binary'].write(
                f"#NEXUS\n\nBEGIN DATA;\n\tDIMENSIONS NTAX={len(sample_names)} NCHAR={snp_biallelic};\n"
                f"\tFORMAT DATATYPE=SNP MISSING=? GAP=- ;\nMATRIX\n"
            )
            self.logger.info(" 二进制NEXUS输出已初始化|Binary NEXUS output initialized")
    
    def get_output_prefix(self):
        """获取输出文件前缀|Get output file prefix"""
        if self.config.output_prefix:
            prefix = self.config.output_prefix
        else:
            parts = Path(self.config.input_file).name.split(".")
            prefix_parts = []
            for p in parts:
                if p.lower() == "vcf":
                    break
                else:
                    prefix_parts.append(p)
            prefix = ".".join(prefix_parts)
        
        return f"{prefix}.min{self.config.min_samples_locus}"
    
    def write_sequences(self, sample_names, temp_file, temp_bin_file=None):
        """写入序列数据|Write sequence data"""
        len_longest_name = max(len(name) for name in sample_names)
        
        # 处理外群|Handle outgroup
        idx_outgroup = None
        outgroup_name = self.config.outgroup.split(",")[0].split(";")[0] if self.config.outgroup else ""
        
        if outgroup_name in sample_names:
            idx_outgroup = sample_names.index(outgroup_name)
            self.logger.info(f" 外群样本已找到|Outgroup sample found: {outgroup_name}")
        
        # 写入外群序列|Write outgroup sequence
        if idx_outgroup is not None:
            self._write_sample_sequence(sample_names[idx_outgroup], idx_outgroup, 
                                      len_longest_name, temp_file, temp_bin_file, is_outgroup=True)
        
        # 写入其他样本序列|Write other sample sequences
        for s in range(len(sample_names)):
            if s != idx_outgroup:
                self._write_sample_sequence(sample_names[s], s, len_longest_name, 
                                          temp_file, temp_bin_file, sample_num=s+1, 
                                          total_samples=len(sample_names))
    
    def _write_sample_sequence(self, sample_name, sample_idx, len_longest_name, 
                             temp_file, temp_bin_file=None, is_outgroup=False, 
                             sample_num=None, total_samples=None):
        """写入单个样本序列|Write individual sample sequence"""
        #  处理核苷酸矩阵|Handle nucleotide matrices
        if temp_file and (not self.config.phylip_disable or self.config.fasta or self.config.nexus):
            with open(temp_file) as tmp_seq:
                seqout = ""
                for line in tmp_seq:
                    seqout += line[sample_idx]
            
            # FASTA格式|FASTA format
            if self.config.fasta:
                self.output_files['fasta'].write(f">{sample_name}\n{seqout}\n")
            
            # PHYLIP和NEXUS格式|PHYLIP and NEXUS formats
            padding = (len_longest_name + 3 - len(sample_name)) * " "
            if not self.config.phylip_disable:
                self.output_files['phylip'].write(f"{sample_name}{padding}{seqout}\n")
            if self.config.nexus:
                self.output_files['nexus'].write(f"{sample_name}{padding}{seqout}\n")
        
        #  处理二进制矩阵|Handle binary matrix
        if temp_bin_file and self.config.nexus_binary:
            with open(temp_bin_file) as bin_tmp_seq:
                seqout = ""
                for line in bin_tmp_seq:
                    seqout += line[sample_idx]
                
                padding = (len_longest_name + 3 - len(sample_name)) * " "
                self.output_files['nexus_binary'].write(f"{sample_name}{padding}{seqout}\n")
        
        # 进度日志|Progress logging
        if is_outgroup:
            self.logger.info(f" 外群序列已添加|Outgroup sequence added: {sample_name}")
        elif sample_num and total_samples:
            self.logger.info(f" 样本 {sample_num}/{total_samples} 已处理|Sample {sample_num}/{total_samples} processed: {sample_name}")
    
    def finalize_output_files(self):
        """完成输出文件写入|Finalize output file writing"""
        outfile_base = str(Path(self.config.output_path, self.get_output_prefix()))
        
        if self.config.nexus and 'nexus' in self.output_files:
            self.output_files['nexus'].write(";\nEND;\n")
            self.logger.info(f" NEXUS矩阵已保存|NEXUS matrix saved: {outfile_base}.nexus")
        
        if self.config.nexus_binary and 'nexus_binary' in self.output_files:
            self.output_files['nexus_binary'].write(";\nEND;\n")
            self.logger.info(f" 二进制NEXUS矩阵已保存|Binary NEXUS matrix saved: {outfile_base}.bin.nexus")
        
        if not self.config.phylip_disable and 'phylip' in self.output_files:
            self.logger.info(f" PHYLIP矩阵已保存|PHYLIP matrix saved: {outfile_base}.phy")
        
        if self.config.fasta and 'fasta' in self.output_files:
            self.logger.info(f" FASTA矩阵已保存|FASTA matrix saved: {outfile_base}.fasta")
        
        # 关闭所有文件|Close all files
        for file_handle in self.output_files.values():
            file_handle.close()
