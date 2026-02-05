"""
VCF文件解析模块|VCF File Parsing Module
"""

import random
from .utils import FileHandler, AMBIG, GEN_BIN

class VCFParser:
    """VCF文件解析器|VCF File Parser"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def extract_sample_names(self):
        """从VCF文件提取样本名称|Extract sample names from VCF file"""
        sample_names = []
        with FileHandler.safe_open(self.config.input_file) as vcf:
            for line in vcf:
                line = line.strip("\n")
                if line.startswith("#CHROM"):
                    record = line.split("\t")
                    sample_names = [record[i].replace("./", "") for i in range(9, len(record))]
                    break
        
        self.logger.info(f" 检测到样本数量|Detected samples: {len(sample_names)}")
        return sample_names
    
    def is_anomalous(self, record, num_samples):
        """检查当前记录的样本数是否与头文件描述一致|Check if record sample count matches header"""
        return bool(len(record) != num_samples + 9)
    
    def is_snp(self, record):
        """判断当前VCF记录是否为SNP|Determine if current VCF record is a SNP"""
        alt = record[4].replace("<NON_REF>", record[3])
        return bool(len(record[3]) == 1 and len(alt) - alt.count(",") == alt.count(",") + 1)
    
    def num_genotypes(self, record, num_samples):
        """获取VCF记录中的基因型数量|Get number of genotypes in VCF record"""
        missing = 0
        for i in range(9, num_samples + 9):
            if record[i].startswith("."):
                missing += 1
        return num_samples - missing
    
    def get_matrix_column(self, record, num_samples):
        """将VCF记录转换为系统发生学矩阵列|Transform VCF record into phylogenetic matrix column"""
        nt_dict = {str(0): record[3].replace("-","*").upper(), ".": "N"}
        alt = record[4].replace("-", "*").replace("<NON_REF>", nt_dict["0"])
        alt = alt.split(",")
        
        for n in range(len(alt)):
            nt_dict[str(n+1)] = alt[n]
        
        column = ""
        for i in range(9, num_samples + 9):
            geno_num = record[i].split(":")[0].replace("/", "").replace("|", "")
            try:
                geno_nuc = "".join(sorted(set([nt_dict[j] for j in geno_num])))
            except KeyError:
                return "malformed"
            
            if self.config.resolve_IUPAC:
                column += AMBIG[nt_dict[random.choice(geno_num)]]
            else:
                column += AMBIG[geno_nuc]
        
        return column
    
    def get_matrix_column_bin(self, record, num_samples):
        """获取二进制NEXUS格式的比对列|Get binary alignment column for NEXUS format"""
        column = ""
        for i in range(9, num_samples + 9):
            genotype = record[i].split(":")[0]
            column += GEN_BIN.get(genotype, "?")
        return column
