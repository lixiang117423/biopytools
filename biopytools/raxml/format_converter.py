
"""
🔄 序列格式转换模块 | Sequence Format Conversion Module
"""

import os
import shutil
import tempfile
from pathlib import Path
from typing import Optional
from .utils import CommandRunner

class FormatConverter:
    """格式转换器 | Format Converter"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def convert_to_phylip(self) -> str:
        """将输入文件转换为PHYLIP格式 | Convert input file to PHYLIP format"""
        input_format = self.config.input_format
        
        self.logger.info(f"🔄 检测到输入格式 | Detected input format: {input_format}")
        
        if input_format == 'phylip':
            self.logger.info("✅ 输入已是PHYLIP格式，无需转换 | Input is already PHYLIP format, no conversion needed")
            return self.config.sequence_file
        
        elif input_format == 'vcf':
            return self._convert_vcf_to_phylip()
        
        elif input_format == 'fasta':
            return self._convert_fasta_to_phylip()
        
        else:
            self.logger.error(f"❌ 不支持的输入格式 | Unsupported input format: {input_format}")
            raise ValueError(f"Unsupported input format: {input_format}")
    
    def _convert_vcf_to_phylip(self) -> str:
        """VCF转PHYLIP | Convert VCF to PHYLIP"""
        self.logger.info("🧬 开始VCF到PHYLIP格式转换 | Starting VCF to PHYLIP conversion")
        
        # 检查vcf2phylip脚本是否存在 | Check if vcf2phylip script exists
        vcf2phylip_script = self.config.output_path / "vcf2phylip.py"
        
        # 如果脚本不存在，创建它 | If script doesn't exist, create it
        if not vcf2phylip_script.exists():
            self._create_vcf2phylip_script(vcf2phylip_script)
        
        # 构建vcf2phylip命令 | Build vcf2phylip command
        phylip_output = self.config.output_path / f"{self.config.output_name}_converted.phy"
        
        cmd_parts = [
            f"python {vcf2phylip_script}",
            f"-i {self.config.sequence_file}",
            f"--output-folder {self.config.output_path}",
            f"--output-prefix {self.config.output_name}_converted",
            f"-m {self.config.min_samples_locus}"
        ]
        
        # 添加可选参数 | Add optional parameters
        if self.config.outgroup_vcf:
            cmd_parts.append(f"-o {self.config.outgroup_vcf}")
        
        if self.config.resolve_iupac:
            cmd_parts.append("-r")
        
        # 禁用其他格式输出，只生成PHYLIP | Disable other formats, only generate PHYLIP
        cmd_parts.extend(["-f", "-n", "-b"])  # disable fasta, nexus, binary nexus
        
        cmd = " ".join(cmd_parts)
        
        success = self.cmd_runner.run(cmd, "VCF到PHYLIP格式转换 | VCF to PHYLIP format conversion")
        
        if success and phylip_output.exists():
            self.logger.info(f"✅ VCF转换成功 | VCF conversion successful: {phylip_output}")
            return str(phylip_output)
        else:
            raise RuntimeError("VCF到PHYLIP转换失败 | VCF to PHYLIP conversion failed")
    
    def _convert_fasta_to_phylip(self) -> str:
        """FASTA转PHYLIP | Convert FASTA to PHYLIP"""
        self.logger.info("📄 开始FASTA到PHYLIP格式转换 | Starting FASTA to PHYLIP conversion")
        
        phylip_output = self.config.output_path / f"{self.config.output_name}_converted.phy"
        
        # 使用seqkit转换FASTA到PHYLIP | Use seqkit to convert FASTA to PHYLIP
        cmd = (
            f"{self.config.seqkit_path} fx2tab {self.config.sequence_file} | "
            f"{self.config.seqkit_path} tab2fx | "
            f"{self.config.seqkit_path} seq -w 0 | "
            f"python -c \""
            f"import sys; "
            f"lines=[l.strip() for l in sys.stdin if l.strip()]; "
            f"seqs=[]; names=[]; i=0; "
            f"while i<len(lines): "
            f"    if lines[i].startswith('>'): "
            f"        names.append(lines[i][1:]); "
            f"        seqs.append(lines[i+1] if i+1<len(lines) else ''); "
            f"        i+=2; "
            f"    else: i+=1; "
            f"print(f'{{len(seqs)}} {{len(seqs[0]) if seqs else 0}}'); "
            f"[print(f'{{name:<10}} {{seq}}') for name,seq in zip(names,seqs)]"
            f"\" > {phylip_output}"
        )
        
        success = self.cmd_runner.run(cmd, "FASTA到PHYLIP格式转换 | FASTA to PHYLIP format conversion")
        
        if success and phylip_output.exists():
            self.logger.info(f"✅ FASTA转换成功 | FASTA conversion successful: {phylip_output}")
            return str(phylip_output)
        else:
            raise RuntimeError("FASTA到PHYLIP转换失败 | FASTA to PHYLIP conversion failed")
    
    def _create_vcf2phylip_script(self, script_path: Path):
        """创建vcf2phylip脚本 | Create vcf2phylip script"""
        self.logger.info("📝 创建vcf2phylip转换脚本 | Creating vcf2phylip conversion script")
        
        vcf2phylip_code = '''#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
VCF to PHYLIP converter for RAxML analysis
Modified version for integration with RAxML pipeline
"""

import argparse
import gzip
import random
import sys
from pathlib import Path

# Dictionary of IUPAC ambiguities for nucleotides
AMBIG = {
    "A"    :"A", "C"    :"C", "G"    :"G", "N"    :"N", "T"    :"T",
    "*A"   :"a", "*C"   :"c", "*G"   :"g", "*N"   :"n", "*T"   :"t",
    "AC"   :"M", "AG"   :"R", "AN"   :"a", "AT"   :"W", "CG"   :"S",
    "CN"   :"c", "CT"   :"Y", "GN"   :"g", "GT"   :"K", "NT"   :"t",
    "*AC"  :"m", "*AG"  :"r", "*AN"  :"a", "*AT"  :"w", "*CG"  :"s",
    "*CN"  :"c", "*CT"  :"y", "*GN"  :"g", "*GT"  :"k", "*NT"  :"t",
    "ACG"  :"V", "ACN"  :"m", "ACT"  :"H", "AGN"  :"r", "AGT"  :"D",
    "ANT"  :"w", "CGN"  :"s", "CGT"  :"B", "CNT"  :"y", "GNT"  :"k",
    "*ACG" :"v", "*ACN" :"m", "*ACT" :"h", "*AGN" :"r", "*AGT" :"d",
    "*ANT" :"w", "*CGN" :"s", "*CGT" :"b", "*CNT" :"y", "*GNT" :"k",
    "ACGN" :"v", "ACGT" :"N", "ACNT" :"h", "AGNT" :"d", "CGNT" :"b",
    "*ACGN":"v", "*ACGT":"N", "*ACNT":"h", "*AGNT":"d", "*CGNT":"b",
    "*"    :"-", "*ACGNT":"N",
}

def extract_sample_names(vcf_file):
    """Extract sample names from VCF file"""
    if vcf_file.lower().endswith(".gz"):
        opener = gzip.open
    else:
        opener = open
    sample_names = []
    with opener(vcf_file, "rt") as vcf:
        for line in vcf:
            line = line.strip("\\n")
            if line.startswith("#CHROM"):
                record = line.split("\\t")
                sample_names = [record[i].replace("./", "") for i in range(9, len(record))]
                break
    return sample_names

def is_anomalous(record, num_samples):
    """Check if record has correct number of samples"""
    return bool(len(record) != num_samples + 9)

def is_snp(record):
    """Check if record is a SNP"""
    alt = record[4].replace("<NON_REF>", record[3])
    return bool(len(record[3]) == 1 and len(alt) - alt.count(",") == alt.count(",") + 1)

def num_genotypes(record, num_samples):
    """Count non-missing genotypes"""
    missing = 0
    for i in range(9, num_samples + 9):
        if record[i].startswith("."):
            missing += 1
    return num_samples - missing

def get_matrix_column(record, num_samples, resolve_IUPAC):
    """Convert VCF record to matrix column"""
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
        if resolve_IUPAC is False:
            column += AMBIG[geno_nuc]
        else:
            column += AMBIG[nt_dict[random.choice(geno_num)]]
    return column

def main():
    parser = argparse.ArgumentParser(description="Convert VCF to PHYLIP")
    parser.add_argument("-i", "--input", required=True, help="Input VCF file")
    parser.add_argument("--output-folder", default="./", help="Output folder")
    parser.add_argument("--output-prefix", help="Output prefix")
    parser.add_argument("-m", "--min-samples-locus", type=int, default=4, help="Minimum samples per locus")
    parser.add_argument("-o", "--outgroup", default="", help="Outgroup name")
    parser.add_argument("-p", "--phylip-disable", action="store_true", help="Disable PHYLIP output")
    parser.add_argument("-f", "--fasta", action="store_true", help="Enable FASTA output")
    parser.add_argument("-n", "--nexus", action="store_true", help="Enable NEXUS output")
    parser.add_argument("-b", "--nexus-binary", action="store_true", help="Enable binary NEXUS output")
    parser.add_argument("-r", "--resolve-IUPAC", action="store_true", help="Resolve IUPAC ambiguities")
    parser.add_argument("-w", "--write-used-sites", action="store_true", help="Write used sites")
    
    args = parser.parse_args()
    
    outgroup = args.outgroup.split(",")[0].split(";")[0]
    
    # Get sample names
    if Path(args.input).exists():
        sample_names = extract_sample_names(args.input)
    else:
        print("\\nInput VCF file not found")
        sys.exit(1)
    
    num_samples = len(sample_names)
    if num_samples == 0:
        print("\\nSample names not found in VCF")
        sys.exit(1)
    
    args.min_samples_locus = min(num_samples, args.min_samples_locus)
    
    if not args.output_prefix:
        parts = Path(args.input).name.split(".")
        args.output_prefix = []
        for p in parts:
            if p.lower() == "vcf":
                break
            else:
                args.output_prefix.append(p)
        args.output_prefix = ".".join(args.output_prefix)
    args.output_prefix += ".min" + str(args.min_samples_locus)
    
    if not Path(args.output_folder).exists():
        Path(args.output_folder).mkdir(parents=True)
    
    outfile = str(Path(args.output_folder, args.output_prefix))
    
    if not args.phylip_disable:
        temporal = open(outfile+".tmp", "w")
    
    # Process VCF
    opener = gzip.open if args.input.lower().endswith(".gz") else open
    
    with opener(args.input, "rt") as vcf:
        snp_num = 0
        snp_accepted = 0
        snp_shallow = 0
        mnp_num = 0
        
        while True:
            vcf_chunk = vcf.readlines(50000)
            if not vcf_chunk:
                break
            
            for line in vcf_chunk:
                line = line.strip()
                if line and not line.startswith("#"):
                    record = line.split("\\t")
                    snp_num += 1
                    
                    if snp_num % 100000 == 0:
                        print(f"{snp_num} genotypes processed.")
                    
                    if is_anomalous(record, num_samples):
                        continue
                    
                    num_samples_locus = num_genotypes(record, num_samples)
                    if num_samples_locus < args.min_samples_locus:
                        snp_shallow += 1
                        continue
                    
                    if is_snp(record):
                        if not args.phylip_disable:
                            site_tmp = get_matrix_column(record, num_samples, args.resolve_IUPAC)
                            if site_tmp == "malformed":
                                continue
                            snp_accepted += 1
                            temporal.write(site_tmp+"\\n")
                    else:
                        mnp_num += 1
    
    print(f"Total genotypes processed: {snp_num}")
    print(f"SNPs that passed filters: {snp_accepted}")
    
    if not args.phylip_disable:
        temporal.close()
        
        # Write PHYLIP output
        output_phy = open(outfile+".phy", "w")
        output_phy.write(f"{len(sample_names)} {snp_accepted}\\n")
        
        # Write outgroup first if specified
        idx_outgroup = None
        if outgroup in sample_names:
            idx_outgroup = sample_names.index(outgroup)
            with open(outfile+".tmp") as tmp_seq:
                seqout = ""
                for line in tmp_seq:
                    seqout += line[idx_outgroup]
                output_phy.write(f"{sample_names[idx_outgroup]:<10} {seqout}\\n")
        
        # Write other samples
        for s in range(len(sample_names)):
            if s != idx_outgroup:
                with open(outfile+".tmp") as tmp_seq:
                    seqout = ""
                    for line in tmp_seq:
                        seqout += line[s]
                    output_phy.write(f"{sample_names[s]:<10} {seqout}\\n")
        
        output_phy.close()
        Path(outfile+".tmp").unlink()
        print(f"PHYLIP matrix saved to: {outfile}.phy")

if __name__ == "__main__":
    main()
'''
        
        with open(script_path, 'w') as f:
            f.write(vcf2phylip_code)
        
        # 设置执行权限 | Set execute permissions
        script_path.chmod(0o755)
        
        self.logger.info(f"📝 vcf2phylip脚本已创建 | vcf2phylip script created: {script_path}")
    
    def _convert_fasta_to_phylip(self) -> str:
        """FASTA转PHYLIP | Convert FASTA to PHYLIP"""
        self.logger.info("📄 开始FASTA到PHYLIP格式转换 | Starting FASTA to PHYLIP conversion")
        
        phylip_output = self.config.output_path / f"{self.config.output_name}_converted.phy"
        
        # 第一步：读取FASTA文件获取序列信息 | Step 1: Read FASTA file to get sequence info
        temp_info = self.config.output_path / "temp_fasta_info.txt"
        
        info_cmd = f"{self.config.seqkit_path} stat {self.config.sequence_file} -T > {temp_info}"
        
        if not self.cmd_runner.run(info_cmd, "获取FASTA文件信息 | Getting FASTA file information"):
            raise RuntimeError("无法获取FASTA文件信息 | Cannot get FASTA file information")
        
        # 第二步：转换为PHYLIP格式 | Step 2: Convert to PHYLIP format
        convert_cmd = (
            f"{self.config.seqkit_path} fx2tab {self.config.sequence_file} | "
            f"awk 'BEGIN{{OFS=\"\\t\"}} {{name=$1; $1=\"\"; seq=$0; gsub(/^\\t/, \"\", seq); "
            f"names[NR]=name; seqs[NR]=seq; if(NR==1) seqlen=length(seq)}} "
            f"END{{print NR, seqlen; for(i=1; i<=NR; i++) printf \"%-10s %s\\n\", names[i], seqs[i]}}' "
            f"> {phylip_output}"
        )
        
        success = self.cmd_runner.run(convert_cmd, "FASTA到PHYLIP格式转换 | FASTA to PHYLIP format conversion")
        
        # 清理临时文件 | Clean up temporary files
        if temp_info.exists():
            temp_info.unlink()
        
        if success and phylip_output.exists():
            self.logger.info(f"✅ FASTA转换成功 | FASTA conversion successful: {phylip_output}")
            return str(phylip_output)
        else:
            raise RuntimeError("FASTA到PHYLIP转换失败 | FASTA to PHYLIP conversion failed")
