"""
VCF数据预处理模块 | VCF Data Preprocessing Module
"""

import os
import pandas as pd
from pathlib import Path
from typing import Dict
from .utils import CommandRunner

class VCFPreprocessor:
    """VCF文件预处理器 | VCF File Preprocessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def quality_control(self):
        """质量控制过滤 | Quality control filtering"""
        self.logger.info("开始VCF质量控制 | Starting VCF quality control")
        
        input_vcf = self.config.vcf_file
        output_vcf = self.config.output_path / f"{self.config.base_name}_qc.vcf.gz"
        
        # 构建过滤命令 | Build filtering command
        filter_params = [
            f"--maf {self.config.maf}",
            f"--max-missing {1 - self.config.missing_rate}",
            f"--hwe {self.config.hwe_pvalue}",
            f"--minDP {self.config.min_dp}",
            f"--maxDP {self.config.max_dp}",
            "--remove-indels",
            "--min-alleles 2",
            "--max-alleles 2",
            "--recode",
            "--recode-INFO-all"
        ]
        
        cmd = f"{self.config.vcftools_path} --vcf {input_vcf} {' '.join(filter_params)} --out {output_vcf.stem}"
        
        success = self.cmd_runner.run(cmd, "VCF质量控制过滤 | VCF quality control filtering")
        
        if success:
            recode_file = self.config.output_path / f"{self.config.base_name}_qc.recode.vcf"
            if recode_file.exists():
                bgzip_cmd = f"bgzip -c {recode_file} > {output_vcf}"
                self.cmd_runner.run(bgzip_cmd, "压缩过滤后的VCF文件 | Compressing filtered VCF file")
                recode_file.unlink()
            
            return str(output_vcf)
        
        return None
    
    def load_group_info(self) -> Dict[str, str]:
        """加载分组信息 | Load group information"""
        if not self.config.group_file:
            return {}
        
        self.logger.info(f"加载分组信息 | Loading group information from: {self.config.group_file}")
        
        try:
            df = pd.read_csv(self.config.group_file, sep='\t', header=None, names=['sample', 'group'])
            group_dict = dict(zip(df['sample'], df['group']))
            
            self.logger.info(f"加载了 {len(group_dict)} 个样本的分组信息 | Loaded group info for {len(group_dict)} samples")
            self.logger.info(f"分组: {set(group_dict.values())} | Groups: {set(group_dict.values())}")
            
            return group_dict
        
        except Exception as e:
            self.logger.error(f"读取分组文件失败 | Failed to read group file: {e}")
            return {}
