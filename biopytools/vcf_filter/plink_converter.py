"""
PLINK格式转换模块|PLINK Format Conversion Module
"""

import subprocess
import os
from typing import List

from ..common.conda_runner import build_conda_command

from ..common.conda_runner import build_conda_command

class PlinkConverter:
    """PLINK格式转换器|PLINK Format Converter"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def convert_with_plink(self, filtered_vcf: str) -> str:
        """使用PLINK进行格式转换|Convert format using PLINK"""
        
        if self.config.verbose:
            self.logger.info("使用PLINK进行格式转换|Converting format with PLINK")
        
        # 构建输出路径|Build output path
        output_prefix = str(self.config.output_path.with_suffix(''))
        
        # 构建PLINK命令|Build PLINK command
        args = [
            "--vcf", filtered_vcf,
            "--recode", "vcf-iid",
            "--out", output_prefix
        ]
        
        if self.config.allow_extra_chr:
            args.append("--allow-extra-chr")
        
        # 添加染色体筛选|Add chromosome filtering
        if self.config.chr_name is not None:
            if isinstance(self.config.chr_name, list):
                chr_list = ",".join(self.config.chr_name)
                args.extend(["--chr", chr_list])
            else:
                args.extend(["--chr", str(self.config.chr_name)])
        
        # 添加位置筛选|Add position filtering
        if self.config.start is not None and self.config.end is not None:
            args.extend(["--from-bp", str(self.config.start),
                       "--to-bp", str(self.config.end)])
        elif self.config.start is not None:
            args.extend(["--from-bp", str(self.config.start)])
        elif self.config.end is not None:
            args.extend(["--to-bp", str(self.config.end)])
        
        # 添加质量控制参数|Add quality control parameters
        if self.config.min_maf is not None:
            args.extend(["--maf", str(self.config.min_maf)])
        
        if self.config.max_missing is not None:
            args.extend(["--geno", str(self.config.max_missing)])

        # 构建完整命令(conda环境自动包装)|Build full command (auto conda wrap)
        cmd = build_conda_command(self.config.plink_path, args)
        
        try:
            # 执行PLINK命令|Execute PLINK command
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            
            result = subprocess.run(cmd, shell=False, capture_output=True, text=True, check=True)
            
            if self.config.verbose:
                self.logger.info("PLINK执行成功|PLINK execution successful")
            
            # PLINK输出文件会自动添加.vcf后缀|PLINK output file will automatically add .vcf suffix
            final_output = output_prefix + ".vcf"
            
            if os.path.exists(final_output):
                return final_output
            else:
                raise FileNotFoundError(f"PLINK输出文件未找到|PLINK output file not found: {final_output}")
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"PLINK执行失败|PLINK execution failed: {e}")
            if e.stderr:
                self.logger.error(f"错误输出|Error output: {e.stderr}")
            raise
        except FileNotFoundError:
            raise FileNotFoundError(f"找不到PLINK可执行文件|PLINK executable not found: {self.config.plink_path}")
