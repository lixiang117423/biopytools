# """
# VCF数据预处理模块 | VCF Data Preprocessing Module
# """

# import os
# import pandas as pd
# from pathlib import Path
# from typing import Dict
# from .utils import CommandRunner

# class VCFPreprocessor:
#     """VCF文件预处理器 | VCF File Preprocessor"""
    
#     def __init__(self, config, logger, cmd_runner: CommandRunner):
#         self.config = config
#         self.logger = logger
#         self.cmd_runner = cmd_runner
    
#     def quality_control(self):
#         """质量控制过滤 | Quality control filtering"""
#         self.logger.info("开始VCF质量控制 | Starting VCF quality control")
        
#         input_vcf = self.config.vcf_file
#         output_vcf = self.config.output_path / f"{self.config.base_name}_qc.vcf.gz"
        
#         # 构建过滤命令 | Build filtering command
#         filter_params = [
#             f"--maf {self.config.maf}",
#             f"--max-missing {1 - self.config.missing_rate}",
#             f"--hwe {self.config.hwe_pvalue}",
#             f"--minDP {self.config.min_dp}",
#             f"--maxDP {self.config.max_dp}",
#             "--remove-indels",
#             "--min-alleles 2",
#             "--max-alleles 2",
#             "--recode",
#             "--recode-INFO-all"
#         ]
        
#         cmd = f"{self.config.vcftools_path} --vcf {input_vcf} {' '.join(filter_params)} --out {output_vcf.stem}"
        
#         success = self.cmd_runner.run(cmd, "VCF质量控制过滤 | VCF quality control filtering")
        
#         if success:
#             recode_file = self.config.output_path / f"{self.config.base_name}_qc.recode.vcf"
#             if recode_file.exists():
#                 bgzip_cmd = f"bgzip -c {recode_file} > {output_vcf}"
#                 self.cmd_runner.run(bgzip_cmd, "压缩过滤后的VCF文件 | Compressing filtered VCF file")
#                 recode_file.unlink()
            
#             return str(output_vcf)
        
#         return None
    
#     def load_group_info(self) -> Dict[str, str]:
#         """加载分组信息 | Load group information"""
#         if not self.config.group_file:
#             return {}
        
#         self.logger.info(f"加载分组信息 | Loading group information from: {self.config.group_file}")
        
#         try:
#             df = pd.read_csv(self.config.group_file, sep='\t', header=None, names=['sample', 'group'])
#             group_dict = dict(zip(df['sample'], df['group']))
            
#             self.logger.info(f"加载了 {len(group_dict)} 个样本的分组信息 | Loaded group info for {len(group_dict)} samples")
#             self.logger.info(f"分组: {set(group_dict.values())} | Groups: {set(group_dict.values())}")
            
#             return group_dict
        
#         except Exception as e:
#             self.logger.error(f"读取分组文件失败 | Failed to read group file: {e}")
#             return {}

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
        # 检查是否跳过质控 | Check if skip QC
        if self.config.skip_qc:
            self.logger.info("跳过VCF质量控制步骤 | Skipping VCF quality control step")
            self.logger.info("使用原始输入文件进行后续分析 | Using original input file for subsequent analysis")
            return self.config.vcf_file
        
        self.logger.info("开始VCF质量控制 | Starting VCF quality control")
        
        input_vcf = self.config.vcf_file
        # 使用绝对路径 | Use absolute path
        output_vcf = self.config.output_path / f"{self.config.base_name}_qc.vcf.gz"
        output_vcf_abs = output_vcf.resolve()
        
        # 获取正确的VCFtools输入参数 | Get correct VCFtools input parameter
        vcf_input_param = self.config.get_vcftools_input_param()
        
        self.logger.info(f"检测到VCF文件格式: {'压缩' if self.config.is_compressed else '未压缩'} | "
                        f"Detected VCF format: {'compressed' if self.config.is_compressed else 'uncompressed'}")
        self.logger.info(f"使用VCFtools参数: {vcf_input_param} | Using VCFtools parameter: {vcf_input_param}")
        self.logger.info(f"输入文件: {input_vcf} | Input file: {input_vcf}")
        self.logger.info(f"输出文件: {output_vcf_abs} | Output file: {output_vcf_abs}")
        
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
        
        # 使用绝对路径构建命令 | Build command with absolute paths
        output_prefix = output_vcf_abs.with_suffix('').with_suffix('')  # 移除.vcf.gz扩展名
        cmd = f"{self.config.vcftools_path} {vcf_input_param} {input_vcf} {' '.join(filter_params)} --out {output_prefix}"
        
        success = self.cmd_runner.run(cmd, "VCF质量控制过滤 | VCF quality control filtering")
        
        if success:
            # 检查生成的.recode.vcf文件 | Check generated .recode.vcf file
            recode_file = Path(f"{output_prefix}.recode.vcf")
            if recode_file.exists():
                # 压缩过滤后的VCF文件 | Compress filtered VCF file
                bgzip_cmd = f"bgzip -c {recode_file} > {output_vcf_abs}"
                compress_success = self.cmd_runner.run(bgzip_cmd, "压缩过滤后的VCF文件 | Compressing filtered VCF file")
                
                if compress_success:
                    # 清理临时文件 | Clean up temporary files
                    try:
                        recode_file.unlink()
                        self.logger.info(f"清理临时文件: {recode_file} | Cleaned up temporary file: {recode_file}")
                    except Exception as e:
                        self.logger.warning(f"清理临时文件失败: {e} | Failed to clean up temporary file: {e}")
                    
                    self.logger.info(f"质控完成，输出文件: {output_vcf_abs} | QC completed, output file: {output_vcf_abs}")
                    return str(output_vcf_abs)
                else:
                    self.logger.error("压缩VCF文件失败 | Failed to compress VCF file")
                    return None
            else:
                self.logger.error(f"VCFtools没有生成预期的文件: {recode_file} | VCFtools did not generate expected file: {recode_file}")
                # 列出输出目录中的文件以帮助调试 | List files in output directory for debugging
                try:
                    output_files = list(self.config.output_path.iterdir())
                    self.logger.info(f"输出目录中的文件: {[f.name for f in output_files]} | Files in output directory: {[f.name for f in output_files]}")
                except:
                    pass
                return None
        
        return None
    
    def check_vcf_format(self, vcf_file: str) -> dict:
        """检查VCF文件格式和基本信息 | Check VCF file format and basic information"""
        self.logger.info(f"检查VCF文件格式 | Checking VCF file format: {vcf_file}")
        
        vcf_info = {
            'is_compressed': vcf_file.endswith('.gz'),
            'exists': os.path.exists(vcf_file),
            'size_mb': 0,
            'sample_count': 0,
            'variant_count': 0
        }
        
        if vcf_info['exists']:
            # 获取文件大小 | Get file size
            vcf_info['size_mb'] = round(os.path.getsize(vcf_file) / (1024 * 1024), 2)
            
            # 使用BCFtools获取基本统计信息 | Use BCFtools to get basic statistics
            try:
                stats_cmd = f"{self.config.bcftools_path} stats {vcf_file}"
                result = self.cmd_runner.run_with_output(stats_cmd, "获取VCF统计信息 | Getting VCF statistics")
                
                if result:
                    for line in result.split('\n'):
                        if line.startswith('SN\t0\tnumber of samples:'):
                            vcf_info['sample_count'] = int(line.split('\t')[3])
                        elif line.startswith('SN\t0\tnumber of records:'):
                            vcf_info['variant_count'] = int(line.split('\t')[3])
                            
            except Exception as e:
                self.logger.warning(f"无法获取VCF统计信息: {e} | Cannot get VCF statistics: {e}")
        
        # 记录VCF文件信息 | Log VCF file information
        self.logger.info(f"VCF文件信息 | VCF file information:")
        self.logger.info(f"  文件路径 | File path: {vcf_file}")
        self.logger.info(f"  文件存在 | File exists: {vcf_info['exists']}")
        self.logger.info(f"  文件大小 | File size: {vcf_info['size_mb']} MB")
        self.logger.info(f"  压缩格式 | Compressed: {vcf_info['is_compressed']}")
        self.logger.info(f"  样本数量 | Sample count: {vcf_info['sample_count']}")
        self.logger.info(f"  变异数量 | Variant count: {vcf_info['variant_count']}")
        
        return vcf_info
    
    def load_group_info(self) -> Dict[str, str]:
        """加载分组信息 | Load group information"""
        if not self.config.group_file:
            self.logger.info("未提供分组文件，将进行全样本分析 | No group file provided, will perform whole-sample analysis")
            return {}
        
        self.logger.info(f"加载分组信息 | Loading group information from: {self.config.group_file}")
        
        try:
            df = pd.read_csv(self.config.group_file, sep='\t', header=None, names=['sample', 'group'])
            group_dict = dict(zip(df['sample'], df['group']))
            
            self.logger.info(f"加载了 {len(group_dict)} 个样本的分组信息 | Loaded group info for {len(group_dict)} samples")
            self.logger.info(f"分组: {set(group_dict.values())} | Groups: {set(group_dict.values())}")
            
            # 验证分组信息 | Validate group information
            unique_groups = set(group_dict.values())
            if len(unique_groups) < 2:
                self.logger.warning(f"警告：只检测到 {len(unique_groups)} 个分组，某些分析可能无法进行 | "
                                  f"Warning: Only {len(unique_groups)} groups detected, some analyses may not be possible")
            
            return group_dict
        
        except Exception as e:
            self.logger.error(f"读取分组文件失败 | Failed to read group file: {e}")
            return {}