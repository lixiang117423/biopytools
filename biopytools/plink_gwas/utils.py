# """
# PLINK GWAS分析工具函数模块 | PLINK GWAS Analysis Utility Functions Module
# """

# import logging
# import subprocess
# import sys
# import shutil
# from pathlib import Path

# class PlinkGWASLogger:
#     """PLINK GWAS分析日志管理器 | PLINK GWAS Analysis Logger Manager"""
    
#     def __init__(self, output_dir: Path, log_name: str = "plink_analysis.log"):
#         self.output_dir = output_dir
#         self.log_file = output_dir / log_name
#         self.setup_logging()
    
#     def setup_logging(self):
#         """设置日志 | Setup logging"""
#         if self.log_file.exists():
#             self.log_file.unlink()
        
#         logging.basicConfig(
#             level=logging.INFO,
#             format='%(asctime)s - %(levelname)s - %(message)s',
#             handlers=[
#                 logging.FileHandler(self.log_file),
#                 logging.StreamHandler(sys.stdout)
#             ]
#         )
#         self.logger = logging.getLogger(__name__)
    
#     def get_logger(self):
#         """获取日志器 | Get logger"""
#         return self.logger

# class CommandRunner:
#     """命令执行器 | Command Runner"""
    
#     def __init__(self, logger, working_dir: Path):
#         self.logger = logger
#         self.working_dir = working_dir
    
#     def run(self, cmd, description: str = "", check: bool = True):
#         """执行命令 | Execute command"""
#         if description:
#             self.logger.info(f"执行步骤 | Executing step: {description}")
        
#         cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
#         self.logger.info(f"命令 | Command: {cmd_str}")
        
#         try:
#             result = subprocess.run(
#                 cmd, 
#                 shell=isinstance(cmd, str),
#                 capture_output=True, 
#                 text=True, 
#                 check=check,
#                 cwd=self.working_dir
#             )
            
#             if result.stdout:
#                 self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
#             if result.stderr:
#                 self.logger.warning(f"标准错误 | Stderr: {result.stderr}")
            
#             if description:
#                 self.logger.info(f"完成 | Completed: {description}")
            
#             return result
            
#         except subprocess.CalledProcessError as e:
#             self.logger.error(f"命令执行失败 | Command failed: {cmd_str}")
#             self.logger.error(f"错误信息 | Error: {e}")
#             if check:
#                 raise
#             return e

# class FileManager:
#     """文件管理器 | File Manager"""
    
#     def __init__(self, logger, working_dir: Path):
#         self.logger = logger
#         self.working_dir = working_dir
    
#     def copy_input_files(self, vcf_file: str, phenotype_file: str):
#         """复制输入文件到工作目录 | Copy input files to working directory"""
#         self.logger.info("复制输入文件到工作目录 | Copying input files to working directory...")
        
#         # 复制VCF文件 | Copy VCF file
#         vcf_dest = self.working_dir / "input.vcf.gz"
#         shutil.copy2(vcf_file, vcf_dest)
#         self.logger.info(f"VCF文件已复制 | VCF file copied: {vcf_dest}")
        
#         # 复制表型文件 | Copy phenotype file
#         pheno_dest = self.working_dir / "phenotype.txt"
#         shutil.copy2(phenotype_file, pheno_dest)
#         self.logger.info(f"表型文件已复制 | Phenotype file copied: {pheno_dest}")
        
#         return vcf_dest, pheno_dest
    
#     def check_file_exists(self, file_path: Path, description: str = "") -> bool:
#         """检查文件是否存在 | Check if file exists"""
#         if file_path.exists():
#             if description:
#                 self.logger.info(f"✓ {description}存在 | exists: {file_path}")
#             return True
#         else:
#             if description:
#                 self.logger.error(f"✗ {description}不存在 | does not exist: {file_path}")
#             return False
"""
PLINK GWAS分析工具函数模块 | PLINK GWAS Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
import shutil
from pathlib import Path

class PlinkGWASLogger:
    """PLINK GWAS分析日志管理器 | PLINK GWAS Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "plink_analysis.log"):
        self.output_dir = Path(output_dir)
        self.log_file = self.output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        # 确保输出目录存在 | Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # 如果日志文件已存在，删除它 | If log file exists, remove it
        if self.log_file.exists():
            self.log_file.unlink()
        
        # 清除可能存在的handler | Clear any existing handlers
        for handler in logging.root.handlers[:]:
            logging.root.removeHandler(handler)
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(str(self.log_file), encoding='utf-8'),
                logging.StreamHandler(sys.stdout)
            ],
            force=True  # 强制重新配置 | Force reconfiguration
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = Path(working_dir).resolve()  # 确保是绝对路径 | Ensure absolute path
        self.plink_param_descriptions = self._init_plink_param_descriptions()
    
    def _init_plink_param_descriptions(self):
        """初始化PLINK参数说明 | Initialize PLINK parameter descriptions"""
        return {
            # 基本参数 | Basic parameters
            "plink": "PLINK v1.9 - 全基因组关联分析工具 | PLINK v1.9 - Genome-wide association analysis tool",
            "--bfile": "指定二进制PLINK文件集合的前缀名 | Specify binary PLINK fileset prefix",
            "--vcf": "指定VCF格式输入文件 | Specify VCF format input file",
            "--out": "指定输出文件前缀 | Specify output file prefix",
            
            # 文件格式转换 | File format conversion
            "--make-bed": "生成二进制PLINK文件格式(.bed/.bim/.fam) | Generate binary PLINK format files (.bed/.bim/.fam)",
            "--recode": "重新编码输出文件格式 | Recode output file format",
            "--keep-allele-order": "保持等位基因顺序不变 | Keep allele order unchanged",
            "--set-missing-var-ids": "为缺失变异ID设置格式(@:#表示chr:pos) | Set format for missing variant IDs (@:# means chr:pos)",
            
            # 表型处理 | Phenotype handling
            "--pheno": "指定表型文件 | Specify phenotype file",
            "--mpheno": "指定多表型文件中的列号 | Specify column number in multi-phenotype file",
            
            # 质量控制 | Quality control
            "--mind": "移除个体缺失率超过阈值的样本 | Remove individuals with missing rate above threshold",
            "--geno": "移除SNP缺失率超过阈值的位点 | Remove SNPs with missing rate above threshold", 
            "--maf": "移除次要等位基因频率低于阈值的SNP | Remove SNPs with minor allele frequency below threshold",
            "--hwe": "移除不符合Hardy-Weinberg平衡的SNP | Remove SNPs failing Hardy-Weinberg equilibrium test",
            "--freq": "计算等位基因频率 | Calculate allele frequencies",
            "--missing": "计算缺失率统计 | Calculate missing rate statistics",
            "--hardy": "进行Hardy-Weinberg平衡检验 | Perform Hardy-Weinberg equilibrium test",
            
            # LD分析 | LD analysis
            "--indep-pairwise": "基于成对LD进行SNP剪枝 | Perform SNP pruning based on pairwise LD",
            "--extract": "仅保留指定列表中的SNP | Keep only SNPs in specified list",
            "--exclude": "排除指定列表中的SNP | Exclude SNPs in specified list",
            
            # 主成分分析 | Principal component analysis
            "--pca": "进行主成分分析，指定输出的主成分数量 | Perform PCA, specify number of components to output",
            
            # 关联分析 | Association analysis
            "--assoc": "进行基本关联分析（卡方检验） | Perform basic association analysis (chi-square test)",
            "--logistic": "进行逻辑回归关联分析（质量性状） | Perform logistic regression association analysis (qualitative traits)",
            "--linear": "进行线性回归关联分析（数量性状） | Perform linear regression association analysis (quantitative traits)",
            "--model": "测试多种遗传模型（加性、显性、隐性） | Test multiple genetic models (additive, dominant, recessive)",
            "--covar": "指定协变量文件 | Specify covariate file",
            "--covar-name": "指定要使用的协变量名称 | Specify covariate names to use",
            "--ci": "计算置信区间（指定置信水平） | Calculate confidence intervals (specify confidence level)",
            "--dominant": "使用显性遗传模型 | Use dominant genetic model",
            "--recessive": "使用隐性遗传模型 | Use recessive genetic model",
            
            # 输出控制 | Output control
            "--allow-extra-chr": "允许非标准染色体编号 | Allow non-standard chromosome codes",
            "--allow-no-sex": "允许性别信息缺失 | Allow missing sex information",
            "--threads": "指定使用的CPU线程数 | Specify number of CPU threads to use",
            
            # LD分析参数 | LD analysis parameters
            "window-size": "LD剪枝窗口大小(kb) | LD pruning window size (kb)",
            "step-size": "LD剪枝步长(SNP数) | LD pruning step size (number of SNPs)", 
            "r2-threshold": "LD剪枝r²阈值 | LD pruning r² threshold"
        }
    
    def _explain_plink_command(self, cmd):
        """解释PLINK命令参数 | Explain PLINK command parameters"""
        if not isinstance(cmd, list) or len(cmd) == 0 or cmd[0] != "plink":
            return
        
        self.logger.info("=" * 60)
        self.logger.info("📋 命令参数详解 | Command Parameter Explanation")
        self.logger.info("=" * 60)
        
        i = 0
        while i < len(cmd):
            param = cmd[i]
            
            if param in self.plink_param_descriptions:
                if i + 1 < len(cmd) and not cmd[i + 1].startswith("--"):
                    # 参数有值 | Parameter has value
                    value = cmd[i + 1]
                    self.logger.info(f"  {param} {value}")
                    self.logger.info(f"    💡 {self.plink_param_descriptions[param]}")
                    
                    # 特殊值解释 | Special value explanations
                    if param == "--indep-pairwise" and i + 3 < len(cmd):
                        window_size = cmd[i + 1]
                        step_size = cmd[i + 2] 
                        r2_threshold = cmd[i + 3]
                        self.logger.info(f"       📊 窗口大小 | Window size: {window_size}kb")
                        self.logger.info(f"       📊 步长 | Step size: {step_size} SNPs")
                        self.logger.info(f"       📊 r²阈值 | r² threshold: {r2_threshold}")
                        i += 3
                    elif param == "--pca":
                        self.logger.info(f"       📊 将计算前{value}个主成分 | Will compute top {value} principal components")
                        i += 1
                    elif param == "--mind":
                        percentage = float(value) * 100
                        self.logger.info(f"       📊 移除个体缺失率>{percentage}%的样本 | Remove samples with missing rate > {percentage}%")
                        i += 1
                    elif param == "--geno":
                        percentage = float(value) * 100
                        self.logger.info(f"       📊 移除SNP缺失率>{percentage}%的位点 | Remove SNPs with missing rate > {percentage}%")
                        i += 1
                    elif param == "--maf":
                        percentage = float(value) * 100
                        self.logger.info(f"       📊 移除次要等位基因频率<{percentage}%的SNP | Remove SNPs with MAF < {percentage}%")
                        i += 1
                    elif param == "--hwe":
                        self.logger.info(f"       📊 移除HWE检验P值<{value}的SNP | Remove SNPs with HWE test p-value < {value}")
                        i += 1
                    elif param == "--threads":
                        self.logger.info(f"       💻 使用{value}个CPU线程并行计算 | Use {value} CPU threads for parallel computing")
                        i += 1
                    elif param == "--ci":
                        confidence = float(value) * 100
                        self.logger.info(f"       📊 计算{confidence}%置信区间 | Calculate {confidence}% confidence intervals")
                        i += 1
                    else:
                        i += 1
                else:
                    # 参数无值（开关参数） | Parameter without value (switch parameter)
                    self.logger.info(f"  {param}")
                    self.logger.info(f"    💡 {self.plink_param_descriptions[param]}")
            else:
                # 未知参数 | Unknown parameter
                if param.startswith("--"):
                    self.logger.info(f"  {param}")
                    self.logger.info(f"    ❓ 未知参数 | Unknown parameter")
            
            i += 1
        
        self.logger.info("=" * 60)
    
    def run(self, cmd, description: str = "", check: bool = True):
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"🚀 执行步骤 | Executing step: {description}")
        
        cmd_str = ' '.join(cmd) if isinstance(cmd, list) else cmd
        self.logger.info(f"⚡ 命令 | Command: {cmd_str}")
        self.logger.info(f"📁 工作目录 | Working directory: {self.working_dir}")
        
        # 如果是PLINK命令，提供详细参数解释 | If it's a PLINK command, provide detailed parameter explanation
        if isinstance(cmd, list) and len(cmd) > 0 and cmd[0] == "plink":
            self._explain_plink_command(cmd)
        
        try:
            result = subprocess.run(
                cmd, 
                shell=isinstance(cmd, str),
                capture_output=True, 
                text=True, 
                check=check,
                cwd=str(self.working_dir)  # 确保使用字符串路径 | Ensure using string path
            )
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"标准错误 | Stderr: {result.stderr}")
            
            if description:
                self.logger.info(f"✅ 完成 | Completed: {description}")
            
            return result
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command failed: {cmd_str}")
            self.logger.error(f"💥 错误信息 | Error: {e}")
            if e.stdout:
                self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            if e.stderr:
                self.logger.error(f"标准错误 | Stderr: {e.stderr}")
            if check:
                raise
            return e

class FileManager:
    """文件管理器 | File Manager"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = Path(working_dir).resolve()  # 确保是绝对路径 | Ensure absolute path
    
    def copy_input_files(self, vcf_file: str, phenotype_file: str):
        """复制输入文件到工作目录 | Copy input files to working directory"""
        self.logger.info("📂 复制输入文件到工作目录 | Copying input files to working directory...")
        
        # 转换为Path对象并检查文件是否存在 | Convert to Path objects and check if files exist
        vcf_path = Path(vcf_file)
        pheno_path = Path(phenotype_file)
        
        if not vcf_path.exists():
            raise FileNotFoundError(f"VCF文件不存在 | VCF file does not exist: {vcf_file}")
        if not pheno_path.exists():
            raise FileNotFoundError(f"表型文件不存在 | Phenotype file does not exist: {phenotype_file}")
        
        # 获取文件大小信息 | Get file size information
        vcf_size = vcf_path.stat().st_size / (1024 * 1024)  # MB
        pheno_size = pheno_path.stat().st_size / 1024  # KB
        
        # 复制VCF文件 | Copy VCF file
        if vcf_path.suffix.lower() == '.gz':
            vcf_dest = self.working_dir / "input.vcf.gz"
        else:
            vcf_dest = self.working_dir / "input.vcf"
        
        try:
            self.logger.info(f"📁 正在复制VCF文件 | Copying VCF file ({vcf_size:.1f} MB)...")
            shutil.copy2(str(vcf_path), str(vcf_dest))
            self.logger.info(f"✅ VCF文件已复制 | VCF file copied: {vcf_dest}")
        except Exception as e:
            raise Exception(f"复制VCF文件失败 | Failed to copy VCF file: {e}")
        
        # 复制表型文件 | Copy phenotype file
        pheno_dest = self.working_dir / "phenotype.txt"
        try:
            self.logger.info(f"📁 正在复制表型文件 | Copying phenotype file ({pheno_size:.1f} KB)...")
            shutil.copy2(str(pheno_path), str(pheno_dest))
            self.logger.info(f"✅ 表型文件已复制 | Phenotype file copied: {pheno_dest}")
        except Exception as e:
            raise Exception(f"复制表型文件失败 | Failed to copy phenotype file: {e}")
        
        return str(vcf_dest), str(pheno_dest)
    
    def check_file_exists(self, file_path: Path, description: str = "") -> bool:
        """检查文件是否存在 | Check if file exists"""
        path_obj = Path(file_path)
        if path_obj.exists():
            if description:
                self.logger.info(f"✅ {description}存在 | exists: {path_obj}")
            return True
        else:
            if description:
                self.logger.error(f"❌ {description}不存在 | does not exist: {path_obj}")
            return False
    
    def get_file_size(self, file_path: Path) -> int:
        """获取文件大小 | Get file size"""
        try:
            return Path(file_path).stat().st_size
        except:
            return 0
    
    def create_directory(self, dir_path: Path):
        """创建目录 | Create directory"""
        try:
            Path(dir_path).mkdir(parents=True, exist_ok=True)
            self.logger.info(f"📁 目录已创建 | Directory created: {dir_path}")
        except Exception as e:
            self.logger.error(f"❌ 创建目录失败 | Failed to create directory {dir_path}: {e}")
            raise
    
    def remove_file(self, file_path: Path):
        """删除文件 | Remove file"""
        try:
            path_obj = Path(file_path)
            if path_obj.exists():
                path_obj.unlink()
                self.logger.info(f"🗑️ 文件已删除 | File removed: {file_path}")
        except Exception as e:
            self.logger.warning(f"⚠️ 删除文件失败 | Failed to remove file {file_path}: {e}")
    
    def list_files(self, directory: Path, pattern: str = "*") -> list:
        """列出目录中的文件 | List files in directory"""
        try:
            dir_path = Path(directory)
            if dir_path.exists() and dir_path.is_dir():
                return list(dir_path.glob(pattern))
            else:
                return []
        except Exception as e:
            self.logger.warning(f"⚠️ 列出文件失败 | Failed to list files in {directory}: {e}")
            return []