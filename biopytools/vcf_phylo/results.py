"""
VCF系统发育分析结果处理模块 | VCF Phylogenetic Analysis Results Module
"""

import os
from pathlib import Path
from datetime import datetime

class ResultsManager:
    """结果管理器 | Results Manager"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self):
        """生成总结报告 | Generate summary report"""
        report_file = f"{self.config.output_prefix}_summary.txt"
        
        self.logger.info("生成分析总结报告 | Generating analysis summary report")
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("="*80 + "\n")
            f.write("VCF系统发育分析总结报告 | VCF Phylogenetic Analysis Summary Report\n")
            f.write("="*80 + "\n\n")
            
            # 基本信息 | Basic information
            f.write("基本信息 | Basic Information:\n")
            f.write(f"  分析时间 | Analysis time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"  输出前缀 | Output prefix: {self.config.output_prefix}\n")
            f.write(f"  工作目录 | Working directory: {self.config.working_dir}\n\n")
            
            # 输入文件 | Input files
            f.write("输入文件 | Input Files:\n")
            if self.config.vcf_file:
                f.write(f"  VCF文件 | VCF file: {self.config.vcf_file}\n")
            if self.config.distance_matrix:
                f.write(f"  距离矩阵 | Distance matrix: {self.config.distance_matrix}\n")
            f.write(f"  跳过VCF2Dis | Skip VCF2Dis: {'是 | Yes' if self.config.skip_vcf2dis else '否 | No'}\n\n")
            
            # 输出文件 | Output files
            f.write("输出文件 | Output Files:\n")
            
            if os.path.exists(self.config.distance_matrix):
                size = os.path.getsize(self.config.distance_matrix)
                f.write(f"  - 距离矩阵 | Distance matrix: {self.config.distance_matrix} ({size} bytes)\n")
            
            if os.path.exists(self.config.tree_output):
                size = os.path.getsize(self.config.tree_output)
                f.write(f"  - 系统发育树 | Phylogenetic tree: {self.config.tree_output} ({size} bytes)\n")
            
            log_file = f"{self.config.output_prefix}.log"
            if os.path.exists(log_file):
                size = os.path.getsize(log_file)
                f.write(f"  - 日志文件 | Log file: {log_file} ({size} bytes)\n")
            
            f.write(f"\n分析方法 | Analysis Methods:\n")
            f.write("  - 距离计算 | Distance calculation: VCF2Dis\n")
            f.write("  - 系统发育树构建 | Tree construction: Neighbor-Joining (NJ)\n")
            f.write("  - 树格式 | Tree format: Newick\n")
            f.write("  - 使用库 | Library used: scikit-bio\n\n")
            
            f.write("使用说明 | Usage Notes:\n")
            f.write("  1. 距离矩阵文件可用于其他系统发育分析软件\n")
            f.write("     Distance matrix can be used with other phylogenetic analysis software\n")
            f.write("  2. Newick格式的树文件可在FigTree、iTOL等软件中查看\n")
            f.write("     Newick tree file can be viewed in FigTree, iTOL and other software\n")
            f.write("  3. 建议使用R的ape包或Python的ete3包进行进一步分析\n")
            f.write("     Recommend using R's ape package or Python's ete3 package for further analysis\n\n")
        
        self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
        return report_file
    
    def log_output_files(self):
        """记录输出文件信息 | Log output file information"""
        self.logger.info("输出文件 | Output files:")
        
        files_to_check = [
            (self.config.distance_matrix, "距离矩阵 | Distance matrix"),
            (self.config.tree_output, "系统发育树 | Phylogenetic tree"),
            (f"{self.config.output_prefix}.log", "日志文件 | Log file"),
            (f"{self.config.output_prefix}_summary.txt", "总结报告 | Summary report")
        ]
        
        for file_path, description in files_to_check:
            if os.path.exists(file_path):
                size = os.path.getsize(file_path)
                self.logger.info(f"  - {description}: {file_path} ({size} bytes)")
    
    def validate_results(self) -> bool:
        """验证结果文件 | Validate result files"""
        self.logger.info("验证输出结果 | Validating output results")
        
        success = True
        
        # 检查距离矩阵 | Check distance matrix
        if not self.config.skip_vcf2dis:
            if not os.path.exists(self.config.distance_matrix):
                self.logger.error(f"距离矩阵文件未生成 | Distance matrix file not generated: {self.config.distance_matrix}")
                success = False
            else:
                self.logger.info("✓ 距离矩阵文件存在 | Distance matrix file exists")
        
        # 检查系统发育树 | Check phylogenetic tree
        if not os.path.exists(self.config.tree_output):
            self.logger.error(f"系统发育树文件未生成 | Phylogenetic tree file not generated: {self.config.tree_output}")
            success = False
        else:
            self.logger.info("✓ 系统发育树文件存在 | Phylogenetic tree file exists")
            
            # 验证树文件内容 | Validate tree file content
            try:
                with open(self.config.tree_output, 'r') as f:
                    tree_content = f.read().strip()
                    if not tree_content:
                        self.logger.error("系统发育树文件为空 | Phylogenetic tree file is empty")
                        success = False
                    elif not tree_content.endswith(';'):
                        self.logger.warning("系统发育树可能格式不正确（不以;结尾）| Tree may have incorrect format (doesn't end with ;)")
                    else:
                        self.logger.info("✓ 系统发育树文件格式正确 | Phylogenetic tree file format is correct")
            except Exception as e:
                self.logger.error(f"无法验证系统发育树文件 | Cannot validate phylogenetic tree file: {e}")
                success = False
        
        return success
