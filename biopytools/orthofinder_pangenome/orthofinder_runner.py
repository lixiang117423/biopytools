# ===== FILE: orthofinder_pangenome/orthofinder_runner.py =====
"""
OrthoFinder执行器模块 | OrthoFinder Runner Module
"""

import os
import shutil
from pathlib import Path
from .utils import CommandRunner

class OrthoFinderRunner:
    """OrthoFinder执行器 | OrthoFinder Runner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.orthofinder_output_dir = None
        self.actual_project_name = None
    
    # def check_existing_results(self) -> bool:
    #     """检查已有OrthoFinder结果 | Check existing OrthoFinder results"""
    #     input_path = Path(self.config.input_dir)
    #     orthofinder_path = input_path / "OrthoFinder"
    #     results_pattern = f"Results_{self.config.project_name}"
        
    #     # 添加调试信息
    #     self.logger.info(f"检查路径 | Checking path:")
    #     self.logger.info(f"   输入目录 | Input directory: {input_path}")
    #     self.logger.info(f"   OrthoFinder目录 | OrthoFinder directory: {orthofinder_path}")
    #     self.logger.info(f"   项目名称 | Project name: {self.config.project_name}")
    #     self.logger.info(f"   结果模式 | Results pattern: {results_pattern}")
    #     self.logger.info(f"   OrthoFinder目录存在 | OrthoFinder directory exists: {orthofinder_path.exists()}")
        
    #     if orthofinder_path.exists():
    #         possible_dirs = list(orthofinder_path.glob(f"{results_pattern}*"))
    #         self.logger.info(f"   找到的目录 | Found directories: {possible_dirs}")
            
    #         if possible_dirs:
    #             results_dir = max(possible_dirs, key=os.path.getctime)
    #             self.logger.info(f"   选择的结果目录 | Selected results directory: {results_dir}")
                
    #             # 检查关键文件是否存在
    #             orthogroups_file = results_dir / "Orthogroups" / "Orthogroups.txt"
    #             gene_count_file = results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
                
    #             self.logger.info(f"   Orthogroups文件 | Orthogroups file: {orthogroups_file}")
    #             self.logger.info(f"   Orthogroups文件存在 | Orthogroups file exists: {orthogroups_file.exists()}")
    #             self.logger.info(f"   GeneCount文件 | GeneCount file: {gene_count_file}")
    #             self.logger.info(f"   GeneCount文件存在 | GeneCount file exists: {gene_count_file.exists()}")
                
    #             if orthogroups_file.exists() and gene_count_file.exists():
    #                 self.logger.info(f"发现已有OrthoFinder结果 | Found existing OrthoFinder results: {results_dir}")
    #                 self.orthofinder_output_dir = results_dir
    #                 return True
        
    #     self.logger.info("未发现已有结果，需要运行OrthoFinder | No existing results found, need to run OrthoFinder")
    #     return False

    def check_existing_results(self) -> bool:
        """检查已有OrthoFinder结果 | Check existing OrthoFinder results"""
        # 检查两个可能的位置：输入目录和输出目录
        locations_to_check = [
            (Path(self.config.input_dir) / "OrthoFinder", "输入目录"),
            (Path(self.config.output_dir) / "Results_OrthoFinder", "输出目录"),
            (Path(self.config.output_dir) / "OrthoFinder_Results", "输出目录（备用名称）")
        ]
        
        results_pattern = f"Results_{self.config.project_name}"
        
        for check_path, location_name in locations_to_check:
            self.logger.info(f"检查{location_name} | Checking {location_name}: {check_path}")
            
            if check_path.exists():
                # 情况1：检查是否为直接的Results目录（在输出目录中）
                if check_path.name.startswith("Results_") or check_path.name.startswith("OrthoFinder_Results"):
                    results_dir = check_path
                    self.logger.info(f"   找到直接结果目录 | Found direct results directory: {results_dir}")
                else:
                    # 情况2：在目录下查找Results_项目名称的子目录（在输入目录中）
                    possible_dirs = list(check_path.glob(f"{results_pattern}*"))
                    if not possible_dirs:
                        self.logger.info(f"   未找到匹配的结果目录 | No matching results directory found")
                        continue
                    
                    results_dir = max(possible_dirs, key=os.path.getctime)
                    self.logger.info(f"   找到结果目录 | Found results directory: {results_dir}")
                
                # 检查关键文件是否存在
                orthogroups_file = results_dir / "Orthogroups" / "Orthogroups.txt"
                gene_count_file = results_dir / "Orthogroups" / "Orthogroups.GeneCount.tsv"
                
                self.logger.info(f"   检查关键文件存在性:")
                self.logger.info(f"     Orthogroups.txt: {orthogroups_file.exists()}")
                self.logger.info(f"     GeneCount.tsv: {gene_count_file.exists()}")
                
                if orthogroups_file.exists() and gene_count_file.exists():
                    self.logger.info(f"发现有效的OrthoFinder结果 | Found valid OrthoFinder results: {results_dir}")
                    self.orthofinder_output_dir = results_dir
                    return True
        
        self.logger.info("未发现有效结果，需要运行OrthoFinder | No valid results found, need to run OrthoFinder")
        return False
    
    def run_orthofinder(self) -> bool:
        """运行OrthoFinder分析 | Run OrthoFinder analysis"""
        
        # 检查是否跳过OrthoFinder步骤
        if self.config.skip_orthofinder:
            self.logger.info("跳过OrthoFinder步骤 | Skipping OrthoFinder step")
            return self.check_existing_results()
        
        # 检查断点续跑
        if self.config.resume_from_existing and not self.config.force_overwrite:
            if self.check_existing_results():
                self.logger.info("使用已有结果续跑 | Resuming from existing results")
                return True
        
        # 如果强制覆盖，清理已有结果
        if self.config.force_overwrite:
            self._cleanup_existing_results()
        
        self.logger.info(f"开始OrthoFinder分析 | Starting OrthoFinder analysis")
        self.logger.info(f"参数设置 | Parameter settings:")
        self.logger.info(f"   线程数 | Threads: {self.config.threads}")
        self.logger.info(f"   搜索程序 | Search program: {self.config.search_program}")
        self.logger.info(f"   MCL inflation: {self.config.mcl_inflation}")
        
        # 构建OrthoFinder命令 | Build OrthoFinder command
        # 让OrthoFinder在默认位置运行，完成后移动结果
        cmd_parts = [
            self.config.orthofinder_path,
            f"-f {self.config.input_dir}",
            f"-t {self.config.threads}",
            f"-S {self.config.search_program}",
            f"-I {self.config.mcl_inflation}",
            f"-n {self.config.project_name}"
        ]
        
        # 记录实际使用的项目名
        self.actual_project_name = self.config.project_name
        
        # 添加序列类型参数 | Add sequence type parameter
        if self.config.sequence_type == 'dna':
            cmd_parts.append("-d")
        
        # 分析模式控制 | Analysis mode control
        if self.config.basic_analysis_only:
            cmd_parts.append("-og")  # 只做到同源基因群推断 | Stop after orthogroups
        else:
            # 完整分析模式 | Full analysis mode
            if self.config.generate_trees:
                # 生成系统发育树 | Generate phylogenetic trees
                cmd_parts.extend([
                    f"-A {self.config.msa_program}",
                    f"-T {self.config.tree_program}"
                ])
            # 如果不生成树但要完整分析，不添加任何停止参数，让OrthoFinder运行完整流程
        
        cmd = " ".join(cmd_parts)
        
        # 执行OrthoFinder | Execute OrthoFinder
        success = self.cmd_runner.run(cmd, "OrthoFinder同源基因群分析 | OrthoFinder orthogroup analysis")
        
        if success:
            self._locate_output_directory()
            self.logger.info(f"OrthoFinder分析完成 | OrthoFinder analysis completed")
            self.logger.info(f"结果目录 | Results directory: {self.orthofinder_output_dir}")
        
        return success
    
    def _locate_output_directory(self):
        """定位OrthoFinder输出目录并移动到指定位置 | Locate OrthoFinder output directory and move to specified location"""
        # OrthoFinder在输入目录下的OrthoFinder子目录中创建Results_项目名目录
        results_pattern = f"Results_{self.actual_project_name}"
        
        input_path = Path(self.config.input_dir)
        orthofinder_path = input_path / "OrthoFinder"
        possible_dirs = list(orthofinder_path.glob(f"{results_pattern}*"))
        
        if possible_dirs:
            source_dir = max(possible_dirs, key=os.path.getctime)
            self.logger.info(f"找到OrthoFinder结果目录 | Found OrthoFinder results directory: {source_dir}")
            
            # 移动到目标输出目录
            target_dir = Path(self.config.output_dir) / "Results_OrthoFinder"
            if target_dir.exists():
                shutil.rmtree(target_dir)
            
            shutil.move(str(source_dir), str(target_dir))
            self.orthofinder_output_dir = target_dir
            self.logger.info(f"结果已移动到 | Results moved to: {target_dir}")
            
            # 同时移动输入目录下的OrthoFinder工作文件夹（如果存在）
            if orthofinder_path.exists():
                target_work_dir = Path(self.config.output_dir) / "OrthoFinder_WorkingDirectory" 
                if target_work_dir.exists():
                    shutil.rmtree(target_work_dir)
                shutil.move(str(orthofinder_path), str(target_work_dir))
                self.logger.info(f"工作目录已移动到 | Working directory moved to: {target_work_dir}")
        else:
            self.logger.error("未找到OrthoFinder结果目录 | OrthoFinder results directory not found")
            self.orthofinder_output_dir = None
    
    # def get_orthogroups_file(self) -> Path:
    #     """获取同源基因群文件路径 | Get orthogroups file path"""
    #     if not self.orthofinder_output_dir:
    #         raise RuntimeError("OrthoFinder结果目录未找到 | OrthoFinder results directory not found")
        
    #     orthogroups_dir = Path(self.orthofinder_output_dir) / "Orthogroups"
        
    #     # 查找Orthogroups.txt文件 | Look for Orthogroups.txt file
    #     orthogroups_file = orthogroups_dir / "Orthogroups.txt"
    #     if orthogroups_file.exists():
    #         return orthogroups_file
        
    #     # 查找其他可能的同源基因群文件 | Look for other possible orthogroups files
    #     possible_files = list(orthogroups_dir.glob("Orthogroups*.txt"))
    #     if possible_files:
    #         return possible_files[0]
        
    #     raise FileNotFoundError(f"未找到同源基因群文件 | Orthogroups file not found in {orthogroups_dir}")
    def get_orthogroups_file(self) -> Path:
        """获取同源基因群文件路径 | Get orthogroups file path"""
        if not self.orthofinder_output_dir:
            raise RuntimeError("OrthoFinder结果目录未找到 | OrthoFinder results directory not found")
        
        orthogroups_dir = Path(self.orthofinder_output_dir) / "Orthogroups"
        
        # 查找Orthogroups.tsv文件 | Look for Orthogroups.tsv file
        orthogroups_file = orthogroups_dir / "Orthogroups.tsv"
        if orthogroups_file.exists():
            return orthogroups_file
        
        raise FileNotFoundError(f"未找到同源基因群详细文件 | Orthogroups detailed file not found: {orthogroups_file}")
    
    def get_gene_count_file(self) -> Path:
        """获取基因计数文件路径 | Get gene count file path"""
        if not self.orthofinder_output_dir:
            raise RuntimeError("OrthoFinder结果目录未找到 | OrthoFinder results directory not found")
        
        orthogroups_dir = Path(self.orthofinder_output_dir) / "Orthogroups"
        gene_count_file = orthogroups_dir / "Orthogroups.GeneCount.tsv"  # 修正为 .tsv
        
        if gene_count_file.exists():
            return gene_count_file
        
        raise FileNotFoundError(f"未找到基因计数文件 | Gene count file not found: {gene_count_file}")
    
    def _cleanup_existing_results(self):
        """清理已有结果 | Cleanup existing results"""
        input_path = Path(self.config.input_dir)
        orthofinder_path = input_path / "OrthoFinder"
        
        if orthofinder_path.exists():
            shutil.rmtree(orthofinder_path)
            self.logger.info("已清理旧的OrthoFinder结果 | Cleaned up old OrthoFinder results")

# ===== END FILE =====