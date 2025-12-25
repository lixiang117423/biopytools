"""
🌳 系统发育树构建模块 | Phylogenetic Tree Building Module
"""

from pathlib import Path

class FastTreeBuilder:
    """FastTree系统发育树构建器 | FastTree Phylogenetic Tree Builder"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def build_tree(self, alignment_file: Path, output_file: Path, seq_type: str) -> bool:
        """构建系统发育树 | Build phylogenetic tree"""
        self.logger.info("=" * 60)
        self.logger.info("🌳 开始构建系统发育树 | Starting phylogenetic tree construction")
        self.logger.info("=" * 60)
        
        # 构建FastTree命令
        fasttree_cmd = f"{self.config.fasttree_path} "
        
        # 根据序列类型添加参数
        if seq_type == 'nucleotide':
            fasttree_cmd += "-nt "
            self.logger.info("🧬 使用核酸模型 | Using nucleotide model")
        else:
            self.logger.info("🧬 使用蛋白质模型 | Using protein model")
        
        # 添加用户自定义参数
        if self.config.fasttree_params:
            fasttree_cmd += f"{self.config.fasttree_params} "
        
        # 添加输入输出文件
        # fasttree_cmd += f"{alignment_file} > {output_file}"
        # 添加输入输出文件
        fasttree_cmd += f"{alignment_file.resolve()} > {output_file.resolve()}"
        
        success = self.cmd_runner.run(
            fasttree_cmd,
            description="FastTree系统发育树构建 | FastTree phylogenetic tree construction"
        )
        
        if success:
            self.logger.info(f"✅ 系统发育树构建完成 | Phylogenetic tree completed: {output_file}")
        else:
            self.logger.error(f"❌ 系统发育树构建失败 | Phylogenetic tree construction failed")
        
        return success
