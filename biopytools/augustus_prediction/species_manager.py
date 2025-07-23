"""
Augustus物种模型管理模块 | Augustus Species Model Management Module
"""

import os

class SpeciesManager:
    """物种模型管理器 | Species Model Manager"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def create_species_model(self):
        """创建新的物种模型 | Create new species model"""
        self.logger.info("="*50)
        self.logger.info("步骤 1: 创建物种模型 | Step 1: Creating species model")
        
        new_species_script = os.path.join(self.config.augustus_path, 'new_species.pl')
        command = f"perl {new_species_script} --species={self.config.species_name}"
        
        try:
            self.cmd_runner.run_command(
                command, 
                f"创建物种模型: {self.config.species_name} | Create species model: {self.config.species_name}"
            )
            self.logger.info(f"物种模型创建成功: {self.config.species_name} | Species model created successfully: {self.config.species_name}")
            return True
            
        except Exception as e:
            self.logger.warning(f"物种模型可能已存在，继续执行... | Species model may already exist, continuing... ({e})")
            return True
