"""
Augustus模型训练模块 | Augustus Model Training Module
"""

import os

class ModelTrainer:
    """模型训练器 | Model Trainer"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def train_model(self):
        """训练模型 | Train model"""
        self.logger.info("="*50)
        self.logger.info("步骤 4: 训练模型 | Step 4: Training model")
        
        # etraining 参数训练 | etraining parameter training
        self._run_etraining()
        
        # 模型参数优化 | Model parameter optimization
        self._run_optimization()
        
        self.logger.info("模型训练完成 | Model training completed")
    
    def _run_etraining(self):
        """运行etraining | Run etraining"""
        self.logger.info("运行etraining参数训练 | Running etraining parameter training")
        
        etraining_bin = os.path.join(self.config.augustus_path, 'etraining')
        command = f"{etraining_bin} --species={self.config.species_name} {self.config.train_file}"
        
        self.cmd_runner.run_command(
            command, 
            "etraining参数训练 | etraining parameter training"
        )
        
        self.logger.info("etraining完成 | etraining completed")
    
    def _run_optimization(self):
        """运行模型优化 | Run model optimization"""
        self.logger.info("运行模型参数优化 | Running model parameter optimization")
        
        optimize_script = os.path.join(self.config.augustus_path, 'optimize_augustus.pl')
        optimize_log = os.path.join(self.config.output_dir, f'optimize_{self.config.species_name}.log')
        
        command = f"perl {optimize_script} --species={self.config.species_name} {self.config.test_file} > {optimize_log} 2>&1"
        
        self.cmd_runner.run_command(
            command, 
            "模型参数优化 | Model parameter optimization",
            timeout=7200  # 2小时超时 | 2-hour timeout
        )
        
        self.logger.info(f"模型优化完成，日志: {optimize_log} | Model optimization completed, log: {optimize_log}")
