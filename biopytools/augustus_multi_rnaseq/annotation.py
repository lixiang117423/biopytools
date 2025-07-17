"""
Augustus预测核心模块 | Augustus Prediction Core Module
"""

import os
import shutil
from pathlib import Path
from .utils import CommandRunner

class AugustusPredictor:
    """Augustus预测器 | Augustus Predictor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def create_extrinsic_config(self) -> str:
        """创建外在信息配置文件 | Create extrinsic configuration file"""
        config_file = self.config.output_path / "extrinsic_multi_rnaseq.cfg"
        
        config_content = """# Extrinsic information configuration file for multiple RNA-seq samples

[SOURCES]
M RM E W b2h

[SOURCE-PARAMETERS]

#   feature type    source  bonus        malus   len_bonus   min_len     max_len     radius
E   exonpart        b2h     1            1       0           1           10000       0
E   exon            b2h     1            1       0           1           10000       0  
I   intron          b2h     10           1       0           1           500000      0
5   dss             b2h     3            1       0           1           1           0
3   ass             b2h     3            1       0           1           1           0

# 默认参数
M   exonpart        M       1            1       0           1           10000       0
M   exon            M       1            1       0           1           10000       0
M   intron          M       1            1       0           1           500000      0
M   CDSpart         M       1            1       0           1           10000       0
M   CDS             M       1            1       0           1           10000       0
"""
        
        with open(config_file, 'w') as f:
            f.write(config_content)
        
        self.logger.info(f"创建外在信息配置文件 | Created extrinsic config file: {config_file}")
        return str(config_file)
    
    def run_augustus_model_only(self) -> str:
        """运行仅使用Augustus模型的预测 | Run Augustus prediction with model only"""
        output_file = self.config.output_path / "augustus_model_only.gff3"
        
        if output_file.exists():
            self.logger.info(f"模型预测结果已存在 | Model-only prediction already exists: {output_file}")
            return str(output_file)
        
        self.logger.info("运行Augustus（仅模型）| Running Augustus (model only)...")
        
        cmd_parts = [
            f"augustus",
            f"--species={self.config.species_model}",
            f"--gff3={'on' if self.config.gff3_output else 'off'}",
            f"--alternatives-from-evidence=false",
            f"{self.config.genome_file}",
            f"> {output_file}"
        ]
        
        cmd = " ".join(cmd_parts)
        success, _ = self.cmd_runner.run(cmd, "Augustus模型预测 | Augustus model prediction")
        
        if success:
            self.logger.info(f"模型预测完成 | Model prediction completed: {output_file}")
            return str(output_file)
        
        return ""
    
    def run_augustus_with_rnaseq(self, hints_file: str) -> str:
        """运行结合转录组数据的Augustus预测 | Run Augustus prediction with RNA-seq data"""
        import os
        
        output_file = self.config.output_path / "augustus_with_multi_rnaseq.gff3"
        
        if output_file.exists():
            self.logger.info(f"转录组预测结果已存在 | RNA-seq prediction already exists: {output_file}")
            return str(output_file)
        
        # 检查hints文件是否为空 | Check if hints file is empty
        if not os.path.exists(hints_file) or os.path.getsize(hints_file) == 0:
            self.logger.warning("Hints文件为空或不存在，将运行仅模型预测 | Hints file is empty or missing, running model-only prediction")
            # 复制模型预测结果 | Copy model-only prediction result
            model_only_file = self.config.output_path / "augustus_model_only.gff3"
            if model_only_file.exists():
                shutil.copy2(model_only_file, output_file)
                self.logger.info(f"已复制模型预测结果 | Copied model-only prediction results: {output_file}")
                return str(output_file)
            else:
                self.logger.error("模型预测结果不存在，无法继续 | Model-only prediction results do not exist, cannot continue")
                return ""
        
        self.logger.info("运行Augustus（结合转录组）| Running Augustus (with RNA-seq)...")
        
        extrinsic_config = self.create_extrinsic_config()
        
        cmd_parts = [
            f"augustus",
            f"--species={self.config.species_model}",
            f"--hintsfile={hints_file}",
            f"--extrinsicCfgFile={extrinsic_config}",
            f"--alternatives-from-evidence={'true' if self.config.alternatives_from_evidence else 'false'}",
            f"--allow_hinted_splicesites={self.config.allow_hinted_splicesites}",
            f"--gff3={'on' if self.config.gff3_output else 'off'}",
            f"{self.config.genome_file}",
            f"> {output_file}"
        ]
        
        cmd = " ".join(cmd_parts)
        success, _ = self.cmd_runner.run(cmd, "Augustus转录组预测 | Augustus RNA-seq prediction")
        
        if success:
            self.logger.info(f"转录组预测完成 | RNA-seq prediction completed: {output_file}")
            return str(output_file)
        
        return ""
