"""
Augustus基因预测模块 | Augustus Gene Prediction Module
"""

import os

class GenePredictor:
    """基因预测器 | Gene Predictor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def predict_test_set(self):
        """预测测试集 | Predict test set"""
        self.logger.info("="*50)
        self.logger.info("步骤 5: 预测测试集 | Step 5: Predicting test set")
        
        augustus_bin = os.path.join(self.config.augustus_path, 'augustus')
        prediction_file = os.path.join(self.config.output_dir, 'prediction_result.gff')
        
        command = f"{augustus_bin} --species={self.config.species_name} {self.config.test_file} > {prediction_file}"
        
        self.cmd_runner.run_command(
            command, 
            "预测测试集 | Predict test set"
        )
        
        self.config.prediction_file = prediction_file
        self.logger.info(f"预测结果已保存: {prediction_file} | Prediction results saved: {prediction_file}")
        
        return prediction_file
    
    def convert_to_gff3(self):
        """转换为GFF3格式 | Convert to GFF3 format"""
        self.logger.info("="*50)
        self.logger.info("步骤 7: 转换为GFF3格式 | Step 7: Converting to GFF3 format")
        
        gff3_file = os.path.join(self.config.output_dir, 'prediction_result.gff3')
        
        # 尝试使用gffread转换 | Try using gffread for conversion
        command = f"gffread {self.config.prediction_file} -o {gff3_file}"
        
        try:
            self.cmd_runner.run_command(command, "转换为GFF3格式 | Convert to GFF3 format")
            self.logger.info(f"GFF3文件已生成: {gff3_file} | GFF3 file generated: {gff3_file}")
        except Exception:
            self.logger.warning("gffread转换失败，尝试简单格式转换... | gffread conversion failed, attempting simple format conversion...")
            self._simple_gff_to_gff3_conversion(gff3_file)
        
        return gff3_file
    
    def _simple_gff_to_gff3_conversion(self, output_file: str):
        """简单GFF到GFF3转换 | Simple GFF to GFF3 conversion"""
        with open(self.config.prediction_file, 'r') as infile, \
             open(output_file, 'w') as outfile:
            
            outfile.write("##gff-version 3\n")
            
            for line in infile:
                if line.startswith('#') or line.strip() == '':
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 9:
                    # 简单处理属性字段确保GFF3格式兼容 | Simple processing of attribute field
                    attributes = fields[8]
                    if 'transcript_id' in attributes and 'gene_id' in attributes:
                        outfile.write(line)
        
        self.logger.info(f"简单格式转换完成: {output_file} | Simple format conversion completed: {output_file}")
