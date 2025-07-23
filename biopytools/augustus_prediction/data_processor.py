"""
Augustus数据处理模块 | Augustus Data Processing Module
"""

import os
import re

class DataProcessor:
    """数据处理器 | Data Processor"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def prepare_training_data(self):
        """准备训练数据 | Prepare training data"""
        self.logger.info("="*50)
        self.logger.info("步骤 2: 准备训练数据 | Step 2: Preparing training data")
        
        gff2gb_script = os.path.join(self.config.augustus_path, 'gff2gbSmallDNA.pl')
        output_file = os.path.join(self.config.output_dir, 'training_set.gb')
        
        command = (f"perl {gff2gb_script} "
                  f"{self.config.gff_file} "
                  f"{self.config.genome_file} "
                  f"{self.config.flank_length} "
                  f"{output_file}")
        
        self.cmd_runner.run_command(
            command, 
            f"生成Augustus训练数据 | Generate Augustus training data"
        )
        
        self.config.training_file = output_file
        self.logger.info(f"训练数据已生成: {output_file} | Training data generated: {output_file}")
        
        return output_file
    
    def split_dataset(self):
        """拆分数据集 | Split dataset"""
        self.logger.info("="*50)
        self.logger.info("步骤 3: 拆分训练和测试集 | Step 3: Splitting training and test sets")
        
        # 统计基因总数 | Count total genes
        with open(self.config.training_file, 'r') as f:
            content = f.read()
            total_genes = content.count('LOCUS')
            
        self.logger.info(f"检测到基因总数: {total_genes} | Detected total genes: {total_genes}")
        
        if total_genes < 100:
            raise ValueError("基因总数少于100个，不足以进行拆分和评估 | Total genes less than 100, insufficient for splitting and evaluation")
        
        # 计算测试集大小 | Calculate test set size
        # 注意：randomSplit.pl 实际上把指定的数量放在 .test 文件中
        test_count = int(total_genes * (1 - self.config.train_ratio))
        train_count = total_genes - test_count
        
        self.logger.info(f"基因总数: {total_genes} | Total genes: {total_genes}")
        self.logger.info(f"训练基因数: {train_count} ({self.config.train_ratio*100:.1f}%) | Training genes: {train_count} ({self.config.train_ratio*100:.1f}%)")
        self.logger.info(f"测试基因数: {test_count} ({(1-self.config.train_ratio)*100:.1f}%) | Testing genes: {test_count} ({(1-self.config.train_ratio)*100:.1f}%)")
        
        # 拆分数据集 | Split dataset
        random_split_script = os.path.join(self.config.augustus_path, 'randomSplit.pl')
        command = f"perl {random_split_script} {self.config.training_file} {test_count}"
        
        self.cmd_runner.run_command(command, "拆分数据集 | Split dataset")
        
        self.config.train_file = self.config.training_file + '.train'
        self.config.test_file = self.config.training_file + '.test'
        
        # 验证拆分结果 | Verify split results
        self._verify_split(total_genes, train_count, test_count)
        
        self.logger.info("数据集拆分完成 | Dataset splitting completed")
        
        return self.config.train_file, self.config.test_file
    
    def _verify_split(self, total_genes: int, expected_train: int, expected_test: int):
        """验证数据集拆分结果 | Verify dataset split results"""
        if os.path.exists(self.config.train_file) and os.path.exists(self.config.test_file):
            # 统计每个文件中的基因数量 | Count genes in each file
            with open(self.config.train_file, 'r') as f:
                actual_train_genes = f.read().count('LOCUS')
            with open(self.config.test_file, 'r') as f:
                actual_test_genes = f.read().count('LOCUS')
            
            train_size = os.path.getsize(self.config.train_file) / 1024 / 1024
            test_size = os.path.getsize(self.config.test_file) / 1024 / 1024
            
            self.logger.info(f"验证 - 训练文件: {actual_train_genes} 个基因 ({train_size:.1f} MB) | Verification - Training file: {actual_train_genes} genes ({train_size:.1f} MB)")
            self.logger.info(f"验证 - 测试文件: {actual_test_genes} 个基因 ({test_size:.1f} MB) | Verification - Testing file: {actual_test_genes} genes ({test_size:.1f} MB)")
            
            # 计算实际比例 | Calculate actual ratios
            actual_train_ratio = actual_train_genes / total_genes
            actual_test_ratio = actual_test_genes / total_genes
            
            self.logger.info(f"实际训练比例: {actual_train_ratio:.3f} | Actual training ratio: {actual_train_ratio:.3f}")
            self.logger.info(f"实际测试比例: {actual_test_ratio:.3f} | Actual testing ratio: {actual_test_ratio:.3f}")
            
            # 检查比例偏差 | Check ratio deviation
            if abs(actual_train_ratio - self.config.train_ratio) > 0.02:
                self.logger.warning(f"训练比例偏离预期: 期望 {self.config.train_ratio:.3f}, 实际 {actual_train_ratio:.3f} | Training ratio deviation: expected {self.config.train_ratio:.3f}, actual {actual_train_ratio:.3f}")
            else:
                self.logger.info("数据集拆分比例正确! | Dataset split ratio is correct!")
