"""
Augustus流水线主程序模块 | Augustus Pipeline Main Module
"""

import argparse
import sys
from .config import PipelineConfig
from .utils import PipelineLogger, CommandRunner, check_dependencies
from .species_manager import SpeciesManager
from .data_processor import DataProcessor
from .model_trainer import ModelTrainer
from .predictor import GenePredictor
from .evaluator import ModelEvaluator
from .results import ResultsSummary

class AugustusPipeline:
    """Augustus流水线主类 | Main Augustus Pipeline Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PipelineConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PipelineLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.species_manager = SpeciesManager(self.config, self.logger, self.cmd_runner)
        self.data_processor = DataProcessor(self.config, self.logger, self.cmd_runner)
        self.model_trainer = ModelTrainer(self.config, self.logger, self.cmd_runner)
        self.gene_predictor = GenePredictor(self.config, self.logger, self.cmd_runner)
        self.model_evaluator = ModelEvaluator(self.config, self.logger)
        self.results_summary = ResultsSummary(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_complete_pipeline(self):
        """运行完整流水线 | Run complete pipeline"""
        try:
            self.logger.info("开始Augustus完整训练和预测流水线 | Starting Augustus complete training and prediction pipeline")
            self.logger.info(f"物种名称: {self.config.species_name} | Species name: {self.config.species_name}")
            
            # 检查依赖 | Check dependencies
            if not self.check_dependencies():
                self.logger.error("依赖检查失败，无法继续 | Dependency check failed, cannot continue")
                return False
            
            # 执行所有步骤 | Execute all steps
            # 步骤1: 创建物种模型 | Step 1: Create species model
            self.species_manager.create_species_model()
            
            # 步骤2: 准备训练数据 | Step 2: Prepare training data
            self.data_processor.prepare_training_data()
            
            # 步骤3: 拆分数据集 | Step 3: Split dataset
            self.data_processor.split_dataset()
            
            # 步骤4: 训练模型 | Step 4: Train model
            self.model_trainer.train_model()
            
            # 步骤5: 预测测试集 | Step 5: Predict test set
            self.gene_predictor.predict_test_set()
            
            # 步骤6: 评估结果 | Step 6: Evaluate results
            evaluation_data = self.model_evaluator.evaluate_predictions()
            
            # 步骤7: 转换格式 | Step 7: Convert format
            self.gene_predictor.convert_to_gff3()
            
            # 生成汇总报告 | Generate summary report
            self.results_summary.generate_summary(evaluation_data)
            
            self.logger.info("="*50)
            self.logger.info("🎉 Augustus流水线执行完成! | Augustus pipeline execution completed!")
            self.logger.info(f"结果文件保存在: {self.config.output_dir} | Result files saved in: {self.config.output_dir}")
            
            # 打印关键评估结果 | Print key evaluation results
            if evaluation_data:
                self.logger.info("\n关键评估结果 | Key evaluation results:")
                if 'nucleotide_sensitivity' in evaluation_data:
                    self.logger.info(f"  核苷酸敏感性: {evaluation_data['nucleotide_sensitivity']:.3f} | Nucleotide sensitivity: {evaluation_data['nucleotide_sensitivity']:.3f}")
                    self.logger.info(f"  核苷酸特异性: {evaluation_data['nucleotide_specificity']:.3f} | Nucleotide specificity: {evaluation_data['nucleotide_specificity']:.3f}")
                if 'gene_sensitivity' in evaluation_data:
                    self.logger.info(f"  基因敏感性: {evaluation_data['gene_sensitivity']:.3f} | Gene sensitivity: {evaluation_data['gene_sensitivity']:.3f}")
                    self.logger.info(f"  基因特异性: {evaluation_data['gene_specificity']:.3f} | Gene specificity: {evaluation_data['gene_specificity']:.3f}")
            
            return True
            
        except Exception as e:
            self.logger.error(f"流水线执行失败: {e} | Pipeline execution failed: {e}")
            return False

def create_argument_parser():
    """创建命令行参数解析器 | Create command line argument parser"""
    parser = argparse.ArgumentParser(
        description='Augustus基因预测完整流水线 | Augustus Gene Prediction Complete Pipeline',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
使用示例 | Usage Examples:
  # 基本用法 | Basic usage
  python run_augustus_pipeline.py -s Rice_NLR_Model -g genome.fa -a annotations.gff3
  
  # 自定义参数 | Custom parameters
  python run_augustus_pipeline.py -s MySpecies -g genome.fa -a genes.gff3 -o results -t 0.8 -f 1000
  
  # 指定Augustus路径 | Specify Augustus path
  python run_augustus_pipeline.py -s MySpecies -g genome.fa -a genes.gff3 --augustus-path /opt/augustus/bin

详细说明 | Detailed Description:
  该脚本自动执行完整的Augustus训练和预测流水线，包括模型训练、参数优化、
  预测评估和结果报告生成。| This script automatically executes the complete 
  Augustus training and prediction pipeline, including model training, parameter 
  optimization, prediction evaluation, and result report generation.
        """
    )
    
    # 必需参数 | Required parameters
    parser.add_argument(
        '--species-name', '-s',
        required=True,
        help='新物种模型名称 | New species model name (e.g., Rice_NLR_Model)'
    )
    
    parser.add_argument(
        '--genome-file', '-g',
        required=True,
        help='基因组FASTA文件路径 | Genome FASTA file path'
    )
    
    parser.add_argument(
        '--gff-file', '-a',
        required=True,
        help='基因注释GFF3文件路径 | Gene annotation GFF3 file path'
    )
    
    # 可选参数 | Optional parameters
    parser.add_argument(
        '--output-dir', '-o',
        default='./augustus_output',
        help='输出目录路径 | Output directory path (default: ./augustus_output)'
    )
    
    parser.add_argument(
        '--augustus-path', '-p',
        default='/share/org/YZWL/yzwl_lixg/miniforge3/envs/Augustus_v.3.5.0/bin',
        help='Augustus安装路径 | Augustus installation path'
    )
    
    parser.add_argument(
        '--train-ratio', '-t',
        type=float,
        default=0.8,
        help='训练集比例 | Training set ratio (default: 0.8)'
    )
    
    parser.add_argument(
        '--flank-length', '-f',
        type=int,
        default=1000,
        help='基因侧翼长度 | Gene flanking length (default: 1000)'
    )
    
    return parser

def main():
    """主函数 | Main function"""
    parser = create_argument_parser()
    args = parser.parse_args()
    
    try:
        # 准备配置参数 | Prepare configuration parameters
        config_kwargs = {
            'species_name': args.species_name,
            'genome_file': args.genome_file,
            'gff_file': args.gff_file,
            'output_dir': args.output_dir,
            'augustus_path': args.augustus_path,
            'train_ratio': args.train_ratio,
            'flank_length': args.flank_length
        }
        
        # 创建流水线 | Create pipeline
        pipeline = AugustusPipeline(**config_kwargs)
        
        # 执行流水线 | Execute pipeline
        success = pipeline.run_complete_pipeline()
        
        if success:
            print(f'\n\033[32m流水线执行完成 | Pipeline execution completed\033[0m')
            print(f'结果目录 | Results directory: {args.output_dir}')
        else:
            print('\n\033[31m流水线执行失败 | Pipeline execution failed\033[0m')
            sys.exit(1)
    
    except KeyboardInterrupt:
        print('\n\033[33m用户中断操作 | User interrupted operation\033[0m')
        sys.exit(1)
    except Exception as e:
        print(f'\n\033[31m发生错误 | Error occurred: {str(e)}\033[0m')
        sys.exit(1)

if __name__ == '__main__':
    main()