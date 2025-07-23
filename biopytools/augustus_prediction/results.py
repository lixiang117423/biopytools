"""
Augustus结果汇总模块 | Augustus Results Summary Module
"""

from pathlib import Path
from datetime import datetime

class ResultsSummary:
    """结果汇总生成器 | Results Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary(self, evaluation_data: dict):
        """生成汇总报告 | Generate summary report"""
        summary_file = self.config.output_path / 'augustus_pipeline_summary.txt'
        
        try:
            with summary_file.open('w', encoding='utf-8') as f:
                f.write("="*70 + "\n")
                f.write("Augustus基因预测流水线汇总报告 | Augustus Gene Prediction Pipeline Summary\n")
                f.write("="*70 + "\n\n")
                
                # 基本信息 | Basic information
                f.write("基本信息 | Basic Information:\n")
                f.write(f"  - 物种名称 | Species Name: {self.config.species_name}\n")
                f.write(f"  - 基因组文件 | Genome File: {self.config.genome_file}\n")
                f.write(f"  - 注释文件 | Annotation File: {self.config.gff_file}\n")
                f.write(f"  - 输出目录 | Output Directory: {self.config.output_dir}\n")
                f.write(f"  - 训练比例 | Training Ratio: {self.config.train_ratio*100}%\n")
                f.write(f"  - 侧翼长度 | Flank Length: {self.config.flank_length} bp\n")
                f.write(f"  - 生成时间 | Generated Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write("\n")
                
                # 关键评估结果 | Key evaluation results
                f.write("关键评估结果 | Key Evaluation Results:\n")
                if evaluation_data:
                    if 'nucleotide_sensitivity' in evaluation_data:
                        f.write(f"  - 核苷酸敏感性 | Nucleotide Sensitivity: {evaluation_data['nucleotide_sensitivity']:.3f}\n")
                        f.write(f"  - 核苷酸特异性 | Nucleotide Specificity: {evaluation_data['nucleotide_specificity']:.3f}\n")
                    
                    if 'gene_sensitivity' in evaluation_data:
                        f.write(f"  - 基因敏感性 | Gene Sensitivity: {evaluation_data['gene_sensitivity']:.3f}\n")
                        f.write(f"  - 基因特异性 | Gene Specificity: {evaluation_data['gene_specificity']:.3f}\n")
                        f.write(f"  - 预测基因数 | Predicted Genes: {evaluation_data['gene_pred']}\n")
                        f.write(f"  - 注释基因数 | Annotated Genes: {evaluation_data['gene_anno']}\n")
                        f.write(f"  - 正确预测 | True Positives: {evaluation_data['gene_tp']}\n")
                        f.write(f"  - 错误预测 | False Positives: {evaluation_data['gene_fp']}\n")
                        f.write(f"  - 漏掉基因 | False Negatives: {evaluation_data['gene_fn']}\n")
                    
                    if 'exon_sensitivity' in evaluation_data:
                        f.write(f"  - 外显子敏感性 | Exon Sensitivity: {evaluation_data['exon_sensitivity']:.3f}\n")
                        f.write(f"  - 外显子特异性 | Exon Specificity: {evaluation_data['exon_specificity']:.3f}\n")
                else:
                    f.write("  - 状态 | Status: 评估数据提取失败 | Evaluation data extraction failed\n")
                
                f.write("\n")
                
                # 输出文件列表 | Output files list
                f.write("输出文件 | Output Files:\n")
                for file_path in sorted(self.config.output_path.iterdir()):
                    if file_path.is_file():
                        size_mb = file_path.stat().st_size / 1024 / 1024
                        f.write(f"  - {file_path.name} ({size_mb:.2f} MB)\n")
                
                f.write("\n")
                
                # 模型质量评估 | Model quality assessment
                f.write("模型质量评估 | Model Quality Assessment:\n")
                if evaluation_data and 'gene_sensitivity' in evaluation_data:
                    gene_sens = evaluation_data['gene_sensitivity']
                    gene_spec = evaluation_data['gene_specificity']
                    
                    if gene_sens >= 0.8 and gene_spec >= 0.8:
                        quality = "优秀 | Excellent"
                    elif gene_sens >= 0.7 and gene_spec >= 0.7:
                        quality = "良好 | Good"
                    elif gene_sens >= 0.6 and gene_spec >= 0.6:
                        quality = "一般 | Fair"
                    else:
                        quality = "需要改进 | Needs Improvement"
                    
                    f.write(f"  - 整体质量 | Overall Quality: {quality}\n")
                    f.write(f"  - 建议 | Recommendation: ")
                    
                    if gene_sens < 0.7:
                        f.write("考虑增加训练数据或调整训练参数 | Consider increasing training data or adjusting parameters")
                    elif gene_spec < 0.7:
                        f.write("考虑优化模型以减少假阳性预测 | Consider optimizing model to reduce false positives")
                    else:
                        f.write("模型性能良好，可用于基因预测 | Model performs well and is suitable for gene prediction")
                    
                    f.write("\n")
                else:
                    f.write("  - 状态 | Status: 无法评估模型质量 | Cannot assess model quality\n")
                
                f.write("\n")
                f.write("="*70 + "\n")
            
            self.logger.info(f"汇总报告已生成 | Summary report generated: {summary_file}")
            
        except Exception as e:
            self.logger.error(f"生成汇总报告失败 | Failed to generate summary report: {str(e)}")
