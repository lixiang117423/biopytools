"""
Augustus模型评估模块 | Augustus Model Evaluation Module
"""

import re
import os
import pandas as pd
from datetime import datetime

class ModelEvaluator:
    """模型评估器 | Model Evaluator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def evaluate_predictions(self):
        """评估预测结果 | Evaluate predictions"""
        self.logger.info("="*50)
        self.logger.info("步骤 6: 解析评估结果 | Step 6: Parsing evaluation results")
        
        with open(self.config.prediction_file, 'r') as f:
            content = f.read()
        
        # 提取评估指标 | Extract evaluation metrics
        evaluation_data = self._extract_evaluation_metrics(content)
        
        # 生成Excel报告 | Generate Excel report
        self._generate_excel_report(evaluation_data)
        
        return evaluation_data
    
    def _extract_evaluation_metrics(self, content: str) -> dict:
        """提取评估指标 | Extract evaluation metrics"""
        evaluation = {}
        
        # 提取核苷酸水平敏感性和特异性 | Extract nucleotide level sensitivity and specificity
        nucleotide_pattern = r'nucleotide level\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|'
        nucleotide_match = re.search(nucleotide_pattern, content)
        if nucleotide_match:
            evaluation['nucleotide_sensitivity'] = float(nucleotide_match.group(1))
            evaluation['nucleotide_specificity'] = float(nucleotide_match.group(2))
        
        # 提取外显子水平数据 | Extract exon level data
        exon_pattern = r'exon level\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|.*?\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|'
        exon_match = re.search(exon_pattern, content, re.DOTALL)
        if exon_match:
            evaluation['exon_pred_total'] = int(exon_match.group(1))
            evaluation['exon_anno_total'] = int(exon_match.group(2))
            evaluation['exon_tp'] = int(exon_match.group(3))
            evaluation['exon_sensitivity'] = float(exon_match.group(4))
            evaluation['exon_specificity'] = float(exon_match.group(5))
        
        # 提取基因水平数据 | Extract gene level data
        gene_pattern = r'gene level\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*(\d+)\s*\|\s*([\d.]+)\s*\|\s*([\d.]+)\s*\|'
        gene_match = re.search(gene_pattern, content)
        if gene_match:
            evaluation['gene_pred'] = int(gene_match.group(1))
            evaluation['gene_anno'] = int(gene_match.group(2))
            evaluation['gene_tp'] = int(gene_match.group(3))
            evaluation['gene_fp'] = int(gene_match.group(4))
            evaluation['gene_fn'] = int(gene_match.group(5))
            evaluation['gene_sensitivity'] = float(gene_match.group(6))
            evaluation['gene_specificity'] = float(gene_match.group(7))
        
        return evaluation
    
    def _generate_excel_report(self, evaluation_data: dict):
        """生成Excel评估报告 | Generate Excel evaluation report"""
        self.logger.info("生成Excel评估报告 | Generating Excel evaluation report")
        
        # 创建英文版评估结果 | Create English evaluation results
        results_data_en = self._create_english_results(evaluation_data)
        results_data_zh = self._create_chinese_results(evaluation_data)
        
        # 创建配置信息 | Create configuration data
        config_data_en = self._create_english_config()
        config_data_zh = self._create_chinese_config()
        
        # 保存Excel文件 | Save Excel files
        self._save_excel_reports(results_data_en, results_data_zh, config_data_en, config_data_zh)
    
    def _create_english_results(self, evaluation_data: dict) -> list:
        """创建英文评估结果 | Create English evaluation results"""
        results_data = []
        
        # 核苷酸水平 | Nucleotide level
        if 'nucleotide_sensitivity' in evaluation_data:
            results_data.extend([
                {
                    'Evaluation Level': 'Nucleotide Level',
                    'Metric': 'Sensitivity',
                    'Value': evaluation_data['nucleotide_sensitivity'],
                    'Description': 'Proportion of correctly predicted nucleotides, reflects model ability to find true genes'
                },
                {
                    'Evaluation Level': 'Nucleotide Level',
                    'Metric': 'Specificity',
                    'Value': evaluation_data['nucleotide_specificity'],
                    'Description': 'Proportion of accurately predicted nucleotides, reflects model prediction precision'
                }
            ])
        
        # 外显子水平 | Exon level
        if 'exon_sensitivity' in evaluation_data:
            results_data.extend([
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Total Predicted Exons',
                    'Value': evaluation_data['exon_pred_total'],
                    'Description': 'Total number of exons predicted by the model'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Total Annotated Exons',
                    'Value': evaluation_data['exon_anno_total'],
                    'Description': 'Total number of exons in reference annotation'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'True Positives',
                    'Value': evaluation_data['exon_tp'],
                    'Description': 'Number of correctly predicted exons (True Positive)'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Sensitivity',
                    'Value': evaluation_data['exon_sensitivity'],
                    'Description': 'Proportion of correctly predicted exons among true exons'
                },
                {
                    'Evaluation Level': 'Exon Level',
                    'Metric': 'Specificity',
                    'Value': evaluation_data['exon_specificity'],
                    'Description': 'Proportion of correct exons among predicted exons'
                }
            ])
        
        # 基因水平 | Gene level
        if 'gene_sensitivity' in evaluation_data:
            results_data.extend([
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Predicted Genes',
                    'Value': evaluation_data['gene_pred'],
                    'Description': 'Total number of genes predicted by the model'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Annotated Genes',
                    'Value': evaluation_data['gene_anno'],
                    'Description': 'Total number of genes in reference annotation'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'True Positives (TP)',
                    'Value': evaluation_data['gene_tp'],
                    'Description': 'Number of completely correctly predicted genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'False Positives (FP)',
                    'Value': evaluation_data['gene_fp'],
                    'Description': 'Number of incorrectly predicted genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'False Negatives (FN)',
                    'Value': evaluation_data['gene_fn'],
                    'Description': 'Number of missed true genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Sensitivity',
                    'Value': evaluation_data['gene_sensitivity'],
                    'Description': 'Proportion of correctly predicted genes among true genes'
                },
                {
                    'Evaluation Level': 'Gene Level',
                    'Metric': 'Specificity',
                    'Value': evaluation_data['gene_specificity'],
                    'Description': 'Proportion of correct genes among predicted genes'
                }
            ])
        
        return results_data
    
    def _create_chinese_results(self, evaluation_data: dict) -> list:
        """创建中文评估结果 | Create Chinese evaluation results"""
        results_data = []
        
        # 核苷酸水平 | Nucleotide level
        if 'nucleotide_sensitivity' in evaluation_data:
            results_data.extend([
                {
                    '评估级别': '核苷酸水平',
                    '评估指标': '敏感性',
                    '数值': evaluation_data['nucleotide_sensitivity'],
                    '说明': '正确预测的核苷酸比例，反映模型找到真实基因的能力'
                },
                {
                    '评估级别': '核苷酸水平',
                    '评估指标': '特异性',
                    '数值': evaluation_data['nucleotide_specificity'],
                    '说明': '预测准确的核苷酸比例，反映模型预测精度'
                }
            ])
        
        # 外显子水平 | Exon level
        if 'exon_sensitivity' in evaluation_data:
            results_data.extend([
                {
                    '评估级别': '外显子水平',
                    '评估指标': '预测外显子总数',
                    '数值': evaluation_data['exon_pred_total'],
                    '说明': '模型预测的外显子总数量'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '注释外显子总数',
                    '数值': evaluation_data['exon_anno_total'],
                    '说明': '参考注释中的外显子总数量'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '正确预测数',
                    '数值': evaluation_data['exon_tp'],
                    '说明': '预测正确的外显子数量(True Positive)'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '敏感性',
                    '数值': evaluation_data['exon_sensitivity'],
                    '说明': '正确预测的外显子占真实外显子的比例'
                },
                {
                    '评估级别': '外显子水平',
                    '评估指标': '特异性',
                    '数值': evaluation_data['exon_specificity'],
                    '说明': '预测的外显子中正确的比例'
                }
            ])
        
        # 基因水平 | Gene level
        if 'gene_sensitivity' in evaluation_data:
            results_data.extend([
                {
                    '评估级别': '基因水平',
                    '评估指标': '预测基因数',
                    '数值': evaluation_data['gene_pred'],
                    '说明': '模型预测的基因总数'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '注释基因数',
                    '数值': evaluation_data['gene_anno'],
                    '说明': '参考注释中的基因总数'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '真阳性(TP)',
                    '数值': evaluation_data['gene_tp'],
                    '说明': '完全正确预测的基因数量'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '假阳性(FP)',
                    '数值': evaluation_data['gene_fp'],
                    '说明': '错误预测的基因数量'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '假阴性(FN)',
                    '数值': evaluation_data['gene_fn'],
                    '说明': '漏掉的真实基因数量'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '敏感性',
                    '数值': evaluation_data['gene_sensitivity'],
                    '说明': '正确预测的基因占真实基因的比例'
                },
                {
                    '评估级别': '基因水平',
                    '评估指标': '特异性',
                    '数值': evaluation_data['gene_specificity'],
                    '说明': '预测基因中正确的比例'
                }
            ])
        
        return results_data
    
    def _create_english_config(self) -> list:
        """创建英文配置信息 | Create English configuration data"""
        return [
            ['Species Name', self.config.species_name],
            ['Genome File', self.config.genome_file],
            ['Annotation File', self.config.gff_file],
            ['Training Ratio', f"{self.config.train_ratio*100}%"],
            ['Flank Length', f"{self.config.flank_length} bp"],
            ['Output Directory', self.config.output_dir],
            ['Augustus Path', self.config.augustus_path],
            ['Generation Time', datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
        ]
    
    def _create_chinese_config(self) -> list:
        """创建中文配置信息 | Create Chinese configuration data"""
        return [
            ['物种名称', self.config.species_name],
            ['基因组文件', self.config.genome_file],
            ['注释文件', self.config.gff_file],
            ['训练集比例', f"{self.config.train_ratio*100}%"],
            ['侧翼长度', f"{self.config.flank_length} bp"],
            ['输出目录', self.config.output_dir],
            ['Augustus路径', self.config.augustus_path],
            ['生成时间', datetime.now().strftime('%Y-%m-%d %H:%M:%S')]
        ]
    
    def _save_excel_reports(self, results_en: list, results_zh: list, config_en: list, config_zh: list):
        """保存Excel报告 | Save Excel reports"""
        # 创建DataFrames | Create DataFrames
        df_results_en = pd.DataFrame(results_en)
        df_results_zh = pd.DataFrame(results_zh)
        df_config_en = pd.DataFrame(config_en, columns=['Parameter', 'Value'])
        df_config_zh = pd.DataFrame(config_zh, columns=['参数', '值'])
        
        # 保存Excel文件 | Save Excel files
        excel_file_en = os.path.join(self.config.output_dir, 'augustus_evaluation_report_EN.xlsx')
        excel_file_zh = os.path.join(self.config.output_dir, 'augustus_evaluation_report_ZH.xlsx')
        
        # 英文版本 | English version
        with pd.ExcelWriter(excel_file_en, engine='openpyxl') as writer:
            df_config_en.to_excel(writer, sheet_name='Configuration', index=False)
            df_results_en.to_excel(writer, sheet_name='Evaluation Results', index=False)
            
            # 添加说明页 | Add explanation sheet
            explanation_data_en = [
                ['Term', 'Explanation'],
                ['Sensitivity', 'Also called recall, represents model ability to correctly identify true genes. Formula: TP/(TP+FN)'],
                ['Specificity', 'Represents model prediction accuracy. Formula: TP/(TP+FP)'],
                ['TP (True Positive)', 'Number of correctly predicted genes'],
                ['FP (False Positive)', 'Number of incorrectly predicted genes'],
                ['FN (False Negative)', 'Number of missed true genes'],
                ['Nucleotide Level', 'Prediction accuracy at DNA sequence base level'],
                ['Exon Level', 'Prediction accuracy at exon structure level'],
                ['Gene Level', 'Prediction accuracy at complete gene level'],
                ['Evaluation Suggestion', 'Generally, models with sensitivity>0.8 and specificity>0.8 are considered excellent']
            ]
            df_explanation_en = pd.DataFrame(explanation_data_en[1:], columns=explanation_data_en[0])
            df_explanation_en.to_excel(writer, sheet_name='Term Explanations', index=False)
        
        # 中文版本 | Chinese version
        with pd.ExcelWriter(excel_file_zh, engine='openpyxl') as writer:
            df_config_zh.to_excel(writer, sheet_name='配置信息', index=False)
            df_results_zh.to_excel(writer, sheet_name='评估结果', index=False)
            
            # 添加说明页 | Add explanation sheet
            explanation_data_zh = [
                ['术语', '解释'],
                ['敏感性(Sensitivity)', '也称召回率，表示模型正确识别真实基因的能力，计算公式: TP/(TP+FN)'],
                ['特异性(Specificity)', '表示模型预测准确度，计算公式: TP/(TP+FP)'],
                ['TP (True Positive)', '真阳性，正确预测的基因数量'],
                ['FP (False Positive)', '假阳性，错误预测的基因数量'],
                ['FN (False Negative)', '假阴性，漏掉的真实基因数量'],
                ['核苷酸水平', '在DNA序列碱基层面的预测准确性'],
                ['外显子水平', '在外显子结构层面的预测准确性'],
                ['基因水平', '在完整基因层面的预测准确性'],
                ['评估建议', '一般认为敏感性>0.8、特异性>0.8的模型较为优秀']
            ]
            df_explanation_zh = pd.DataFrame(explanation_data_zh[1:], columns=explanation_data_zh[0])
            df_explanation_zh.to_excel(writer, sheet_name='术语解释', index=False)
        
        self.logger.info(f"Excel评估报告已生成: | Excel evaluation reports generated:")
        self.logger.info(f"  英文版本: {excel_file_en} | English version: {excel_file_en}")
        self.logger.info(f"  中文版本: {excel_file_zh} | Chinese version: {excel_file_zh}")
