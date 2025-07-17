"""
Augustus预测结果处理模块 | Augustus Prediction Results Processing Module
"""

import os
from pathlib import Path
from datetime import datetime

class ResultsComparator:
    """结果比较器 | Results Comparator"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def compare_predictions(self, model_only_file: str, rnaseq_file: str, hints_file: str) -> dict:
        """比较预测结果 | Compare prediction results"""
        self.logger.info("比较预测结果 | Comparing prediction results...")
        
        comparison_results = {}
        
        # 统计方案A (仅模型)
        model_stats = self._count_features(model_only_file)
        comparison_results['model_only'] = model_stats
        
        # 统计方案B (结合转录组)
        rnaseq_stats = self._count_features(rnaseq_file)
        comparison_results['with_rnaseq'] = rnaseq_stats
        
        # 转录组支持情况
        support_stats = self._calculate_support(rnaseq_file, hints_file)
        comparison_results['support'] = support_stats
        
        # 输出比较结果
        self._log_comparison_results(comparison_results)
        
        return comparison_results
    
    def _count_features(self, gff_file: str) -> dict:
        """统计GFF文件中的特征 | Count features in GFF file"""
        import os
        
        if not os.path.exists(gff_file):
            return {'genes': 0, 'mRNAs': 0, 'CDSs': 0, 'exons': 0}
        
        feature_counts = {'genes': 0, 'mRNAs': 0, 'CDSs': 0, 'exons': 0}
        
        try:
            with open(gff_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        feature_type = parts[2]
                        if feature_type == 'gene':
                            feature_counts['genes'] += 1
                        elif feature_type == 'mRNA':
                            feature_counts['mRNAs'] += 1
                        elif feature_type == 'CDS':
                            feature_counts['CDSs'] += 1
                        elif feature_type == 'exon':
                            feature_counts['exons'] += 1
        
        except Exception as e:
            self.logger.warning(f"统计特征失败 | Failed to count features: {e}")
        
        return feature_counts
    
    def _calculate_support(self, gff_file: str, hints_file: str) -> dict:
        """计算转录组支持情况 | Calculate RNA-seq support"""
        support_stats = {'supported_genes': 0, 'total_genes': 0, 'support_ratio': 0.0}
        
        try:
            # 检查bedtools是否可用
            if not self._check_bedtools():
                self.logger.warning("需要安装bedtools来计算转录组支持情况 | Need bedtools to calculate RNA-seq support")
                return support_stats
            
            # 使用bedtools计算支持情况
            genes_with_support_file = self.config.output_path / "genes_with_support.gff"
            bedtools_cmd = f"bedtools intersect -a {gff_file} -b {hints_file} -c > {genes_with_support_file}"
            
            success, _ = self.cmd_runner.run(bedtools_cmd, "计算基因支持情况 | Calculate gene support")
            
            if success:
                supported_genes = 0
                total_genes = 0
                
                with open(genes_with_support_file, 'r') as f:
                    for line in f:
                        if '\tgene\t' in line:
                            total_genes += 1
                            # 最后一列是支持数
                            support_count = int(line.strip().split('\t')[-1])
                            if support_count > 0:
                                supported_genes += 1
                
                support_ratio = (supported_genes / total_genes * 100) if total_genes > 0 else 0
                
                support_stats = {
                    'supported_genes': supported_genes,
                    'total_genes': total_genes,
                    'support_ratio': support_ratio
                }
        
        except Exception as e:
            self.logger.warning(f"计算支持情况失败 | Failed to calculate support: {e}")
        
        return support_stats
    
    def _check_bedtools(self) -> bool:
        """检查bedtools是否可用 | Check if bedtools is available"""
        try:
            self.cmd_runner.run("bedtools --version", "检查bedtools | Check bedtools", capture_output=True)
            return True
        except:
            return False
    
    def _log_comparison_results(self, results: dict):
        """记录比较结果 | Log comparison results"""
        self.logger.info("=== Augustus预测结果比较 | Augustus Prediction Results Comparison ===")
        
        model_stats = results.get('model_only', {})
        rnaseq_stats = results.get('with_rnaseq', {})
        support_stats = results.get('support', {})
        
        self.logger.info("\n方案A (仅模型) | Method A (Model only):")
        self.logger.info(f"  基因数 | Genes: {model_stats.get('genes', 0)}")
        self.logger.info(f"  mRNA数 | mRNAs: {model_stats.get('mRNAs', 0)}")
        self.logger.info(f"  CDS数 | CDSs: {model_stats.get('CDSs', 0)}")
        self.logger.info(f"  外显子数 | Exons: {model_stats.get('exons', 0)}")
        
        self.logger.info("\n方案B (结合转录组) | Method B (With RNA-seq):")
        self.logger.info(f"  基因数 | Genes: {rnaseq_stats.get('genes', 0)}")
        self.logger.info(f"  mRNA数 | mRNAs: {rnaseq_stats.get('mRNAs', 0)}")
        self.logger.info(f"  CDS数 | CDSs: {rnaseq_stats.get('CDSs', 0)}")
        self.logger.info(f"  外显子数 | Exons: {rnaseq_stats.get('exons', 0)}")
        
        if support_stats.get('total_genes', 0) > 0:
            self.logger.info("\n转录组支持情况 | RNA-seq Support:")
            self.logger.info(f"  有转录组支持的基因 | Genes with RNA-seq support: {support_stats.get('supported_genes', 0)} / {support_stats.get('total_genes', 0)}")
            self.logger.info(f"  支持比例 | Support ratio: {support_stats.get('support_ratio', 0):.2f}%")

class ReportGenerator:
    """报告生成器 | Report Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_final_report(self, comparison_results: dict, bam_files: list):
        """生成最终报告 | Generate final report"""
        report_file = self.config.output_path / "final_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("多转录组Augustus基因预测报告\n")
            f.write("Multiple RNA-seq Augustus Gene Prediction Report\n")
            f.write("=" * 50 + "\n\n")
            
            # 基本信息
            f.write(f"处理时间 | Processing time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"基因组文件 | Genome file: {self.config.genome_file}\n")
            f.write(f"Augustus模型 | Augustus model: {self.config.species_model}\n")
            f.write(f"转录组样本数 | RNA-seq samples: {len(self.config.samples)}\n\n")
            
            # 样本信息
            f.write("样本信息 | Sample Information:\n")
            for i, sample in enumerate(self.config.samples, 1):
                f.write(f"  {i}. {sample['name']}: {sample['r1_file']}, {sample['r2_file']}\n")
            f.write("\n")
            
            # 预测结果比较
            f.write("预测结果比较 | Prediction Results Comparison:\n")
            model_stats = comparison_results.get('model_only', {})
            rnaseq_stats = comparison_results.get('with_rnaseq', {})
            support_stats = comparison_results.get('support', {})
            
            f.write(f"方案A (仅模型) | Method A (Model only):\n")
            f.write(f"  基因数 | Genes: {model_stats.get('genes', 0)}\n")
            f.write(f"  mRNA数 | mRNAs: {model_stats.get('mRNAs', 0)}\n")
            f.write(f"  CDS数 | CDSs: {model_stats.get('CDSs', 0)}\n")
            f.write(f"  外显子数 | Exons: {model_stats.get('exons', 0)}\n\n")
            
            f.write(f"方案B (结合转录组) | Method B (With RNA-seq):\n")
            f.write(f"  基因数 | Genes: {rnaseq_stats.get('genes', 0)}\n")
            f.write(f"  mRNA数 | mRNAs: {rnaseq_stats.get('mRNAs', 0)}\n")
            f.write(f"  CDS数 | CDSs: {rnaseq_stats.get('CDSs', 0)}\n")
            f.write(f"  外显子数 | Exons: {rnaseq_stats.get('exons', 0)}\n\n")
            
            if support_stats.get('total_genes', 0) > 0:
                f.write(f"转录组支持情况 | RNA-seq Support:\n")
                f.write(f"  有转录组支持的基因 | Genes with RNA-seq support: {support_stats.get('supported_genes', 0)} / {support_stats.get('total_genes', 0)}\n")
                f.write(f"  支持比例 | Support ratio: {support_stats.get('support_ratio', 0):.2f}%\n\n")
            
            # 输出文件
            f.write("结果文件 | Result Files:\n")
            f.write("- augustus_model_only.gff3: 仅使用Augustus模型的预测结果 | Augustus model-only prediction\n")
            f.write("- augustus_with_multi_rnaseq.gff3: 结合多转录组的预测结果 | Augustus prediction with multiple RNA-seq\n")
            f.write("- filtered_hints.gff: 合并过滤后的转录组hints | Merged and filtered RNA-seq hints\n")
            f.write("- */: 各样本的比对和处理结果 | Alignment and processing results for each sample\n\n")
            
            f.write("建议后续分析 | Suggested Further Analysis:\n")
            f.write("- 使用BUSCO评估基因集完整性 | Use BUSCO to assess gene set completeness\n")
            f.write("- 使用gffcompare与参考注释比较 | Use gffcompare to compare with reference annotation\n")
            f.write("- 进行功能注释分析 | Perform functional annotation analysis\n")
        
        self.logger.info(f"最终报告已生成 | Final report generated: {report_file}")
