"""
重复序列分析结果处理模块 | Repeat Sequence Analysis Results Processing Module
"""

import os
import pandas as pd
from pathlib import Path
from collections import defaultdict, Counter
from .utils import RepeatResultParser

class StatisticsGenerator:
    """统计信息生成器 | Statistics Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.parser = RepeatResultParser(logger)
    
    def generate_repeatmasker_stats(self) -> dict:
        """生成RepeatMasker统计信息 | Generate RepeatMasker statistics"""
        stats = {
            'total_repeats': 0,
            'total_masked_bp': 0,
            'repeat_classes': Counter(),
            'repeat_families': Counter(),
            'divergence_distribution': []
        }
        
        if hasattr(self.config, 'repeatmasker_outputs'):
            out_files = [f for f in self.config.repeatmasker_outputs if f.endswith('.out')]
            
            for out_file in out_files:
                repeats = self.parser.parse_repeatmasker_out(out_file)
                
                for repeat in repeats:
                    stats['total_repeats'] += 1
                    repeat_length = repeat['query_end'] - repeat['query_start'] + 1
                    stats['total_masked_bp'] += repeat_length
                    
                    # 分类统计 | Classification statistics
                    repeat_class = repeat['repeat_class'].split('/')[0]  # 取主要分类
                    stats['repeat_classes'][repeat_class] += 1
                    stats['repeat_families'][repeat['repeat_name']] += 1
                    
                    # 分化度分布 | Divergence distribution
                    stats['divergence_distribution'].append(repeat['divergence'])
        
        return stats
    
    def generate_trf_stats(self) -> dict:
        """生成TRF统计信息 | Generate TRF statistics"""
        stats = {
            'total_tandem_repeats': 0,
            'total_tr_bp': 0,
            'period_size_distribution': Counter(),
            'copy_number_distribution': []
        }
        
        if hasattr(self.config, 'trf_outputs'):
            for trf_file in self.config.trf_outputs:
                repeats = self.parser.parse_trf_output(trf_file)
                
                for repeat in repeats:
                    stats['total_tandem_repeats'] += 1
                    repeat_length = repeat['end'] - repeat['start'] + 1
                    stats['total_tr_bp'] += repeat_length
                    
                    stats['period_size_distribution'][repeat['period_size']] += 1
                    stats['copy_number_distribution'].append(repeat['copy_number'])
        
        return stats
    
    def generate_edta_stats(self) -> dict:
        """生成EDTA统计信息 | Generate EDTA statistics"""
        stats = {
            'total_tes': 0,
            'total_te_bp': 0,
            'te_classifications': Counter(),
            'te_families': Counter(),
            'intact_tes': 0
        }
        
        if hasattr(self.config, 'edta_outputs'):
            # 解析EDTA GFF3注释文件 | Parse EDTA GFF3 annotation files
            gff_files = [f for f in self.config.edta_outputs if f.endswith('.gff3')]
            
            for gff_file in gff_files:
                tes = self.parser.parse_edta_output(os.path.dirname(gff_file))
                
                for te in tes:
                    stats['total_tes'] += 1
                    te_length = te['end'] - te['start'] + 1
                    stats['total_te_bp'] += te_length
                    
                    # TE分类统计 | TE classification statistics
                    te_class = te['te_classification']
                    stats['te_classifications'][te_class] += 1
                    
                    te_name = te['te_name']
                    stats['te_families'][te_name] += 1
                    
                    # 完整TE统计 | Intact TE statistics
                    if 'intact' in gff_file.lower():
                        stats['intact_tes'] += 1
        
        return stats
    
    def calculate_genome_coverage(self) -> dict:
        """计算基因组覆盖率 | Calculate genome coverage"""
        # 获取基因组总长度 | Get total genome length
        from .utils import SequenceValidator
        validator = SequenceValidator(self.logger)
        genome_stats = validator.get_sequence_stats(self.config.genome_file)
        total_genome_bp = genome_stats['total_length']
        
        # RepeatMasker覆盖率 | RepeatMasker coverage
        rm_stats = self.generate_repeatmasker_stats()
        rm_coverage = (rm_stats['total_masked_bp'] / total_genome_bp * 100) if total_genome_bp > 0 else 0
        
        # TRF覆盖率 | TRF coverage
        trf_stats = self.generate_trf_stats()
        trf_coverage = (trf_stats['total_tr_bp'] / total_genome_bp * 100) if total_genome_bp > 0 else 0
        
        # EDTA覆盖率 | EDTA coverage
        edta_stats = self.generate_edta_stats()
        edta_coverage = (edta_stats['total_te_bp'] / total_genome_bp * 100) if total_genome_bp > 0 else 0
        
        coverage_stats = {
            'total_genome_bp': total_genome_bp,
            'repeatmasker_coverage_percent': rm_coverage,
            'trf_coverage_percent': trf_coverage,
            'edta_coverage_percent': edta_coverage,
            'repeatmasker_stats': rm_stats,
            'trf_stats': trf_stats,
            'edta_stats': edta_stats
        }
        
        return coverage_stats

class SummaryGenerator:
    """总结生成器 | Summary Generator"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.stats_generator = StatisticsGenerator(config, logger)
    
    def generate_summary_report(self):
        """生成总结报告 | Generate summary report"""
        report_file = os.path.join(self.config.output_dir, "repeat_analysis_summary.txt")
        
        # 生成统计信息 | Generate statistics
        coverage_stats = self.stats_generator.calculate_genome_coverage()
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("基因组重复序列分析总结报告 | Genome Repeat Sequence Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            # 输入文件信息 | Input file information
            f.write("输入文件信息 | Input File Information:\n")
            f.write(f"  基因组文件 | Genome file: {self.config.genome_file}\n")
            f.write(f"  物种 | Species: {self.config.species}\n")
            f.write(f"  基因组总长度 | Total genome length: {coverage_stats['total_genome_bp']:,} bp\n")
            if self.config.custom_lib:
                f.write(f"  自定义重复库 | Custom repeat library: {self.config.custom_lib}\n")
            f.write("\n")
            
            # 分析参数 | Analysis parameters
            f.write("分析参数 | Analysis Parameters:\n")
            f.write(f"  RepeatMasker线程数 | RepeatMasker threads: {self.config.rm_threads}\n")
            f.write(f"  软屏蔽 | Soft masking: {'是 | Yes' if self.config.rm_soft_mask else '否 | No'}\n")
            f.write(f"  生成GFF文件 | Generate GFF: {'是 | Yes' if self.config.rm_generate_gff else '否 | No'}\n")
            f.write(f"  运行RepeatModeler | Run RepeatModeler: {'是 | Yes' if self.config.rm_run_modeler else '否 | No'}\n")
            f.write(f"  运行TRF分析 | Run TRF analysis: {'是 | Yes' if self.config.trf_run_analysis else '否 | No'}\n")
            f.write(f"  运行EDTA分析 | Run EDTA analysis: {'是 | Yes' if self.config.edta_run_analysis else '否 | No'}\n")
            if self.config.edta_run_analysis:
                f.write(f"  EDTA物种类型 | EDTA species type: {self.config.edta_species}\n")
                f.write(f"  EDTA分析步骤 | EDTA analysis step: {self.config.edta_step}\n")
            f.write("\n")
            
            # RepeatMasker结果 | RepeatMasker results
            rm_stats = coverage_stats['repeatmasker_stats']
            f.write("RepeatMasker分析结果 | RepeatMasker Analysis Results:\n")
            f.write(f"  重复序列总数 | Total repeats: {rm_stats['total_repeats']:,}\n")
            f.write(f"  屏蔽碱基数 | Masked base pairs: {rm_stats['total_masked_bp']:,} bp\n")
            f.write(f"  基因组覆盖率 | Genome coverage: {coverage_stats['repeatmasker_coverage_percent']:.2f}%\n")
            
            # 重复序列分类 | Repeat classification
            if rm_stats['repeat_classes']:
                f.write("  主要重复序列类型 | Major repeat classes:\n")
                for class_name, count in rm_stats['repeat_classes'].most_common(10):
                    percentage = (count / rm_stats['total_repeats'] * 100) if rm_stats['total_repeats'] > 0 else 0
                    f.write(f"    {class_name}: {count:,} ({percentage:.1f}%)\n")
            f.write("\n")
            
            # TRF结果 | TRF results
            trf_stats = coverage_stats['trf_stats']
            f.write("TRF分析结果 | TRF Analysis Results:\n")
            f.write(f"  串联重复总数 | Total tandem repeats: {trf_stats['total_tandem_repeats']:,}\n")
            f.write(f"  串联重复碱基数 | Tandem repeat base pairs: {trf_stats['total_tr_bp']:,} bp\n")
            f.write(f"  基因组覆盖率 | Genome coverage: {coverage_stats['trf_coverage_percent']:.2f}%\n")
            
            # 周期长度分布 | Period size distribution
            if trf_stats['period_size_distribution']:
                f.write("  周期长度分布 | Period size distribution:\n")
                for period_size, count in sorted(trf_stats['period_size_distribution'].items()):
                    f.write(f"    {period_size}bp: {count:,}\n")
            f.write("\n")
            
            # EDTA结果 | EDTA results
            edta_stats = coverage_stats['edta_stats']
            f.write("EDTA分析结果 | EDTA Analysis Results:\n")
            f.write(f"  转座元件总数 | Total transposable elements: {edta_stats['total_tes']:,}\n")
            f.write(f"  转座元件碱基数 | Transposable element base pairs: {edta_stats['total_te_bp']:,} bp\n")
            f.write(f"  基因组覆盖率 | Genome coverage: {coverage_stats['edta_coverage_percent']:.2f}%\n")
            f.write(f"  完整转座元件数 | Intact transposable elements: {edta_stats['intact_tes']:,}\n")
            
            # TE分类分布 | TE classification distribution
            if edta_stats['te_classifications']:
                f.write("  转座元件分类 | Transposable element classifications:\n")
                for te_class, count in edta_stats['te_classifications'].most_common(10):
                    percentage = (count / edta_stats['total_tes'] * 100) if edta_stats['total_tes'] > 0 else 0
                    f.write(f"    {te_class}: {count:,} ({percentage:.1f}%)\n")
            f.write("\n")
            
            # 输出文件列表 | Output file list
            f.write("输出文件 | Output Files:\n")
            if hasattr(self.config, 'repeatmasker_outputs'):
                f.write("  RepeatMasker输出文件 | RepeatMasker output files:\n")
                for file in self.config.repeatmasker_outputs:
                    f.write(f"    {file}\n")
            
            if hasattr(self.config, 'trf_outputs'):
                f.write("  TRF输出文件 | TRF output files:\n")
                for file in self.config.trf_outputs:
                    f.write(f"    {file}\n")
            
            if hasattr(self.config, 'edta_outputs'):
                f.write("  EDTA输出文件 | EDTA output files:\n")
                for file in self.config.edta_outputs:
                    f.write(f"    {os.path.basename(file)}\n")
            
            if hasattr(self.config, 'custom_repeat_lib'):
                f.write(f"  RepeatModeler重复库 | RepeatModeler library: {self.config.custom_repeat_lib}\n")
            
            if hasattr(self.config, 'edta_lib'):
                f.write(f"  EDTA转座元件库 | EDTA TE library: {os.path.basename(self.config.edta_lib)}\n")
        
        self.logger.info(f"分析总结报告已生成 | Analysis summary report generated: {report_file}")
    
    def export_statistics_table(self):
        """导出统计表格 | Export statistics table"""
        coverage_stats = self.stats_generator.calculate_genome_coverage()
        
        # 创建RepeatMasker统计表 | Create RepeatMasker statistics table
        rm_stats = coverage_stats['repeatmasker_stats']
        if rm_stats['repeat_classes']:
            rm_df = pd.DataFrame([
                {'Repeat_Class': class_name, 'Count': count, 
                 'Percentage': count / rm_stats['total_repeats'] * 100}
                for class_name, count in rm_stats['repeat_classes'].items()
            ])
            rm_table_file = os.path.join(self.config.output_dir, "repeatmasker_classification.csv")
            rm_df.to_csv(rm_table_file, index=False)
            self.logger.info(f"RepeatMasker分类统计表已导出 | RepeatMasker classification table exported: {rm_table_file}")
        
        # 创建TRF统计表 | Create TRF statistics table
        trf_stats = coverage_stats['trf_stats']
        if trf_stats['period_size_distribution']:
            trf_df = pd.DataFrame([
                {'Period_Size': period_size, 'Count': count}
                for period_size, count in trf_stats['period_size_distribution'].items()
            ])
            trf_table_file = os.path.join(self.config.output_dir, "trf_period_distribution.csv")
            trf_df.to_csv(trf_table_file, index=False)
            self.logger.info(f"TRF周期分布统计表已导出 | TRF period distribution table exported: {trf_table_file}")
        
        # 创建EDTA统计表 | Create EDTA statistics table
        edta_stats = coverage_stats['edta_stats']
        if edta_stats['te_classifications']:
            edta_df = pd.DataFrame([
                {'TE_Classification': te_class, 'Count': count, 
                 'Percentage': count / edta_stats['total_tes'] * 100}
                for te_class, count in edta_stats['te_classifications'].items()
            ])
            edta_table_file = os.path.join(self.config.output_dir, "edta_te_classification.csv")
            edta_df.to_csv(edta_table_file, index=False)
            self.logger.info(f"EDTA转座元件分类统计表已导出 | EDTA TE classification table exported: {edta_table_file}")
