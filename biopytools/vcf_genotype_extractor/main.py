"""
VCF基因型提取主程序模块 | VCF Genotype Extraction Main Module
"""

import argparse
import sys
import os
from .config import VCFConfig
from .utils import VCFLogger, check_dependencies
from .vcf_parser import VCFParserFactory
from .genotype_processor import GenotypeProcessor
from .output_formatter import OutputFormatter
from pathlib import Path

class VCFGenotypeExtractor:
    """VCF基因型提取主类 | Main VCF Genotype Extractor Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = VCFConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = VCFLogger(Path(self.config.output_dir))
        self.logger = self.logger_manager.get_logger()
        
        # 检查依赖 | Check dependencies
        self.dependencies = check_dependencies()
        self._log_dependencies()
        
        # 初始化各个处理器 | Initialize processors
        self.parser = VCFParserFactory.create_parser(
            self.config.vcf_file, 
            self.logger, 
            use_fast=self.dependencies['cyvcf2']
        )
        self.processor = GenotypeProcessor(self.config, self.logger)
        self.formatter = OutputFormatter(self.config, self.logger)
    
    def _log_dependencies(self):
        """记录依赖信息 | Log dependencies info"""
        self.logger.info("依赖检查 | Dependency check:")
        for dep, available in self.dependencies.items():
            status = "可用 | Available" if available else "不可用 | Not available"
            self.logger.info(f"  {dep}: {status}")
    
    def _get_target_samples(self):
        """获取目标样本列表 | Get target sample list"""
        all_samples = self.parser.get_samples()
        
        if self.config.samples == "all":
            return all_samples
        elif isinstance(self.config.samples, list):
            # 验证样本是否存在 | Validate sample existence
            valid_samples = []
            for sample in self.config.samples:
                if sample in all_samples:
                    valid_samples.append(sample)
                else:
                    self.logger.warning(f"样本 '{sample}' 在VCF文件中未找到 | Sample '{sample}' not found in VCF file")
            return valid_samples
        else:
            return all_samples
    
    def run_extraction(self):
        """运行基因型提取 | Run genotype extraction"""
        try:
            self.logger.info("开始VCF基因型提取 | Starting VCF genotype extraction")
            self.logger.info(f"输入文件 | Input file: {self.config.vcf_file}")
            self.logger.info(f"输出前缀 | Output prefix: {self.config.output_prefix}")
            self.logger.info(f"输出格式 | Output format: {self.config.output_type}")
            
            # 获取目标样本 | Get target samples
            target_samples = self._get_target_samples()
            if not target_samples:
                self.logger.error("没有有效的目标样本 | No valid target samples")
                return False
            
            self.logger.info(f"目标样本数 | Number of target samples: {len(target_samples)}")
            
            # 处理VCF记录 | Process VCF records
            processed_count = 0
            for record in self.parser.parse_records(target_samples):
                row = self.processor.process_record(record)
                if row:
                    self.processor.add_record(row)
                    processed_count += 1
                
                if processed_count % 10000 == 0:
                    self.logger.info(f"已处理 | Processed {processed_count} variants")
            
            self.logger.info(f"处理完成，共处理 | Processing completed, total processed: {processed_count} variants")
            
            # 输出结果 | Output results
            self._write_results()
            
            # 生成汇总 | Generate summary
            stats = self.processor.get_summary_stats()
            self.formatter.write_summary(stats)
            
            self.logger.info("VCF基因型提取完成 | VCF genotype extraction completed")
            return True
            
        except Exception as e:
            self.logger.error(f"提取过程中发生错误 | Error during extraction: {e}")
            return False
    
    # def _write_results(self):
    #     """写入结果文件 | Write result files"""
    #     if self.config.split_by_chromosome:
    #         # 按染色体分别输出 | Output by chromosome
    #         results_by_chrom = self.processor.get_results_by_chromosome()
    #         for chrom, data in results_by_chrom.items():
    #             self.formatter.write_output(data, suffix=chrom)
    #     else:
    #         # 输出到单个文件 | Output to single file
    #         all_results = self.processor.get_all_results()
    #         self.formatter.write_output(all_results)

    def _write_results(self):
        """写入结果文件 | Write result files"""
        
        # 将 'y'/'yes' 转换为 True，其他（包括'n'/'no'）转换为 False
        should_split = str(self.config.split_by_chromosome).lower() in ['y', 'yes']

        if should_split:
            # 按染色体分别输出 | Output by chromosome
            results_by_chrom = self.processor.get_results_by_chromosome()
            for chrom, data in results_by_chrom.items():
                self.formatter.write_output(data, suffix=chrom)
        else:
            # 输出到单个文件 | Output to single file
            all_results = self.processor.get_all_results()
            self.formatter.write_output(all_results)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='VCF基因型提取工具 | VCF Genotype Extraction Tool',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter
    )
    
    # # 必需参数 | Required arguments
    # parser.add_argument('vcf_file', 
    #                    help='VCF文件路径（支持.gz压缩格式） | VCF file path (supports .gz compressed format)')

    parser.add_argument('-i', '--input', 
                   dest='vcf_file',  # 保持变量名不变，避免修改其他代码
                   required=True,    # 设为必需参数
                   help='VCF文件路径(支持.gz压缩格式) | VCF file path (supports .gz compressed format)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='vcf_genotype',
                       help='输出文件前缀 | Output file prefix')
    parser.add_argument('-s', '--samples', default='all',
                       help='样本选择：all（所有样本）或逗号分隔的样本名称 | Sample selection: all (all samples) or comma-separated sample names')
    parser.add_argument('--biallelic-only', action='store_true',
                       help='只保留双等位位点 | Keep only biallelic sites')
    # parser.add_argument('-e', '--each', action='store_true',
    #                    help='按染色体拆分输出文件 | Split output files by chromosome')
    parser.add_argument('-e', '--each', choices=['yes', 'y', 'no', 'n'], default='n',
                   help='按染色体拆分输出文件：yes/y（是）或no/n（否） | Split output files by chromosome: yes/y or no/n')
    parser.add_argument('-t', '--output-type', choices=['txt', 'csv', 'excel'], default='txt',
                       help='输出文件格式 | Output file format')
    parser.add_argument('--output-dir', default='./',
                       help='输出目录 | Output directory')
    
    args = parser.parse_args()
    
    # 创建提取器并运行 | Create extractor and run
    extractor = VCFGenotypeExtractor(
        vcf_file=args.vcf_file,
        output_prefix=args.output,
        samples=args.samples,
        biallelic_only=args.biallelic_only,
        split_by_chromosome=args.each,
        output_type=args.output_type,
        output_dir=args.output_dir
    )
    
    success = extractor.run_extraction()
    if not success:
        sys.exit(1)

if __name__ == "__main__":
    main()
