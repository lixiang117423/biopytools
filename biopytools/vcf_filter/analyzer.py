"""
主分析器模块 | Main Analyzer Module
"""

import sys
import os
from .vcf_reader import VCFReader
from .python_filter import PythonFilter
from .plink_converter import PlinkConverter
from .utils import check_plink_availability

class VCFFilterAnalyzer:
    """VCF筛选分析器 | VCF Filter Analyzer"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        
        # 初始化各个组件 | Initialize components
        self.vcf_reader = VCFReader(config, logger)
        self.python_filter = PythonFilter(config, logger)
        self.plink_converter = PlinkConverter(config, logger)
    
    def run_analysis(self):
        """运行VCF筛选流程 | Run VCF filtering pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始VCF文件筛选 | Starting VCF file filtering")
            self.logger.info("=" * 60)
            self.logger.info(f"输入文件 | Input file: {self.config.vcf_file}")
            self.logger.info(f"输出文件 | Output file: {self.config.output_file}")
            
            # 步骤1: 验证输入文件 | Step 1: Validate input file
            self.logger.info("\n步骤1: 验证输入文件 | Step 1: Validating input file")
            total_variants = self.vcf_reader.count_variants()
            self.logger.info(f"输入文件包含 {total_variants} 个变异位点 | "
                           f"Input file contains {total_variants} variants")
            
            if total_variants == 0:
                raise ValueError("错误: 输入文件中没有找到有效的变异位点 | "
                               "Error: No valid variants found in input file")
            
            # 步骤2: 选择处理方式 | Step 2: Choose processing method
            if self.config.convert_format:
                # 检查PLINK可用性 | Check PLINK availability
                if not check_plink_availability(self.config.plink_path, self.logger):
                    self.logger.warning("PLINK不可用，转为使用Python筛选 | "
                                      "PLINK not available, switching to Python filtering")
                    self.config.convert_format = False
                else:
                    self.logger.info("\n步骤2: 使用PLINK进行筛选和格式转换 | "
                                   "Step 2: Using PLINK for filtering and format conversion")
                    output_file = self.plink_converter.convert_with_plink(self.config.vcf_file)
            
            if not self.config.convert_format:
                self.logger.info("\n步骤2: 使用Python进行筛选 | "
                               "Step 2: Using Python for filtering")
                output_file = self.python_filter.filter_vcf(self.vcf_reader)
            
            # 步骤3: 验证输出文件 | Step 3: Validate output file
            self.logger.info("\n步骤3: 验证输出文件 | Step 3: Validating output file")
            if not os.path.exists(output_file):
                raise FileNotFoundError(f"输出文件未生成 | Output file not generated: {output_file}")
            
            # 统计输出文件 | Count output file
            self.config.vcf_file = output_file  # 临时更改以统计输出
            output_reader = VCFReader(self.config, self.logger)
            final_variants = output_reader.count_variants()
            
            # 完成信息 | Completion information
            self.logger.info("\n" + "=" * 60)
            self.logger.info("VCF文件筛选完成！| VCF file filtering completed!")
            self.logger.info("=" * 60)
            self.logger.info(f"筛选结果: {total_variants} -> {final_variants} 个变异位点 | "
                           f"Filtering result: {total_variants} -> {final_variants} variants")
            self.logger.info(f"筛选率: {final_variants/total_variants*100:.2f}% | "
                           f"Retention rate: {final_variants/total_variants*100:.2f}%")
            self.logger.info(f"输出文件: {output_file} | Output file: {output_file}")
            
            return output_file
            
        except Exception as e:
            self.logger.error(f"筛选流程在执行过程中意外终止 | "
                            f"Filtering pipeline terminated unexpectedly: {e}")
            sys.exit(1)
