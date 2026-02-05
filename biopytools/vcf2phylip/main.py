"""
VCF转换工具主程序模块|VCF Converter Main Module
"""

import argparse
import sys
from .config import ConverterConfig
from .utils import ConverterLogger
from .vcf_parser import VCFParser
from .matrix_writer import MatrixWriter
from .processor import VCFProcessor

class VCFConverter:
    """VCF转换器主类|Main VCF Converter Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = ConverterConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志|Initialize logging
        self.logger_manager = ConverterLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个组件|Initialize components
        self.parser = VCFParser(self.config, self.logger)
        self.writer = MatrixWriter(self.config, self.logger)
        self.processor = VCFProcessor(self.config, self.logger, self.parser, self.writer)
        
        self.logger.info("VCF转换工具已初始化|VCF Converter initialized")
    
    def run_conversion(self):
        """运行完整的VCF转换流程|Run complete VCF conversion pipeline"""
        try:
            self.logger.info("开始VCF转换流程|Starting VCF conversion pipeline")
            self.logger.info(f"输入文件|Input file: {self.config.input_file}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"使用线程数|Using threads: {self.config.threads}")
            
            # Step 1: 提取样本名称|Extract sample names
            sample_names = self.parser.extract_sample_names()
            if not sample_names:
                raise ValueError("VCF文件中未找到样本名称，文件可能损坏或缺少头部信息|Sample names not found in VCF")
            
            num_samples = len(sample_names)
            self.logger.info(f"VCF中的样本数量|Number of samples in VCF: {num_samples}")
            
            # 调整最小样本数|Adjust minimum sample count
            if self.config.min_samples_locus > num_samples:
                self.config.min_samples_locus = num_samples
                self.logger.info(f" 最小样本数已调整为|Min samples adjusted to: {num_samples}")
            
            # Step 2: 处理VCF基因型|Process VCF genotypes
            stats, temp_file, temp_bin_file = self.processor.process_vcf(sample_names, num_samples)
            
            if stats['snp_accepted'] == 0:
                raise ValueError("没有SNP通过筛选条件|No SNPs passed the filtering criteria")
            
            # Step 3: 初始化输出文件|Initialize output files
            self.writer.initialize_output_files(sample_names, stats['snp_accepted'], stats['snp_biallelic'])
            
            # Step 4: 写入序列矩阵|Write sequence matrices
            self.logger.info("开始写入输出矩阵|Starting to write output matrices")
            self.writer.write_sequences(sample_names, temp_file, temp_bin_file)
            
            # Step 5: 完成文件写入|Finalize file writing
            self.writer.finalize_output_files()
            
            # Step 6: 清理临时文件|Cleanup temporary files
            self.processor.cleanup_temp_files()
            
            self.logger.info("VCF转换完成！| VCF conversion completed successfully!")
            self.logger.info(f"结果保存在|Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"转换过程中发生错误|Error during conversion: {e}")
            self.processor.cleanup_temp_files()  # 确保清理临时文件|Ensure cleanup of temp files
            sys.exit(1)

def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='VCF转换工具|VCF Converter',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -i variants.vcf -o converted_results
        """
    )
    
    # 必需参数|Required arguments
    parser.add_argument('-i', '--input',
                       required=True,
                       help='输入VCF文件路径|Input VCF file path')

    # 输出参数|Output arguments
    parser.add_argument('-o', '--output',
                       default='./converted_output',
                       help='输出目录|Output directory')
    parser.add_argument('--output-prefix',
                       help='输出文件名前缀|Output filename prefix')

    # 转换参数|Conversion parameters
    parser.add_argument('-m', '--min-samples-locus',
                       type=int, default=4,
                       help='位点最少样本数|Minimum samples per locus')
    parser.add_argument('-g', '--outgroup',
                       default="",
                       help='外群样本名称|Outgroup sample name')

    # 输出格式控制|Output format control
    parser.add_argument('-p', '--phylip-disable',
                       action='store_true',
                       help='禁用PHYLIP输出|Disable PHYLIP output')
    parser.add_argument('-f', '--fasta',
                       action='store_true',
                       help='启用FASTA输出|Enable FASTA output')
    parser.add_argument('-n', '--nexus',
                       action='store_true',
                       help='启用NEXUS输出|Enable NEXUS output')
    parser.add_argument('-b', '--nexus-binary',
                       action='store_true',
                       help='启用二进制NEXUS输出|Enable binary NEXUS output')

    # 处理选项|Processing options
    parser.add_argument('-r', '--resolve-IUPAC',
                       action='store_true',
                       help='随机解析杂合子基因型|Resolve heterozygous genotypes')
    parser.add_argument('-w', '--write-used-sites',
                       action='store_true',
                       help='保存筛选通过的位点坐标|Save used sites coordinates')
    parser.add_argument('-t', '--threads',
                       type=int, default=88,
                       help='线程数|Number of threads')
    
    # 版本信息|Version info
    parser.add_argument('-v', '--version',
                       action='version',
                       version='%(prog)s v2.9.1 - 模块化版本|Modular Version')
    
    args = parser.parse_args()
    
    # 创建转换器并运行|Create converter and run
    converter = VCFConverter(
        input_file=args.input,
        output_dir=args.output,
        output_prefix=args.output_prefix,
        min_samples_locus=args.min_samples_locus,
        outgroup=args.outgroup,
        phylip_disable=args.phylip_disable,
        fasta=args.fasta,
        nexus=args.nexus,
        nexus_binary=args.nexus_binary,
        resolve_IUPAC=args.resolve_IUPAC,
        write_used_sites=args.write_used_sites,
        threads=args.threads
    )
    
    converter.run_conversion()

if __name__ == "__main__":
    main()
