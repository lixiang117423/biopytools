"""
GFF格式转换主程序模块 | GFF Format Conversion Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import GFFConfig
from .utils import GFFLogger, get_file_stats
from .gff_parser import GFFParser
from .id_generator import IDGenerator
from .formatter import GFFFormatter

class GFFConverter:
    """GFF格式转换器主类 | Main GFF Format Converter Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = GFFConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        output_dir = Path(self.config.output_file).parent
        self.logger_manager = GFFLogger(output_dir)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化各个处理器 | Initialize processors
        self.parser = GFFParser(self.config, self.logger)
        self.id_generator = IDGenerator(self.config, self.logger)
        self.formatter = GFFFormatter(self.config, self.logger)
    
    def run_conversion(self):
        """运行完整的格式转换流程 | Run complete format conversion pipeline"""
        try:
            self.logger.info("🚀 开始GFF格式转换流程 | Starting GFF format conversion pipeline")
            self.logger.info(f"📄 输入文件 | Input file: {self.config.input_file}")
            self.logger.info(f"📝 输出文件 | Output file: {self.config.output_file}")
            self.logger.info(f"🧬 物种名称 | Species name: {self.config.species_name}")
            self.logger.info(f"🏷️  物种缩写 | Species prefix: {self.config.species_prefix}")
            
            # 获取输入文件统计 | Get input file statistics
            input_stats = get_file_stats(self.config.input_file)
            self.logger.info("📊 输入文件统计 | Input file statistics:")
            for key, value in input_stats.items():
                self.logger.info(f"   {key}: {value}")
            
            # 步骤1: 解析GFF文件 | Step 1: Parse GFF file
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤1: 解析GFF文件 | Step 1: Parse GFF file")
            self.logger.info(f"{'=' * 60}")
            
            features = self.parser.parse_file()
            header_lines = self.parser.header_lines
            
            if not features:
                raise RuntimeError("未找到有效的GFF特征 | No valid GFF features found")
            
            # 步骤2: 生成新的基因ID | Step 2: Generate new gene IDs
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤2: 生成新的基因ID | Step 2: Generate new gene IDs")
            self.logger.info(f"{'=' * 60}")
            
            id_mapping = self.id_generator.generate_new_ids(features)
            
            if not id_mapping:
                self.logger.warning("⚠️  未生成任何ID映射 | No ID mappings generated")
            
            # 步骤3: 格式化特征 | Step 3: Format features
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤3: 格式化GFF特征 | Step 3: Format GFF features")
            self.logger.info(f"{'=' * 60}")
            
            formatted_features = self.formatter.format_features(features, id_mapping)
            
            # 步骤4: 写入输出文件 | Step 4: Write output file
            self.logger.info(f"{'=' * 60}")
            self.logger.info("步骤4: 写入输出文件 | Step 4: Write output file")
            self.logger.info(f"{'=' * 60}")
            
            self.formatter.write_output(header_lines, formatted_features)
            
            # 获取输出文件统计 | Get output file statistics
            output_stats = get_file_stats(self.config.output_file)
            self.logger.info("📈 输出文件统计 | Output file statistics:")
            for key, value in output_stats.items():
                self.logger.info(f"   {key}: {value}")
            
            self.logger.info(f"{'=' * 60}")
            self.logger.info("🎉 GFF格式转换完成! | GFF format conversion completed!")
            self.logger.info(f"{'=' * 60}")
            self.logger.info(f"✅ 处理基因数量 | Processed genes: {len(id_mapping)}")
            self.logger.info(f"📁 输出文件保存在 | Output file saved at: {self.config.output_file}")
            
        except Exception as e:
            self.logger.error(f"❌ 格式转换流程意外终止 | Format conversion pipeline terminated unexpectedly: {e}")
            sys.exit(1)
    
    def show_sample_conversion(self, num_samples: int = 5):
        """显示转换示例 | Show conversion samples"""
        self.logger.info(f"📋 转换示例 | Conversion samples (showing first {num_samples}):")
        
        # 简化版解析来获取示例 | Simplified parsing for samples
        sample_count = 0
        with open(self.config.input_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) == 9 and fields[2] == 'gene':
                    # 提取原始基因ID | Extract original gene ID
                    attrs = {}
                    for pair in fields[8].split(';'):
                        if '=' in pair:
                            key, value = pair.split('=', 1)
                            attrs[key.strip()] = value.strip()
                    
                    original_id = attrs.get('ID', '').replace('gene-', '')
                    if original_id:
                        self.logger.info(f"   原始 | Original: {original_id}")
                        # 这里只是示例，实际的新ID需要通过完整流程生成
                        self.logger.info(f"   转换后 | Converted: [需要运行完整转换] | [Run full conversion needed]")
                        sample_count += 1
                        if sample_count >= num_samples:
                            break

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='GFF格式转换脚本 (模块化版本) | GFF Format Conversion Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -i input.gff -o output.gff -s OV53 -p Ov
  %(prog)s --input data.gff3 --output converted.gff --species-name "Species01" --species-prefix "Sp"
  %(prog)s -i annotation.gff -o result.gff -s OV53 -p Ov --start-num 20 --step 5
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input', required=True, 
                       help='🧬 输入GFF文件路径 | Input GFF file path')
    parser.add_argument('-o', '--output', required=True,
                       help='📝 输出GFF文件路径 | Output GFF file path')
    parser.add_argument('-s', '--species-name', required=True,
                       help='🦠 物种名称 (如: OV53) | Species name (e.g., OV53)')
    parser.add_argument('-p', '--species-prefix', required=True,
                       help='🏷️  物种缩写 (如: Ov) | Species prefix (e.g., Ov)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('--start-num', type=int, default=10,
                       help='🔢 起始编号 | Starting number for gene numbering')
    parser.add_argument('--step', type=int, default=10,
                       help='📏 编号步长 | Step size for gene numbering')
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🧵 线程数 | Number of threads')
    
    # 输出选项 | Output options
    parser.add_argument('-v', '--verbose', action='store_true',
                       help='📢 详细输出模式 | Verbose output mode')
    parser.add_argument('--keep-intermediate', action='store_true',
                       help='💾 保留中间文件 | Keep intermediate files')
    parser.add_argument('--show-sample', type=int, metavar='N',
                       help='🔍 显示N个转换示例然后退出 | Show N conversion samples and exit')
    
    args = parser.parse_args()
    
    # 创建转换器并运行 | Create converter and run
    converter = GFFConverter(
        input_file=args.input,
        output_file=args.output,
        species_name=args.species_name,
        species_prefix=args.species_prefix,
        start_number=args.start_num,
        step=args.step,
        threads=args.threads,
        verbose=args.verbose,
        keep_intermediate=args.keep_intermediate
    )
    
    # 如果只是显示示例 | If just showing samples
    if args.show_sample:
        converter.show_sample_conversion(args.show_sample)
        return
    
    # 运行转换 | Run conversion
    converter.run_conversion()

if __name__ == "__main__":
    main()
