"""
GATK Joint Genotyping 主程序模块 | GATK Joint Genotyping Main Module
"""

import argparse
import sys
from .config import JointConfig
from .utils import JointLogger, CommandRunner, check_dependencies
from .file_processor import FileTypeDetector
from .joint_genotyping import JointGenotyper
from .filtering import VariantFilter
from .results import ResultsSummary

class GATKJointGenotyper:
    """GATK Joint Genotyping 主类 | Main GATK Joint Genotyping Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = JointConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = JointLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.file_detector = FileTypeDetector(self.config, self.logger)
        self.joint_genotyper = JointGenotyper(self.config, self.logger, self.cmd_runner)
        self.variant_filter = VariantFilter(self.config, self.logger, self.cmd_runner)
        self.results_summary = ResultsSummary(self.config, self.logger)
    
    def run_pipeline(self):
        """运行完整的Joint Genotyping流程 | Run complete joint genotyping pipeline"""
        try:
            self.logger.info("🚀 开始GATK Joint Genotyping流程 | Starting GATK Joint Genotyping pipeline")
            self.logger.info(f"📂 输入目录 | Input directory: {self.config.input_dir}")
            self.logger.info(f"🧬 参考基因组 | Reference genome: {self.config.reference}")
            self.logger.info(f"📁 输出目录 | Output directory: {self.config.output_dir}")
            
            # 检查依赖 | Check dependencies
            check_dependencies(self.config, self.logger)
            
            # 步骤1: 检测文件类型 | Step 1: Detect file type
            self.logger.info(f"{'=' * 80}")
            self.logger.info("步骤1: 检测输入文件类型 | Step 1: Detect input file type")
            self.logger.info(f"{'=' * 80}")
            
            vcf_files = self.file_detector.detect_file_type()
            sample_map_file = self.file_detector.create_sample_map(vcf_files)
            
            # 步骤2: GenomicsDB导入 | Step 2: GenomicsDB import
            self.logger.info(f"{'=' * 80}")
            self.logger.info("步骤2: GenomicsDB导入 | Step 2: GenomicsDB import")
            self.logger.info(f"{'=' * 80}")
            
            if not self.joint_genotyper.run_genomicsdb_import(sample_map_file):
                raise RuntimeError("❌ GenomicsDB导入失败 | GenomicsDB import failed")
            
            # 步骤3: 联合分型 | Step 3: Joint genotyping
            self.logger.info(f"{'=' * 80}")
            self.logger.info("步骤3: 联合分型 | Step 3: Joint genotyping")
            self.logger.info(f"{'=' * 80}")
            
            if not self.joint_genotyper.run_genotype_gvcfs():
                raise RuntimeError("❌ 联合分型失败 | Joint genotyping failed")
            
            # 步骤4: 提取和过滤SNP | Step 4: Extract and filter SNPs
            self.logger.info(f"{'=' * 80}")
            self.logger.info("步骤4: 提取和过滤SNP | Step 4: Extract and filter SNPs")
            self.logger.info(f"{'=' * 80}")
            
            if not self.variant_filter.select_snps():
                raise RuntimeError("❌ SNP提取失败 | SNP extraction failed")
            
            if not self.variant_filter.filter_snps():
                raise RuntimeError("❌ SNP过滤失败 | SNP filtering failed")
            
            # 步骤5: 提取和过滤INDEL | Step 5: Extract and filter INDELs
            self.logger.info(f"{'=' * 80}")
            self.logger.info("步骤5: 提取和过滤INDEL | Step 5: Extract and filter INDELs")
            self.logger.info(f"{'=' * 80}")
            
            if not self.variant_filter.select_indels():
                raise RuntimeError("❌ INDEL提取失败 | INDEL extraction failed")
            
            if not self.variant_filter.filter_indels():
                raise RuntimeError("❌ INDEL过滤失败 | INDEL filtering failed")
            
            # 步骤6: 合并SNP和INDEL | Step 6: Merge SNPs and INDELs
            self.logger.info(f"{'=' * 80}")
            self.logger.info("步骤6: 合并SNP和INDEL | Step 6: Merge SNPs and INDELs")
            self.logger.info(f"{'=' * 80}")
            
            if not self.variant_filter.merge_variants():
                raise RuntimeError("❌ 变异合并失败 | Variant merging failed")
            
            # 生成汇总报告 | Generate summary report
            self.results_summary.generate_summary()
            
            self.logger.info(f"{'=' * 80}")
            self.logger.info("🎉 Joint Genotyping分析完成! | Joint Genotyping analysis completed!")
            self.logger.info(f"{'=' * 80}")
            self.logger.info(f"📂 结果保存在 | Results saved in: {self.config.output_dir}")
            
        except Exception as e:
            self.logger.error(f"❌ 分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 GATK Joint Genotyping 脚本 (模块化版本) | GATK Joint Genotyping Script (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  %(prog)s -i gvcf_folder/ -r ref.fasta -o results/
  %(prog)s -i vcf_dir/ -r genome.fa -o output/ -t 88 -m 100g
  %(prog)s -i data/ -r ref.fa -o out/ -L chr1 --snp-qd 3.0
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-i', '--input-dir', required=True,
                       help='📂 输入目录 (包含VCF/GVCF文件) | Input directory (containing VCF/GVCF files)')
    parser.add_argument('-r', '--reference', required=True,
                       help='🧬 参考基因组文件 | Reference genome file (.fasta/.fa)')
    
    # 可选参数 | Optional arguments
    parser.add_argument('-o', '--output', default='./joint_genotyping_output',
                       help='📁 输出目录 | Output directory')
    
    # 计算资源 | Computing resources
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='⚙️ 线程数 | Number of threads')
    parser.add_argument('-m', '--memory', default='100g',
                       help='💾 Java内存设置 (如: 100g, 50g) | Java memory setting (e.g., 100g, 50g)')
    
    # 区间设置 | Interval settings
    parser.add_argument('-L', '--intervals',
                       help='📍 分析区间 (染色体或区间文件) | Analysis intervals (chromosome or interval file)')
    
    # SNP过滤参数 | SNP filtering parameters
    snp_group = parser.add_argument_group('SNP过滤参数 | SNP Filtering Parameters')
    snp_group.add_argument('--snp-qd', type=float, default=2.0,
                          help='SNP QD阈值 | SNP QD threshold')
    snp_group.add_argument('--snp-fs', type=float, default=60.0,
                          help='SNP FS阈值 | SNP FS threshold')
    snp_group.add_argument('--snp-mq', type=float, default=40.0,
                          help='SNP MQ阈值 | SNP MQ threshold')
    snp_group.add_argument('--snp-mqrs', type=float, default=-12.5,
                          help='SNP MQRankSum阈值 | SNP MQRankSum threshold')
    snp_group.add_argument('--snp-rprs', type=float, default=-8.0,
                          help='SNP ReadPosRankSum阈值 | SNP ReadPosRankSum threshold')
    snp_group.add_argument('--snp-sor', type=float, default=3.0,
                          help='SNP SOR阈值 | SNP SOR threshold')
    
    # INDEL过滤参数 | INDEL filtering parameters
    indel_group = parser.add_argument_group('INDEL过滤参数 | INDEL Filtering Parameters')
    indel_group.add_argument('--indel-qd', type=float, default=2.0,
                            help='INDEL QD阈值 | INDEL QD threshold')
    indel_group.add_argument('--indel-fs', type=float, default=200.0,
                            help='INDEL FS阈值 | INDEL FS threshold')
    indel_group.add_argument('--indel-rprs', type=float, default=-20.0,
                            help='INDEL ReadPosRankSum阈值 | INDEL ReadPosRankSum threshold')
    indel_group.add_argument('--indel-sor', type=float, default=10.0,
                            help='INDEL SOR阈值 | INDEL SOR threshold')
    
    # 工具路径 | Tool paths
    tool_group = parser.add_argument_group('工具路径 | Tool Paths')
    tool_group.add_argument('--gatk-path', default='gatk',
                           help='GATK软件路径 | GATK software path')
    tool_group.add_argument('--bcftools-path', default='bcftools',
                           help='BCFtools软件路径 | BCFtools software path')
    
    args = parser.parse_args()
    
    # 创建分析器并运行 | Create genotyper and run
    genotyper = GATKJointGenotyper(
        input_dir=args.input_dir,
        reference=args.reference,
        output_dir=args.output,
        threads=args.threads,
        memory=args.memory,
        intervals=args.intervals,
        snp_qd=args.snp_qd,
        snp_fs=args.snp_fs,
        snp_mq=args.snp_mq,
        snp_mqrs=args.snp_mqrs,
        snp_rprs=args.snp_rprs,
        snp_sor=args.snp_sor,
        indel_qd=args.indel_qd,
        indel_fs=args.indel_fs,
        indel_rprs=args.indel_rprs,
        indel_sor=args.indel_sor,
        gatk_path=args.gatk_path,
        bcftools_path=args.bcftools_path
    )
    
    genotyper.run_pipeline()

if __name__ == "__main__":
    main()
