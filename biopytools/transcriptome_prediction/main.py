"""
🧬 转录组预测分析主模块 | Transcriptome Prediction Analysis Main Module
"""

import argparse
import sys
from pathlib import Path
from .config import TranscriptomeConfig
from .utils import TranscriptomeLogger, CommandRunner, check_dependencies
from .alignment import HISAT2Aligner
from .assembly import StringTieAssembler, TrinityAssembler
from .annotation import PASAAnnotator
from .coding_prediction import TransDecoderPredictor

class TranscriptomeAnalyzer:
    """🚀 转录组预测分析器 | Transcriptome Prediction Analyzer"""
    
    def __init__(self, genome_file: str, rna_seq_files: list, output_dir: str = './transcriptome_output', **kwargs):
        """🔧 初始化转录组分析器 | Initialize transcriptome analyzer"""
        
        # 创建配置对象 | Create configuration object
        self.config = TranscriptomeConfig(
            genome_file=genome_file,
            rna_seq_files=rna_seq_files,
            output_dir=output_dir,
            **kwargs
        )
        
        # 验证配置 | Validate configuration
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        logger_manager = TranscriptomeLogger(self.config.output_path)
        self.logger = logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个组件 | Initialize components
        self.aligner = HISAT2Aligner(self.config, self.logger, self.cmd_runner)
        self.stringtie_assembler = StringTieAssembler(self.config, self.logger, self.cmd_runner)
        self.trinity_assembler = TrinityAssembler(self.config, self.logger, self.cmd_runner)
        self.annotator = PASAAnnotator(self.config, self.logger, self.cmd_runner)
        self.coding_predictor = TransDecoderPredictor(self.config, self.logger, self.cmd_runner)
        
        self.logger.info("🎉 转录组预测分析器初始化完成 | Transcriptome prediction analyzer initialized")
        self.logger.info(f"🧵 使用线程数 | Thread count: {self.config.threads}")

    def _check_step_completed(self, step_name: str) -> bool:
        """🔍 检查步骤是否已完成 | Check if step is completed"""
        
        if step_name == "hisat2_alignment":
            # 检查BAM文件是否存在
            bam_dir = self.config.output_path / "hisat2_alignment"
            if bam_dir.exists():
                bam_files = list(bam_dir.glob("*_sorted.bam"))
                if bam_files:
                    self.config.bam_files = [str(f) for f in bam_files]
                    self.logger.info(f"✅ 发现已存在的BAM文件 | Found existing BAM files: {len(bam_files)} files")
                    return True
        
        elif step_name == "stringtie_assembly":
            # 检查StringTie输出是否存在
            assembly_dir = self.config.output_path / "stringtie_assembly"
            merged_gtf = assembly_dir / "merged.gtf"
            if merged_gtf.exists():
                self.config.merged_gtf = str(merged_gtf)
                self.logger.info(f"✅ 发现已存在的StringTie合并GTF | Found existing StringTie merged GTF: {merged_gtf}")
                return True
            # 检查单样本GTF文件
            gtf_files = list(assembly_dir.glob("sample_*.gtf")) if assembly_dir.exists() else []
            if gtf_files:
                self.config.gtf_files = [str(f) for f in gtf_files]
                if len(gtf_files) == 1:
                    self.config.merged_gtf = str(gtf_files[0])
                self.logger.info(f"✅ 发现已存在的StringTie GTF文件 | Found existing StringTie GTF files: {len(gtf_files)} files")
                return True
        
        elif step_name == "trinity_assembly":
            # 检查Trinity组装输出是否存在
            assembly_dir = self.config.output_path / "trinity_assembly"
            
            # 多个可能的Trinity输出文件位置
            possible_trinity_files = [
                assembly_dir / "Trinity.fasta",
                self.config.output_path / f"{assembly_dir.name}.Trinity.fasta",
                assembly_dir / f"{assembly_dir.name}.Trinity.fasta",
                self.config.output_path / "Trinity.fasta"
            ]
            
            for trinity_fasta in possible_trinity_files:
                if trinity_fasta.exists():
                    self.config.trinity_assembly = str(trinity_fasta)
                    self.logger.info(f"✅ 发现已存在的Trinity组装 | Found existing Trinity assembly: {trinity_fasta}")
                    return True
        
        elif step_name == "pasa_annotation":
            # 检查PASA注释输出是否存在
            annotation_dir = self.config.output_path / "pasa_annotation"
            if annotation_dir.exists():
                possible_outputs = [
                    annotation_dir / "pasa_assemblies.gff3",
                    annotation_dir / "pasa_assemblies.gtf",
                    annotation_dir / f"{self.config.base_name}.pasa_assemblies.gff3",
                    annotation_dir / f"{self.config.base_name}.pasa_assemblies.gtf"
                ]
                
                for output_file in possible_outputs:
                    if output_file.exists():
                        self.config.pasa_annotation = str(output_file)
                        self.logger.info(f"✅ 发现已存在的PASA注释 | Found existing PASA annotation: {output_file}")
                        return True
        
        elif step_name == "transdecoder_prediction":
            # 检查TransDecoder输出是否存在
            prediction_dir = self.config.output_path / "transdecoder_prediction"
            if prediction_dir.exists():
                expected_files = [
                    prediction_dir / "transcripts.fasta.transdecoder.cds",
                    prediction_dir / "transcripts.fasta.transdecoder.pep",
                    prediction_dir / "transcripts.fasta.transdecoder.gff3"
                ]
                
                existing_files = [f for f in expected_files if f.exists()]
                if existing_files:
                    # 更完整的文件类型检查
                    outputs = {}
                    for f in existing_files:
                        if f.name.endswith('.cds'):
                            outputs['cds'] = str(f)
                        elif f.name.endswith('.pep'):
                            outputs['proteins'] = str(f)
                        elif f.name.endswith('.gff3'):
                            outputs['gff'] = str(f)
                    
                    self.config.transdecoder_outputs = outputs
                    self.logger.info(f"✅ 发现已存在的TransDecoder输出 | Found existing TransDecoder outputs: {len(existing_files)} files")
                    return True
        
        return False
    
    def _find_existing_trinity_assembly(self) -> bool:
        """🔍 查找已存在的Trinity组装文件 | Find existing Trinity assembly file"""
        # 可能的Trinity文件位置
        possible_locations = [
            # 在主输出目录中
            self.config.output_path / "trinity_assembly.Trinity.fasta",
            self.config.output_path / "Trinity.fasta",
            # 在trinity_assembly子目录中
            self.config.output_path / "trinity_assembly" / "Trinity.fasta",
            self.config.output_path / "trinity_assembly" / "trinity_assembly.Trinity.fasta",
        ]
        
        for trinity_file in possible_locations:
            if trinity_file.exists():
                self.config.trinity_assembly = str(trinity_file)
                self.logger.info(f"✅ 找到已存在的Trinity组装文件 | Found existing Trinity assembly: {trinity_file}")
                return True
        
        return False
    
    def run_analysis(self, resume=True):
        """🚀 运行完整的转录组预测分析 | Run complete transcriptome prediction analysis
        
        Args:
            resume (bool): 是否跳过已完成的步骤 | Whether to skip completed steps
        """
        self.logger.info("🎬 开始转录组预测分析 | Starting transcriptome prediction analysis")
        if resume:
            self.logger.info("🔄 启用断点续传模式 | Resume mode enabled - will skip completed steps")
        
        try:
            # 检查依赖软件 | Check dependencies
            check_dependencies(self.config, self.logger)
            
            # Step 1: HISAT2比对 | HISAT2 alignment
            self.logger.info("🎯 === 步骤1: HISAT2比对 | Step 1: HISAT2 Alignment ===")
            if resume and self._check_step_completed("hisat2_alignment"):
                self.logger.info("⭐ 跳过HISAT2比对步骤 - 已完成 | Skipping HISAT2 alignment - already completed")
            else:
                if not self.aligner.align_reads():
                    raise RuntimeError("❌ HISAT2比对失败 | HISAT2 alignment failed")
            
            # Step 2: StringTie转录本重构 | StringTie transcript reconstruction
            self.logger.info("🧩 === 步骤2: StringTie转录本重构 | Step 2: StringTie Transcript Reconstruction ===")
            if resume and self._check_step_completed("stringtie_assembly"):
                self.logger.info("⭐ 跳过StringTie转录本重构步骤 - 已完成 | Skipping StringTie reconstruction - already completed")
            else:
                if not self.stringtie_assembler.reconstruct_transcripts():
                    raise RuntimeError("❌ StringTie转录本重构失败 | StringTie transcript reconstruction failed")
            
            # Step 3: Trinity de novo组装 | Trinity de novo assembly
            if self.config.skip_trinity:
                self.logger.info("🔗 === 步骤3: Trinity de novo组装 - 跳过 | Step 3: Trinity de novo Assembly - SKIPPED ===")
                # 尝试查找已存在的Trinity组装文件
                if not self._find_existing_trinity_assembly():
                    self.logger.warning("⚠️ 跳过了Trinity组装但未找到已存在的Trinity文件，PASA步骤可能会失败")
                    self.logger.warning("⚠️ Skipped Trinity assembly but no existing Trinity file found, PASA step may fail")
            else:
                self.logger.info("🔗 === 步骤3: Trinity de novo组装 | Step 3: Trinity de novo Assembly ===")
                if resume and self._check_step_completed("trinity_assembly"):
                    self.logger.info("⭐ 跳过Trinity de novo组装步骤 - 已完成 | Skipping Trinity assembly - already completed")
                else:
                    if not self.trinity_assembler.de_novo_assembly():
                        raise RuntimeError("❌ Trinity de novo组装失败 | Trinity de novo assembly failed")
            
            # Step 4: PASA基因结构注释 | PASA gene structure annotation
            self.logger.info("📍 === 步骤4: PASA基因结构注释 | Step 4: PASA Gene Structure Annotation ===")
            
            # 检查是否有Trinity组装文件用于PASA
            if not hasattr(self.config, 'trinity_assembly') or not self.config.trinity_assembly:
                if not self._find_existing_trinity_assembly():
                    self.logger.error("❌ PASA步骤需要Trinity组装文件，但未找到 | PASA step requires Trinity assembly file, but not found")
                    if self.config.skip_trinity:
                        self.logger.error("💡 提示：如果跳过了Trinity步骤，请确保之前已经运行过Trinity并且输出文件存在")
                        self.logger.error("💡 Hint: If Trinity step was skipped, ensure Trinity was run previously and output file exists")
                    raise RuntimeError("❌ 缺少Trinity组装文件，无法运行PASA | Missing Trinity assembly file, cannot run PASA")
            
            if resume and self._check_step_completed("pasa_annotation"):
                self.logger.info("⭐ 跳过PASA基因结构注释步骤 - 已完成 | Skipping PASA annotation - already completed")
            else:
                if not self.annotator.map_transcripts_to_genome():
                    raise RuntimeError("❌ PASA基因结构注释失败 | PASA gene structure annotation failed")
            
            # Step 5: TransDecoder编码区预测 | TransDecoder coding region prediction
            self.logger.info("🔍 === 步骤5: TransDecoder编码区预测 | Step 5: TransDecoder Coding Region Prediction ===")
            if resume and self._check_step_completed("transdecoder_prediction"):
                self.logger.info("⭐ 跳过TransDecoder编码区预测步骤 - 已完成 | Skipping TransDecoder prediction - already completed")
            else:
                if not self.coding_predictor.identify_coding_regions():
                    raise RuntimeError("❌ TransDecoder编码区预测失败 | TransDecoder coding region prediction failed")
            
            # 生成分析报告 | Generate analysis report
            self._generate_report()
            
            self.logger.info("🎉🎊 转录组预测分析完成！ | Transcriptome prediction analysis completed! 🎊🎉")
            
        except Exception as e:
            self.logger.error(f"💥 分析过程中发生错误 | Error occurred during analysis: {str(e)}")
            sys.exit(1)


    def run_specific_step(self, step_name: str):
        """🎯 运行特定步骤 | Run specific step"""
        self.logger.info(f"🎯 运行特定步骤 | Running specific step: {step_name}")
        
        try:
            if step_name == "alignment":
                if not self.aligner.align_reads():
                    raise RuntimeError("❌ HISAT2比对失败 | HISAT2 alignment failed")
            
            elif step_name == "stringtie":
                if not self.stringtie_assembler.reconstruct_transcripts():
                    raise RuntimeError("❌ StringTie转录本重构失败 | StringTie transcript reconstruction failed")
            
            elif step_name == "trinity":
                if not self.trinity_assembler.de_novo_assembly():
                    raise RuntimeError("❌ Trinity de novo组装失败 | Trinity de novo assembly failed")
            
            elif step_name == "pasa":
                if not self.annotator.map_transcripts_to_genome():
                    raise RuntimeError("❌ PASA基因结构注释失败 | PASA gene structure annotation failed")
            
            elif step_name == "transdecoder":
                if not self.coding_predictor.identify_coding_regions():
                    raise RuntimeError("❌ TransDecoder编码区预测失败 | TransDecoder coding region prediction failed")
            
            else:
                raise ValueError(f"❌ 未知步骤 | Unknown step: {step_name}")
            
            self.logger.info(f"✅ 步骤 {step_name} 完成 | Step {step_name} completed")
            
        except Exception as e:
            self.logger.error(f"💥 步骤 {step_name} 失败 | Step {step_name} failed: {str(e)}")
            sys.exit(1)
    
    def _generate_report(self):
        """📊 生成分析报告 | Generate analysis report"""
        report_file = self.config.output_path / "transcriptome_analysis_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("🧬 转录组预测分析报告 | Transcriptome Prediction Analysis Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write("📍 输入文件 | Input Files:\n")
            f.write(f"  🧬 基因组文件 | Genome file: {self.config.genome_file}\n")
            f.write(f"  📊 RNA-seq文件 | RNA-seq files: {len(self.config.rna_seq_files)} 个文件 | files\n")
            for i, rna_file in enumerate(self.config.rna_seq_files, 1):
                f.write(f"    {i}. {rna_file}\n")
            
            if self.config.samples_file:
                f.write(f"  📋 样本文件 | Samples file: {self.config.samples_file}\n")
            f.write("\n")
            
            f.write("⚙️ 分析参数 | Analysis Parameters:\n")
            f.write(f"  🧵 线程数 | Threads: {self.config.threads}\n")
            f.write(f"  🎯 HISAT2最小内含子长度 | HISAT2 min intron length: {self.config.hisat2_min_intron}\n")
            f.write(f"  🎯 HISAT2最大内含子长度 | HISAT2 max intron length: {self.config.hisat2_max_intron}\n")
            f.write(f"  🧩 StringTie最小转录本长度 | StringTie min transcript length: {self.config.stringtie_min_length}\n")
            f.write(f"  🔗 Trinity最小contig长度 | Trinity min contig length: {self.config.trinity_min_contig_length}\n")
            f.write(f"  🔍 TransDecoder最小蛋白质长度 | TransDecoder min protein length: {self.config.transdecoder_min_protein_len}\n")
            f.write("\n")
            
            f.write("📄 输出文件 | Output Files:\n")
            if hasattr(self.config, 'bam_files'):
                f.write(f"  🎯 比对文件 | Alignment files: {len(self.config.bam_files)} 个BAM文件 | BAM files\n")
                for i, bam_file in enumerate(self.config.bam_files, 1):
                    f.write(f"    {i}. {Path(bam_file).name}\n")
            
            if hasattr(self.config, 'merged_gtf'):
                f.write(f"  🧩 StringTie转录本 | StringTie transcripts: {Path(self.config.merged_gtf).name}\n")
            
            if hasattr(self.config, 'trinity_assembly'):
                f.write(f"  🔗 Trinity组装 | Trinity assembly: {Path(self.config.trinity_assembly).name}\n")
            
            if hasattr(self.config, 'pasa_annotation'):
                f.write(f"  📍 PASA注释 | PASA annotation: {Path(self.config.pasa_annotation).name}\n")
            
            if hasattr(self.config, 'transdecoder_outputs'):
                f.write("  🔍 TransDecoder输出 | TransDecoder outputs:\n")
                for key, file_path in self.config.transdecoder_outputs.items():
                    if file_path:
                        f.write(f"    {key}: {Path(file_path).name}\n")
            
            f.write("\n")
            f.write("🎉 分析完成 | Analysis completed\n")
        
        self.logger.info(f"📊 分析报告已生成 | Analysis report generated: {report_file}")

# def main():
#     """🚀 主函数 | Main function"""
#     parser = argparse.ArgumentParser(
#         description="🧬 转录组预测分析工具 | Transcriptome-based prediction analysis tool",
#         formatter_class=argparse.RawDescriptionHelpFormatter,
#         epilog="""
# 🌟 示例 | Examples:
#   基本分析 | Basic analysis:
#     python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results
  
#   断点续传 | Resume analysis:
#     python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results --resume
  
#   重新开始分析 | Restart analysis:
#     python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results --no-resume
  
#   运行特定步骤 | Run specific step:
#     python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results --step alignment
#         """
#     )
    
#     # 必需参数 | Required parameters
#     parser.add_argument('-g', '--genome', required=True,
#                        help='🧬 基因组FASTA文件 | Genome FASTA file')
    
#     # RNA-seq文件或样本文件 | RNA-seq files or samples file
#     input_group = parser.add_mutually_exclusive_group(required=True)
#     input_group.add_argument('-r', '--rna-seq', nargs='+',
#                             help='📊 RNA-seq FASTQ文件 (支持单端或配对末端) | RNA-seq FASTQ files (supports single-end or paired-end)')
#     input_group.add_argument('--samples-file',
#                             help='📋 Trinity样本文件格式 | Trinity samples file format')
    
#     parser.add_argument('-o', '--output', required=True,
#                        help='📁 输出目录 | Output directory')

def main():
    """🚀 主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description="🧬 转录组预测分析工具 | Transcriptome-based prediction analysis tool",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
🌟 示例 | Examples:
  基本分析 | Basic analysis:
    python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results
  
  跳过Trinity步骤 | Skip Trinity step:
    python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results --skip-trinity
  
  断点续传 | Resume analysis:
    python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results --resume
  
  重新开始分析 | Restart analysis:
    python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results --no-resume
  
  运行特定步骤 | Run specific step:
    python transcriptome_prediction/main.py -g genome.fa -r sample_R1.fq sample_R2.fq -o results --step alignment
        """
    )
    
    # 必需参数 | Required parameters
    parser.add_argument('-g', '--genome', required=True,
                       help='🧬 基因组FASTA文件 | Genome FASTA file')
    
    # RNA-seq文件或样本文件 | RNA-seq files or samples file
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument('-r', '--rna-seq', nargs='+',
                            help='📊 RNA-seq FASTQ文件 (支持单端或配对末端) | RNA-seq FASTQ files (supports single-end or paired-end)')
    input_group.add_argument('--samples-file',
                            help='📋 Trinity样本文件格式 | Trinity samples file format')
    
    parser.add_argument('-o', '--output', required=True,
                       help='📁 输出目录 | Output directory')
    
    # 流程控制参数 | Workflow control parameters
    # workflow_group = parser.add_argument_group("流程控制 | Workflow Control")
    # # workflow_group.add_argument('--resume', action='store_true', default=True,
    # #                            help='🔄 启用断点续传，跳过已完成的步骤（默认启用） | Enable resume mode, skip completed steps (enabled by default)')
    # workflow_group.add_argument('--no-resume', action='store_true',
    #                            help='🚫 禁用断点续传，重新运行所有步骤 | Disable resume mode, rerun all steps')
    # workflow_group.add_argument('--skip-trinity', action='store_true',
    #                            help='🔗 跳过Trinity de novo组装步骤 | Skip Trinity de novo assembly step')
    # workflow_group.add_argument('--step', choices=['alignment', 'stringtie', 'trinity', 'pasa', 'transdecoder'],
    #                            help='🎯 只运行特定步骤 | Run only specific step')
    
    # 流程控制参数 | Workflow control parameters
    workflow_group = parser.add_argument_group("流程控制 | Workflow Control")
    workflow_group.add_argument('--resume', action='store_true', default=True,
                               help='🔄 启用断点续传，跳过已完成的步骤（默认启用） | Enable resume mode, skip completed steps (enabled by default)')
    workflow_group.add_argument('--no-resume', action='store_true',
                               help='🚫 禁用断点续传，重新运行所有步骤 | Disable resume mode, rerun all steps')
    workflow_group.add_argument('--step', choices=['alignment', 'stringtie', 'trinity', 'pasa', 'transdecoder'],
                               help='🎯 只运行特定步骤 | Run only specific step')
    
    # 可选参数 | Optional parameters
    parser.add_argument('-t', '--threads', type=int, default=88,
                       help='🧵 线程数 (默认: 88) | Number of threads (default: 88)')
    
    # HISAT2参数 | HISAT2 parameters
    parser.add_argument('--hisat2-min-intron', type=int, default=20,
                       help='🎯 HISAT2最小内含子长度 (默认: 20) | HISAT2 minimum intron length (default: 20)')
    parser.add_argument('--hisat2-max-intron', type=int, default=500000,
                       help='🎯 HISAT2最大内含子长度 (默认: 500000) | HISAT2 maximum intron length (default: 500000)')
    parser.add_argument('--hisat2-novel-splicesite', action='store_true',
                       help='🎯 HISAT2输出新剪接位点 | HISAT2 output novel splice sites')
    parser.add_argument('--no-dta', action='store_true',
                       help='🎯 禁用HISAT2的--dta选项 | Disable HISAT2 --dta option')
    
    # StringTie参数 | StringTie parameters
    parser.add_argument('--stringtie-min-length', type=int, default=200,
                       help='🧩 StringTie最小转录本长度 (默认: 200) | StringTie minimum transcript length (default: 200)')
    parser.add_argument('--stringtie-min-coverage', type=float, default=1.0,
                       help='🧩 StringTie最小覆盖度 (默认: 1.0) | StringTie minimum coverage (default: 1.0)')
    parser.add_argument('--stringtie-min-fpkm', type=float, default=1.0,
                       help='🧩 StringTie最小FPKM (默认: 1.0) | StringTie minimum FPKM (default: 1.0)')
    parser.add_argument('--stringtie-min-iso', type=float, default=0.01,
                       help='🧩 StringTie最小isoform (默认: 1.0) | StringTie isoform fraction (default: 0.01)')
    parser.add_argument('--stringtie-conservative', action='store_true',
                       help='🧩 StringTie保守组装模式 | StringTie conservative assembly mode')
    
    # Trinity参数 | Trinity parameters
    parser.add_argument('--trinity-min-contig-length', type=int, default=200,
                       help='🔗 Trinity最小contig长度 (默认: 200) | Trinity minimum contig length (default: 200)')
    parser.add_argument('--trinity-max-memory', default='200G',
                       help='🔗 Trinity最大内存使用 (默认: 200G) | Trinity maximum memory usage (default: 200G)')
    parser.add_argument('--trinity-cpu', type=int, default=88,
                       help='🔗 Trinity CPU数量 (默认: 88) | Trinity CPU count (default: 88)')
    parser.add_argument('--trinity-ss-lib-type', choices=['FR', 'RF', 'F', 'R'],
                       help='🔗 Trinity链特异性类型 | Trinity strand-specific library type')
    
    # PASA参数 | PASA parameters
    parser.add_argument('--pasa-max-intron-length', type=int, default=100000,
                       help='📍 PASA最大内含子长度 (默认: 100000) | PASA maximum intron length (default: 100000)')
    parser.add_argument('--pasa-min-percent-aligned', type=int, default=90,
                       help='📍 PASA最小比对百分比 (默认: 90) | PASA minimum percent aligned (default: 90)')
    parser.add_argument('--pasa-min-avg-per-id', type=int, default=95,
                       help='📍 PASA最小平均身份百分比 (默认: 95) | PASA minimum average percent identity (default: 95)')
    parser.add_argument('--pasa-aligners', default='gmap,blat',
                       help='📍 PASA比对器 (默认: gmap,blat) | PASA aligners (default: gmap,blat)')
    parser.add_argument('--pasa-cpu', type=int, default=88,
                       help='📍 PASA CPU数量 (默认: 88) | PASA CPU count (default: 88)')
    
    # TransDecoder参数 | TransDecoder parameters
    parser.add_argument('--transdecoder-min-protein-len', type=int, default=100,
                       help='🔍 TransDecoder最小蛋白质长度 (默认: 100) | TransDecoder minimum protein length (default: 100)')
    parser.add_argument('--transdecoder-genetic-code', default='universal',
                       help='🔍 TransDecoder遗传密码 (默认: universal) | TransDecoder genetic code (default: universal)')
    parser.add_argument('--transdecoder-complete-orfs-only', action='store_true',
                       help='🔍 TransDecoder只保留完整ORF | TransDecoder keep only complete ORFs')
    
    # 工具路径 | Tool paths
    parser.add_argument('--hisat2-path', default='hisat2',
                       help='🛠️ HISAT2可执行文件路径 (默认: hisat2) | HISAT2 executable path (default: hisat2)')
    parser.add_argument('--stringtie-path', default='stringtie',
                       help='🛠️ StringTie可执行文件路径 (默认: stringtie) | StringTie executable path (default: stringtie)')
    parser.add_argument('--trinity-path', default='Trinity',
                       help='🛠️ Trinity可执行文件路径 (默认: Trinity) | Trinity executable path (default: Trinity)')
    parser.add_argument('--pasa-path', default='Launch_PASA_pipeline.pl',
                       help='🛠️ PASA可执行文件路径 (默认: Launch_PASA_pipeline.pl) | PASA executable path (default: Launch_PASA_pipeline.pl)')
    parser.add_argument('--transdecoder-longorfs-path', default='TransDecoder.LongOrfs',
                       help='🛠️ TransDecoder.LongOrfs可执行文件路径 | TransDecoder.LongOrfs executable path')
    parser.add_argument('--transdecoder-predict-path', default='TransDecoder.Predict',
                       help='🛠️ TransDecoder.Predict可执行文件路径 | TransDecoder.Predict executable path')
    parser.add_argument('--samtools-path', default='samtools',
                       help='🛠️ SAMtools可执行文件路径 (默认: samtools) | SAMtools executable path (default: samtools)')
    
    args = parser.parse_args()
    
    # 处理输入文件 | Process input files
    if args.rna_seq:
        rna_seq_files = args.rna_seq
        samples_file = None
    else:
        rna_seq_files = []  # 样本文件模式下留空 | Empty for samples file mode
        samples_file = args.samples_file
    
    # 处理流程控制参数 | Handle workflow control parameters
    # resume_mode = args.resume and not args.no_resume
    # 处理流程控制参数 | Handle workflow control parameters  
    resume_mode = not getattr(args, 'no_resume', False)  # 默认启用断点续传
    skip_trinity = getattr(args, 'skip_trinity', False)  # 默认不跳过Trinity
    step = getattr(args, 'step', None)  # 默认运行所有步骤
        
    # 创建分析器并运行 | Create analyzer and run
    analyzer = TranscriptomeAnalyzer(
        genome_file=args.genome,
        rna_seq_files=rna_seq_files,
        output_dir=args.output,
        samples_file=samples_file,
        # skip_trinity=args.skip_trinity,
        skip_trinity=skip_trinity,
        threads=args.threads,
        hisat2_min_intron=args.hisat2_min_intron,
        hisat2_max_intron=args.hisat2_max_intron,
        hisat2_novel_splicesite=args.hisat2_novel_splicesite,
        hisat2_dta=not args.no_dta,
        stringtie_min_length=args.stringtie_min_length,
        stringtie_min_coverage=args.stringtie_min_coverage,
        stringtie_min_fpkm=args.stringtie_min_fpkm,
        stringtie_min_iso=args.stringtie_min_iso,
        stringtie_conservative=args.stringtie_conservative,
        trinity_min_contig_length=args.trinity_min_contig_length,
        trinity_max_memory=args.trinity_max_memory,
        trinity_cpu=args.trinity_cpu,
        trinity_ss_lib_type=args.trinity_ss_lib_type,
        pasa_max_intron_length=args.pasa_max_intron_length,
        pasa_min_percent_aligned=args.pasa_min_percent_aligned,
        pasa_min_avg_per_id=args.pasa_min_avg_per_id,
        pasa_aligners=args.pasa_aligners,
        pasa_cpu=args.pasa_cpu,
        transdecoder_min_protein_len=args.transdecoder_min_protein_len,
        transdecoder_genetic_code=args.transdecoder_genetic_code,
        transdecoder_complete_orfs_only=args.transdecoder_complete_orfs_only,
        hisat2_path=args.hisat2_path,
        stringtie_path=args.stringtie_path,
        trinity_path=args.trinity_path,
        pasa_path=args.pasa_path,
        transdecoder_longorfs_path=args.transdecoder_longorfs_path,
        transdecoder_predict_path=args.transdecoder_predict_path,
        samtools_path=args.samtools_path
    )
    
    step = getattr(args, 'step', None)

    if step:
        # 运行特定步骤 | Run specific step
        if step == 'trinity' and skip_trinity:
            analyzer.logger.error("❌ 无法运行Trinity步骤：指定了--skip-trinity参数")
            sys.exit(1)
        analyzer.run_specific_step(step)
    else:
        # 运行完整分析 | Run full analysis
        analyzer.run_analysis(resume=resume_mode)

if __name__ == "__main__":
    main()