"""
基因组组装主程序模块 | Genome Assembly Main Module 🧬
"""

import argparse
import sys
import time
import traceback
from pathlib import Path
from typing import Dict, List, Optional, Any

from .config import AssemblyConfig
from .utils import (
    AssemblyLogger, CommandRunner, check_dependencies, 
    create_directory_structure, format_time
)
from .data_processing import DataQualityController
from .assembly import HifiasmAssembler
from .hic_processing import HiCProcessor
from .quality_control import AssemblyQualityAssessor
from .results import ResultsManager

class GenomeAssembler:
    """基因组组装器主类 | Main Genome Assembler Class"""
    
    def __init__(self, **kwargs):
        """初始化基因组组装器 | Initialize genome assembler"""
        # 初始化配置
        self.config = AssemblyConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志系统
        self.logger_manager = AssemblyLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 记录开始信息
        self.logger.info("🧬 基因组组装流程启动 | Genome assembly pipeline started")
        self.logger.info(f"📋 项目名称: {self.config.project_name}")
        self.logger.info(f"🗂️ 输出目录: {self.config.output_dir}")
        self.logger.info(f"🔍 检测到的数据类型: {', '.join(self.config.detected_data_types)}")
        self.logger.info(f"🎯 组装策略: {self.config.assembly_strategy}")
        self.logger.info(f"🔗 Hi-C策略: {self.config.hic_strategy}")
        
        # 初始化命令执行器
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各功能模块
        self.data_qc = DataQualityController(self.config, self.logger, self.cmd_runner)
        self.assembler = HifiasmAssembler(self.config, self.logger, self.cmd_runner)
        self.hic_processor = HiCProcessor(self.config, self.logger, self.cmd_runner)
        self.quality_assessor = AssemblyQualityAssessor(self.config, self.logger, self.cmd_runner)
        self.results_manager = ResultsManager(self.config, self.logger, self.cmd_runner)
        
        # 记录运行时间
        self.start_time = time.time()
        
    def run_assembly(self) -> Dict[str, Any]:
        """运行完整的基因组组装流程 | Run complete genome assembly pipeline"""
        try:
            self.logger.info("🚀 开始基因组组装完整流程 | Starting complete genome assembly pipeline")
            
            # 创建目录结构
            create_directory_structure(self.config.output_path, self.logger)
            
            # Phase 1: 环境检查和依赖验证
            self._phase1_environment_check()
            
            # Phase 2: 数据质量控制
            qc_results = self._phase2_data_quality_control()
            
            # Phase 3: Hifiasm基因组组装
            assembly_files = self._phase3_hifiasm_assembly()
            
            # Phase 4: Hi-C处理 (如果有Hi-C数据)
            scaffold_results = self._phase4_hic_processing(assembly_files)
            
            # Phase 5: 质量评估
            quality_results = self._phase5_quality_assessment(assembly_files, scaffold_results)
            
            # Phase 6: 结果整理和报告
            final_results = self._phase6_results_organization(assembly_files, scaffold_results, quality_results)
            
            # 流程完成
            self._pipeline_completion(final_results)
            
            return final_results
            
        except Exception as e:
            self._handle_pipeline_error(e)
            raise
        
    def _phase1_environment_check(self):
        """Phase 1: 环境检查 | Environment check"""
        self.logger.info("🔧 Phase 1: 环境检查和依赖验证 | Environment check and dependency verification")
        
        try:
            # 如果跳过依赖检查
            if self.config.skip_dependency_check:
                self.logger.warning("⚠️ 跳过依赖软件检查 | Skipping dependency check")
            else:
                # 检查依赖软件
                check_dependencies(self.config, self.logger)
            
            # 记录系统信息
            self._log_system_info()
            
            # 验证输入文件完整性
            self._verify_input_files()
            
            self.logger.info("✅ Phase 1 完成: 环境检查通过 | Phase 1 completed: Environment check passed")
            
        except Exception as e:
            self.logger.error(f"❌ Phase 1 失败: 环境检查失败 | Phase 1 failed: {e}")
            raise
    
    def _phase2_data_quality_control(self) -> Dict[str, Any]:
        """Phase 2: 数据质量控制 | Data quality control"""
        self.logger.info("🔬 Phase 2: 数据质量控制 | Data quality control")
        
        try:
            qc_results = self.data_qc.run_quality_control()
            
            # 检查是否通过质量控制
            self._evaluate_qc_results(qc_results)
            
            self.logger.info("✅ Phase 2 完成: 数据质量控制完成 | Phase 2 completed: Data QC completed")
            return qc_results
            
        except Exception as e:
            self.logger.error(f"❌ Phase 2 失败: 数据质量控制失败 | Phase 2 failed: {e}")
            raise
    
    def _phase3_hifiasm_assembly(self) -> Dict[str, List[str]]:
        """Phase 3: Hifiasm组装 | Hifiasm assembly"""
        self.logger.info("🧬 Phase 3: Hifiasm基因组组装 | Hifiasm genome assembly")
        
        try:
            assembly_files = self.assembler.run_assembly()
            
            # 验证组装结果
            self._validate_assembly_outputs(assembly_files)
            
            self.logger.info("✅ Phase 3 完成: Hifiasm组装完成 | Phase 3 completed: Hifiasm assembly completed")
            return assembly_files
            
        except Exception as e:
            self.logger.error(f"❌ Phase 3 失败: Hifiasm组装失败 | Phase 3 failed: {e}")
            raise
    
    def _phase4_hic_processing(self, assembly_files: Dict[str, List[str]]) -> Optional[Dict[str, Dict[str, str]]]:
        """Phase 4: Hi-C处理 | Hi-C processing"""
        if 'Hi-C' not in self.config.detected_data_types:
            self.logger.info("⏭️ Phase 4: 跳过Hi-C处理 (无Hi-C数据) | Phase 4: Skipping Hi-C processing (no Hi-C data)")
            return None
        
        self.logger.info("🔗 Phase 4: Hi-C染色体挂载处理 | Hi-C chromosome scaffolding processing")
        
        try:
            scaffold_results = self.hic_processor.process_all_assemblies(assembly_files)
            
            # 验证挂载结果
            self._validate_scaffold_outputs(scaffold_results)
            
            self.logger.info("✅ Phase 4 完成: Hi-C处理完成 | Phase 4 completed: Hi-C processing completed")
            return scaffold_results
            
        except Exception as e:
            self.logger.error(f"❌ Phase 4 失败: Hi-C处理失败 | Phase 4 failed: {e}")
            # Hi-C处理失败不应该终止整个流程
            self.logger.warning("⚠️ 继续进行质量评估，但不包含挂载结果 | Continuing with quality assessment without scaffolding results")
            return None
    
    def _phase5_quality_assessment(self, assembly_files: Dict[str, List[str]], 
                                 scaffold_results: Optional[Dict[str, Dict[str, str]]]) -> Dict[str, Dict[str, Any]]:
        """Phase 5: 质量评估 | Quality assessment"""
        self.logger.info("🏆 Phase 5: 组装质量评估 | Assembly quality assessment")
        
        try:
            quality_results = self.quality_assessor.assess_all_results(assembly_files, scaffold_results)
            
            # 记录质量评估摘要
            self._log_quality_summary(quality_results)
            
            self.logger.info("✅ Phase 5 完成: 质量评估完成 | Phase 5 completed: Quality assessment completed")
            return quality_results
            
        except Exception as e:
            self.logger.error(f"❌ Phase 5 失败: 质量评估失败 | Phase 5 failed: {e}")
            # 质量评估失败不应终止流程，返回空结果
            self.logger.warning("⚠️ 继续进行结果整理 | Continuing with results organization")
            return {}
    
    def _phase6_results_organization(self, assembly_files: Dict[str, List[str]], 
                                   scaffold_results: Optional[Dict[str, Dict[str, str]]],
                                   quality_results: Dict[str, Dict[str, Any]]) -> Dict[str, Any]:
        """Phase 6: 结果整理 | Results organization"""
        self.logger.info("📋 Phase 6: 结果整理和报告生成 | Results organization and report generation")
        
        try:
            final_results = self.results_manager.organize_final_results(
                assembly_files, scaffold_results, quality_results
            )
            
            self.logger.info("✅ Phase 6 完成: 结果整理完成 | Phase 6 completed: Results organization completed")
            return final_results
            
        except Exception as e:
            self.logger.error(f"❌ Phase 6 失败: 结果整理失败 | Phase 6 failed: {e}")
            raise
    
    def _pipeline_completion(self, final_results: Dict[str, Any]):
        """流程完成处理 | Pipeline completion handling"""
        end_time = time.time()
        total_time = end_time - self.start_time
        
        self.logger.info("🎉 基因组组装流程完成！| Genome assembly pipeline completed!")
        self.logger.info(f"⏱️ 总耗时 | Total time: {format_time(total_time)}")
        self.logger.info(f"📁 结果目录 | Results directory: {self.config.output_dir}")
        
        # 输出主要结果摘要
        self._print_results_summary(final_results)
        
        # 输出下一步建议
        self._print_next_steps(final_results)
    
    def _handle_pipeline_error(self, error: Exception):
        """处理流程错误 | Handle pipeline error"""
        self.logger.error("💥 基因组组装流程意外终止！| Genome assembly pipeline terminated unexpectedly!")
        self.logger.error(f"❌ 错误信息 | Error message: {str(error)}")
        self.logger.error(f"📍 错误详情 | Error details: {traceback.format_exc()}")
        
        end_time = time.time()
        total_time = end_time - self.start_time
        self.logger.error(f"⏱️ 运行时间 | Runtime before failure: {format_time(total_time)}")
    
    def _log_system_info(self):
        """记录系统信息 | Log system information"""
        import platform
        import psutil
        
        self.logger.info("💻 系统信息 | System Information:")
        self.logger.info(f"  操作系统 | OS: {platform.system()} {platform.release()}")
        self.logger.info(f"  Python版本 | Python: {platform.python_version()}")
        self.logger.info(f"  CPU核心数 | CPU cores: {psutil.cpu_count()}")
        self.logger.info(f"  配置线程数 | Configured threads: {self.config.threads}")
        
        memory = psutil.virtual_memory()
        self.logger.info(f"  系统内存 | System memory: {memory.total // (1024**3)} GB")
        self.logger.info(f"  可用内存 | Available memory: {memory.available // (1024**3)} GB")
    
    def _verify_input_files(self):
        """验证输入文件 | Verify input files"""
        self.logger.info("📁 验证输入文件 | Verifying input files")
        
        files_to_check = []
        
        # HiFi数据
        if self.config.hifi_reads:
            files_to_check.append(("HiFi reads", self.config.hifi_reads))
        
        # Hi-C数据
        if self.config.hic_r1:
            files_to_check.append(("Hi-C R1", self.config.hic_r1))
        if self.config.hic_r2:
            files_to_check.append(("Hi-C R2", self.config.hic_r2))
        
        # ONT数据
        if self.config.ont_reads:
            files_to_check.append(("ONT reads", self.config.ont_reads))
        
        # NGS数据
        if self.config.ngs_r1:
            files_to_check.append(("NGS R1", self.config.ngs_r1))
        if self.config.ngs_r2:
            files_to_check.append(("NGS R2", self.config.ngs_r2))
        
        for name, file_path in files_to_check:
            if not Path(file_path).exists():
                raise FileNotFoundError(f"输入文件不存在 | Input file not found: {name} - {file_path}")
            self.logger.info(f"  ✅ {name}: {Path(file_path).name}")
    
    def _evaluate_qc_results(self, qc_results: Dict[str, Any]):
        """评估质控结果 | Evaluate QC results"""
        critical_issues = []
        
        # 检查HiFi数据
        if 'hifi' in qc_results:
            hifi = qc_results['hifi']
            if not hifi.get('pass_qc', True):
                critical_issues.append(f"HiFi数据质量不达标: 覆盖度 {hifi.get('coverage', 0):.1f}X")
        
        # 检查Hi-C数据
        if 'hic' in qc_results:
            hic = qc_results['hic']
            if 'pass_qc' in hic and not hic['pass_qc']:
                critical_issues.append(f"Hi-C数据质量不达标: 覆盖度 {hic.get('coverage', 0):.1f}X")
        
        if critical_issues:
            self.logger.warning("⚠️ 数据质量警告 | Data quality warnings:")
            for issue in critical_issues:
                self.logger.warning(f"  - {issue}")
            self.logger.warning("📝 建议检查数据质量，但流程将继续执行 | Recommend checking data quality, but pipeline will continue")
    
    def _validate_assembly_outputs(self, assembly_files: Dict[str, List[str]]):
        """验证组装输出 | Validate assembly outputs"""
        if not assembly_files:
            raise RuntimeError("未生成任何组装文件 | No assembly files generated")
        
        for assembly_type, files in assembly_files.items():
            self.logger.info(f"✅ {assembly_type} 组装: {len(files)} 个文件")
            for file_path in files:
                if not Path(file_path).exists():
                    raise FileNotFoundError(f"组装输出文件缺失 | Assembly output file missing: {file_path}")
    
    def _validate_scaffold_outputs(self, scaffold_results: Optional[Dict[str, Dict[str, str]]]):
        """验证挂载输出 | Validate scaffold outputs"""
        if not scaffold_results:
            self.logger.warning("⚠️ 未生成挂载结果 | No scaffolding results generated")
            return
        
        success_count = 0
        total_count = 0
        
        for assembly_type, results in scaffold_results.items():
            for assembly_name, files in results.items():
                total_count += 1
                if files:  # 至少有一个输出文件
                    success_count += 1
                    self.logger.info(f"✅ {assembly_type}_{assembly_name} 挂载成功")
        
        self.logger.info(f"📊 挂载成功率: {success_count}/{total_count} ({success_count/total_count*100:.1f}%)")
    
    def _log_quality_summary(self, quality_results: Dict[str, Dict[str, Any]]):
        """记录质量评估摘要 | Log quality assessment summary"""
        if 'grades' in quality_results:
            self.logger.info("🏆 质量评级摘要 | Quality grades summary:")
            for name, grade in quality_results['grades'].items():
                self.logger.info(f"  📋 {name}: {grade}")
    
    def _print_results_summary(self, final_results: Dict[str, Any]):
        """打印结果摘要 | Print results summary"""
        print("\n" + "="*70)
        print("🎉 基因组组装完成！| Genome Assembly Completed!")
        print("="*70)
        
        # 项目信息
        project_info = final_results.get('project_info', {})
        print(f"📋 项目: {project_info.get('project_name', 'N/A')}")
        print(f"🧬 策略: {project_info.get('assembly_strategy', 'N/A')}")
        
        # 主要结果
        assemblies = final_results.get('assemblies', {})
        scaffolds = final_results.get('scaffolds', {})
        
        print(f"📊 组装文件: {sum(len(files) for files in assemblies.values())} 个")
        if scaffolds:
            print(f"🔗 挂载文件: {sum(len(results) for results in scaffolds.values())} 个")
        
        # 质量摘要
        quality_summary = final_results.get('quality_summary', {})
        if quality_summary:
            print(f"🏆 总体评级: {quality_summary.get('overall_grade', 'N/A')}")
            
            key_metrics = quality_summary.get('key_metrics', {})
            if 'best_assembly_n50' in key_metrics:
                print(f"📏 最佳N50: {key_metrics['best_assembly_n50']:,} bp")
        
        print(f"📁 结果目录: {self.config.output_dir}")
        print("="*70)
    
    def _print_next_steps(self, final_results: Dict[str, Any]):
        """打印下一步建议 | Print next steps"""
        print("\n💡 下一步建议 | Next Steps:")
        
        # 基于结果给出建议
        juicebox_files = final_results.get('juicebox_files', {})
        if juicebox_files:
            print("1. 🔗 使用Juicebox进行手动校正 (参考 for_juicebox/JUICEBOX_INSTRUCTIONS.txt)")
        
        print("2. 🧬 进行基因组注释 (基因预测、重复序列注释)")
        print("3. 📊 比较基因组学分析")
        print("4. 🔍 变异检测和群体分析")
        
        quality_summary = final_results.get('quality_summary', {})
        recommendations = quality_summary.get('recommendations', [])
        if recommendations:
            print("\n🔧 改进建议 | Improvement Recommendations:")
            for i, rec in enumerate(recommendations[:3], 1):  # 只显示前3个
                print(f"{i}. {rec}")
        
        print(f"\n📋 详细报告: {self.config.output_dir}/results/FINAL_ASSEMBLY_REPORT.txt")


# def main():
#     """主函数 | Main function"""
#     parser = argparse.ArgumentParser(
#         description='🧬 基因组组装工具 (多数据整合版本) | Genome Assembly Tool (Multi-data Integration Version)',
#         formatter_class=argparse.ArgumentDefaultsHelpFormatter,
#         epilog="""
# 🌟 示例用法 | Examples:
#   # 基础HiFi组装 | Basic HiFi assembly:
#   python -m genomeasm -i raw_data/ -o assembly_results/
  
#   # HiFi + Hi-C染色体级组装 | HiFi + Hi-C chromosome-level assembly:
#   python -m genomeasm -i data/ -o results/ --hic-strategy complete_juicer
  
#   # 指定项目参数 | Specify project parameters:
#   python -m genomeasm -i input/ -o output/ -n my_genome --genome-size 3g --threads 64
  
#   # 使用简化Hi-C流程 | Use simplified Hi-C pipeline:
#   python -m genomeasm -i data/ -o results/ --hic-strategy simplified_salsa2

# 📚 支持的数据类型 | Supported Data Types:
#   - HiFi: 高准确长读长数据 (必需) | High-accuracy long reads (required)
#   - Hi-C: 染色体构象捕获数据 | Chromosome conformation capture data
#   - ONT: Oxford Nanopore长读长数据 | Oxford Nanopore long reads
#   - NGS: Illumina短读长数据 | Illumina short reads

# 🔗 Hi-C处理策略 | Hi-C Processing Strategies:
#   - complete_juicer: 完整Juicer + 3D-DNA流程 (最高质量)
#   - standard_3ddna: 简化3D-DNA流程 (平衡质量与复杂度)
#   - simplified_salsa2: SALSA2流程 (最简化)

# 💡 更多信息请访问项目文档 | For more information, visit project documentation
#         """
#     )
    
#     # 必需参数 | Required arguments
#     required = parser.add_argument_group('必需参数 | Required Arguments')
#     required.add_argument('-i', '--input-dir', required=True,
#                          help='输入数据目录 (自动检测文件类型) | Input data directory (auto-detect file types)')
    
#     # 基本参数 | Basic arguments
#     basic = parser.add_argument_group('基本参数 | Basic Arguments')
#     basic.add_argument('-o', '--output-dir', default='./assembly_output',
#                       help='输出目录 | Output directory')
#     basic.add_argument('-n', '--project-name', default='genome_assembly',
#                       help='项目名称 | Project name')
#     basic.add_argument('-t', '--threads', type=int, default=88,
#                       help='线程数 | Number of threads')
    
#     # 质量控制参数 | Quality control arguments
#     qc_group = parser.add_argument_group('质量控制参数 | Quality Control Arguments')
#     qc_group.add_argument('--skip-fastqc', action='store_true', default=True,
#                          help='跳过FastQC质量检查 (默认跳过，节省时间) | Skip FastQC quality check (default: skip to save time)')
#     qc_group.add_argument('--run-fastqc', action='store_true', 
#                          help='运行FastQC质量检查 | Run FastQC quality check')
#     qc_group.add_argument('--min-hifi-coverage', type=int, default=30,
#                          help='最小HiFi覆盖度 | Minimum HiFi coverage')
#     qc_group.add_argument('--min-hic-coverage', type=int, default=50,
#                          help='最小Hi-C覆盖度 | Minimum Hi-C coverage')
    
#     #     # 质量控制参数 | Quality control arguments
#     # qc_group = parser.add_argument_group('质量控制参数 | Quality Control Arguments')
#     # qc_group.add_argument('--min-hifi-coverage', type=int, default=30,
#     #                      help='最小HiFi覆盖度 | Minimum HiFi coverage')
#     # qc_group.add_argument('--min-hic-coverage', type=int, default=50,
#     #                      help='最小Hi-C覆盖度 | Minimum Hi-C coverage')
#     # qc_group.add_argument('--min-mapping-rate', type=float, default=0.7,
#     #                      help='最小映射率 | Minimum mapping rate')
#     # qc_group.add_argument('--busco-lineage', default='auto',
#     #                      help='BUSCO谱系数据库 | BUSCO lineage database')
    
#     # Hi-C参数 | Hi-C arguments
#     hic_group = parser.add_argument_group('Hi-C参数 | Hi-C Arguments')
#     hic_group.add_argument('--hic-strategy', 
#                           choices=['complete_juicer', 'standard_3ddna', 'simplified_salsa2'],
#                           default='complete_juicer',
#                           help='Hi-C处理策略 | Hi-C processing strategy')
#     hic_group.add_argument('--restriction-enzyme', default='MboI',
#                           choices=['MboI', 'DpnII', 'HindIII', 'EcoRI'],
#                           help='限制性酶类型 | Restriction enzyme type')
#     hic_group.add_argument('--min-contig-size', type=int, default=15000,
#                           help='最小contig大小阈值 | Minimum contig size threshold')
#     hic_group.add_argument('--edit-rounds', type=int, default=2,
#                           help='3D-DNA编辑轮数 | 3D-DNA editing rounds')
    
#     # 组装参数 | Assembly arguments  
#     assembly = parser.add_argument_group('组装参数 | Assembly Arguments')
#     assembly.add_argument('--genome-size', default='3g',
#                          help='预估基因组大小 | Estimated genome size (e.g., 3g, 500m)')
#     assembly.add_argument('--species-type', choices=['diploid', 'haploid', 'polyploid'],
#                          default='diploid', help='物种倍性 | Species ploidy')
#     assembly.add_argument('--telomere-motif', default='CCCTAA',
#                          help='端粒序列motif | Telomere sequence motif')
#     assembly.add_argument('--purge-level', type=int, choices=[0, 1, 2, 3], default=1,
#                          help='Purging级别 | Purging level')
#     assembly.add_argument('--purge-max', type=int, default=80,
#                          help='Purging覆盖度上限 | Purging coverage upper limit')
#     assembly.add_argument('--similarity-threshold', type=float, default=0.75,
#                          help='相似度阈值 | Similarity threshold')
#     assembly.add_argument('--n-haplotypes', type=int, default=2,
#                          help='单倍型数量 | Number of haplotypes')
    
#     # 工具路径参数 | Tool path arguments
#     tools = parser.add_argument_group('工具路径 | Tool Paths')
#     tools.add_argument('--hifiasm-path', default='hifiasm',
#                       help='Hifiasm程序路径 | Hifiasm program path')
#     tools.add_argument('--bwa-path', default='bwa',
#                       help='BWA程序路径 | BWA program path')
#     tools.add_argument('--samtools-path', default='samtools',
#                       help='Samtools程序路径 | Samtools program path')
#     tools.add_argument('--juicer-path', default='juicer.sh',
#                       help='Juicer脚本路径 | Juicer script path')
#     tools.add_argument('--pipeline-3ddna', default='3d-dna/run-asm-pipeline.sh',
#                       help='3D-DNA pipeline路径 | 3D-DNA pipeline path')
#     tools.add_argument('--juicer-tools', default='juicer_tools.jar',
#                       help='Juicer tools JAR路径 | Juicer tools JAR path')
#     tools.add_argument('--salsa2-path', default='run_pipeline.py',
#                       help='SALSA2脚本路径 | SALSA2 script path')
#     # 在main函数的parser.add_argument部分添加：
#     tools.add_argument('--skip-dependency-check', action='store_true',
#                       help='跳过依赖软件检查 | Skip dependency check')
    
#     # 解析参数
#     args = parser.parse_args()

#     # 处理FastQC参数逻辑
#     if args.run_fastqc:
#         skip_fastqc = False  # 如果明确指定运行FastQC
#     else:
#         skip_fastqc = args.skip_fastqc  # 否则按skip_fastqc参数决定
    
#     try:
#         # 创建组装器实例
#         assembler = GenomeAssembler(
#             input_dir=args.input_dir,
#             output_dir=args.output_dir,
#             project_name=args.project_name,
#             threads=args.threads,
#             hic_strategy=args.hic_strategy,
#             restriction_enzyme=args.restriction_enzyme,
#             min_contig_size=args.min_contig_size,
#             edit_rounds=args.edit_rounds,
#             genome_size=args.genome_size,
#             species_type=args.species_type,
#             telomere_motif=args.telomere_motif,
#             purge_level=args.purge_level,
#             purge_max=args.purge_max,
#             similarity_threshold=args.similarity_threshold,
#             n_haplotypes=args.n_haplotypes,
#             skip_fastqc=skip_fastqc,  # 添加这个参数
#             min_hifi_coverage=args.min_hifi_coverage,
#             min_hic_coverage=args.min_hic_coverage,
#             min_mapping_rate=args.min_mapping_rate,
#             busco_lineage=args.busco_lineage,
#             hifiasm_path=args.hifiasm_path,
#             bwa_path=args.bwa_path,
#             samtools_path=args.samtools_path,
#             juicer_path=args.juicer_path,
#             pipeline_3ddna=args.pipeline_3ddna,
#             juicer_tools=args.juicer_tools,
#             salsa2_path=args.salsa2_path
#         )
        
#         # 运行组装流程
#         final_results = assembler.run_assembly()
        
#         # 成功退出
#         sys.exit(0)
        
#     except KeyboardInterrupt:
#         print("\n⏹️ 用户中断程序执行 | User interrupted program execution")
#         sys.exit(1)
#     except Exception as e:
#         print(f"\n💥 程序执行失败 | Program execution failed: {e}")
#         sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 基因组组装工具 (多数据整合版本) | Genome Assembly Tool (Multi-data Integration Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
🌟 示例用法 | Examples:
  # 基础HiFi组装 | Basic HiFi assembly:
  python -m genomeasm -i raw_data/ -o assembly_results/
  
  # HiFi + Hi-C染色体级组装 | HiFi + Hi-C chromosome-level assembly:
  python -m genomeasm -i data/ -o results/ --hic-strategy complete_juicer
  
  # 指定项目参数 | Specify project parameters:
  python -m genomeasm -i input/ -o output/ -n my_genome --genome-size 3g --threads 64
  
  # 使用简化Hi-C流程 | Use simplified Hi-C pipeline:
  python -m genomeasm -i data/ -o results/ --hic-strategy simplified_salsa2

📚 支持的数据类型 | Supported Data Types:
  - HiFi: 高准确长读长数据 (必需) | High-accuracy long reads (required)
  - Hi-C: 染色体构象捕获数据 | Chromosome conformation capture data
  - ONT: Oxford Nanopore长读长数据 | Oxford Nanopore long reads
  - NGS: Illumina短读长数据 | Illumina short reads

🔗 Hi-C处理策略 | Hi-C Processing Strategies:
  - complete_juicer: 完整Juicer + 3D-DNA流程 (最高质量)
  - standard_3ddna: 简化3D-DNA流程 (平衡质量与复杂度)
  - simplified_salsa2: SALSA2流程 (最简化)

💡 更多信息请访问项目文档 | For more information, visit project documentation
        """
    )
    
    # 必需参数 | Required arguments
    required = parser.add_argument_group('必需参数 | Required Arguments')
    required.add_argument('-i', '--input-dir', required=True,
                         help='输入数据目录 (自动检测文件类型) | Input data directory (auto-detect file types)')
    
    # 基本参数 | Basic arguments
    basic = parser.add_argument_group('基本参数 | Basic Arguments')
    basic.add_argument('-o', '--output-dir', default='./assembly_output',
                      help='输出目录 | Output directory')
    basic.add_argument('-n', '--project-name', default='genome_assembly',
                      help='项目名称 | Project name')
    basic.add_argument('-t', '--threads', type=int, default=88,
                      help='线程数 | Number of threads')
    
    # Hi-C参数 | Hi-C arguments
    hic_group = parser.add_argument_group('Hi-C参数 | Hi-C Arguments')
    hic_group.add_argument('--hic-strategy', 
                          choices=['complete_juicer', 'standard_3ddna', 'simplified_salsa2'],
                          default='complete_juicer',
                          help='Hi-C处理策略 | Hi-C processing strategy')
    hic_group.add_argument('--restriction-enzyme', default='MboI',
                          choices=['MboI', 'DpnII', 'HindIII', 'EcoRI'],
                          help='限制性酶类型 | Restriction enzyme type')
    hic_group.add_argument('--min-contig-size', type=int, default=15000,
                          help='最小contig大小阈值 | Minimum contig size threshold')
    hic_group.add_argument('--edit-rounds', type=int, default=2,
                          help='3D-DNA编辑轮数 | 3D-DNA editing rounds')
    
    # 组装参数 | Assembly arguments  
    assembly = parser.add_argument_group('组装参数 | Assembly Arguments')
    assembly.add_argument('--genome-size', default='3g',
                         help='预估基因组大小 | Estimated genome size (e.g., 3g, 500m)')
    assembly.add_argument('--species-type', choices=['diploid', 'haploid', 'polyploid'],
                         default='diploid', help='物种倍性 | Species ploidy')
    assembly.add_argument('--telomere-motif', default='CCCTAA',
                         help='端粒序列motif | Telomere sequence motif')
    assembly.add_argument('--purge-level', type=int, choices=[0, 1, 2, 3], default=1,
                         help='Purging级别 | Purging level')
    assembly.add_argument('--purge-max', type=int, default=80,
                         help='Purging覆盖度上限 | Purging coverage upper limit')
    assembly.add_argument('--similarity-threshold', type=float, default=0.75,
                         help='相似度阈值 | Similarity threshold')
    assembly.add_argument('--n-haplotypes', type=int, default=2,
                         help='单倍型数量 | Number of haplotypes')
    
    # 质量控制参数 | Quality control arguments
    qc_group = parser.add_argument_group('质量控制参数 | Quality Control Arguments')
    qc_group.add_argument('--skip-fastqc', action='store_true', default=True,
                         help='跳过FastQC质量检查 (默认跳过，节省时间) | Skip FastQC quality check (default: skip to save time)')
    qc_group.add_argument('--run-fastqc', action='store_true', 
                         help='运行FastQC质量检查 | Run FastQC quality check')
    qc_group.add_argument('--min-hifi-coverage', type=int, default=30,
                         help='最小HiFi覆盖度 | Minimum HiFi coverage')
    qc_group.add_argument('--min-hic-coverage', type=int, default=50,
                         help='最小Hi-C覆盖度 | Minimum Hi-C coverage')
    qc_group.add_argument('--min-mapping-rate', type=float, default=0.7,
                         help='最小映射率 | Minimum mapping rate')
    qc_group.add_argument('--busco-lineage', default='auto',
                         help='BUSCO谱系数据库 | BUSCO lineage database')
    
    # 工具路径参数 | Tool path arguments
    tools = parser.add_argument_group('工具路径 | Tool Paths')
    tools.add_argument('--hifiasm-path', default='hifiasm',
                      help='Hifiasm程序路径 | Hifiasm program path')
    tools.add_argument('--bwa-path', default='bwa',
                      help='BWA程序路径 | BWA program path')
    tools.add_argument('--samtools-path', default='samtools',
                      help='Samtools程序路径 | Samtools program path')
    tools.add_argument('--juicer-path', default='juicer.sh',
                      help='Juicer脚本路径 | Juicer script path')
    tools.add_argument('--pipeline-3ddna', default='3d-dna/run-asm-pipeline.sh',
                      help='3D-DNA pipeline路径 | 3D-DNA pipeline path')
    tools.add_argument('--juicer-tools', default='juicer_tools.jar',
                      help='Juicer tools JAR路径 | Juicer tools JAR path')
    tools.add_argument('--salsa2-path', default='run_pipeline.py',
                      help='SALSA2脚本路径 | SALSA2 script path')
    tools.add_argument('--skip-dependency-check', action='store_true',
                      help='跳过依赖软件检查 | Skip dependency check')
    
    # 解析参数
    args = parser.parse_args()
    
    # 处理FastQC参数逻辑
    if args.run_fastqc:
        skip_fastqc = False  # 如果明确指定运行FastQC
    else:
        skip_fastqc = args.skip_fastqc  # 否则按skip_fastqc参数决定
    
    try:
        # 创建组装器实例
        assembler = GenomeAssembler(
            # 基本参数
            input_dir=args.input_dir,
            output_dir=args.output_dir,
            project_name=args.project_name,
            threads=args.threads,
            
            # Hi-C相关参数
            hic_strategy=args.hic_strategy,
            restriction_enzyme=args.restriction_enzyme,
            min_contig_size=args.min_contig_size,
            edit_rounds=args.edit_rounds,
            
            # 组装参数
            genome_size=args.genome_size,
            species_type=args.species_type,
            telomere_motif=args.telomere_motif,
            purge_level=args.purge_level,
            purge_max=args.purge_max,
            similarity_threshold=args.similarity_threshold,
            n_haplotypes=args.n_haplotypes,
            
            # 质量控制参数
            skip_fastqc=skip_fastqc,
            min_hifi_coverage=args.min_hifi_coverage,
            min_hic_coverage=args.min_hic_coverage,
            min_mapping_rate=args.min_mapping_rate,
            busco_lineage=args.busco_lineage,
            
            # 工具路径参数
            hifiasm_path=args.hifiasm_path,
            bwa_path=args.bwa_path,
            samtools_path=args.samtools_path,
            juicer_path=args.juicer_path,
            pipeline_3ddna=args.pipeline_3ddna,
            juicer_tools=args.juicer_tools,
            salsa2_path=args.salsa2_path,
            skip_dependency_check=getattr(args, 'skip_dependency_check', False)
        )
        
        # 运行组装流程
        final_results = assembler.run_assembly()
        
        # 成功退出
        sys.exit(0)
        
    except KeyboardInterrupt:
        print("\n⏹️ 用户中断程序执行 | User interrupted program execution")
        sys.exit(1)
    except Exception as e:
        print(f"\n💥 程序执行失败 | Program execution failed: {e}")
        sys.exit(1)

if __name__ == "__main__":
    main()
