"""
ALLHiC流水线主程序模块 🎯 | ALLHiC Pipeline Main Module
"""

import os
import sys
import time
import argparse
from pathlib import Path

# 添加当前目录到sys.path，解决相对导入问题
current_dir = Path(__file__).parent
if str(current_dir) not in sys.path:
    sys.path.insert(0, str(current_dir))

try:
    from config import PipelineConfig
    from logger import PipelineLogger
    from utils import check_dependencies, setup_conda_environment
    from steps import ALLHiCSteps
    from asmkit import AsmkitProcessor
except ImportError:
    # 如果相对导入失败，尝试绝对导入
    sys.path.append(str(current_dir.parent))
    from allhic.config import PipelineConfig
    from allhic.logger import PipelineLogger
    from allhic.utils import check_dependencies, setup_conda_environment
    from allhic.steps import ALLHiCSteps
    from allhic.asmkit import AsmkitProcessor

class ALLHiCPipeline:
    """ALLHiC流水线主类 | Main ALLHiC Pipeline Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = PipelineConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = PipelineLogger(self.config.work_dir_abs)

        # 初始化组件 | Initialize components
        self.steps = ALLHiCSteps(self.config, self.logger_manager)
        self.asmkit_processor = AsmkitProcessor(self.config, self.logger_manager)
    
    def run_pipeline(self):
        """运行完整的流水线 | Run complete pipeline"""
        start_time = time.time()
        
        try:
            self.logger_manager.log_section("🚀 启动ALLHiC流水线 v5.4 (Asmkit版本) | Starting ALLHiC Pipeline v5.4 (Asmkit Edition)")
            
            # 0. 环境设置
            self._setup_environment()
            
            # 1. 依赖检查
            check_dependencies(self.config, self.logger_manager)
            
            # 2. 运行诊断模式（如果启用）
            if self.config.diagnose_mode:
                self._run_diagnostic()
                return
            
            # 3. 执行ALLHiC步骤
            self._run_allhic_steps()
            
            # 4. 执行asmkit步骤
            self.asmkit_processor.run_asmkit_jbat()
            
            # 5. 输出完成信息
            end_time = time.time()
            elapsed = int(end_time - start_time)
            self._output_completion_summary(elapsed)
            
        except Exception as e:
            self.logger_manager.log_error(f"流水线执行过程中意外终止 | Pipeline execution terminated unexpectedly: {e}")
            sys.exit(1)
    
    def _setup_environment(self):
        """设置环境 | Setup environment"""
        self.logger_manager.log_section("🔧 设置环境 | Setting up environment")
        
        # 设置conda环境 | Setup conda environment
        setup_conda_environment(self.logger_manager)
        
        # 设置PATH | Setup PATH
        bwa_bin = self.config.bwa_path
        allhic_dir = self.config.allhic_bin_path  # ALLHiC工具在根目录
        allhic_scripts = os.path.join(self.config.allhic_software_path, "scripts")

        current_path = os.environ.get('PATH', '')
        new_path = f"{bwa_bin}:{allhic_dir}:{allhic_scripts}:{current_path}"
        os.environ['PATH'] = new_path
        
        # 进入工作目录 | Enter working directory
        os.chdir(self.config.work_dir_abs)

        self.logger_manager.log(f"工作目录 | Working directory: {self.config.work_dir_abs}")
        self.logger_manager.log(f"线程数 | Threads: {self.config.threads}")
        self.logger_manager.log_step_time()
    
    def _run_allhic_steps(self):
        """运行ALLHiC步骤 | Run ALLHiC steps"""
        # Step 0-7
        self.steps.step0_prepare_data()
        self.steps.step1_mapping()
        self.steps.step1_5_allele_detection()
        self.steps.step2_pruning()
        self.steps.step3_partition()
        self.steps.step3_5_extract()
        self.steps.step4_rescue()
        self.steps.step5_optimize()
        self.steps.step6_build()
        self.steps.step7_plot()
    
    def _run_diagnostic(self):
        """运行诊断模式 | Run diagnostic mode"""
        self.logger_manager.log_section("🔍 流水线状态诊断 | Diagnosing Pipeline State")
        
        critical_files = [
            (os.path.join(self.config.directories["extract"], "sample.clean.clm"), "CLM矩阵"),
            (self.config.reference, "参考基因组"),
            (os.path.join(self.config.directories["mapping"], "sample.clean.bam"), "比对BAM文件"),
            (os.path.join(self.config.directories["build"], "groups.agp"), "AGP文件")
        ]
        
        for file_path, description in critical_files:
            if os.path.exists(file_path):
                size = self.logger_manager.format_file_size(file_path)
                self.logger_manager.log(f"  ✅ 找到 | Found: {description} ({size})")
            else:
                self.logger_manager.log(f"  ❌ 缺失 | Missing: {description}")

        self.logger_manager.log("诊断完成 | Diagnosis complete")
    
    def _output_completion_summary(self, elapsed: int):
        """输出完成摘要 | Output completion summary"""
        self.logger_manager.log_section(f"🎉 流水线完成 | Pipeline Completed in {self.logger_manager.format_time(elapsed)}")

        self.logger_manager.log("")
        self.logger_manager.log("📂 输出摘要 | Output Summary:")
        for dir_name in sorted(self.config.directories.values()):
            self.logger_manager.log(f"   📂 {dir_name}")

        self.logger_manager.log("")
        self.logger_manager.log(f"📝 完整日志 | Full log: {self.logger_manager.log_file}")

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='🧬 ALLHiC Pipeline v5.4 (Asmkit版本) - Hi-C基因组支架构建工具',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
示例 | Examples:
  # 基本用法 | Basic usage
  %(prog)s -r draft.asm.fasta -1 reads_R1.fastq.gz -2 reads_R2.fastq.gz -k 12
  
  # 自定义参数 | Custom parameters
  %(prog)s --reference genome.fa --read1 hic_R1.fq.gz --read2 hic_R2.fq.gz \\
      --chr-num 24 --threads 88 --enzyme GATC
  
  # 跳过某些步骤 | Skip certain steps
  %(prog)s -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -k 12 \\
      --skip-prune --skip-rescue
  
  # 诊断模式 | Diagnostic mode
  %(prog)s -r genome.fa -1 hic_R1.fq.gz -2 hic_R2.fq.gz -k 12 --diagnose
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('-r', '--reference', required=True,
                       help='📄 参考基因组文件 | Reference genome file')
    parser.add_argument('-1', '--read1', required=True,
                       help='📄 Hi-C读段1文件 | Hi-C read 1 file')
    parser.add_argument('-2', '--read2', required=True,
                       help='📄 Hi-C读段2文件 | Hi-C read 2 file')
    parser.add_argument('-k', '--chr-num', required=True, type=int,
                       help='🔢 染色体数量 | Number of chromosomes')
    
    # 基本参数 | Basic parameters
    parser.add_argument('-e', '--enzyme', default="GATC",
                       help='✂️ ALLHiC酶切位点 | ALLHiC enzyme motif')
    parser.add_argument('-t', '--threads', default=88, type=int,
                       help='⚡ CPU线程数 | CPU threads')
    parser.add_argument('-w', '--workdir', default="./allhic_output",
                       help='📂 工作目录 | Working directory')
    
    # MapQ过滤参数 | MapQ filtering parameters
    parser.add_argument('--mapq-step1', default=1, type=int,
                       help='📏 Step 1 MapQ阈值 | Step 1 MapQ threshold')
    
    # 绘图参数 | Plotting parameters
    parser.add_argument('--bin-size', default="500k",
                       help='📊 二进制大小 | Bin size')
    parser.add_argument('--min-bin-size', default="50k",
                       help='📊 最小二进制大小 | Minimum bin size')
    
    # 跳过步骤选项 | Skip step options
    parser.add_argument('--skip-mapping', action='store_true',
                       help='⏭️ 跳过步骤1: 比对 | Skip Step 1: Mapping')
    parser.add_argument('--skip-allele', action='store_true',
                       help='⏭️ 跳过步骤1.5: 等位基因检测 | Skip Step 1.5: Allele Detection')
    parser.add_argument('--skip-prune', action='store_true',
                       help='⏭️ 跳过步骤2: 修剪 | Skip Step 2: Pruning')
    parser.add_argument('--skip-partition', action='store_true',
                       help='⏭️ 跳过步骤3: 分区 | Skip Step 3: Partition')
    parser.add_argument('--skip-extract', action='store_true',
                       help='⏭️ 跳过步骤3.5: 矩阵提取 | Skip Step 3.5: Extract Matrix')
    parser.add_argument('--skip-rescue', action='store_true',
                       help='⏭️ 跳过步骤4: 拯救 | Skip Step 4: Rescue')
    parser.add_argument('--skip-optimize', action='store_true',
                       help='⏭️ 跳过步骤5: 优化 | Skip Step 5: Optimization')
    parser.add_argument('--skip-build', action='store_true',
                       help='⏭️ 跳过步骤6: 构建 | Skip Step 6: Build FASTA')
    parser.add_argument('--skip-plot', action='store_true',
                       help='⏭️ 跳过步骤7: 绘图 | Skip Step 7: Plot Heatmap')
    parser.add_argument('--skip-asmkit', action='store_true',
                       help='⏭️ 跳过步骤8: JBAT生成 | Skip Step 8: JBAT Generation')
    
    # 其他选项 | Other options
    parser.add_argument('--diagnose', action='store_true',
                       help='🔍 运行诊断模式 | Run diagnostic mode')
    
    args = parser.parse_args()
    
    # 构建跳过步骤字典 | Build skip steps dictionary
    skip_steps = {
        "mapping": args.skip_mapping,
        "allele": args.skip_allele,
        "prune": args.skip_prune,
        "partition": args.skip_partition,
        "extract": args.skip_extract,
        "rescue": args.skip_rescue,
        "optimize": args.skip_optimize,
        "build": args.skip_build,
        "plot": args.skip_plot,
        "asmkit": args.skip_asmkit
    }
    
    # 创建流水线并运行 | Create pipeline and run
    pipeline = ALLHiCPipeline(
        reference=args.reference,
        read1=args.read1,
        read2=args.read2,
        chr_num=args.chr_num,
        enzyme=args.enzyme,
        threads=args.threads,
        workdir=args.workdir,
        mapq_step1=args.mapq_step1,
        bin_size=args.bin_size,
        min_bin_size=args.min_bin_size,
        skip_steps=skip_steps,
        diagnose_mode=args.diagnose
    )
    
    pipeline.run_pipeline()

if __name__ == "__main__":
    main()