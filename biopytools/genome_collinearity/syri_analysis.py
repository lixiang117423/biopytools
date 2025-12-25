"""
SyRI结构变异分析模块 | SyRI Structural Variation Analysis Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class SyRIAnalyzer:
    """SyRI分析器 | SyRI Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        
        # 创建SyRI结果目录 | Create SyRI results directory
        self.syri_dir = Path(self.config.output_dir) / "syri_results"
        self.syri_dir.mkdir(parents=True, exist_ok=True)
    
    def run_syri_analysis(self, alignment_files: dict) -> dict:
        """运行SyRI结构变异分析 | Run SyRI structural variation analysis"""
        self.logger.info("🔬 开始SyRI结构变异分析 | Starting SyRI structural variation analysis")
        
        syri_files = {}
        
        for pair_name, bam_file in alignment_files.items():
            # 解析样本名称 | Parse sample names
            ref_sample, query_sample = pair_name.split('_')
            
            ref_genome = self.config.genome_paths[ref_sample]
            query_genome = self.config.genome_paths[query_sample]
            
            # 生成输出文件前缀 | Generate output file prefix
            output_prefix = self.syri_dir / pair_name
            
            # 运行SyRI | Run SyRI
            success = self._run_syri(
                bam_file, ref_genome, query_genome, output_prefix, pair_name
            )
            
            if success:
                # SyRI使用--dir和--prefix参数后，输出文件名格式为: {prefix}syri.out
                syri_output = self.syri_dir / f"{pair_name}syri.out"
                
                # 验证输出文件是否存在 | Verify output file exists
                if syri_output.exists():
                    syri_files[pair_name] = str(syri_output)
                    self.logger.info(f"✅ SyRI输出文件已生成 | SyRI output file generated: {syri_output}")
                else:
                    self.logger.error(f"❌ SyRI输出文件不存在 | SyRI output file not found: {syri_output}")
                    return {}
            else:
                self.logger.error(f"❌ SyRI分析失败 | SyRI analysis failed: {pair_name}")
                return {}
        
        self.logger.info(f"🎉 完成 {len(syri_files)} 个SyRI分析 | Completed {len(syri_files)} SyRI analyses")
        return syri_files
    
    def _run_syri(self, bam_file: str, ref_genome: str, query_genome: str, 
                  output_prefix: Path, pair_name: str) -> bool:
        """运行SyRI分析单个比对 | Run SyRI analysis for single alignment"""
        
        self.logger.info(f"🔍 SyRI分析 | SyRI analysis: {pair_name}")
        
        # 构建SyRI命令 | Build SyRI command
        # SyRI需要使用--dir指定输出目录，--prefix只用于文件名前缀
        output_dir = output_prefix.parent
        prefix_name = output_prefix.name
        
        self.logger.info(f"📁 SyRI输出目录 | SyRI output directory: {output_dir}")
        self.logger.info(f"🏷️ SyRI文件前缀 | SyRI file prefix: {prefix_name}")
        
        cmd = (
            f"{self.config.syri_path} -c {bam_file} -r {ref_genome} "
            f"-q {query_genome} -F B --dir {output_dir} --prefix {prefix_name}"
        )
        
        description = f"🔬 SyRI analysis: {pair_name}"
        return self.cmd_runner.run(cmd, description)
    
    def check_syri_completed(self) -> dict:
        """检查SyRI分析是否已完成 | Check if SyRI analysis is completed"""
        self.logger.info("🔍 检查现有SyRI结果文件 | Checking existing SyRI result files")
        
        existing_files = {}
        
        for i in range(len(self.config.sample_list) - 1):
            ref_sample = self.config.sample_list[i]
            query_sample = self.config.sample_list[i + 1]
            
            pair_name = f"{ref_sample}_{query_sample}"
            syri_output = self.syri_dir / f"{pair_name}syri.out"
            
            if syri_output.exists():
                existing_files[pair_name] = str(syri_output)
                self.logger.info(f"✅ 找到已完成的SyRI分析 | Found completed SyRI analysis: {pair_name}")
            else:
                self.logger.info(f"❌ SyRI结果文件不存在 | SyRI result file missing: {pair_name}")
        
        return existing_files
