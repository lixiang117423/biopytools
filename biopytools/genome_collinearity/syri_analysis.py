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
                syri_output = f"{output_prefix}syri.out"
                syri_files[pair_name] = syri_output
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
        cmd = (
            f"{self.config.syri_path} -c {bam_file} -r {ref_genome} "
            f"-q {query_genome} -F B --prefix {output_prefix}"
        )
        
        description = f"🔬 SyRI analysis: {pair_name}"
        return self.cmd_runner.run(cmd, description)
