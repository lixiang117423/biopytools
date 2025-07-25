"""
基因组组装质量控制模块 | Genome Assembly Quality Control Module
"""

from pathlib import Path
from typing import Dict
from .utils import CommandRunner

class ContaminationScreener:
    """外源污染筛查器 | Contamination Screener"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def screen_contamination(self, assembly_file: str) -> str:
        """筛查外源污染 | Screen for contamination"""
        if not self.config.run_contamination_screen:
            self.logger.info("跳过外源污染筛查 | Skipping contamination screening")
            return assembly_file
        
        self.logger.info("开始外源污染筛查 | Starting contamination screening")
        
        output_dir = self.config.output_path / "contamination_screen"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建FCS命令 | Build FCS command
        cmd = f"{self.config.fcs_path} --fasta {assembly_file} --output-dir {output_dir} --euk"
        
        # 执行筛查 | Execute screening
        success = self.cmd_runner.run(cmd, "NCBI FCS外源污染筛查 | NCBI FCS contamination screening", timeout=7200)
        
        if success:
            clean_file = output_dir / "cleaned_assembly.fasta"
            if clean_file.exists():
                self.logger.info(f"外源污染筛查完成 | Contamination screening completed: {clean_file}")
                return str(clean_file)
        
        self.logger.warning("外源污染筛查失败，使用原始组装 | Contamination screening failed, using original assembly")
        return assembly_file

class AssemblyAnnotator:
    """组装注释器 | Assembly Annotator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_flagger(self, assembly_file: str) -> Dict:
        """运行Flagger错误注释 | Run Flagger error annotation"""
        if not self.config.run_flagger:
            self.logger.info("跳过Flagger错误注释 | Skipping Flagger error annotation")
            return {}
        
        self.logger.info("开始Flagger错误注释 | Starting Flagger error annotation")
        
        output_dir = self.config.output_path / "flagger_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建Flagger命令 | Build Flagger command
        cmd = f"{self.config.flagger_path} --assembly {assembly_file} --reads {self.config.hifi_reads} --output {output_dir}"
        
        # 执行注释 | Execute annotation
        success = self.cmd_runner.run(cmd, "Flagger组装错误注释 | Flagger assembly error annotation", timeout=3600)
        
        results = {}
        if success:
            results_file = output_dir / "flagger_results.tsv"
            if results_file.exists():
                self.logger.info(f"Flagger错误注释完成 | Flagger error annotation completed: {results_file}")
                results['flagger_results'] = str(results_file)
        
        return results
    
    def run_nucfreq(self, assembly_file: str) -> Dict:
        """运行NucFreq分析 | Run NucFreq analysis"""
        self.logger.info("开始NucFreq核苷酸频率分析 | Starting NucFreq nucleotide frequency analysis")
        
        output_dir = self.config.output_path / "nucfreq_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建NucFreq命令 | Build NucFreq command
        cmd = f"{self.config.nucfreq_path} {assembly_file} > {output_dir}/nucfreq_results.txt"
        
        # 执行分析 | Execute analysis
        success = self.cmd_runner.run(cmd, "NucFreq核苷酸频率分析 | NucFreq nucleotide frequency analysis")
        
        results = {}
        if success:
            results_file = output_dir / "nucfreq_results.txt"
            if results_file.exists():
                self.logger.info(f"NucFreq分析完成 | NucFreq analysis completed: {results_file}")
                results['nucfreq_results'] = str(results_file)
        
        return results
