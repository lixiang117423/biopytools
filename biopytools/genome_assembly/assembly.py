"""
基因组组装核心模块 | Genome Assembly Core Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class VerkkoAssembler:
    """Verkko组装器 | Verkko Assembler"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_verkko_assembly(self) -> str:
        """运行Verkko组装 | Run Verkko assembly"""
        self.logger.info("开始Verkko组装 | Starting Verkko assembly")
        
        output_prefix = self.config.output_path / "verkko_assembly"
        
        # 构建Verkko命令 | Build Verkko command
        cmd_parts = [
            self.config.verkko_path,
            f"--hifi {self.config.hifi_reads}"
        ]
        
        if self.config.ont_reads:
            cmd_parts.append(f"--nano {self.config.ont_reads}")
        
        cmd_parts.extend([
            f"--output {output_prefix}",
            f"--local-cpus {self.config.verkko_local_cpus}",
            f"--local-memory {self.config.verkko_local_memory}g"
        ])
        
        if self.config.verkko_screen:
            cmd_parts.append("--screen")
        
        if self.config.trio_mode:
            cmd_parts.extend([
                f"--hifi-hap1 {self.config.parent1_reads}",
                f"--hifi-hap2 {self.config.parent2_reads}",
                f"--hifi-trio {self.config.child_reads}"
            ])
        
        cmd = " ".join(cmd_parts)
        
        # 执行组装 | Execute assembly
        success = self.cmd_runner.run(cmd, "Verkko基因组组装 | Verkko genome assembly", timeout=86400)  # 24小时超时
        
        if success:
            assembly_file = output_prefix / "assembly.fasta"
            if assembly_file.exists():
                self.logger.info(f"Verkko组装完成 | Verkko assembly completed: {assembly_file}")
                return str(assembly_file)
        
        self.logger.error("Verkko组装失败 | Verkko assembly failed")
        return None

class HifiasmAssembler:
    """hifiasm组装器 | hifiasm Assembler"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_hifiasm_assembly(self) -> str:
        """运行hifiasm组装 | Run hifiasm assembly"""
        self.logger.info("开始hifiasm组装 | Starting hifiasm assembly")
        
        output_prefix = self.config.output_path / "hifiasm_assembly" / "assembly"
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        
        # 构建hifiasm命令 | Build hifiasm command
        cmd_parts = [
            self.config.hifiasm_path,
            f"-o {output_prefix}",
            f"-t {self.config.threads}"
        ]
        
        if self.config.hifiasm_ultra_long:
            cmd_parts.append("--ul")
        
        cmd_parts.extend([
            f"--purge-level {self.config.hifiasm_purge_level}",
            f"-s {self.config.hifiasm_similarity}",
            self.config.hifi_reads
        ])
        
        if self.config.ont_reads:
            cmd_parts.append(self.config.ont_reads)
        
        cmd = " ".join(cmd_parts)
        
        # 执行组装 | Execute assembly
        success = self.cmd_runner.run(cmd, "hifiasm基因组组装 | hifiasm genome assembly", timeout=86400)
        
        if success:
            # 转换GFA到FASTA | Convert GFA to FASTA
            gfa_file = f"{output_prefix}.bp.p_ctg.gfa"
            fasta_file = f"{output_prefix}.fasta"
            
            if os.path.exists(gfa_file):
                convert_cmd = f"awk '/^S/{{print \">\"$2; print $3}}' {gfa_file} > {fasta_file}"
                if self.cmd_runner.run(convert_cmd, "转换GFA到FASTA | Convert GFA to FASTA"):
                    self.logger.info(f"hifiasm组装完成 | hifiasm assembly completed: {fasta_file}")
                    return fasta_file
        
        self.logger.error("hifiasm组装失败 | hifiasm assembly failed")
        return None

class GraphasingProcessor:
    """Graphasing分期处理器 | Graphasing Phasing Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_graphasing(self, assembly_file: str) -> str:
        """运行Graphasing分期 | Run Graphasing phasing"""
        self.logger.info("开始Graphasing分期分析 | Starting Graphasing phasing analysis")
        
        output_prefix = self.config.output_path / "graphasing_output"
        output_prefix.mkdir(parents=True, exist_ok=True)
        
        # 构建Graphasing命令 | Build Graphasing command
        cmd_parts = [
            self.config.graphasing_path,
            f"--assembly {assembly_file}",
            f"--reads {self.config.hifi_reads}",
            f"--output {output_prefix}",
            f"--kmer-size {self.config.graphasing_kmer_size}",
            f"--threads {self.config.threads}"
        ]
        
        cmd = " ".join(cmd_parts)
        
        # 执行分期 | Execute phasing
        success = self.cmd_runner.run(cmd, "Graphasing分期分析 | Graphasing phasing analysis", timeout=21600)  # 6小时超时
        
        if success:
            phased_file = output_prefix / "phased_assembly.fasta"
            if phased_file.exists():
                self.logger.info(f"Graphasing分期完成 | Graphasing phasing completed: {phased_file}")
                return str(phased_file)
        
        self.logger.warning("Graphasing分期失败，使用原始组装 | Graphasing phasing failed, using original assembly")
        return assembly_file
