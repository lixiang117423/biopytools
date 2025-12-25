"""
基因组组装核心模块 | Genome Assembly Core Module 🧬
"""

import os
from pathlib import Path
from typing import Dict, List, Tuple, Optional
from .utils import CommandRunner, get_file_stats

class HifiasmAssembler:
    """Hifiasm组装器 | Hifiasm Assembler"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.assembly_dir = Path(config.output_dir) / "assembly"
        self.assembly_dir.mkdir(parents=True, exist_ok=True)
    
    def run_assembly(self) -> Dict[str, List[str]]:
        """运行Hifiasm组装 | Run Hifiasm assembly"""
        self.logger.info("🧬 开始Hifiasm基因组组装 | Starting Hifiasm genome assembly")
        
        # 构建hifiasm命令
        hifiasm_cmd = self._build_hifiasm_command()  # 移除attempt参数
        
        # 执行组装
        success = self.cmd_runner.run(
            hifiasm_cmd, 
            f"Hifiasm组装 - {self.config.assembly_strategy}策略",
            timeout=86400  # 24小时超时
        )
        
        if not success:
            raise RuntimeError("❌ Hifiasm组装失败 | Hifiasm assembly failed")
        
        # 检查输出文件
        output_files = self._check_assembly_outputs()
        
        # 转换GFA到FASTA
        fasta_files = self._convert_gfa_to_fasta(output_files)
        
        # 组装统计
        self._generate_assembly_stats(fasta_files)
        
        self.logger.info("✅ Hifiasm组装完成 | Hifiasm assembly completed")
        return fasta_files
    
    def _build_hifiasm_command(self) -> str:
        """构建hifiasm命令 | Build hifiasm command"""
        cmd_parts = [
            self.config.hifiasm_path,
            f"-o {self.assembly_dir / self.config.project_name}",
            f"-t {self.config.threads}"
        ]
        
        # 基本参数
        cmd_parts.extend([
            "--primary",  # 输出primary assembly
            f"-l {self.config.purge_level}",  # purging级别
            f"--purge-max {self.config.purge_max}",  # purging覆盖度上限
            f"--n-hap {self.config.n_haplotypes}"  # 单倍型数量
        ])
        
        # Hi-C数据整合
        if 'Hi-C' in self.config.detected_data_types:
            cmd_parts.extend([
                f"--h1 {self.config.hic_r1}",
                f"--h2 {self.config.hic_r2}",
                "--dual-scaf"  # 双重脚手架
            ])
        
        # ONT数据整合
        if 'ONT' in self.config.detected_data_types:
            cmd_parts.extend([
                f"--ul {self.config.ont_reads}",
                "--ul-rate 0.15"  # ONT错误率
            ])
        
        # 端粒序列保护
        if self.config.telomere_motif:
            cmd_parts.append(f"--telo-m {self.config.telomere_motif}")
        
        # 物种特异参数
        if self.config.species_type == "haploid":
            cmd_parts.append("--hom-cov auto")
        
        # 相似度阈值调整
        cmd_parts.append(f"-s {self.config.similarity_threshold}")
        
        # 添加HiFi reads文件
        cmd_parts.append(self.config.hifi_reads)
        
        return " ".join(cmd_parts)
    
    def _check_assembly_outputs(self) -> Dict[str, List[str]]:
        """检查组装输出文件 | Check assembly output files"""
        self.logger.info("🔍 检查组装输出文件 | Checking assembly output files")
        
        expected_outputs = {}
        base_name = self.config.project_name
        
        # Primary contigs (基础版本)
        primary_bp = self.assembly_dir / f"{base_name}.bp.p_ctg.gfa"
        if primary_bp.exists():
            expected_outputs.setdefault('primary', []).append(str(primary_bp))
            self.logger.info(f"✅ 找到基础primary contigs: {primary_bp}")
        
        # Hi-C优化版本 (如果有Hi-C数据)
        if 'Hi-C' in self.config.detected_data_types:
            primary_hic = self.assembly_dir / f"{base_name}.hic.p_ctg.gfa"
            if primary_hic.exists():
                expected_outputs.setdefault('primary', []).append(str(primary_hic))
                self.logger.info(f"✅ 找到Hi-C优化primary contigs: {primary_hic}")
        
        # Alternate contigs
        alternate_bp = self.assembly_dir / f"{base_name}.bp.a_ctg.gfa"
        if alternate_bp.exists():
            expected_outputs.setdefault('alternate', []).append(str(alternate_bp))
            self.logger.info(f"✅ 找到基础alternate contigs: {alternate_bp}")
        
        if 'Hi-C' in self.config.detected_data_types:
            alternate_hic = self.assembly_dir / f"{base_name}.hic.a_ctg.gfa"
            if alternate_hic.exists():
                expected_outputs.setdefault('alternate', []).append(str(alternate_hic))
                self.logger.info(f"✅ 找到Hi-C优化alternate contigs: {alternate_hic}")
        
        # Haplotype-resolved contigs (如果是杂合基因组)
        for hap in [1, 2]:
            hap_file = self.assembly_dir / f"{base_name}.dip.hap{hap}.p_ctg.gfa"
            if hap_file.exists():
                expected_outputs.setdefault(f'haplotype_{hap}', []).append(str(hap_file))
                self.logger.info(f"✅ 找到单倍型{hap}组装: {hap_file}")
        
        if not expected_outputs:
            raise RuntimeError("❌ 未找到预期的组装输出文件 | No expected assembly output files found")
        
        return expected_outputs
    
    def _convert_gfa_to_fasta(self, gfa_files: Dict[str, List[str]]) -> Dict[str, List[str]]:
        """转换GFA文件到FASTA格式 | Convert GFA files to FASTA format"""
        self.logger.info("🔄 转换GFA文件到FASTA格式 | Converting GFA files to FASTA format")
        
        fasta_files = {}
        
        for assembly_type, gfa_list in gfa_files.items():
            fasta_files[assembly_type] = []
            
            for i, gfa_file in enumerate(gfa_list):
                gfa_path = Path(gfa_file)
                
                # 构造FASTA文件名，区分不同版本
                if len(gfa_list) > 1:
                    # 多个版本时添加后缀区分
                    version_suffix = "hic" if "hic" in gfa_path.name else "bp"
                    fasta_name = f"{self.config.project_name}_{assembly_type}_{version_suffix}.fa"
                else:
                    fasta_name = f"{self.config.project_name}_{assembly_type}.fa"
                
                fasta_path = self.assembly_dir / fasta_name
                
                # 转换GFA到FASTA
                cmd = f"awk '/^S/{{print \">\"$2; print $3}}' {gfa_file} > {fasta_path}"
                success = self.cmd_runner.run(cmd, f"转换{assembly_type}组装到FASTA")
                
                if success and fasta_path.exists():
                    fasta_files[assembly_type].append(str(fasta_path))
                    self.logger.info(f"✅ 转换完成: {fasta_path}")
                else:
                    self.logger.error(f"❌ 转换失败: {gfa_file} -> {fasta_path}")
        
        return fasta_files
    
    def _generate_assembly_stats(self, fasta_files: Dict[str, List[str]]):
        """生成组装统计信息 | Generate assembly statistics"""
        self.logger.info("📊 生成组装统计信息 | Generating assembly statistics")
        
        stats_file = self.assembly_dir / "assembly_statistics.txt"
        
        with open(stats_file, 'w', encoding='utf-8') as f:
            f.write("🧬 基因组组装统计报告 | Genome Assembly Statistics Report\n")
            f.write("=" * 70 + "\n\n")
            f.write(f"项目名称 | Project: {self.config.project_name}\n")
            f.write(f"组装策略 | Assembly strategy: {self.config.assembly_strategy}\n")
            f.write(f"数据类型 | Data types: {', '.join(self.config.detected_data_types)}\n\n")
            
            for assembly_type, fasta_list in fasta_files.items():
                f.write(f"📋 {assembly_type.upper()} 组装统计:\n")
                f.write("-" * 50 + "\n")
                
                for fasta_file in fasta_list:
                    fasta_name = Path(fasta_file).name
                    f.write(f"文件 | File: {fasta_name}\n")
                    
                    # 获取统计信息
                    stats = get_file_stats(fasta_file, self.logger)
                    if stats:
                        f.write(f"  序列数 | Sequences: {stats.get('num', 'N/A')}\n")
                        f.write(f"  总长度 | Total length: {stats.get('sum_len', 'N/A')}\n")
                        f.write(f"  最大长度 | Max length: {stats.get('max_len', 'N/A')}\n")
                        f.write(f"  平均长度 | Average length: {stats.get('avg_len', 'N/A')}\n")
                        if 'N50' in stats:
                            f.write(f"  N50: {stats.get('N50', 'N/A')}\n")
                    f.write("\n")
        
        self.logger.info(f"📋 组装统计已保存: {stats_file}")
