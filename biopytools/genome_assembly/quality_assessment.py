"""
基因组组装质量评估模块 | Genome Assembly Quality Assessment Module
"""

from pathlib import Path
from typing import Dict
from .utils import CommandRunner

class QualityAssessor:
    """组装质量评估器 | Assembly Quality Assessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_merqury(self, assembly_file: str) -> Dict:
        """运行Merqury质量评估 | Run Merqury quality assessment"""
        if not self.config.run_merqury:
            self.logger.info("跳过Merqury质量评估 | Skipping Merqury quality assessment")
            return {}
        
        self.logger.info("开始Merqury质量评估 | Starting Merqury quality assessment")
        
        output_dir = self.config.output_path / "merqury_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 首先生成k-mer数据库 | First generate k-mer database
        kmer_db = output_dir / "reads.meryl"
        cmd_meryl = f"meryl count k=21 {self.config.hifi_reads} output {kmer_db}"
        
        if not self.cmd_runner.run(cmd_meryl, "生成k-mer数据库 | Generate k-mer database", timeout=3600):
            return {}
        
        # 运行Merqury | Run Merqury
        cmd_merqury = f"cd {output_dir} && {self.config.merqury_path} {kmer_db} {assembly_file} assembly"
        
        success = self.cmd_runner.run(cmd_merqury, "Merqury质量评估 | Merqury quality assessment", timeout=3600)
        
        results = {}
        if success:
            qv_file = output_dir / "assembly.qv"
            completeness_file = output_dir / "assembly.completeness.stats"
            
            if qv_file.exists():
                self.logger.info(f"Merqury质量评估完成 | Merqury quality assessment completed")
                results['merqury_qv'] = str(qv_file)
                results['merqury_completeness'] = str(completeness_file)
        
        return results
    
    def run_inspector(self, assembly_file: str) -> Dict:
        """运行Inspector组装检查 | Run Inspector assembly inspection"""
        if not self.config.run_inspector:
            self.logger.info("跳过Inspector组装检查 | Skipping Inspector assembly inspection")
            return {}
        
        self.logger.info("开始Inspector组装检查 | Starting Inspector assembly inspection")
        
        output_dir = self.config.output_path / "inspector_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建Inspector命令 | Build Inspector command
        cmd = f"{self.config.inspector_path} -c {assembly_file} -r {self.config.hifi_reads} -o {output_dir} --thread {self.config.threads}"
        
        # 执行检查 | Execute inspection
        success = self.cmd_runner.run(cmd, "Inspector组装质量检查 | Inspector assembly quality inspection", timeout=3600)
        
        results = {}
        if success:
            results_file = output_dir / "inspector_summary.txt"
            if results_file.exists():
                self.logger.info(f"Inspector检查完成 | Inspector inspection completed: {results_file}")
                results['inspector_results'] = str(results_file)
        
        return results
    
    def run_deepvariant_qv(self, assembly_file: str) -> Dict:
        """运行DeepVariant质量值估计 | Run DeepVariant quality value estimation"""
        if not self.config.run_deepvariant:
            self.logger.info("跳过DeepVariant质量值估计 | Skipping DeepVariant quality value estimation")
            return {}
        
        self.logger.info("开始DeepVariant质量值估计 | Starting DeepVariant quality value estimation")
        
        output_dir = self.config.output_path / "deepvariant_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 首先比对reads到组装 | First align reads to assembly
        bam_file = output_dir / "aligned_reads.bam"
        cmd_align = f"{self.config.minimap2_path} -ax map-hifi {assembly_file} {self.config.hifi_reads} | samtools sort -o {bam_file}"
        
        if not self.cmd_runner.run(cmd_align, "比对reads到组装 | Align reads to assembly", timeout=3600):
            return {}
        
        # 索引BAM文件 | Index BAM file
        cmd_index = f"samtools index {bam_file}"
        if not self.cmd_runner.run(cmd_index, "索引BAM文件 | Index BAM file"):
            return {}
        
        # 运行DeepVariant | Run DeepVariant
        vcf_file = output_dir / "variants.vcf.gz"
        cmd_dv = f"{self.config.deepvariant_path} --model_type=PACBIO --ref={assembly_file} --reads={bam_file} --output_vcf={vcf_file} --num_shards={self.config.threads}"
        
        success = self.cmd_runner.run(cmd_dv, "DeepVariant变异检测 | DeepVariant variant calling", timeout=7200)
        
        results = {}
        if success and vcf_file.exists():
            self.logger.info(f"DeepVariant质量值估计完成 | DeepVariant quality value estimation completed: {vcf_file}")
            results['deepvariant_vcf'] = str(vcf_file)
        
        return results
    
    def run_compleasm(self, assembly_file: str) -> Dict:
        """运行compleasm基因完整性评估 | Run compleasm gene completeness assessment"""
        if not self.config.run_compleasm:
            self.logger.info("跳过compleasm基因完整性评估 | Skipping compleasm gene completeness assessment")
            return {}
        
        self.logger.info("开始compleasm基因完整性评估 | Starting compleasm gene completeness assessment")
        
        output_dir = self.config.output_path / "compleasm_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建compleasm命令 | Build compleasm command
        cmd = f"{self.config.compleasm_path} run -a {assembly_file} -o {output_dir} -l primates_odb{self.config.orthodb_version} -t {self.config.threads}"
        
        # 执行评估 | Execute assessment
        success = self.cmd_runner.run(cmd, "compleasm基因完整性评估 | compleasm gene completeness assessment", timeout=3600)
        
        results = {}
        if success:
            results_file = output_dir / "summary.txt"
            if results_file.exists():
                self.logger.info(f"compleasm评估完成 | compleasm assessment completed: {results_file}")
                results['compleasm_results'] = str(results_file)
        
        return results
