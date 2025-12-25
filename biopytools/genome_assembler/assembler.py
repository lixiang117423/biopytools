"""
基因组组装执行模块 🧬 | Genome Assembly Execution Module
"""

import os
from pathlib import Path

try:
    from utils import run_command
except ImportError:
    from biopytools.genome_assembler.utils import run_command

class HifiasmAssembler:
    """Hifiasm组装器 | Hifiasm Assembler"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def run_assembly(self) -> bool:
        """运行hifiasm组装 | Run hifiasm assembly"""
        self.logger.info("🚀 运行hifiasm组装 (HiFi + Hi-C模式) | Running hifiasm assembly (HiFi + Hi-C mode)")
        
        # 构建hifiasm命令 | Build hifiasm command
        cmd_parts = [
            "hifiasm",
            f"-o {self.config.prefix}",
            f"-t {self.config.threads}",
            f"--hg-size {self.config.genome_size}",
            f"--n-hap {self.config.n_hap}",
            f"--h1 {self.config.hic_r1}",
            f"--h2 {self.config.hic_r2}",
            "--primary",
            "-l 3",
            self.config.hifi_data
        ]
        
        cmd = " ".join(cmd_parts)
        
        # 执行命令并记录日志 | Execute command and log
        log_file = os.path.join(self.config.log_dir, "hifiasm_assembly.log")
        full_cmd = f"{cmd} 2>&1 | tee {log_file}"
        
        success = run_command(
            full_cmd, 
            self.logger, 
            work_dir=self.config.raw_dir,
            capture_output=False
        )
        
        if success:
            self.logger.info("✅ 组装完成！| Assembly completed!")
            return True
        else:
            self.logger.error("❌ 组装失败，请检查错误信息| Assembly failed, check error messages")
            return False
    
    def convert_gfa_to_fasta(self) -> dict:
        """转换GFA为FASTA格式 | Convert GFA to FASTA format"""
        self.logger.info("🚀 转换GFA为FASTA格式| Converting GFA to FASTA format")
        
        # 定义需要转换的GFA文件和对应的FASTA文件名 (去除前缀)
        gfa_files = {
            f"{self.config.prefix}.hic.hap1.p_ctg.gfa": "hap1.fa",
            f"{self.config.prefix}.hic.hap2.p_ctg.gfa": "hap2.fa",
            f"{self.config.prefix}.hic.p_ctg.gfa": "primary.fa",
            f"{self.config.prefix}.hic.a_ctg.gfa": "alternate.fa"
        }
        
        results = {}
        
        for gfa_file, fasta_name in gfa_files.items():
            gfa_path = os.path.join(self.config.raw_dir, gfa_file)
            fasta_path = os.path.join(self.config.fasta_dir, fasta_name)
            
            if os.path.exists(gfa_path):
                self.logger.info(f"转换 | Converting: {gfa_file} -> {fasta_name}")

                # 使用awk转换 | Convert using awk
                cmd = f"""awk '/^S/{{print ">"$2; print $3}}' {gfa_path} > {fasta_path}"""

                if run_command(cmd, self.logger):
                    stats = self._get_fasta_summary(fasta_path)
                    results[fasta_name] = stats
                    self.logger.info(f"✅ 转换成功 | Successfully converted: {fasta_name}")
                    self.logger.info(f"   序列数 | Sequences: {stats['num_seqs']}")
                    self.logger.info(f"   总长度 | Total length: {stats['total_len']:,} bp")
                else:
                    self.logger.error(f"❌ 转换失败 | Failed to convert: {gfa_file}")
            else:
                self.logger.warning(f"⚠️ 文件不存在 | File not found: {gfa_file}")
        
        return results
    
    def _get_fasta_summary(self, fasta_file: str) -> dict:
        """获取FASTA文件摘要 | Get FASTA file summary"""
        if not os.path.exists(fasta_file):
            return {}
        
        # 使用shell命令快速统计 | Use shell commands for quick statistics
        cmd = f"""grep -c "^>" {fasta_file}"""
        num_seqs = os.popen(cmd).read().strip()
        
        cmd = f"""awk '/^>/{{next}}{{sum+=length($0)}}END{{print sum}}' {fasta_file}"""
        total_len = os.popen(cmd).read().strip()
        
        return {
            'num_seqs': int(num_seqs) if num_seqs else 0,
            'total_len': int(total_len) if total_len else 0
        }
