"""
多序列比对核心模块 | MSA Core Alignment Module
"""

from pathlib import Path
from .utils import CommandRunner

class MSAAligner:
    """MSA比对器 | MSA Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_file = f"{self.config.output_prefix}.{self.config.output_format}"
    
    def run_mafft(self):
        """运行MAFFT比对 | Run MAFFT alignment"""
        self.logger.info("🧬 开始MAFFT比对 | Starting MAFFT alignment")
        
        # 构建命令 | Build command
        cmd_parts = [self.config.mafft_path]
        
        # 添加策略参数 | Add strategy parameters
        if self.config.mafft_strategy == 'linsi':
            cmd_parts.append('--localpair --maxiterate 1000')
        elif self.config.mafft_strategy == 'ginsi':
            cmd_parts.append('--globalpair --maxiterate 1000')
        elif self.config.mafft_strategy == 'einsi':
            cmd_parts.append('--ep 0 --genafpair --maxiterate 1000')
        elif self.config.mafft_strategy == 'fftns':
            cmd_parts.append('--retree 2')
        elif self.config.mafft_strategy == 'fftnsi':
            cmd_parts.append('--retree 2 --maxiterate 2')
        else:  # auto
            cmd_parts.append(f'--auto --maxiterate {self.config.mafft_maxiterate}')
        
        # 添加线程数 | Add thread count
        cmd_parts.append(f'--thread {self.config.threads}')
        
        # 添加输入输出 | Add input/output
        cmd_parts.append(f'{self.config.input_file} > {self.output_file}')
        
        cmd = ' '.join(cmd_parts)
        
        return self.cmd_runner.run(cmd, "MAFFT序列比对 | MAFFT alignment")
    
    def run_clustalo(self):
        """运行Clustal Omega比对 | Run Clustal Omega alignment"""
        self.logger.info("🧬 开始Clustal Omega比对 | Starting Clustal Omega alignment")
        
        cmd_parts = [
            self.config.clustalo_path,
            f'-i {self.config.input_file}',
            f'-o {self.output_file}',
            f'--threads={self.config.threads}',
            '--force'  # 强制覆盖输出文件 | Force overwrite
        ]
        
        # 添加迭代参数 | Add iteration parameter
        if self.config.clustalo_iterations > 0:
            cmd_parts.append(f'--iter={self.config.clustalo_iterations}')
        
        # 输出格式 | Output format
        if self.config.output_format == 'clustal':
            cmd_parts.append('--outfmt=clustal')
        elif self.config.output_format == 'phylip':
            cmd_parts.append('--outfmt=phylip')
        else:  # fasta
            cmd_parts.append('--outfmt=fasta')
        
        cmd = ' '.join(cmd_parts)
        
        return self.cmd_runner.run(cmd, "Clustal Omega序列比对 | Clustal Omega alignment")
    
    def run_muscle(self):
        """运行MUSCLE比对 | Run MUSCLE alignment"""
        self.logger.info("🧬 开始MUSCLE比对 | Starting MUSCLE alignment")
        
        cmd_parts = [
            self.config.muscle_path,
            f'-in {self.config.input_file}',
            f'-out {self.output_file}',
            f'-maxiters {self.config.muscle_maxiters}'
        ]
        
        # MUSCLE v5的参数格式不同 | MUSCLE v5 has different parameter format
        # 这里使用v3/v4的格式 | Using v3/v4 format here
        
        cmd = ' '.join(cmd_parts)
        
        return self.cmd_runner.run(cmd, "MUSCLE序列比对 | MUSCLE alignment")
    
    def run_tcoffee(self):
        """运行T-Coffee比对 | Run T-Coffee alignment"""
        self.logger.info("🧬 开始T-Coffee比对 | Starting T-Coffee alignment")
        
        cmd_parts = [
            self.config.tcoffee_path,
            f'-seq {self.config.input_file}',
            f'-outfile {self.output_file}',
            f'-n_core {self.config.threads}',
            '-quiet'
        ]
        
        # 输出格式 | Output format
        if self.config.output_format == 'clustal':
            cmd_parts.append('-output clustalw')
        elif self.config.output_format == 'phylip':
            cmd_parts.append('-output phylip')
        else:  # fasta
            cmd_parts.append('-output fasta_aln')
        
        cmd = ' '.join(cmd_parts)
        
        return self.cmd_runner.run(cmd, "T-Coffee序列比对 | T-Coffee alignment")
    
    def align(self):
        """执行比对 | Perform alignment"""
        method_map = {
            'mafft': self.run_mafft,
            'clustalo': self.run_clustalo,
            'muscle': self.run_muscle,
            't_coffee': self.run_tcoffee
        }
        
        align_func = method_map.get(self.config.method)
        if align_func:
            return align_func()
        else:
            self.logger.error(f"❌ 未知的比对方法 | Unknown alignment method: {self.config.method}")
            return False
