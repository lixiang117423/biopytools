"""
多序列比对核心模块|MSA Core Alignment Module
"""

from pathlib import Path
from .utils import CommandRunner

class MSAAligner:
    """MSA比对器|MSA Aligner"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_file = f"{self.config.output_prefix}.{self.config.output_format}"
    
    def run_mafft(self):
        """运行MAFFT比对|Run MAFFT alignment"""
        self.logger.info(" 开始MAFFT比对|Starting MAFFT alignment")
        
        # 构建命令参数|Build command arguments
        args = []
        
        # 添加策略参数|Add strategy parameters
        if self.config.mafft_strategy == 'linsi':
            args.extend(['--localpair', '--maxiterate', '1000'])
        elif self.config.mafft_strategy == 'ginsi':
            args.extend(['--globalpair', '--maxiterate', '1000'])
        elif self.config.mafft_strategy == 'einsi':
            args.extend(['--ep', '0', '--genafpair', '--maxiterate', '1000'])
        elif self.config.mafft_strategy == 'fftns':
            args.extend(['--retree', '2'])
        elif self.config.mafft_strategy == 'fftnsi':
            args.extend(['--retree', '2', '--maxiterate', '2'])
        else:  # auto
            args.extend(['--auto', '--maxiterate', str(self.config.mafft_maxiterate)])
        
        # 添加线程数|Add thread count
        args.extend(['--thread', str(self.config.threads)])
        
        # 添加输入|Add input
        args.append(self.config.input_file)
        
        return self.cmd_runner.run(
            [self.config.mafft_path] + args,
            "MAFFT序列比对|MAFFT alignment",
            output_file=self.output_file
        )
    
    def run_clustalo(self):
        """运行Clustal Omega比对|Run Clustal Omega alignment"""
        self.logger.info(" 开始Clustal Omega比对|Starting Clustal Omega alignment")
        
        args = [
            '-i', self.config.input_file,
            '-o', self.output_file,
            f'--threads={self.config.threads}',
            '--force'  # 强制覆盖输出文件|Force overwrite
        ]
        
        # 添加迭代参数|Add iteration parameter
        if self.config.clustalo_iterations > 0:
            args.append(f'--iter={self.config.clustalo_iterations}')
        
        # 输出格式|Output format
        if self.config.output_format == 'clustal':
            args.append('--outfmt=clustal')
        elif self.config.output_format == 'phylip':
            args.append('--outfmt=phylip')
        else:  # fasta
            args.append('--outfmt=fasta')
        
        return self.cmd_runner.run(
            [self.config.clustalo_path] + args, "Clustal Omega序列比对|Clustal Omega alignment")
    
    def run_muscle(self):
        """运行MUSCLE比对|Run MUSCLE alignment"""
        self.logger.info(" 开始MUSCLE比对|Starting MUSCLE alignment")
        
        args = [
            '-in', self.config.input_file,
            '-out', self.output_file,
            '-maxiters', str(self.config.muscle_maxiters)
        ]
        
        # MUSCLE v5的参数格式不同|MUSCLE v5 has different parameter format
        # 这里使用v3/v4的格式|Using v3/v4 format here
        
        return self.cmd_runner.run(
            [self.config.muscle_path] + args, "MUSCLE序列比对|MUSCLE alignment")
    
    def run_tcoffee(self):
        """运行T-Coffee比对|Run T-Coffee alignment"""
        self.logger.info(" 开始T-Coffee比对|Starting T-Coffee alignment")
        
        args = [
            '-seq', self.config.input_file,
            '-outfile', self.output_file,
            '-n_core', str(self.config.threads),
            '-quiet'
        ]
        
        # 输出格式|Output format
        if self.config.output_format == 'clustal':
            args.extend(['-output', 'clustalw'])
        elif self.config.output_format == 'phylip':
            args.extend(['-output', 'phylip'])
        else:  # fasta
            args.extend(['-output', 'fasta_aln'])
        
        return self.cmd_runner.run(
            [self.config.tcoffee_path] + args, "T-Coffee序列比对|T-Coffee alignment")
    
    def align(self):
        """执行比对|Perform alignment"""
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
            self.logger.error(f" 未知的比对方法|Unknown alignment method: {self.config.method}")
            return False
