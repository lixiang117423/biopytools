"""
GenomeSyn2核心分析逻辑模块|GenomeSyn2 Core Analysis Logic Module
"""

import os
from pathlib import Path
from .utils import GenomeSyn2CommandRunner


class GenomeSyn2Calculator:
    """GenomeSyn2核心计算类|GenomeSyn2 Core Calculator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.cmd_runner = GenomeSyn2CommandRunner(logger, config.output_path if hasattr(config, 'output_path') else None)

    def run_alignment(self):
        """运行序列比对|Run sequence alignment"""
        self.logger.info(f"开始序列比对|Starting sequence alignment: {self.config.align}")

        # 构建参数列表|Build argument list
        args = []

        # 必需参数|Required parameters
        args.extend(['--align', self.config.align])
        args.extend(['--genome', self.config.genome_dir])
        args.extend(['--outdir', self.config.outdir])
        args.extend(['--thread', str(self.config.threads)])

        # 蛋白质比对需要基因目录|Protein alignment requires gene directory
        if self.config.align in ['blastp', 'mmseqs', 'diamond']:
            if not self.config.gene_dir:
                self.logger.error(
                    f"{self.config.align}模式需要--gene参数|"
                    f"{self.config.align} mode requires --gene parameter"
                )
                return False
            args.extend(['--gene', self.config.gene_dir])

        # 执行命令|Execute command
        success = self.cmd_runner.run_perl_script_via_shell(
            perl_path=self.config.perl_path,
            script_path=self.config.genomesyn2_pl,
            args=args,
            description=f"{self.config.align}序列比对|{self.config.align} sequence alignment"
        )

        if success:
            self.logger.info(f"序列比对完成|Sequence alignment completed: {self.config.align}")
        else:
            self.logger.error(f"序列比对失败|Sequence alignment failed: {self.config.align}")

        return success

    def calculate_snp_from_vcf(self):
        """从VCF计算SNP密度和一致性|Calculate SNP density and identity from VCF"""
        self.logger.info("开始计算SNP密度和一致性|Starting SNP density and identity calculation")

        # 构建参数列表|Build argument list
        args = []
        args.extend(['--vcf', self.config.vcf])
        args.extend(['--bin', str(self.config.bin)])
        args.extend(['--type', 'identity'])

        # 执行命令|Execute command
        success = self.cmd_runner.run_perl_script_via_shell(
            perl_path=self.config.perl_path,
            script_path=self.config.genomesyn2_pl,
            args=args,
            description="计算SNP密度和一致性|Calculate SNP density and identity"
        )

        if success:
            self.logger.info("SNP计算完成|SNP calculation completed")
        else:
            self.logger.error("SNP计算失败|SNP calculation failed")

        return success

    def plot_ancestry_deconvolution(self):
        """绘制祖先血统解析图|Plot ancestry deconvolution"""
        self.logger.info("开始绘制祖先血统解析图|Starting ancestry deconvolution plotting")

        # 构建参数列表|Build argument list
        args = []
        args.extend(['--type', 'identity'])
        args.extend(['--identity', self.config.identity])
        args.extend(['--density', self.config.density])

        # 执行命令|Execute command
        success = self.cmd_runner.run_perl_script_via_shell(
            perl_path=self.config.perl_path,
            script_path=self.config.genomesyn2_pl,
            args=args,
            description="绘制祖先血统解析图|Plot ancestry deconvolution"
        )

        if success:
            self.logger.info("血统解析图绘制完成|Ancestry deconvolution plotting completed")
        else:
            self.logger.error("血统解析图绘制失败|Ancestry deconvolution plotting failed")

        return success

    def run_plotting(self):
        """运行共线性绘图|Run synteny plotting"""
        self.logger.info("开始绘制共线性图|Starting synteny plotting")

        # 构建参数列表|Build argument list
        args = ['--conf', self.config.conf]

        # 如果指定了anno参数，添加到命令|Add anno parameter if specified
        if self.config.anno:
            args.append('--anno')
            # 将anno视为True时，添加问号生成配置模板
            # 实际使用中，用户应该先运行 --conf ? 生成模板
            # 然后再运行 --conf total.conf 进行绘图
            # 这里我们假设配置文件已经准备好

        # 执行命令|Execute command
        success = self.cmd_runner.run_perl_script_via_shell(
            perl_path=self.config.perl_path,
            script_path=self.config.genomesyn2_pl,
            args=args,
            description="绘制共线性图|Plot synteny diagram"
        )

        if success:
            self.logger.info("共线性图绘制完成|Synteny plotting completed")
        else:
            self.logger.error("共线性图绘制失败|Synteny plotting failed")

        return success

    def generate_file_list(self):
        """生成文件列表|Generate file list"""
        self.logger.info(f"开始生成{self.config.type}文件列表|Starting {self.config.type} file list generation")

        # 构建参数列表|Build argument list
        args = []
        args.extend(['--type', self.config.type])
        args.extend(['--path', self.config.path])
        args.extend(['--out', self.config.out])

        # 执行命令|Execute command
        success = self.cmd_runner.run_perl_script_via_shell(
            perl_path=self.config.perl_path,
            script_path=self.config.genomesyn2_pl,
            args=args,
            description=f"生成{self.config.type}文件列表|Generate {self.config.type} file list"
        )

        if success:
            self.logger.info(f"{self.config.type}文件列表生成完成|{self.config.type} file list generation completed")
        else:
            self.logger.error(f"{self.config.type}文件列表生成失败|{self.config.type} file list generation failed")

        return success

    def run_analysis(self):
        """运行分析|Run analysis"""
        # 判断运行模式|Determine run mode
        if self.config.align:
            # 比对模式|Alignment mode
            return self.run_alignment()
        elif self.config.vcf:
            # VCF计算模式|VCF calculation mode
            return self.calculate_snp_from_vcf()
        elif self.config.identity and self.config.density:
            # 绘图模式|Plotting mode
            return self.plot_ancestry_deconvolution()
        elif self.config.conf:
            # 绘图模式（配置文件）|Plotting mode (config file)
            return self.run_plotting()
        elif self.config.type:
            # 文件生成模式|File generation mode
            return self.generate_file_list()
        else:
            self.logger.error("未指定有效的运行模式|No valid run mode specified")
            return False
