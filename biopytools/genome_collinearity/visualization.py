"""
PlotSR可视化模块 | PlotSR Visualization Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class CollinearityVisualizer:
    """共线性可视化器 | Collinearity Visualizer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        
        # 创建可视化结果目录 | Create visualization results directory
        self.plot_dir = Path(self.config.output_dir) / "plots"
        self.plot_dir.mkdir(parents=True, exist_ok=True)
    
    def create_plotsr_visualization(self, syri_files: dict) -> bool:
        """创建PlotSR可视化 | Create PlotSR visualization"""
        self.logger.info("🎨 开始PlotSR可视化 | Starting PlotSR visualization")
        
        # 生成genomes.txt文件 | Generate genomes.txt file
        genomes_file = self._create_genomes_file()
        
        # 生成输出文件名 | Generate output filename
        if self.config.chromosome:
            output_file = self.plot_dir / f"collinearity_{self.config.chromosome}.{self.config.plotsr_format}"
        else:
            output_file = self.plot_dir / f"collinearity_all.{self.config.plotsr_format}"
        
        # 构建plotsr命令 | Build plotsr command
        cmd_parts = [f"{self.config.plotsr_path}"]
        
        # 添加SyRI结果文件 | Add SyRI result files
        for pair_name in sorted(syri_files.keys()):
            cmd_parts.append(f"--sr {syri_files[pair_name]}")
        
        # 添加基因组文件 | Add genome file
        cmd_parts.append(f"--genomes {genomes_file}")
        
        # 添加输出文件 | Add output file
        cmd_parts.append(f"-o {output_file}")
        
        # 添加染色体筛选 | Add chromosome filtering
        if self.config.chromosome:
            cmd_parts.append(f"--chr {self.config.chromosome}")
        
        # 添加其他参数 | Add other parameters
        if self.config.skip_synteny:
            cmd_parts.append("--nosyn")
        
        # 添加图形参数 | Add figure parameters
        cmd_parts.append(f"--figsize {self.config.figure_width} {self.config.figure_height}")
        
        cmd = " ".join(cmd_parts)
        description = "🎨 PlotSR可视化 | PlotSR visualization"
        
        success = self.cmd_runner.run(cmd, description)
        
        if success:
            self.logger.info(f"🖼️ 可视化文件已生成 | Visualization file generated: {output_file}")
        
        return success
    
    def _create_genomes_file(self) -> str:
        """创建genomes.txt文件 | Create genomes.txt file"""
        genomes_file = self.plot_dir / "genomes.txt"
        
        with open(genomes_file, 'w') as f:
            f.write("#file\tname\ttags\n")
            for sample in self.config.sample_list:
                genome_path = self.config.genome_paths[sample]
                f.write(f"{genome_path}\t{sample}\tlw:{self.config.line_width}\n")
        
        self.logger.info(f"📝 基因组配置文件已创建 | Genome configuration file created: {genomes_file}")
        return str(genomes_file)
