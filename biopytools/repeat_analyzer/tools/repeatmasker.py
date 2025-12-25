"""
RepeatMasker工具包装器 | RepeatMasker Tool Wrapper
"""

import os
from pathlib import Path

class RepeatMaskerRunner:
    """🎭 RepeatMasker运行器 | RepeatMasker Runner"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_files = {}
    
    def run_repeatmasker(self, library_path: Path) -> bool:
        """🚀 运行RepeatMasker | Run RepeatMasker"""
        if not library_path or not library_path.exists():
            self.logger.error(f"❌ 重复库文件不存在 | Repeat library file does not exist: {library_path}")
            return False
        
        # RepeatMasker输出目录 | RepeatMasker output directory
        output_dir = self.config.output_path / "repeatmasker_output"
        output_dir.mkdir(exist_ok=True)
        
        # 复制基因组文件到输出目录 | Copy genome file to output directory
        genome_copy = output_dir / Path(self.config.genome_file).name
        if not genome_copy.exists():
            import shutil
            shutil.copy2(self.config.genome_file, genome_copy)
        
        cmd = (
            f"cd {output_dir} && "
            f"{self.config.repeatmasker_path} -lib {library_path} "
            # 新增 -xsmall 参数以启用软屏蔽（小写字母）
            # Add -xsmall parameter to enable soft masking (lowercase letters)
            f"-xsmall "
            f"-pa {self.config.threads} -gff -dir {output_dir} {genome_copy.name}"
        )
        
        success = self.cmd_runner.run(
            cmd,
            "🎭 运行RepeatMasker识别重复序列 | Run RepeatMasker to identify repeat sequences"
        )
        
        if success:
            # 收集输出文件 | Collect output files
            genome_name = Path(self.config.genome_file).name
            possible_outputs = {
                'out': output_dir / f"{genome_name}.out",
                'masked': output_dir / f"{genome_name}.masked",
                'gff': output_dir / f"{genome_name}.out.gff",
                'tbl': output_dir / f"{genome_name}.tbl"
            }
            
            for file_type, file_path in possible_outputs.items():
                if file_path.exists():
                    self.output_files[file_type] = file_path
                    self.logger.info(f"✅ RepeatMasker {file_type} 文件生成 | RepeatMasker {file_type} file generated: {file_path}")
        
        return success
    
    def get_output_files(self) -> dict:
        """📁 获取输出文件字典 | Get output files dictionary"""
        return self.output_files
    
    def get_masked_genome(self) -> Path:
        """🎭 获取屏蔽基因组文件 | Get masked genome file"""
        return self.output_files.get('masked')
    
    def get_annotation_file(self) -> Path:
        """📋 获取注释文件 | Get annotation file"""
        return self.output_files.get('out')
