"""
TEsorter工具包装器 | TEsorter Tool Wrapper
"""

import os
from pathlib import Path

class TESorterRunner:
    """🏷️ TEsorter运行器 | TEsorter Runner"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_files = {}
    
    def run_tesorter(self, repeat_library: Path) -> bool:
        """🚀 运行TEsorter | Run TEsorter"""
        if not repeat_library or not repeat_library.exists():
            self.logger.error(f"❌ 重复库文件不存在 | Repeat library file does not exist: {repeat_library}")
            return False
        
        # TEsorter输出目录 | TEsorter output directory  
        output_dir = self.config.output_path / "tesorter_output"
        output_dir.mkdir(exist_ok=True)
        
        output_prefix = output_dir / f"{self.config.base_name}_tesorter"
        
        cmd = (
            f"{self.config.tesorter_path} {repeat_library} "
            f"-db rexdb-plant "
            f"-pre {output_prefix} "
            f"-p {self.config.threads}"
        )
        
        success = self.cmd_runner.run(
            cmd,
            "🏷️ 运行TEsorter对转座元件进行分类 | Run TEsorter to classify transposable elements"
        )
        
        if success:
            # 收集输出文件 | Collect output files
            possible_outputs = {
                'cls': Path(f"{output_prefix}.cls.lib"),
                'cls_tsv': Path(f"{output_prefix}.cls.tsv"),
                'rexdb': Path(f"{output_prefix}.rexdb.cls.lib"),
                'rexdb_tsv': Path(f"{output_prefix}.rexdb.cls.tsv")
            }
            
            for file_type, file_path in possible_outputs.items():
                if file_path.exists():
                    self.output_files[file_type] = file_path
                    self.logger.info(f"✅ TEsorter {file_type} 文件生成 | TEsorter {file_type} file generated: {file_path}")
        
        return success
    
    def get_output_files(self) -> dict:
        """📁 获取输出文件字典 | Get output files dictionary"""
        return self.output_files
    
    def get_classified_library(self) -> Path:
        """📚 获取分类后的库文件 | Get classified library file"""
        return self.output_files.get('cls') or self.output_files.get('rexdb')
    
    def get_classification_table(self) -> Path:
        """📊 获取分类结果表格 | Get classification result table"""
        return self.output_files.get('cls_tsv') or self.output_files.get('rexdb_tsv')
