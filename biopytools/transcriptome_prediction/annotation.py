"""
📍 基因结构注释模块 | Gene Structure Annotation Module
"""

from pathlib import Path
from .utils import CommandRunner

class PASAAnnotator:
    """📍 PASA基因结构注释器 | PASA Gene Structure Annotator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def map_transcripts_to_genome(self):
        """🗺️ 将转录本映射到基因组获得完整基因结构 | Map transcripts to genome to obtain complete gene structures"""
        if not hasattr(self.config, 'trinity_assembly'):
            self.logger.error("❌ 未找到Trinity组装文件，请先执行de novo组装步骤 | No Trinity assembly found, please run de novo assembly step first")
            return False
        
        output_dir = self.config.output_path / "pasa_annotation"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        self.logger.info("🗺️ 使用PASA将转录本映射到基因组 | Using PASA to map transcripts to genome")
        
        # 准备PASA配置文件 | Prepare PASA configuration file
        pasa_config = self._create_pasa_config(output_dir)
        
        # 运行PASA | Run PASA
        cmd = (
            f"{self.config.pasa_path} "
            f"--config {pasa_config} "
            f"--create --run "
            f"--genome {self.config.genome_file} "
            f"--transcripts {self.config.trinity_assembly} "
            f"--ALIGNERS {self.config.pasa_aligners} "
            f"--MAX_INTRON_LENGTH {self.config.pasa_max_intron_length} "
            f"--CPU {self.config.pasa_cpu}"
        )
        
        if not self.cmd_runner.run(cmd, "🗺️ PASA转录本映射 | PASA transcript mapping"):
            return False
        
        # 查找输出文件 | Find output files
        # PASA输出文件命名模式可能因版本而异 | PASA output file naming may vary by version
        possible_outputs = [
            output_dir / "pasa_assemblies.gff3",
            output_dir / "pasa_assemblies.gtf",
            output_dir / f"{self.config.base_name}.pasa_assemblies.gff3",
            output_dir / f"{self.config.base_name}.pasa_assemblies.gtf"
        ]
        
        pasa_output = None
        for output_file in possible_outputs:
            if output_file.exists():
                pasa_output = output_file
                break
        
        if pasa_output:
            self.config.pasa_annotation = str(pasa_output)
            self.logger.info(f"🎉 PASA注释完成 | PASA annotation completed: {pasa_output}")
            return True
        else:
            self.logger.error("❌ PASA注释失败，未找到输出文件 | PASA annotation failed, output file not found")
            # 列出输出目录内容以便调试 | List output directory contents for debugging
            try:
                files = list(output_dir.iterdir())
                self.logger.info(f"📁 输出目录内容 | Output directory contents: {[f.name for f in files]}")
            except:
                pass
            return False
    
    def _create_pasa_config(self, output_dir: Path) -> str:
        """📝 创建PASA配置文件 | Create PASA configuration file"""
        config_file = output_dir / "pasa.config"
        db_name = f"pasa_db_{self.config.base_name}"
        
        # PASA配置文件内容 | PASA configuration file content
        config_content = f"""
# PASA configuration file
# Database configuration
MYSQLDB={db_name}

# Host configuration (use localhost for local MySQL)
MYSQL_RW_HOST=localhost
MYSQL_RO_HOST=localhost

# User configuration (adjust according to your MySQL setup)
MYSQL_RW_USER=root
MYSQL_RO_USER=root

# Password (leave empty for no password, or set according to your MySQL setup)
MYSQL_RW_PASSWORD=
MYSQL_RO_PASSWORD=

# Other parameters
MAX_INTRON_LENGTH={self.config.pasa_max_intron_length}
MIN_PERCENT_ALIGNED={self.config.pasa_min_percent_aligned}
MIN_AVG_PER_ID={self.config.pasa_min_avg_per_id}
"""
        
        with open(config_file, 'w') as f:
            f.write(config_content.strip())
        
        return str(config_file)
