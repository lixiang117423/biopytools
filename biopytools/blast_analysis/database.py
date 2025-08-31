"""
🗄️ BLAST数据库管理模块
"""

from pathlib import Path
from .utils import CommandRunner

class DatabaseManager:
    """🗄️ BLAST数据库管理器"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.target_db_path = None
    
    def create_target_database(self) -> str:
        """创建目标基因数据库"""
        self.logger.info("🗄️ 创建BLAST数据库...")
        
        target_file = self.config.target_file
        target_path = Path(target_file)
        
        # 检查是否已存在数据库文件
        db_extensions = ['.nhr', '.nin', '.nsq'] if self.config.get_db_type() == 'nucl' else ['.phr', '.pin', '.psq']
        db_exists = all((target_path.parent / f"{target_path.stem}{ext}").exists() for ext in db_extensions)
        
        if db_exists:
            self.logger.info("✅ 数据库文件已存在，跳过创建")
            self.target_db_path = str(target_path.parent / target_path.stem)
            return self.target_db_path
        
        # 创建数据库
        db_name = target_path.parent / f"{target_path.stem}_db"
        self.target_db_path = str(db_name)
        
        cmd = (
            f"{self.config.makeblastdb_path} "
            f"-in {target_file} "
            f"-dbtype {self.config.get_db_type()} "
            f"-out {self.target_db_path}"
        )
        
        success = self.cmd_runner.run(cmd, "创建BLAST数据库")
        
        if success:
            self.logger.info(f"✅ BLAST数据库创建成功: {self.target_db_path}")
        else:
            raise RuntimeError("BLAST数据库创建失败")
        
        return self.target_db_path
    
    def create_target_database(self) -> str:
        """创建目标基因数据库"""
        self.logger.info("🗄️ 创建BLAST数据库...")
        
        target_file = self.config.target_file
        target_path = Path(target_file)
        
        # 在输出目录下创建blast_db子目录
        blast_db_dir = self.config.output_path / "blast_db"
        blast_db_dir.mkdir(exist_ok=True)
        self.logger.info(f"📁 数据库目录: {blast_db_dir}")
        
        # 设置数据库文件路径到blast_db目录
        db_name = blast_db_dir / f"{target_path.stem}_db"
        self.target_db_path = str(db_name)
        
        # 检查是否已存在数据库文件
        db_extensions = ['.nhr', '.nin', '.nsq'] if self.config.get_db_type() == 'nucl' else ['.phr', '.pin', '.psq']
        db_exists = all((Path(f"{self.target_db_path}{ext}").exists() for ext in db_extensions))
        
        if db_exists:
            self.logger.info("✅ 数据库文件已存在，跳过创建")
            return self.target_db_path
        
        # 创建数据库
        cmd = (
            f"{self.config.makeblastdb_path} "
            f"-in {target_file} "
            f"-dbtype {self.config.get_db_type()} "
            f"-out {self.target_db_path}"
        )
        
        success = self.cmd_runner.run(cmd, "创建BLAST数据库")
        
        if success:
            self.logger.info(f"✅ BLAST数据库创建成功: {self.target_db_path}")
        else:
            raise RuntimeError("BLAST数据库创建失败")
        
        return self.target_db_path
