"""
RepeatModeler工具包装器 | RepeatModeler Tool Wrapper
"""

import os
from pathlib import Path

class RepeatModelerRunner:
    """🔄 RepeatModeler运行器 | RepeatModeler Runner"""
    
    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.output_lib = None
    
    # def build_database(self) -> bool:
    #     """🏗️ 构建RepeatModeler数据库 | Build RepeatModeler database"""
    #     db_name = f"{self.config.base_name}_db"
    #     db_path = self.config.output_path / db_name
        
    #     cmd = f"{self.config.repeatmodeler_path} -database {db_name} -genomefile {self.config.genome_file}"
        
    #     success = self.cmd_runner.run(
    #         cmd, 
    #         f"🗄️ 构建RepeatModeler数据库 | Build RepeatModeler database: {db_name}"
    #     )
        
    #     if success:
    #         self.db_name = db_name
    #         self.logger.info(f"✅ RepeatModeler数据库创建成功 | RepeatModeler database created: {db_name}")
        
    #     return success
    def build_database(self) -> bool:
        """🏗️ 构建RepeatModeler数据库 | Build RepeatModeler database"""
        db_name = f"{self.config.base_name}_db"
        
        cmd = f"BuildDatabase -name {db_name} {self.config.genome_file}"
        
        success = self.cmd_runner.run(
            cmd, 
            f"🗄️ 构建RepeatModeler数据库 | Build RepeatModeler database: {db_name}"
        )
        
        # 添加这部分成功处理逻辑
        if success:
            self.db_name = db_name
            self.logger.info(f"✅ RepeatModeler数据库创建成功 | RepeatModeler database created: {db_name}")
        
        return success
    
    # def run_repeatmodeler(self) -> bool:
    #     """🚀 运行RepeatModeler | Run RepeatModeler"""
    #     if not hasattr(self, 'db_name'):
    #         self.logger.error("❌ 数据库未创建，无法运行RepeatModeler | Database not created, cannot run RepeatModeler")
    #         return False
        
    #     cmd = f"{self.config.repeatmodeler_path} -database {self.db_name} -threads {self.config.threads}"
        
    #     success = self.cmd_runner.run(
    #         cmd,
    #         f"🧬 运行RepeatModeler从头构建重复库 | Run RepeatModeler de novo repeat library construction",
    #         timeout=None  # RepeatModeler可能需要很长时间
    #     )
        
    #     if success:
    #         # 查找输出库文件 | Find output library file
    #         expected_lib = self.config.output_path / f"{self.db_name}-families.fa"
    #         if expected_lib.exists():
    #             self.output_lib = expected_lib
    #             self.logger.info(f"✅ RepeatModeler库文件生成 | RepeatModeler library file generated: {expected_lib}")
    #         else:
    #             self.logger.warning(f"⚠️ 未找到预期的库文件 | Expected library file not found: {expected_lib}")
        
    #     return success
    # def run_repeatmodeler(self) -> bool:
    #     """🚀 运行RepeatModeler | Run RepeatModeler"""
    #     if not hasattr(self, 'db_name'):
    #         self.logger.error("❌ 数据库未创建，无法运行RepeatModeler | Database not created, cannot run RepeatModeler")
    #         return False
        
    #     # RepeatModeler只需要database参数，不需要genomefile
    #     cmd = f"{self.config.repeatmodeler_path} -database {self.db_name} -threads {self.config.threads}"
        
    #     success = self.cmd_runner.run(
    #         cmd,
    #         f"🧬 运行RepeatModeler从头构建重复库 | Run RepeatModeler de novo repeat library construction",
    #         timeout=None
    #     )

    #     if success:
    #         # 查找输出库文件 | Find output library file
    #         expected_lib = self.config.output_path / f"{self.db_name}-families.fa"
    #         if expected_lib.exists():
    #             self.output_lib = expected_lib
    #             self.logger.info(f"✅ RepeatModeler库文件生成 | RepeatModeler library file generated: {expected_lib}")
    #         else:
    #             self.logger.warning(f"⚠️ 未找到预期的库文件 | Expected library file not found: {expected_lib}")

    def run_repeatmodeler(self) -> bool:
        """🚀 运行RepeatModeler | Run RepeatModeler"""
        if not hasattr(self, 'db_name'):
            self.logger.error("❌ 数据库未创建，无法运行RepeatModeler | Database not created, cannot run RepeatModeler")
            return False
        
        cmd = f"{self.config.repeatmodeler_path} -database {self.db_name} -threads {self.config.threads}"
        
        success = self.cmd_runner.run(
            cmd,
            f"🧬 运行RepeatModeler从头构建重复库 | Run RepeatModeler de novo repeat library construction",
            timeout=None
        )
        
        # 不管命令是否"成功"，都检查是否生成了库文件
        # RepeatModeler有时会生成文件但返回非零退出码
        expected_lib = self.config.output_path / f"{self.db_name}-families.fa"
        if expected_lib.exists():
            self.output_lib = expected_lib
            self.logger.info(f"✅ RepeatModeler库文件生成 | RepeatModeler library file generated: {expected_lib}")
            return True  # 强制返回True，因为库文件存在
        else:
            self.logger.warning(f"⚠️ 未找到预期的库文件 | Expected library file not found: {expected_lib}")
            return False
    
    def get_library_path(self) -> Path:
        """📁 获取生成的库文件路径 | Get generated library file path"""
        return self.output_lib
