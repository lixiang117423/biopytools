"""
FTP客户端模块 | FTP Client Module
"""

import ftplib
from typing import List, Optional

# 注意: 在模块化模式下，提取脚本会添加以下导入语句:
# from .utils import Logger

class FTPClient:
    """FTP客户端 | FTP Client"""
    
    def __init__(self, host: str, timeout: int = 60, logger=None):
        self.host = host
        self.timeout = timeout
        self.logger = logger
        self.ftp = None
        self.connected = False
    
    def connect(self) -> bool:
        """连接到FTP服务器 | Connect to FTP server"""
        try:
            if self.logger:
                self.logger.info(f"🌐 正在连接到FTP服务器: {self.host}...")
            else:
                print(f"🌐 正在连接到FTP服务器: {self.host}...")
            
            self.ftp = ftplib.FTP(self.host, timeout=self.timeout)
            self.ftp.login()
            self.connected = True
            
            if self.logger:
                self.logger.success("✅ FTP连接成功！")
            else:
                print("✅ FTP连接成功！")
            
            return True
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"❌ FTP连接失败: {e}")
            else:
                print(f"❌ FTP连接失败: {e}")
            
            self.connected = False
            return False
    
    def disconnect(self):
        """断开FTP连接 | Disconnect from FTP server"""
        if self.ftp:
            try:
                self.ftp.quit()
                if self.logger:
                    self.logger.info("🔌 FTP连接已断开")
                else:
                    print("🔌 FTP连接已断开")
            except:
                pass
            finally:
                self.connected = False
    
    def list_directory(self, path: str) -> List[str]:
        """列出目录内容 | List directory contents"""
        if not self.connected:
            return []
        
        try:
            return self.ftp.nlst(path)
        except ftplib.error_perm:
            return []
    
    def change_directory(self, path: str) -> bool:
        """切换目录 | Change directory"""
        if not self.connected:
            return False
        
        try:
            self.ftp.cwd(path)
            return True
        except ftplib.error_perm:
            return False
    
    def get_file_size(self, file_path: str) -> Optional[int]:
        """获取文件大小 | Get file size"""
        if not self.connected:
            return None
        
        try:
            return self.ftp.size(file_path)
        except ftplib.error_perm:
            return None
    
    def __enter__(self):
        """上下文管理器入口 | Context manager entry"""
        self.connect()
        return self
    
    def __exit__(self, exc_type, exc_val, exc_tb):
        """上下文管理器出口 | Context manager exit"""
        self.disconnect()
