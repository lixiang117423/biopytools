"""
CNCB FTP搜索器模块 | CNCB FTP Searcher Module
"""

import ftplib
import time
from typing import Optional, List, Tuple
import logging


class CNCBFTPSearcher:
    """CNCB FTP搜索器 | CNCB FTP Searcher"""

    def __init__(self, ftp_client, base_dirs: List[str], path_cache):
        self.ftp = ftp_client
        self.base_dirs = base_dirs
        self.path_cache = path_cache

    def find_base_path_for_id(self, archive: str, run_id: str) -> Optional[str]:
        """
        为Run ID查找基础路径 | Find base path for Run ID

        Args:
            archive: 数据库类型 (SRA/ERA/DRA) | Database type
            run_id: Run ID | Run ID

        Returns:
            基础路径 | Base path, or None if not found
        """
        # 生成候选父目录 | Generate candidate parent directories
        candidate_parents = []
        for i in range(len(run_id) - 2, 5, -1):
            candidate_parents.append(run_id[:i])

        # 首先检查缓存 | Check cache first
        for parent in candidate_parents:
            cached_path = self.path_cache.get(parent)
            if cached_path:
                logging.debug(f"🎯 缓存命中 | Cache hit for {parent}: {cached_path}")
                return cached_path

        logging.info(f"  🔍 正在为 {run_id} 进行全面寻址...")

        # 搜索FTP服务器 | Search FTP server
        for parent_pattern in candidate_parents:
            for base_dir in self.base_dirs:
                level1_path_base = f"/{base_dir}/{archive}"

                try:
                    # 列出一级目录
                    level1_dirs_from_ftp = self.ftp.nlst(level1_path_base)
                except ftplib.error_perm:
                    logging.debug(f"  ⏭️ 无法访问目录 | Cannot access directory: {level1_path_base}")
                    continue

                # 提取目录名 | Extract directory names
                level1_dirs_to_try = [d.split('/')[-1] for d in level1_dirs_from_ftp]

                # 搜索匹配的路径 | Search for matching paths
                for l1_dir_name in level1_dirs_to_try:
                    path_to_verify = f"{level1_path_base}/{l1_dir_name}/{parent_pattern}/{run_id}"

                    try:
                        # 验证路径是否存在 | Verify if path exists
                        self.ftp.cwd(path_to_verify)
                        correct_base_path = f"{level1_path_base}/{l1_dir_name}"

                        logging.info(f"  👍 路径找到并缓存! 父目录 '{parent_pattern}' 位于: {correct_base_path}")
                        self.path_cache.set(parent_pattern, correct_base_path)
                        return correct_base_path

                    except ftplib.error_perm:
                        continue

        logging.warning(f"  😥 未能为 {run_id} 找到路径")
        return None

    def find_files_for_run(self, run_id: str, base_path: str,
                          filename_templates: List[str]) -> List[str]:
        """
        为Run ID查找文件 | Find files for Run ID

        Args:
            run_id: Run ID | Run ID
            base_path: 基础路径 | Base path
            filename_templates: 文件名模板列表 | Filename template list

        Returns:
            找到的文件URL列表 | List of found file URLs
        """
        found_urls = []

        # 确定父目录名 | Determine parent directory name
        parent_dir_name = ""
        for i in range(len(run_id) - 2, 5, -1):
            key = run_id[:i]
            cached_path = self.path_cache.get(key)
            if cached_path and cached_path == base_path:
                parent_dir_name = key
                break

        if not parent_dir_name:
            logging.error(f"  ❌ 无法确定 {run_id} 的父目录")
            return found_urls

        final_dir_path = f"{base_path}/{parent_dir_name}/{run_id}"
        logging.debug(f"  🔍 搜索目录: {final_dir_path}")

        # 查找文件 | Search for files
        for template in filename_templates:
            filename = template.format(run_id=run_id)
            file_path = f"{final_dir_path}/{filename}"

            try:
                # 检查文件是否存在 | Check if file exists
                file_size = self.ftp.size(file_path)
                ftp_url = f"ftp://download2.cncb.ac.cn{file_path}"
                found_urls.append(ftp_url)
                logging.debug(f"  📄 找到文件: {filename} (大小: {file_size:,} bytes)")

            except ftplib.error_perm:
                logging.debug(f"  ⏭️ 文件不存在: {filename}")
                continue

        return found_urls

    def test_ftp_access(self) -> bool:
        """测试FTP访问权限 | Test FTP access permissions"""
        try:
            # 尝试列出根目录
            self.ftp.nlst()
            return True
        except ftplib.error_perm as e:
            logging.error(f"❌ FTP访问权限测试失败 | FTP access test failed: {e}")
            return False

    def get_server_info(self) -> dict:
        """获取FTP服务器信息 | Get FTP server information"""
        try:
            welcome_msg = self.ftp.getwelcome()
            return {
                "welcome_message": welcome_msg,
                "pwd": self.ftp.pwd(),
                "base_dirs_accessible": self._test_base_dirs()
            }
        except Exception as e:
            logging.error(f"❌ 获取FTP服务器信息失败 | Failed to get FTP server info: {e}")
            return {}

    def _test_base_dirs(self) -> List[str]:
        """测试基础目录可访问性 | Test base directories accessibility"""
        accessible_dirs = []
        for base_dir in self.base_dirs:
            try:
                self.ftp.nlst(f"/{base_dir}")
                accessible_dirs.append(base_dir)
            except ftplib.error_perm:
                logging.debug(f"  ⏭️ 基础目录不可访问 | Base directory not accessible: {base_dir}")
        return accessible_dirs