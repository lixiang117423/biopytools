"""
链接查找模块 | Link Finder Module
"""

from typing import Dict, List, Optional, Tuple
from collections import defaultdict

# 注意: 在模块化模式下，提取脚本会添加以下导入语句:
# from .utils import Logger
# from .ftp_client import FTPClient
# from .config import Config

class LinkFinder:
    """链接查找器 | Link Finder"""
    
    def __init__(self, config=None, logger=None):
        self.config = config or Config()
        self.logger = logger
        self.path_cache = {}  # 路径缓存
        self.ftp_client = None
    
    def find_base_path_for_id(self, ftp, archive: str, run_id: str) -> Optional[str]:
        """为ID查找基础路径 | Find base path for ID"""
        candidate_parents = []
        for i in range(len(run_id) - 2, 5, -1):
            candidate_parents.append(run_id[:i])
        
        # 检查缓存
        for parent in candidate_parents:
            if parent in self.path_cache:
                return self.path_cache[parent]
        
        if self.logger:
            self.logger.debug(f"🔍 正在为 {run_id} 进行全面寻址...")
        else:
            print(f"🔍 正在为 {run_id} 进行全面寻址...")
        
        for parent_pattern in candidate_parents:
            for base_dir in self.config.base_dirs:
                level1_path_base = f"/{base_dir}/{archive}"
                
                # 获取一级目录列表
                level1_dirs_from_ftp = ftp.list_directory(level1_path_base)
                level1_dirs_to_try = [d.split('/')[-1] for d in level1_dirs_from_ftp]
                
                # 尝试每个一级目录
                for l1_dir_name in level1_dirs_to_try:
                    path_to_verify = f"{level1_path_base}/{l1_dir_name}/{parent_pattern}/{run_id}"
                    
                    if ftp.change_directory(path_to_verify):
                        correct_base_path = f"{level1_path_base}/{l1_dir_name}"
                        
                        if self.logger:
                            self.logger.debug(f"👍 路径找到并缓存! 父目录 '{parent_pattern}' 位于: {correct_base_path}")
                        else:
                            print(f"👍 路径找到并缓存! 父目录 '{parent_pattern}' 位于: {correct_base_path}")
                        
                        self.path_cache[parent_pattern] = correct_base_path
                        return correct_base_path
        
        if self.logger:
            self.logger.debug(f"😥 未能为 {run_id} 找到路径")
        else:
            print(f"😥 未能为 {run_id} 找到路径")
        
        return None
    
    def find_download_links(self, run_id: str, base_path: str) -> List[str]:
        """查找下载链接 | Find download links"""
        # 获取父目录名
        parent_dir_name = ""
        for i in range(len(run_id) - 2, 5, -1):
            key = run_id[:i]
            if key in self.path_cache and self.path_cache[key] == base_path:
                parent_dir_name = key
                break
        
        if not parent_dir_name:
            return []
        
        final_dir_path = f"{base_path}/{parent_dir_name}/{run_id}"
        found_urls = []
        
        # 尝试不同的文件名模板
        for template in self.config.filename_templates:
            filename = template.format(run_id=run_id)
            file_path = f"{final_dir_path}/{filename}"
            
            if self.ftp_client.get_file_size(file_path) is not None:
                full_url = f"ftp://{self.config.ftp_host}{file_path}"
                found_urls.append(full_url)
        
        return found_urls
    
    def search_project_ids(self, ftp, projects: Dict[str, List[str]]) -> Tuple[List[str], List[Tuple[str, str]]]:
        """搜索项目ID | Search project IDs"""
        all_found_urls = []
        all_failed_ids = []
        
        project_count = 0
        for project_id, id_list in projects.items():
            project_count += 1
            
            if self.logger:
                self.logger.info(f"📁 --- [{project_count}/{len(projects)}] 正在处理项目: {project_id} ({len(id_list)} 个ID) ---")
            else:
                print(f"📁 --- [{project_count}/{len(projects)}] 正在处理项目: {project_id} ({len(id_list)} 个ID) ---")
            
            for run_id in id_list:
                prefix = run_id[:3]
                archive = self.config.archive_map.get(prefix)
                
                if not archive:
                    if self.logger:
                        self.logger.warning(f"ID {run_id} 的前缀未知，跳过")
                    else:
                        print(f"⚠️ ID {run_id} 的前缀未知，跳过")
                    
                    all_failed_ids.append((project_id, run_id))
                    continue
                
                # 查找基础路径
                base_path = self.find_base_path_for_id(ftp, archive, run_id)
                
                if not base_path:
                    all_failed_ids.append((project_id, run_id))
                    continue
                
                # 查找下载链接
                found_urls = self.find_download_links(run_id, base_path)
                
                if found_urls:
                    all_found_urls.extend(found_urls)
                else:
                    all_failed_ids.append((project_id, run_id))
        
        return all_found_urls, all_failed_ids
