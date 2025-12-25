"""
文件处理模块 | File Handler Module
"""

from collections import defaultdict
from pathlib import Path
from typing import Dict, List, Tuple, Optional

# 注意: 在模块化模式下，提取脚本会添加以下导入语句:
# from .utils import Logger

class FileHandler:
    """文件处理器 | File Handler"""
    
    def __init__(self, logger=None):
        self.logger = logger
    
    def read_and_group_by_project(self, input_file: str) -> Optional[Dict[str, List[str]]]:
        """读取并按项目分组 | Read and group by project"""
        projects = defaultdict(list)
        
        try:
            with open(input_file, 'r', encoding='utf-8') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line or line.startswith('#'):
                        continue
                    
                    parts = line.split('\t')
                    if len(parts) == 2:
                        project_id, run_id = parts
                        projects[project_id].append(run_id)
                    else:
                        if self.logger:
                            self.logger.warning(f"第{line_num}行格式不正确，已忽略: '{line}'")
                        else:
                            print(f"⚠️ 第{line_num}行格式不正确，已忽略: '{line}'")
            
            # 对每个项目的ID列表进行排序
            for project_id in projects:
                projects[project_id].sort()
            
            if self.logger:
                self.logger.info(f"📄 成功读取 {len(projects)} 个项目")
            else:
                print(f"📄 成功读取 {len(projects)} 个项目")
            
            return dict(projects)
            
        except FileNotFoundError:
            if self.logger:
                self.logger.error(f"输入文件未找到: {input_file}")
            else:
                print(f"❌ 输入文件未找到: {input_file}")
            return None
        except Exception as e:
            if self.logger:
                self.logger.error(f"读取文件时发生错误: {e}")
            else:
                print(f"❌ 读取文件时发生错误: {e}")
            return None
    
    def write_links(self, output_file: str, urls: List[str]) -> bool:
        """写入链接文件 | Write links file"""
        try:
            output_path = Path(output_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_path, 'w', encoding='utf-8') as f:
                for url in sorted(urls):
                    f.write(url + '\n')
            
            if self.logger:
                self.logger.success(f"💾 已保存 {len(urls)} 个链接到: {output_file}")
            else:
                print(f"✅ 已保存 {len(urls)} 个链接到: {output_file}")
            
            return True
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"写入链接文件失败: {e}")
            else:
                print(f"❌ 写入链接文件失败: {e}")
            return False
    
    def write_failed_ids(self, failed_file: str, failed_ids: List[Tuple[str, str]]) -> bool:
        """写入失败的ID文件 | Write failed IDs file"""
        try:
            output_path = Path(failed_file)
            output_path.parent.mkdir(parents=True, exist_ok=True)
            
            with open(output_path, 'w', encoding='utf-8') as f:
                # 按项目排序，再按run_id排序
                for project, run_id in sorted(failed_ids):
                    f.write(f"{project}\t{run_id}\n")
            
            if self.logger:
                self.logger.success(f"📝 已保存 {len(failed_ids)} 个失败ID到: {failed_file}")
            else:
                print(f"✅ 已保存 {len(failed_ids)} 个失败ID到: {failed_file}")
            
            return True
            
        except Exception as e:
            if self.logger:
                self.logger.error(f"写入失败ID文件失败: {e}")
            else:
                print(f"❌ 写入失败ID文件失败: {e}")
            return False
