"""
进度管理模块|Progress Manager Module

用于管理和跟踪分析任务的进度，支持断点续传功能
Manage and track analysis task progress, support checkpoint resumption
"""

import json
from pathlib import Path
from typing import Dict, Optional
from datetime import datetime


class ProgressManager:
    """进度管理器|Progress Manager

    管理分析任务的执行进度，支持断点续传
    Manage analysis task execution progress, support checkpoint resumption
    """

    # 步骤定义|Step definitions
    STEPS = {
        'check_dependencies': '步骤1: 检查依赖软件|Step 1: Checking dependencies',
        'setup_temp_dir': '步骤2: 设置临时目录|Step 2: Setting up temporary directory',
        'find_input_files': '步骤3: 查找输入文件|Step 3: Finding input files',
        'prepare_kmer_lib': '步骤4: 准备k-mer库|Step 4: Preparing k-mer library',
        'parse_annotations': '步骤5: 解析注释信息|Step 5: Parsing annotation information',
        'process_samples': '步骤6: 处理样本|Step 6: Processing samples',
        'merge_results': '步骤7: 合并结果|Step 7: Merging results',
        'cleanup': '步骤8: 清理临时文件|Step 8: Cleaning up temporary files'
    }

    def __init__(self, progress_file: Path, logger):
        """初始化进度管理器|Initialize progress manager

        Args:
            progress_file: 进度文件路径|Progress file path
            logger: 日志记录器|Logger instance
        """
        self.progress_file = progress_file
        self.logger = logger
        self.progress_data = self._load_progress()

    def _load_progress(self) -> Dict:
        """加载进度数据|Load progress data"""
        if self.progress_file.exists():
            try:
                with open(self.progress_file, 'r', encoding='utf-8') as f:
                    data = json.load(f)
                self.logger.info(f"发现进度文件，将支持断点续传|Progress file found, checkpoint resumption enabled")
                return data
            except Exception as e:
                self.logger.warning(f"读取进度文件失败，将创建新的进度记录|Failed to read progress file, creating new: {e}")
                return self._init_progress()
        else:
            return self._init_progress()

    def _init_progress(self) -> Dict:
        """初始化进度数据|Initialize progress data"""
        return {
            'start_time': datetime.now().isoformat(),
            'last_update': datetime.now().isoformat(),
            'steps': {},
            'completed_steps': [],
            'current_step': None
        }

    def _save_progress(self):
        """保存进度数据|Save progress data"""
        self.progress_data['last_update'] = datetime.now().isoformat()
        with open(self.progress_file, 'w', encoding='utf-8') as f:
            json.dump(self.progress_data, f, indent=2, ensure_ascii=False)

    def is_step_completed(self, step_name: str) -> bool:
        """检查步骤是否已完成|Check if step is completed

        Args:
            step_name: 步骤名称|Step name

        Returns:
            bool: 是否已完成|Whether completed
        """
        return step_name in self.progress_data.get('completed_steps', [])

    def mark_step_completed(self, step_name: str, metadata: Optional[Dict] = None):
        """标记步骤为已完成|Mark step as completed

        Args:
            step_name: 步骤名称|Step name
            metadata: 额外的元数据|Additional metadata
        """
        if step_name not in self.progress_data.get('completed_steps', []):
            self.progress_data.setdefault('completed_steps', []).append(step_name)

        self.progress_data.setdefault('steps', {})[step_name] = {
            'status': 'completed',
            'completed_at': datetime.now().isoformat(),
            'metadata': metadata or {}
        }
        self.progress_data['current_step'] = None
        self._save_progress()

        step_desc = self.STEPS.get(step_name, step_name)
        self.logger.info(f"{step_desc} 已完成|{step_desc} completed")
        self.logger.info(f"进度已保存|Progress saved to: {self.progress_file}")

    def mark_step_in_progress(self, step_name: str):
        """标记步骤为进行中|Mark step as in progress

        Args:
            step_name: 步骤名称|Step name
        """
        self.progress_data['current_step'] = step_name
        self.progress_data.setdefault('steps', {})[step_name] = {
            'status': 'in_progress',
            'started_at': datetime.now().isoformat()
        }
        self._save_progress()

    def get_step_metadata(self, step_name: str) -> Optional[Dict]:
        """获取步骤元数据|Get step metadata

        Args:
            step_name: 步骤名称|Step name

        Returns:
            元数据字典|Metadata dictionary
        """
        step_info = self.progress_data.get('steps', {}).get(step_name, {})
        return step_info.get('metadata')

    def print_progress_summary(self):
        """打印进度摘要|Print progress summary"""
        completed = self.progress_data.get('completed_steps', [])
        total = len(self.STEPS)
        completed_count = len(completed)

        self.logger.info("=" * 60)
        self.logger.info(f"进度摘要|Progress Summary:")
        self.logger.info(f"  已完成步骤|Completed steps: {completed_count}/{total}")
        self.logger.info(f"  开始时间|Start time: {self.progress_data.get('start_time')}")
        self.logger.info(f"  最后更新|Last update: {self.progress_data.get('last_update')}")

        if completed_count > 0:
            self.logger.info(f"  已完成的步骤|Completed step names:")
            for step in completed:
                step_desc = self.STEPS.get(step, step)
                self.logger.info(f"     {step_desc}")

        remaining = total - completed_count
        if remaining > 0:
            self.logger.info(f"  剩余步骤|Remaining steps: {remaining}")
        else:
            self.logger.info(f"  所有步骤已完成|All steps completed")

        self.logger.info("=" * 60)

    def reset_progress(self):
        """重置进度|Reset progress"""
        if self.progress_file.exists():
            self.progress_file.unlink()
            self.logger.info(f"进度文件已删除|Progress file deleted: {self.progress_file}")
        self.progress_data = self._init_progress()
