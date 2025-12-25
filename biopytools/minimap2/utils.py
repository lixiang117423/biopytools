"""
Minimap2分析工具函数模块 | Minimap2 Analysis Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class Minimap2Logger:
    """Minimap2分析日志管理器 | Minimap2 Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "minimap2_analysis.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False

class IntervalProcessor:
    """区间处理器 | Interval Processor"""
    
    def __init__(self, logger, config=None):
        self.logger = logger
        self.config = config
    
    def subtract_interval_from_list(self, unmapped_intervals, mapped_start, mapped_end):
        """从未比对区间列表中砍掉一个已比对区间 | Remove a mapped interval from unmapped interval list"""
        new_intervals = []
        
        for start, end in unmapped_intervals:
            # 检查是否有重叠 | Check for overlap
            if mapped_end < start or mapped_start > end:
                # 无重叠，保留原区间 | No overlap, keep original interval
                new_intervals.append((start, end))
            else:
                # 有重叠，需要"砍掉"重叠部分 | Has overlap, need to "cut out" the overlapping part
                
                # 保留左侧未重叠部分 | Keep left non-overlapping part
                if start < mapped_start:
                    new_intervals.append((start, mapped_start - 1))
                
                # 保留右侧未重叠部分 | Keep right non-overlapping part  
                if end > mapped_end:
                    new_intervals.append((mapped_end + 1, end))
                
                # 中间重叠的部分被"砍掉"了，不保留 | The overlapping middle part is "cut out", not kept
        
        return new_intervals
    
    def find_unmapped_intervals_by_subtraction(self, df_filtered):
        """使用减法找到未比对区间 | Find unmapped intervals using subtraction method"""
        self.logger.info("使用减法策略查找未比对区间 | Finding unmapped regions using subtraction strategy")
        
        # 检查必要的列是否存在 | Check if necessary columns exist
        required_columns = ['query_name', 'query_len', 'query_start', 'query_end', 'number_match']
        missing_columns = [col for col in required_columns if col not in df_filtered.columns]
        
        if missing_columns:
            self.logger.error(f"缺少必要的列 | Missing required columns: {missing_columns}")
            self.logger.error(f"可用的列 | Available columns: {list(df_filtered.columns)}")
            return []
        
        # 按query_name和number_match排序（number_match从大到小） | Sort by query_name and number_match (descending)
        df_sorted = df_filtered.sort_values(['query_name', 'number_match'], ascending=[True, False])
        self.logger.info(f"排序后的记录数 | Records after sorting: {len(df_sorted)}")
        
        # 按query_name分组处理 | Process by query_name groups
        all_unmapped_regions = []
        
        for query_name, group in df_sorted.groupby('query_name'):
            query_len = int(group['query_len'].iloc[0])
            self.logger.info(f"\n处理序列 | Processing sequence: {query_name} (长度 | length: {query_len:,})")
            
            # 初始化未比对区间为整个序列 | Initialize unmapped intervals as full sequence
            unmapped_intervals = [(1, query_len)]
            self.logger.info(f"初始未比对区间 | Initial unmapped intervals: {unmapped_intervals}")
            
            # 依次砍掉每个比对区间 | Cut out each mapped interval
            for idx, (_, row) in enumerate(group.iterrows()):
                # 从PAF(0-base, half-open)转换为内部使用的1-base, closed
                # [start, end) -> [start + 1, end]
                mapped_start = int(row['query_start']) + 1  # 转换为1-based
                mapped_end = int(row['query_end'])
                number_match = int(row['number_match'])
                
                # 确保start < end | Ensure start < end
                if mapped_start > mapped_end:
                    mapped_start, mapped_end = mapped_end, mapped_start
                
                self.logger.info(f"  砍掉区间 #{idx+1} | Cutting out interval #{idx+1}: {mapped_start:,}-{mapped_end:,} (match: {number_match:,})")
                
                # 从未比对区间中砍掉当前比对区间 | Cut out current mapped interval from unmapped intervals
                old_intervals = unmapped_intervals[:]
                unmapped_intervals = self.subtract_interval_from_list(
                    unmapped_intervals, mapped_start, mapped_end
                )
                
                self.logger.info(f"    砍掉前 | Before cutting: {old_intervals}")
                self.logger.info(f"    砍掉后 | After cutting: {unmapped_intervals}")
                
                # 如果没有剩余区间，提前结束 | If no remaining intervals, break early
                if not unmapped_intervals:
                    self.logger.info(f"    所有区间都被砍掉了 | All intervals have been cut out")
                    break
            
            # 筛选长度足够的区间 | Filter intervals with sufficient length
            for start, end in unmapped_intervals:
                length = end - start + 1
                if length >= self.config.min_unmapped_length:
                    all_unmapped_regions.append({
                        'query_name': query_name,
                        'start': start,
                        'end': end,
                        'length': length
                    })
                    self.logger.info(f"  保留未比对区间 | Keeping unmapped interval: {start:,}-{end:,} (长度 | length: {length:,})")
                else:
                    self.logger.info(f"  过滤短区间 | Filtering short interval: {start:,}-{end:,} (长度 | length: {length:,} < {self.config.min_unmapped_length:,})")
        
        self.logger.info(f"\n总计找到 {len(all_unmapped_regions)} 个符合条件的未比对区间")
        self.logger.info(f"Total found {len(all_unmapped_regions)} qualified unmapped regions")
        
        return all_unmapped_regions
