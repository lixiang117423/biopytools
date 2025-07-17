"""
Augustus多转录组预测工具函数模块 | Augustus Multiple RNA-seq Prediction Utility Functions Module
"""

import os
import logging
import subprocess
import sys
from pathlib import Path

class AugustusLogger:
    """Augustus预测日志管理器 | Augustus Prediction Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "augustus_multi_rnaseq.log"):
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
    
    def run(self, cmd: str, description: str = "", capture_output: bool = True) -> tuple:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=capture_output, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            if capture_output and result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True, result.stdout if capture_output else ""
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            if capture_output:
                self.logger.error(f"错误信息 | Error message: {e.stderr}")
            return False, e.stderr if capture_output else ""

class FileValidator:
    """文件验证器 | File Validator"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def check_file_exists(self, file_path: str, description: str = "") -> bool:
        """检查文件是否存在 | Check if file exists"""
        if os.path.exists(file_path):
            if description:
                self.logger.info(f"✓ {description}已存在 | already exists: {file_path}")
            return True
        return False
    
    def check_tool_availability(self, tool_name: str) -> bool:
        """检查工具是否可用 | Check if tool is available"""
        try:
            subprocess.run([tool_name, "--help"], 
                         capture_output=True, 
                         check=True)
            return True
        except (subprocess.CalledProcessError, FileNotFoundError):
            return False

def check_dependencies():
    """检查依赖工具 | Check dependencies"""
    # 核心必需工具 | Core required tools
    core_required_tools = [
        'hisat2', 'hisat2-build', 'samtools', 'augustus', 'bedtools'
    ]
    
    # Augustus相关工具 (可选，但推荐) | Augustus tools (optional but recommended)
    augustus_tools = [
        'filterBam', 'bam2hints'
    ]
    
    missing_core = []
    missing_augustus = []
    
    def find_tool_in_paths(tool_name):
        """在多个路径中查找工具"""
        import os
        
        print(f"正在搜索工具 | Searching for tool: {tool_name}")
        
        # 扩展搜索路径 | Extended search paths
        search_paths = []
        
        # 1. PATH环境变量中的路径
        if 'PATH' in os.environ:
            path_dirs = os.environ['PATH'].split(os.pathsep)
            search_paths.extend(path_dirs)
            print(f"  从PATH添加了 {len(path_dirs)} 个路径 | Added {len(path_dirs)} paths from PATH")
        
        # 2. 常见的安装路径
        common_paths = [
            os.path.expanduser('~/.local/bin'),
            '/usr/local/bin',
            '/usr/bin',
            '/opt/conda/bin',
            '/opt/miniconda3/bin'
        ]
        search_paths.extend(common_paths)
        print(f"  添加常见路径 | Added common paths: {common_paths}")
        
        # 3. conda环境路径
        if 'CONDA_PREFIX' in os.environ:
            conda_bin = os.path.join(os.environ['CONDA_PREFIX'], 'bin')
            search_paths.insert(0, conda_bin)  # 优先搜索当前conda环境
            print(f"  添加conda环境路径 | Added conda env path: {conda_bin}")
        
        # 4. Augustus特定路径
        if 'AUGUSTUS_CONFIG_PATH' in os.environ:
            augustus_dir = os.path.dirname(os.environ['AUGUSTUS_CONFIG_PATH'])
            augustus_bin = os.path.join(augustus_dir, 'bin')
            augustus_scripts = os.path.join(augustus_dir, 'scripts')
            search_paths.extend([augustus_bin, augustus_scripts])
            print(f"  添加Augustus路径 | Added Augustus paths: {augustus_bin}, {augustus_scripts}")
        
        # 5. 确保~/.local/bin被正确处理
        home_local_bin = os.path.expanduser('~/.local/bin')
        if home_local_bin not in search_paths:
            search_paths.append(home_local_bin)
        print(f"  确保包含~/.local/bin | Ensuring ~/.local/bin is included: {home_local_bin}")
        
        # 去重并保持顺序，但保留所有路径用于调试
        seen = set()
        unique_paths = []
        for path in search_paths:
            if path and path not in seen:
                seen.add(path)
                if os.path.isdir(path):
                    unique_paths.append(path)
                    print(f"    有效路径 | Valid path: {path}")
                else:
                    print(f"    无效路径 | Invalid path: {path}")
        
        print(f"  总共搜索 {len(unique_paths)} 个有效路径 | Total {len(unique_paths)} valid paths to search")
        
        # 在各个路径中查找工具
        for i, path in enumerate(unique_paths, 1):
            tool_path = os.path.join(path, tool_name)
            print(f"    [{i}/{len(unique_paths)}] 检查: {tool_path}")
            
            if os.path.isfile(tool_path):
                if os.access(tool_path, os.X_OK):
                    print(f"    ✅ 找到可执行文件 | Found executable: {tool_path}")
                    return tool_path
                else:
                    print(f"    ⚠️  文件存在但不可执行 | File exists but not executable: {tool_path}")
            else:
                print(f"    ❌ 文件不存在 | File not found: {tool_path}")
        
        # 特别检查~/.local/bin（如果还没找到）
        home_local_bin = os.path.expanduser('~/.local/bin')
        if home_local_bin not in unique_paths:
            print(f"  额外检查~/.local/bin | Extra check for ~/.local/bin: {home_local_bin}")
            if os.path.isdir(home_local_bin):
                tool_path = os.path.join(home_local_bin, tool_name)
                print(f"    检查: {tool_path}")
                if os.path.isfile(tool_path) and os.access(tool_path, os.X_OK):
                    print(f"    ✅ 在~/.local/bin中找到 | Found in ~/.local/bin: {tool_path}")
                    return tool_path
        
        print(f"  ❌ 未找到工具 {tool_name} | Tool {tool_name} not found in any path")
        return None
    
    def check_tool(tool_name):
        """检查单个工具是否可用"""
        print(f"\n🔍 检查工具 | Checking tool: {tool_name}")
        
        # 首先尝试直接运行（已在PATH中）
        def test_tool_execution(cmd_list, test_name):
            try:
                print(f"  测试 | Testing: {' '.join(cmd_list)}")
                result = subprocess.run(cmd_list, 
                                      capture_output=True, 
                                      text=True,
                                      timeout=10)
                print(f"    返回码 | Return code: {result.returncode}")
                print(f"    输出长度 | Output length - stdout: {len(result.stdout)}, stderr: {len(result.stderr)}")
                
                # 对于filterBam和bam2hints，只要有输出就认为工具可用，不管退出码
                if tool_name in ['filterBam', 'bam2hints']:
                    if result.stdout or result.stderr:
                        if 'Usage:' in result.stdout or 'Usage:' in result.stderr or 'usage:' in result.stdout.lower() or 'usage:' in result.stderr.lower():
                            print(f"    ✅ {test_name}: 找到使用说明，工具可用")
                            return True
                        elif result.stdout or result.stderr:
                            print(f"    ✅ {test_name}: 有输出，工具可用")
                            return True
                else:
                    # 其他工具的正常检查
                    if result.returncode == 0 or result.stdout or result.stderr:
                        print(f"    ✅ {test_name}: 工具可用")
                        return True
                
                print(f"    ❌ {test_name}: 无有效输出")
                return False
                
            except FileNotFoundError:
                print(f"    ❌ {test_name}: 命令未找到")
                return False
            except subprocess.TimeoutExpired:
                print(f"    ⚠️  {test_name}: 超时，但工具可能存在")
                return True
            except Exception as e:
                print(f"    ⚠️  {test_name}: 异常但工具可能可用 - {e}")
                return True
        
        # 测试不同的帮助参数
        help_commands = []
        if tool_name == 'filterBam':
            help_commands = [
                [tool_name, '-h'],
                [tool_name, '--help']
            ]
        elif tool_name == 'bam2hints':
            help_commands = [
                [tool_name, '--help'],
                [tool_name, '-h']
            ]
        else:
            help_commands = [
                [tool_name, '--help'],
                [tool_name, '-h'],
                [tool_name, '--version'],
                [tool_name, '-v']
            ]
        
        # 尝试直接运行
        for cmd in help_commands:
            if test_tool_execution(cmd, "直接运行"):
                return True
        
        # 如果直接运行失败，尝试在路径中查找
        print(f"  直接运行失败，搜索路径...")
        tool_path = find_tool_in_paths(tool_name)
        if tool_path:
            print(f"  找到工具路径: {tool_path}")
            for cmd in help_commands:
                cmd_with_path = [tool_path] + cmd[1:]  # 替换命令名为完整路径
                if test_tool_execution(cmd_with_path, f"完整路径运行"):
                    return True
        
        print(f"  ❌ 所有测试失败，工具 {tool_name} 不可用")
        return False
    
    # 检查核心工具 | Check core tools
    for tool in core_required_tools:
        if not check_tool(tool):
            missing_core.append(tool)
    
    # 检查Augustus工具 | Check Augustus tools
    for tool in augustus_tools:
        if not check_tool(tool):
            missing_augustus.append(tool)
    
    if missing_core:
        print(f"错误 | Error: 缺少核心工具 | Missing core tools: {', '.join(missing_core)}")
        print("请安装缺少的工具后再运行 | Please install missing tools before running")
        sys.exit(1)
    
    if missing_augustus:
        print(f"警告 | Warning: 缺少Augustus工具 | Missing Augustus tools: {', '.join(missing_augustus)}")
        print("这可能影响分析质量，建议安装：| This may affect analysis quality, recommend installing:")
        print("conda install -c bioconda bamtools")
        print("\n或者重新安装完整的Augustus包 | Or reinstall complete Augustus package:")
        print("conda install -c bioconda augustus")
        print("\n脚本将以降级模式运行（跳过BAM过滤和hints生成）| Script will run in degraded mode (skip BAM filtering and hints generation)")
        return False  # 返回False表示有工具缺失但可以继续
    
    return True  # 返回True表示所有工具都可用
