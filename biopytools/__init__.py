"""
biopytools - 生物信息学工具包
Bioinformatics tools package

这是一个用于生物信息学分析的Python工具包，包含多个常用的生物信息学工具模块。
This is a Python toolkit for bioinformatics analysis, containing multiple commonly used bioinformatics tool modules.

主要模块 | Main modules:
    - fastp: 高速FASTQ质控工具 | High-speed FASTQ quality control tool
    - 更多模块正在开发中... | More modules are under development...

作者 | Author: biopytools team
许可证 | License: MIT
"""

# 导入版本信息
from ._version import __version__, __version_info__, get_version_string, get_build_info

# 导入主要模块
try:
    from . import fastp
    _available_modules = ['fastp']
except ImportError:
    _available_modules = []

# 公共接口
__all__ = ['__version__', '__version_info__', 'get_version_string', 'get_build_info', 
           'show_version', 'show_info', 'list_modules'] + _available_modules

def show_version(detailed: bool = False):
    """
    显示版本信息
    Show version information
    
    Args:
        detailed: 是否显示详细信息 | Whether to show detailed information
    """
    if detailed:
        build_info = get_build_info()
        print(f"biopytools {build_info['version']}")
        print(f"状态 | Status: {build_info['status']}")
        print(f"Git版本 | Git revision: {build_info['git_revision']}")
        print(f"构建日期 | Build date: {build_info['build_date']}")
        print(f"Python版本 | Python version: {build_info['python_version']}")
        print(f"Python实现 | Python implementation: {build_info['python_implementation']}")
        print(f"平台 | Platform: {build_info['platform']}")
    else:
        print(f"biopytools {get_version_string()}")

def show_info():
    """显示包的基本信息 | Show basic package information"""
    print("biopytools - 生物信息学工具包 | Bioinformatics Tools Package")
    print(f"版本 | Version: {get_version_string()}")
    print(f"可用模块 | Available modules: {', '.join(_available_modules) if _available_modules else 'None'}")
    print("许可证 | License: MIT")
    print("项目主页 | Homepage: https://github.com/yourusername/biopytools")

def list_modules():
    """列出所有可用的模块 | List all available modules"""
    print("可用模块 | Available modules:")
    if not _available_modules:
        print("  (暂无可用模块 | No modules available)")
        return
    
    for module_name in _available_modules:
        try:
            module = globals()[module_name]
            if hasattr(module, '__doc__') and module.__doc__:
                description = module.__doc__.split('\n')[0].strip()
            else:
                description = "无描述 | No description"
            print(f"  - {module_name}: {description}")
        except (KeyError, AttributeError):
            print(f"  - {module_name}: 无法加载模块信息 | Cannot load module info")

def check_dependencies():
    """检查依赖包是否安装 | Check if dependencies are installed"""
    dependencies = {
        'pandas': 'pandas',
        'numpy': 'numpy'
    }
    
    missing_deps = []
    installed_deps = []
    
    for dep_name, import_name in dependencies.items():
        try:
            __import__(import_name)
            installed_deps.append(dep_name)
        except ImportError:
            missing_deps.append(dep_name)
    
    print("依赖检查 | Dependency Check:")
    if installed_deps:
        print(f"  ✓ 已安装 | Installed: {', '.join(installed_deps)}")
    
    if missing_deps:
        print(f"  ✗ 缺失 | Missing: {', '.join(missing_deps)}")
        print("  请运行 | Please run: pip install " + " ".join(missing_deps))
        return False
    else:
        print("  ✓ 所有依赖都已安装 | All dependencies are installed")
        return True

# 兼容性检查
def check_python_version():
    """检查Python版本兼容性 | Check Python version compatibility"""
    import sys
    
    min_version = (3, 8)
    current_version = sys.version_info[:2]
    
    if current_version < min_version:
        print(f"警告：Python版本过低 | Warning: Python version too old")
        print(f"当前版本 | Current: {sys.version}")
        print(f"最低要求 | Minimum required: {'.'.join(map(str, min_version))}")
        return False
    
    return True

# 包级别的配置和初始化
def _initialize_package():
    """初始化包 | Initialize package"""
    # 检查Python版本
    if not check_python_version():
        import warnings
        warnings.warn(
            "biopytools 可能无法在当前Python版本下正常工作 | "
            "biopytools may not work properly with current Python version",
            UserWarning
        )

# 执行包初始化
_initialize_package()