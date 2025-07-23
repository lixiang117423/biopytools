"""
VCF LD热图分析工具函数模块 | VCF LD Heatmap Analysis Utility Functions Module
"""

import logging
import sys
from pathlib import Path

class LDLogger:
    """LD分析日志管理器 | LD Analysis Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "ld_analysis.log"):
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

def check_dependencies(logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    missing_deps = []
    
    # 检查Python包 | Check Python packages
    required_packages = [
        ('numpy', 'NumPy'),
        ('pandas', 'Pandas'),
        ('matplotlib', 'Matplotlib'),
        ('seaborn', 'Seaborn'),
        ('scipy', 'SciPy'),
        ('allel', 'scikit-allel')
    ]
    
    for package, name in required_packages:
        try:
            __import__(package)
            logger.info(f"✓ {name} 可用 | available")
        except ImportError:
            missing_deps.append(name)
            logger.error(f"✗ {name} 不可用 | not available")
    
    if missing_deps:
        error_msg = f"缺少依赖包 | Missing dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        logger.error("请安装缺少的包 | Please install missing packages:")
        logger.error("pip install numpy pandas matplotlib seaborn scipy scikit-allel")
        
        # 特别提示 scikit-allel 的安装
        if 'scikit-allel' in missing_deps:
            logger.error("\n特别注意 | Special note:")
            logger.error("scikit-allel 是VCF文件处理的核心依赖 | scikit-allel is a core dependency for VCF file processing")
            logger.error("如果安装失败，请尝试 | If installation fails, try:")
            logger.error("conda install -c conda-forge scikit-allel")
            logger.error("或者 | or:")
            logger.error("pip install --upgrade pip setuptools wheel")
            logger.error("pip install scikit-allel")
        
        raise RuntimeError(error_msg)
    
    return True
