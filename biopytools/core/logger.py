# biopytools/core/logger.py

import logging
import sys
from pathlib import Path

# 定义一个全局图标字典
LOG_ICONS = {
    "start": "🚀", "finish": "✅", "process": "🔄", "data": "📁",
    "plot": "📊", "stats": "📈", "config": "⚙️", "run": "⚡",
    "success": "✅", "warning": "⚠️", "error": "❌", "debug": "🐞",
    "info": "ℹ️", "param": "📋", "step": "➡️", "save": "📄",
    "find": "🔍", "clean": "🧹", "star": "⭐", "pin": "📍", "summary": "📋"
}

class IconLogger(logging.Logger):
    """一个可以添加图标的自定义Logger类"""
    def __init__(self, name, level=logging.NOTSET):
        super().__init__(name, level)

    def log_with_icon(self, level, icon_key, msg, *args, **kwargs):
        """记录带图标的日志"""
        icon = LOG_ICONS.get(icon_key, "➡️")
        super().log(level, f"{icon} {msg}", *args, **kwargs)
    
    # 为常用日志级别创建便捷方法
    def start(self, msg, *args, **kwargs): self.log_with_icon(logging.INFO, "start", msg, *args, **kwargs)
    def success(self, msg, *args, **kwargs): self.log_with_icon(logging.INFO, "success", msg, *args, **kwargs)
    def icon_info(self, icon_key, msg, *args, **kwargs): self.log_with_icon(logging.INFO, icon_key, msg, *args, **kwargs)
    def warning(self, msg, *args, **kwargs): self.log_with_icon(logging.WARNING, "warning", msg, *args, **kwargs)
    def error(self, msg, *args, **kwargs): self.log_with_icon(logging.ERROR, "error", msg, *args, **kwargs)


def setup_logger(output_dir: Path, log_name: str = "analysis.log") -> IconLogger:
    """设置并返回一个IconLogger实例"""
    logging.setLoggerClass(IconLogger)
    
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    log_file = output_dir / log_name

    # 清理旧的handlers，防止重复打印日志
    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)

    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(levelname)s - %(message)s',
        handlers=[
            logging.FileHandler(str(log_file), encoding='utf-8'),
            logging.StreamHandler(sys.stdout)
        ],
        force=True
    )
    
    logger = logging.getLogger(__name__)
    return logger