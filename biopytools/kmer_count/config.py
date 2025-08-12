"""
⚙️ 配置管理 | Configuration Management
"""

import tempfile
from pathlib import Path

class KmerCountConfig:
    """🔧 K-mer计数配置类 | K-mer counting configuration class"""
    
    def __init__(self, args):
        # 输入文件 | Input files
        self.input_dir = Path(args.input)
        self.pattern = args.pattern
        self.kmer_lib = Path(args.kmer_lib)
        self.bed_file = Path(args.bed_file) if args.bed_file else None
        
        # 输出设置 | Output settings
        self.output_dir = Path(args.output)
        self.output_dir.mkdir(parents=True, exist_ok=True)
        
        # Jellyfish参数 | Jellyfish parameters
        self.kmer_size = args.kmer_size
        self.hash_size = args.hash_size
        self.threads = args.threads
        self.canonical = args.canonical
        
        # 滑动窗口参数 | Sliding window parameters
        self.window_size = args.window_size
        self.step_size = args.step_size if args.step_size else args.window_size // 5
        
        # 其他参数 | Other parameters
        # self.keep_temp = args.keep_temp
        self.keep_temp = args.keep_temp if hasattr(args, 'keep_temp') else True  # 默认保留临时文件
        self.keep_binary = args.keep_binary
        self.verbose = args.verbose
        
        # 工具路径 | Tool paths
        self.jellyfish_path = args.jellyfish_path
        
        # 内部变量 | Internal variables
        self.temp_dir = None
        self.samples = []
        
    def setup_temp_dir(self):
        """📁 设置临时目录 | Setup temporary directory"""
        # self.temp_dir = Path(tempfile.mkdtemp(prefix="kmer_count_"))
        import time
        import random
        temp_name = f"kmer_count_{int(time.time())}_{random.randint(1000,9999)}"
        self.temp_dir = Path.cwd() / "tmp" / temp_name
        self.temp_dir.mkdir(parents=True, exist_ok=True)
        return self.temp_dir
