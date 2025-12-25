"""
HiFiasm配置管理模块 | HiFiasm Configuration Management Module
"""

import os
import shutil
from pathlib import Path
from typing import Optional, List, Union
from dataclasses import dataclass, field

class ConfigurationError(Exception):
    """配置错误异常 | Configuration error exception"""
    pass

@dataclass
class HifiasmConfig:
    """HiFiasm分析配置类 | HiFiasm Analysis Configuration Class"""
    
    # ===== 必需参数 | Required parameters =====
    input_reads: str
    
    # ===== 基本参数 | Basic parameters =====
    output_dir: str = './hifiasm_output'
    prefix: str = 'sample'
    threads: int = 32
    
    # ===== HiFiasm组装参数 | HiFiasm assembly parameters =====
    hg_size: str = 'auto'
    purge_level: int = 3
    purge_max: int = 65
    similarity_threshold: float = 0.75
    ont_reads: Optional[str] = None
    hi_c_1: Optional[str] = None
    hi_c_2: Optional[str] = None
    extra_hifiasm_args: str = ''
    
    # ===== 质量评估参数 | Quality assessment parameters =====
    skip_busco: bool = False
    busco_lineage: str = 'auto'
    busco_mode: str = 'genome'
    skip_quast: bool = False
    reference_genome: Optional[str] = None
    
    # ===== 分析参数 | Analysis parameters =====
    analyze_haplotypes: bool = False
    min_contig_length: int = 1000
    generate_plots: bool = False
    assembly_type: str = 'auto'
    
    # ===== 输出控制参数 | Output control parameters =====
    keep_intermediate: bool = False
    compress_output: bool = False
    output_formats: List[str] = field(default_factory=lambda: ['both'])
    
    # ===== 系统参数 | System parameters =====
    memory: int = 64
    tmp_dir: str = '/tmp'
    max_runtime: int = 48
    resume: bool = False
    
    # ===== 工具路径参数 | Tool paths parameters =====
    hifiasm_path: str = 'hifiasm'
    busco_path: str = 'busco'
    quast_path: str = 'quast'
    python_path: str = 'python3'
    samtools_path: str = 'samtools'
    
    # ===== 数据库路径参数 | Database paths parameters =====
    busco_db_path: Optional[str] = None
    busco_download_path: Optional[str] = None
    
    # ===== 高级参数 | Advanced parameters =====
    debug: bool = False
    verbose: int = 0
    log_level: str = 'INFO'
    config_file: Optional[str] = None
    dry_run: bool = False
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        # 转换路径为绝对路径
        self.input_reads = os.path.abspath(self.input_reads)
        self.output_dir = os.path.abspath(self.output_dir)
        
        if self.ont_reads:
            self.ont_reads = os.path.abspath(self.ont_reads)
        if self.hi_c_1:
            self.hi_c_1 = os.path.abspath(self.hi_c_1)
        if self.hi_c_2:
            self.hi_c_2 = os.path.abspath(self.hi_c_2)
        if self.reference_genome:
            self.reference_genome = os.path.abspath(self.reference_genome)
        
        # 创建输出目录
        os.makedirs(self.output_dir, exist_ok=True)
        
        # 设置工作目录路径
        self.working_dir = self.output_dir
        self.log_dir = os.path.join(self.output_dir, 'logs')
        self.tmp_work_dir = os.path.join(self.output_dir, 'tmp')
        
        # 创建子目录
        for dir_path in [self.log_dir, self.tmp_work_dir]:
            os.makedirs(dir_path, exist_ok=True)
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # ===== 验证输入文件 | Validate input files =====
        if not os.path.exists(self.input_reads):
            errors.append(f"输入HiFi数据文件不存在 | Input HiFi data file does not exist: {self.input_reads}")
        
        if self.ont_reads and not os.path.exists(self.ont_reads):
            errors.append(f"ONT数据文件不存在 | ONT data file does not exist: {self.ont_reads}")
        
        if self.hi_c_1 and not os.path.exists(self.hi_c_1):
            errors.append(f"Hi-C第一端文件不存在 | Hi-C first-end file does not exist: {self.hi_c_1}")
        
        if self.hi_c_2 and not os.path.exists(self.hi_c_2):
            errors.append(f"Hi-C第二端文件不存在 | Hi-C second-end file does not exist: {self.hi_c_2}")
        
        if self.reference_genome and not os.path.exists(self.reference_genome):
            errors.append(f"参考基因组文件不存在 | Reference genome file does not exist: {self.reference_genome}")
        
        # ===== 验证Hi-C数据配对 | Validate Hi-C data pairing =====
        if bool(self.hi_c_1) != bool(self.hi_c_2):
            errors.append("Hi-C数据需要同时提供两端文件 | Hi-C data requires both end files")
        
        # ===== 验证数值参数 | Validate numeric parameters =====
        if self.threads <= 0:
            errors.append("线程数必须大于0 | Number of threads must be greater than 0")
        
        if self.memory <= 0:
            errors.append("内存大小必须大于0 | Memory size must be greater than 0")
        
        if self.purge_level not in range(0, 4):
            errors.append("purge级别必须在0-3之间 | Purge level must be between 0-3")
        
        if self.purge_max <= 0:
            errors.append("最大purge覆盖度必须大于0 | Maximum purge coverage must be greater than 0")
        
        if self.similarity_threshold <= 0 or self.similarity_threshold > 1:
            errors.append("相似性阈值必须在0-1之间 | Similarity threshold must be between 0-1")
        
        if self.min_contig_length <= 0:
            errors.append("最小contig长度必须大于0 | Minimum contig length must be greater than 0")
        
        if self.max_runtime <= 0:
            errors.append("最大运行时间必须大于0 | Maximum runtime must be greater than 0")
        
        # ===== 验证字符串参数 | Validate string parameters =====
        valid_busco_modes = ['genome', 'proteins', 'transcriptome']
        if self.busco_mode not in valid_busco_modes:
            errors.append(f"BUSCO模式必须是以下之一 | BUSCO mode must be one of: {valid_busco_modes}")
        
        valid_assembly_types = ['auto', 'diploid', 'triploid', 'polyploid']
        if self.assembly_type not in valid_assembly_types:
            errors.append(f"组装类型必须是以下之一 | Assembly type must be one of: {valid_assembly_types}")
        
        valid_output_formats = ['fasta', 'gfa', 'both']
        for fmt in self.output_formats:
            if fmt not in valid_output_formats:
                errors.append(f"输出格式必须是以下之一 | Output format must be one of: {valid_output_formats}")
        
        valid_log_levels = ['DEBUG', 'INFO', 'WARNING', 'ERROR']
        if self.log_level not in valid_log_levels:
            errors.append(f"日志级别必须是以下之一 | Log level must be one of: {valid_log_levels}")
        
        # ===== 验证基因组大小格式 | Validate genome size format =====
        if self.hg_size != 'auto':
            if not self._validate_genome_size_format(self.hg_size):
                errors.append("基因组大小格式错误，应为如 '1.4g', '2100m', '2.1G' | Invalid genome size format, should be like '1.4g', '2100m', '2.1G'")
        
        # ===== 验证工具可用性 | Validate tool availability =====
        required_tools = [
            (self.hifiasm_path, "HiFiasm"),
            (self.python_path, "Python3")
        ]
        
        if not self.skip_busco:
            required_tools.append((self.busco_path, "BUSCO"))
        
        if not self.skip_quast:
            required_tools.append((self.quast_path, "QUAST"))
        
        for tool_path, tool_name in required_tools:
            if not shutil.which(tool_path):
                errors.append(f"{tool_name}工具未找到或不可执行 | {tool_name} tool not found or not executable: {tool_path}")
        
        # ===== 验证目录权限 | Validate directory permissions =====
        try:
            test_file = os.path.join(self.output_dir, '.write_test')
            with open(test_file, 'w') as f:
                f.write('test')
            os.remove(test_file)
        except (OSError, IOError):
            errors.append(f"输出目录不可写 | Output directory is not writable: {self.output_dir}")
        
        # ===== 抛出错误 | Raise errors =====
        if errors:
            error_message = "\n配置验证失败 | Configuration validation failed:\n" + "\n".join(f"  - {error}" for error in errors)
            raise ConfigurationError(error_message)
    
    def _validate_genome_size_format(self, size_str: str) -> bool:
        """验证基因组大小格式 | Validate genome size format"""
        import re
        pattern = r'^\d+(\.\d+)?[kmgtKMGT]?b?$'
        return bool(re.match(pattern, size_str))
    
    def estimate_genome_size(self) -> str:
        """自动估计基因组大小 | Automatically estimate genome size"""
        if self.hg_size != 'auto':
            return self.hg_size
        
        # 基于组装类型的默认估计
        default_sizes = {
            'diploid': '1.4g',
            'triploid': '2.1g', 
            'polyploid': '2.8g'
        }
        
        if self.assembly_type in default_sizes:
            estimated_size = default_sizes[self.assembly_type]
        else:
            # 尝试从输入文件大小估计
            try:
                file_size = os.path.getsize(self.input_reads)
                # 假设覆盖度约为50x，转换为基因组大小
                estimated_bp = file_size // 50
                if estimated_bp > 2e9:
                    estimated_size = f"{estimated_bp/1e9:.1f}g"
                elif estimated_bp > 2e6:
                    estimated_size = f"{estimated_bp/1e6:.0f}m"
                else:
                    estimated_size = "1.4g"  # 默认值
            except:
                estimated_size = "1.4g"  # 默认值
        
        return estimated_size
    
    def get_busco_lineage(self) -> str:
        """获取BUSCO谱系 | Get BUSCO lineage"""
        if self.busco_lineage != 'auto':
            return self.busco_lineage
        
        # 自动选择合适的BUSCO谱系
        recommended_lineages = [
            'brassicales_odb10',  # 十字花科
            'eudicots_odb10',     # 真双子叶植物
            'embryophyta_odb10',  # 陆地植物
            'eukaryota_odb10'     # 真核生物
        ]
        
        return recommended_lineages[0]  # 默认返回第一个
    
    def get_hifiasm_command_args(self) -> List[str]:
        """获取HiFiasm命令参数 | Get HiFiasm command arguments"""
        args = [
            self.hifiasm_path,
            '-o', self.prefix,
            '-t', str(self.threads),
            '--hg-size', self.estimate_genome_size(),
            '-l', str(self.purge_level),
            '--purge-max', str(self.purge_max),
            '-s', str(self.similarity_threshold)
        ]
        
        # 添加ONT数据
        if self.ont_reads:
            args.extend(['--ul', self.ont_reads])
        
        # 添加Hi-C数据
        if self.hi_c_1 and self.hi_c_2:
            args.extend(['--h1', self.hi_c_1, '--h2', self.hi_c_2])
        
        # 添加额外参数
        if self.extra_hifiasm_args:
            args.extend(self.extra_hifiasm_args.split())
        
        # 添加输入文件
        args.append(self.input_reads)
        
        return args
    
    def to_dict(self) -> dict:
        """转换为字典格式 | Convert to dictionary format"""
        return {
            field.name: getattr(self, field.name)
            for field in self.__dataclass_fields__.values()
        }
    
    def save_config(self, config_file: Optional[str] = None) -> str:
        """保存配置到文件 | Save configuration to file"""
        import json
        
        if config_file is None:
            config_file = os.path.join(self.output_dir, 'hifiasm_config.json')
        
        config_data = self.to_dict()
        
        with open(config_file, 'w', encoding='utf-8') as f:
            json.dump(config_data, f, indent=2, ensure_ascii=False)
        
        return config_file
    
    @classmethod
    def from_config_file(cls, config_file: str) -> 'HifiasmConfig':
        """从配置文件加载 | Load from configuration file"""
        import json
        
        with open(config_file, 'r', encoding='utf-8') as f:
            config_data = json.load(f)
        
        return cls(**config_data)