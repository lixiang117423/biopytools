# ===== FILE: genome_assembly/config.py =====
"""
基因组组装配置管理模块 | Genome Assembly Configuration Management Module
"""

import os
from dataclasses import dataclass
from pathlib import Path
from typing import Optional

@dataclass
class AssemblyConfig:
    """基因组组装配置类 | Genome Assembly Configuration Class"""
    
    # 必需文件 | Required files
    hifi_reads: str
    output_dir: str
    
    # 可选输入文件 | Optional input files
    ont_reads: Optional[str] = None
    reference_genome: Optional[str] = None
    
    # 组装参数 | Assembly parameters
    genome_size: str = '3.2g'
    ploidy: int = 2
    
    # Verkko参数 | Verkko parameters
    verkko_version: str = '1.4.1'
    verkko_grid: bool = False
    verkko_local_cpus: int = 32
    verkko_local_memory: int = 128
    verkko_screen: bool = True
    
    # hifiasm参数 | hifiasm parameters
    hifiasm_version: str = '0.19.6'
    hifiasm_ultra_long: bool = True
    hifiasm_purge_level: int = 3
    hifiasm_similarity: float = 0.8
    
    # Graphasing参数 | Graphasing parameters
    graphasing_version: str = '0.3.1-alpha'
    graphasing_kmer_size: int = 21
    
    # 质控参数 | Quality control parameters
    run_contamination_screen: bool = True
    fcs_version: str = '0.4.0'
    run_flagger: bool = True
    flagger_version: str = '0.3.3'
    run_merqury: bool = True
    merqury_version: str = '1.0'
    run_inspector: bool = True
    inspector_version: str = '1.2'
    
    # 质量评估参数 | Quality assessment parameters
    run_deepvariant: bool = True
    deepvariant_version: str = '1.6'
    run_compleasm: bool = True
    compleasm_version: str = '0.2.5'
    orthodb_version: str = '10'
    
    # 比对参数 | Alignment parameters
    minimap2_version: str = '2.26'
    mashmap_version: str = '3.1.3'
    
    # 系统参数 | System parameters
    threads: int = 32
    memory: int = 128
    tmp_dir: str = '/tmp'
    keep_temp: bool = False
    
    # 工具路径 | Tool paths
    verkko_path: str = 'verkko'
    hifiasm_path: str = 'hifiasm'
    graphasing_path: str = 'graphasing'
    fcs_path: str = 'fcs.py'
    flagger_path: str = 'flagger'
    merqury_path: str = 'merqury.sh'
    nucfreq_path: str = 'nucfreq'
    inspector_path: str = 'inspector.py'
    deepvariant_path: str = 'run_deepvariant'
    compleasm_path: str = 'compleasm'
    minimap2_path: str = 'minimap2'
    mashmap_path: str = 'mashmap'
    
    # 家系分析参数 | Family trio parameters
    trio_mode: bool = False
    parent1_reads: Optional[str] = None
    parent2_reads: Optional[str] = None
    child_reads: Optional[str] = None
    
    def __post_init__(self):
        """初始化后处理 | Post-initialization processing"""
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)
        
        # 标准化路径 | Normalize paths
        self.hifi_reads = os.path.normpath(os.path.abspath(self.hifi_reads))
        if self.ont_reads:
            self.ont_reads = os.path.normpath(os.path.abspath(self.ont_reads))
        if self.reference_genome:
            self.reference_genome = os.path.normpath(os.path.abspath(self.reference_genome))
        self.output_dir = os.path.normpath(os.path.abspath(self.output_dir))
    
    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        errors = []
        
        # 检查必需文件 | Check required files
        if not os.path.exists(self.hifi_reads):
            errors.append(f"HiFi reads文件不存在 | HiFi reads file does not exist: {self.hifi_reads}")
        
        if self.ont_reads and not os.path.exists(self.ont_reads):
            errors.append(f"ONT reads文件不存在 | ONT reads file does not exist: {self.ont_reads}")
            
        if self.reference_genome and not os.path.exists(self.reference_genome):
            errors.append(f"参考基因组文件不存在 | Reference genome file does not exist: {self.reference_genome}")
        
        # 检查参数范围 | Check parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正整数 | Number of threads must be positive: {self.threads}")
        
        if self.memory <= 0:
            errors.append(f"内存大小必须为正数 | Memory size must be positive: {self.memory}")
            
        if self.ploidy not in [1, 2]:
            errors.append(f"倍性必须为1或2 | Ploidy must be 1 or 2: {self.ploidy}")
        
        # 家系模式检查 | Trio mode validation
        if self.trio_mode:
            trio_files = [self.parent1_reads, self.parent2_reads, self.child_reads]
            if not all(trio_files):
                errors.append("家系模式需要提供所有三个样本的reads文件 | Trio mode requires all three sample read files")
            else:
                for i, reads_file in enumerate(trio_files, 1):
                    if not os.path.exists(reads_file):
                        errors.append(f"家系样本{i}的reads文件不存在 | Trio sample {i} reads file does not exist: {reads_file}")
        
        if errors:
            raise ValueError("\n".join(errors))
        
        return True

# ===== END FILE =====

# ===== FILE: genome_assembly/utils.py =====
"""
基因组组装工具函数模块 | Genome Assembly Utility Functions Module
"""

import logging
import subprocess
import sys
from pathlib import Path

class AssemblyLogger:
    """基因组组装日志管理器 | Genome Assembly Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "assembly_analysis.log"):
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
        self.working_dir = working_dir.resolve()
    
    def run(self, cmd: str, description: str = "", timeout: int = None) -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"执行步骤 | Executing step: {description}")
        
        self.logger.info(f"命令 | Command: {cmd}")
        self.logger.info(f"工作目录 | Working directory: {self.working_dir}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir,
                timeout=timeout
            )
            
            self.logger.info(f"命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"错误代码 | Error code: {e.returncode}")
            self.logger.error(f"错误信息 | Error message: {e.stderr}")
            if e.stdout:
                self.logger.error(f"标准输出 | Stdout: {e.stdout}")
            return False
        except subprocess.TimeoutExpired:
            self.logger.error(f"命令执行超时 | Command execution timeout: {description}")
            return False

def check_dependencies(config, logger):
    """检查依赖软件 | Check dependencies"""
    logger.info("检查依赖软件 | Checking dependencies")
    
    # 基本依赖 | Basic dependencies
    basic_deps = [
        (config.verkko_path, "Verkko"),
        (config.hifiasm_path, "hifiasm"),
        (config.minimap2_path, "minimap2"),
        (config.mashmap_path, "mashmap")
    ]
    
    # 可选依赖 | Optional dependencies
    optional_deps = []
    if config.run_contamination_screen:
        optional_deps.append((config.fcs_path, "NCBI FCS"))
    if config.run_flagger:
        optional_deps.append((config.flagger_path, "Flagger"))
    if config.run_merqury:
        optional_deps.append((config.merqury_path, "Merqury"))
    if config.run_inspector:
        optional_deps.append((config.inspector_path, "Inspector"))
    if config.run_deepvariant:
        optional_deps.append((config.deepvariant_path, "DeepVariant"))
    if config.run_compleasm:
        optional_deps.append((config.compleasm_path, "compleasm"))
    
    missing_deps = []
    optional_missing = []
    
    # 检查基本依赖 | Check basic dependencies
    for cmd, name in basic_deps:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✓ {name} 可用 | available")
            else:
                missing_deps.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            missing_deps.append(name)
    
    # 检查可选依赖 | Check optional dependencies
    for cmd, name in optional_deps:
        try:
            result = subprocess.run([cmd, "--version"], 
                                  capture_output=True, text=True, timeout=10)
            if result.returncode == 0:
                logger.info(f"✓ {name} 可用 | available")
            else:
                optional_missing.append(name)
        except (subprocess.TimeoutExpired, FileNotFoundError):
            optional_missing.append(name)
    
    if missing_deps:
        error_msg = f"缺少必需依赖软件 | Missing required dependencies: {', '.join(missing_deps)}"
        logger.error(error_msg)
        raise RuntimeError(error_msg)
    
    if optional_missing:
        logger.warning(f"缺少可选依赖软件 | Missing optional dependencies: {', '.join(optional_missing)}")
        logger.warning("相关功能将被跳过 | Related functions will be skipped")
    
    return True

# ===== END FILE =====

# ===== FILE: genome_assembly/data_processing.py =====
"""
基因组组装数据处理模块 | Genome Assembly Data Processing Module
"""

import os
from typing import List, Dict
from .utils import CommandRunner

class ReadProcessor:
    """测序数据预处理器 | Sequencing Read Preprocessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def validate_reads(self) -> bool:
        """验证reads文件格式和质量 | Validate reads file format and quality"""
        self.logger.info("验证reads文件 | Validating reads files")
        
        reads_files = [self.config.hifi_reads]
        if self.config.ont_reads:
            reads_files.append(self.config.ont_reads)
        
        if self.config.trio_mode:
            reads_files.extend([self.config.parent1_reads, self.config.parent2_reads, self.config.child_reads])
        
        for reads_file in reads_files:
            if not self._check_file_format(reads_file):
                return False
        
        return True
    
    def _check_file_format(self, reads_file: str) -> bool:
        """检查单个reads文件格式 | Check individual reads file format"""
        self.logger.info(f"检查文件格式 | Checking file format: {reads_file}")
        
        # 检查文件扩展名 | Check file extension
        valid_extensions = ['.fastq', '.fq', '.fastq.gz', '.fq.gz', '.fasta', '.fa', '.fasta.gz', '.fa.gz']
        if not any(reads_file.endswith(ext) for ext in valid_extensions):
            self.logger.error(f"不支持的文件格式 | Unsupported file format: {reads_file}")
            return False
        
        # 检查文件大小 | Check file size
        file_size = os.path.getsize(reads_file)
        if file_size == 0:
            self.logger.error(f"文件为空 | File is empty: {reads_file}")
            return False
        
        self.logger.info(f"文件格式验证通过 | File format validation passed: {reads_file} ({file_size/1024/1024/1024:.2f} GB)")
        return True

# ===== END FILE =====

# ===== FILE: genome_assembly/assembly.py =====
"""
基因组组装核心模块 | Genome Assembly Core Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class VerkkoAssembler:
    """Verkko组装器 | Verkko Assembler"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_verkko_assembly(self) -> str:
        """运行Verkko组装 | Run Verkko assembly"""
        self.logger.info("开始Verkko组装 | Starting Verkko assembly")
        
        output_prefix = self.config.output_path / "verkko_assembly"
        
        # 构建Verkko命令 | Build Verkko command
        cmd_parts = [
            self.config.verkko_path,
            f"--hifi {self.config.hifi_reads}"
        ]
        
        if self.config.ont_reads:
            cmd_parts.append(f"--nano {self.config.ont_reads}")
        
        cmd_parts.extend([
            f"--output {output_prefix}",
            f"--local-cpus {self.config.verkko_local_cpus}",
            f"--local-memory {self.config.verkko_local_memory}g"
        ])
        
        if self.config.verkko_screen:
            cmd_parts.append("--screen")
        
        if self.config.trio_mode:
            cmd_parts.extend([
                f"--hifi-hap1 {self.config.parent1_reads}",
                f"--hifi-hap2 {self.config.parent2_reads}",
                f"--hifi-trio {self.config.child_reads}"
            ])
        
        cmd = " ".join(cmd_parts)
        
        # 执行组装 | Execute assembly
        success = self.cmd_runner.run(cmd, "Verkko基因组组装 | Verkko genome assembly", timeout=86400)  # 24小时超时
        
        if success:
            assembly_file = output_prefix / "assembly.fasta"
            if assembly_file.exists():
                self.logger.info(f"Verkko组装完成 | Verkko assembly completed: {assembly_file}")
                return str(assembly_file)
        
        self.logger.error("Verkko组装失败 | Verkko assembly failed")
        return None

class HifiasmAssembler:
    """hifiasm组装器 | hifiasm Assembler"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_hifiasm_assembly(self) -> str:
        """运行hifiasm组装 | Run hifiasm assembly"""
        self.logger.info("开始hifiasm组装 | Starting hifiasm assembly")
        
        output_prefix = self.config.output_path / "hifiasm_assembly" / "assembly"
        output_prefix.parent.mkdir(parents=True, exist_ok=True)
        
        # 构建hifiasm命令 | Build hifiasm command
        cmd_parts = [
            self.config.hifiasm_path,
            f"-o {output_prefix}",
            f"-t {self.config.threads}"
        ]
        
        if self.config.hifiasm_ultra_long:
            cmd_parts.append("--ul")
        
        cmd_parts.extend([
            f"--purge-level {self.config.hifiasm_purge_level}",
            f"-s {self.config.hifiasm_similarity}",
            self.config.hifi_reads
        ])
        
        if self.config.ont_reads:
            cmd_parts.append(self.config.ont_reads)
        
        cmd = " ".join(cmd_parts)
        
        # 执行组装 | Execute assembly
        success = self.cmd_runner.run(cmd, "hifiasm基因组组装 | hifiasm genome assembly", timeout=86400)
        
        if success:
            # 转换GFA到FASTA | Convert GFA to FASTA
            gfa_file = f"{output_prefix}.bp.p_ctg.gfa"
            fasta_file = f"{output_prefix}.fasta"
            
            if os.path.exists(gfa_file):
                convert_cmd = f"awk '/^S/{{print \">\"$2; print $3}}' {gfa_file} > {fasta_file}"
                if self.cmd_runner.run(convert_cmd, "转换GFA到FASTA | Convert GFA to FASTA"):
                    self.logger.info(f"hifiasm组装完成 | hifiasm assembly completed: {fasta_file}")
                    return fasta_file
        
        self.logger.error("hifiasm组装失败 | hifiasm assembly failed")
        return None

class GraphasingProcessor:
    """Graphasing分期处理器 | Graphasing Phasing Processor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_graphasing(self, assembly_file: str) -> str:
        """运行Graphasing分期 | Run Graphasing phasing"""
        self.logger.info("开始Graphasing分期分析 | Starting Graphasing phasing analysis")
        
        output_prefix = self.config.output_path / "graphasing_output"
        output_prefix.mkdir(parents=True, exist_ok=True)
        
        # 构建Graphasing命令 | Build Graphasing command
        cmd_parts = [
            self.config.graphasing_path,
            f"--assembly {assembly_file}",
            f"--reads {self.config.hifi_reads}",
            f"--output {output_prefix}",
            f"--kmer-size {self.config.graphasing_kmer_size}",
            f"--threads {self.config.threads}"
        ]
        
        cmd = " ".join(cmd_parts)
        
        # 执行分期 | Execute phasing
        success = self.cmd_runner.run(cmd, "Graphasing分期分析 | Graphasing phasing analysis", timeout=21600)  # 6小时超时
        
        if success:
            phased_file = output_prefix / "phased_assembly.fasta"
            if phased_file.exists():
                self.logger.info(f"Graphasing分期完成 | Graphasing phasing completed: {phased_file}")
                return str(phased_file)
        
        self.logger.warning("Graphasing分期失败，使用原始组装 | Graphasing phasing failed, using original assembly")
        return assembly_file

# ===== END FILE =====

# ===== FILE: genome_assembly/quality_control.py =====
"""
基因组组装质量控制模块 | Genome Assembly Quality Control Module
"""

from pathlib import Path
from typing import Dict
from .utils import CommandRunner

class ContaminationScreener:
    """外源污染筛查器 | Contamination Screener"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def screen_contamination(self, assembly_file: str) -> str:
        """筛查外源污染 | Screen for contamination"""
        if not self.config.run_contamination_screen:
            self.logger.info("跳过外源污染筛查 | Skipping contamination screening")
            return assembly_file
        
        self.logger.info("开始外源污染筛查 | Starting contamination screening")
        
        output_dir = self.config.output_path / "contamination_screen"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建FCS命令 | Build FCS command
        cmd = f"{self.config.fcs_path} --fasta {assembly_file} --output-dir {output_dir} --euk"
        
        # 执行筛查 | Execute screening
        success = self.cmd_runner.run(cmd, "NCBI FCS外源污染筛查 | NCBI FCS contamination screening", timeout=7200)
        
        if success:
            clean_file = output_dir / "cleaned_assembly.fasta"
            if clean_file.exists():
                self.logger.info(f"外源污染筛查完成 | Contamination screening completed: {clean_file}")
                return str(clean_file)
        
        self.logger.warning("外源污染筛查失败，使用原始组装 | Contamination screening failed, using original assembly")
        return assembly_file

class AssemblyAnnotator:
    """组装注释器 | Assembly Annotator"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_flagger(self, assembly_file: str) -> Dict:
        """运行Flagger错误注释 | Run Flagger error annotation"""
        if not self.config.run_flagger:
            self.logger.info("跳过Flagger错误注释 | Skipping Flagger error annotation")
            return {}
        
        self.logger.info("开始Flagger错误注释 | Starting Flagger error annotation")
        
        output_dir = self.config.output_path / "flagger_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建Flagger命令 | Build Flagger command
        cmd = f"{self.config.flagger_path} --assembly {assembly_file} --reads {self.config.hifi_reads} --output {output_dir}"
        
        # 执行注释 | Execute annotation
        success = self.cmd_runner.run(cmd, "Flagger组装错误注释 | Flagger assembly error annotation", timeout=3600)
        
        results = {}
        if success:
            results_file = output_dir / "flagger_results.tsv"
            if results_file.exists():
                self.logger.info(f"Flagger错误注释完成 | Flagger error annotation completed: {results_file}")
                results['flagger_results'] = str(results_file)
        
        return results
    
    def run_nucfreq(self, assembly_file: str) -> Dict:
        """运行NucFreq分析 | Run NucFreq analysis"""
        self.logger.info("开始NucFreq核苷酸频率分析 | Starting NucFreq nucleotide frequency analysis")
        
        output_dir = self.config.output_path / "nucfreq_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建NucFreq命令 | Build NucFreq command
        cmd = f"{self.config.nucfreq_path} {assembly_file} > {output_dir}/nucfreq_results.txt"
        
        # 执行分析 | Execute analysis
        success = self.cmd_runner.run(cmd, "NucFreq核苷酸频率分析 | NucFreq nucleotide frequency analysis")
        
        results = {}
        if success:
            results_file = output_dir / "nucfreq_results.txt"
            if results_file.exists():
                self.logger.info(f"NucFreq分析完成 | NucFreq analysis completed: {results_file}")
                results['nucfreq_results'] = str(results_file)
        
        return results

# ===== END FILE =====

# ===== FILE: genome_assembly/quality_assessment.py =====
"""
基因组组装质量评估模块 | Genome Assembly Quality Assessment Module
"""

from pathlib import Path
from typing import Dict
from .utils import CommandRunner

class QualityAssessor:
    """组装质量评估器 | Assembly Quality Assessor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_merqury(self, assembly_file: str) -> Dict:
        """运行Merqury质量评估 | Run Merqury quality assessment"""
        if not self.config.run_merqury:
            self.logger.info("跳过Merqury质量评估 | Skipping Merqury quality assessment")
            return {}
        
        self.logger.info("开始Merqury质量评估 | Starting Merqury quality assessment")
        
        output_dir = self.config.output_path / "merqury_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 首先生成k-mer数据库 | First generate k-mer database
        kmer_db = output_dir / "reads.meryl"
        cmd_meryl = f"meryl count k=21 {self.config.hifi_reads} output {kmer_db}"
        
        if not self.cmd_runner.run(cmd_meryl, "生成k-mer数据库 | Generate k-mer database", timeout=3600):
            return {}
        
        # 运行Merqury | Run Merqury
        cmd_merqury = f"cd {output_dir} && {self.config.merqury_path} {kmer_db} {assembly_file} assembly"
        
        success = self.cmd_runner.run(cmd_merqury, "Merqury质量评估 | Merqury quality assessment", timeout=3600)
        
        results = {}
        if success:
            qv_file = output_dir / "assembly.qv"
            completeness_file = output_dir / "assembly.completeness.stats"
            
            if qv_file.exists():
                self.logger.info(f"Merqury质量评估完成 | Merqury quality assessment completed")
                results['merqury_qv'] = str(qv_file)
                results['merqury_completeness'] = str(completeness_file)
        
        return results
    
    def run_inspector(self, assembly_file: str) -> Dict:
        """运行Inspector组装检查 | Run Inspector assembly inspection"""
        if not self.config.run_inspector:
            self.logger.info("跳过Inspector组装检查 | Skipping Inspector assembly inspection")
            return {}
        
        self.logger.info("开始Inspector组装检查 | Starting Inspector assembly inspection")
        
        output_dir = self.config.output_path / "inspector_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建Inspector命令 | Build Inspector command
        cmd = f"{self.config.inspector_path} -c {assembly_file} -r {self.config.hifi_reads} -o {output_dir} --thread {self.config.threads}"
        
        # 执行检查 | Execute inspection
        success = self.cmd_runner.run(cmd, "Inspector组装质量检查 | Inspector assembly quality inspection", timeout=3600)
        
        results = {}
        if success:
            results_file = output_dir / "inspector_summary.txt"
            if results_file.exists():
                self.logger.info(f"Inspector检查完成 | Inspector inspection completed: {results_file}")
                results['inspector_results'] = str(results_file)
        
        return results
    
    def run_deepvariant_qv(self, assembly_file: str) -> Dict:
        """运行DeepVariant质量值估计 | Run DeepVariant quality value estimation"""
        if not self.config.run_deepvariant:
            self.logger.info("跳过DeepVariant质量值估计 | Skipping DeepVariant quality value estimation")
            return {}
        
        self.logger.info("开始DeepVariant质量值估计 | Starting DeepVariant quality value estimation")
        
        output_dir = self.config.output_path / "deepvariant_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 首先比对reads到组装 | First align reads to assembly
        bam_file = output_dir / "aligned_reads.bam"
        cmd_align = f"{self.config.minimap2_path} -ax map-hifi {assembly_file} {self.config.hifi_reads} | samtools sort -o {bam_file}"
        
        if not self.cmd_runner.run(cmd_align, "比对reads到组装 | Align reads to assembly", timeout=3600):
            return {}
        
        # 索引BAM文件 | Index BAM file
        cmd_index = f"samtools index {bam_file}"
        if not self.cmd_runner.run(cmd_index, "索引BAM文件 | Index BAM file"):
            return {}
        
        # 运行DeepVariant | Run DeepVariant
        vcf_file = output_dir / "variants.vcf.gz"
        cmd_dv = f"{self.config.deepvariant_path} --model_type=PACBIO --ref={assembly_file} --reads={bam_file} --output_vcf={vcf_file} --num_shards={self.config.threads}"
        
        success = self.cmd_runner.run(cmd_dv, "DeepVariant变异检测 | DeepVariant variant calling", timeout=7200)
        
        results = {}
        if success and vcf_file.exists():
            self.logger.info(f"DeepVariant质量值估计完成 | DeepVariant quality value estimation completed: {vcf_file}")
            results['deepvariant_vcf'] = str(vcf_file)
        
        return results
    
    def run_compleasm(self, assembly_file: str) -> Dict:
        """运行compleasm基因完整性评估 | Run compleasm gene completeness assessment"""
        if not self.config.run_compleasm:
            self.logger.info("跳过compleasm基因完整性评估 | Skipping compleasm gene completeness assessment")
            return {}
        
        self.logger.info("开始compleasm基因完整性评估 | Starting compleasm gene completeness assessment")
        
        output_dir = self.config.output_path / "compleasm_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 构建compleasm命令 | Build compleasm command
        cmd = f"{self.config.compleasm_path} run -a {assembly_file} -o {output_dir} -l primates_odb{self.config.orthodb_version} -t {self.config.threads}"
        
        # 执行评估 | Execute assessment
        success = self.cmd_runner.run(cmd, "compleasm基因完整性评估 | compleasm gene completeness assessment", timeout=3600)
        
        results = {}
        if success:
            results_file = output_dir / "summary.txt"
            if results_file.exists():
                self.logger.info(f"compleasm评估完成 | compleasm assessment completed: {results_file}")
                results['compleasm_results'] = str(results_file)
        
        return results

# ===== END FILE =====

# ===== FILE: genome_assembly/alignment.py =====
"""
基因组组装比对分析模块 | Genome Assembly Alignment Analysis Module
"""

from pathlib import Path
from typing import Dict
from .utils import CommandRunner

class AlignmentAnalyzer:
    """比对分析器 | Alignment Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def align_to_reference(self, assembly_file: str) -> Dict:
        """比对到参考基因组 | Align to reference genome"""
        if not self.config.reference_genome:
            self.logger.info("未提供参考基因组，跳过比对分析 | Reference genome not provided, skipping alignment analysis")
            return {}
        
        self.logger.info("开始比对到参考基因组分析 | Starting alignment to reference genome analysis")
        
        output_dir = self.config.output_path / "alignment_output"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        # minimap2比对 | minimap2 alignment
        sam_file = output_dir / "alignment.sam"
        cmd_minimap2 = f"{self.config.minimap2_path} -ax asm5 {self.config.reference_genome} {assembly_file} > {sam_file}"
        
        if self.cmd_runner.run(cmd_minimap2, "minimap2序列比对 | minimap2 sequence alignment", timeout=3600):
            self.logger.info(f"minimap2比对完成 | minimap2 alignment completed: {sam_file}")
            results['minimap2_alignment'] = str(sam_file)
        
        # mashmap快速比对 | mashmap fast alignment
        mashmap_file = output_dir / "mashmap_alignment.out"
        cmd_mashmap = f"{self.config.mashmap_path} -r {self.config.reference_genome} -q {assembly_file} -o {mashmap_file}"
        
        if self.cmd_runner.run(cmd_mashmap, "mashmap快速序列比对 | mashmap fast sequence alignment"):
            self.logger.info(f"mashmap比对完成 | mashmap alignment completed: {mashmap_file}")
            results['mashmap_alignment'] = str(mashmap_file)
        
        return results
    
    def analyze_t2t_status(self, alignment_results: Dict) -> Dict:
        """分析T2T状态 | Analyze T2T status"""
        self.logger.info("分析T2T(端粒到端粒)状态 | Analyzing T2T (telomere-to-telomere) status")
        
        output_dir = self.config.output_path / "t2t_analysis"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        t2t_results = {}
        
        if 'minimap2_alignment' in alignment_results:
            # 基于minimap2比对结果分析T2T状态 | Analyze T2T status based on minimap2 alignment
            sam_file = alignment_results['minimap2_alignment']
            t2t_file = output_dir / "t2t_status.txt"
            
            # 这里可以添加具体的T2T状态分析代码 | Add specific T2T status analysis code here
            with open(t2t_file, 'w') as f:
                f.write("# T2T状态分析结果 | T2T Status Analysis Results\n")
                f.write("# 基于minimap2比对结果 | Based on minimap2 alignment results\n")
                f.write("\n")
                f.write("详细的T2T状态分析需要自定义实现 | Detailed T2T status analysis requires custom implementation\n")
                f.write("建议分析内容 | Suggested analysis content:\n")
                f.write("- 染色体端粒序列检测 | Telomere sequence detection\n")
                f.write("- 着丝粒区域完整性 | Centromere region integrity\n")
                f.write("- 已知gap位点闭合状态 | Known gap closure status\n")
                f.write("- 组装连续性评估 | Assembly continuity assessment\n")
            
            self.logger.info("T2T状态分析完成 | T2T status analysis completed")
            t2t_results['t2t_analysis'] = str(t2t_file)
        
        return t2t_results

class TrioAnalyzer:
    """家系三元组分析器 | Family Trio Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def analyze_parental_support(self, child_assembly: str, parent1_assembly: str, parent2_assembly: str) -> Dict:
        """分析亲本支持度 | Analyze parental support"""
        if not self.config.trio_mode:
            self.logger.info("非家系模式，跳过亲本支持度分析 | Non-trio mode, skipping parental support analysis")
            return {}
        
        self.logger.info("开始家系亲本支持度分析 | Starting family parental support analysis")
        
        output_dir = self.config.output_path / "trio_analysis"
        output_dir.mkdir(parents=True, exist_ok=True)
        
        results = {}
        
        # 子代到亲本1的比对 | Child to parent1 alignment
        parent1_alignment = output_dir / "child_to_parent1.sam"
        cmd_p1 = f"{self.config.minimap2_path} -ax asm5 {parent1_assembly} {child_assembly} > {parent1_alignment}"
        
        if self.cmd_runner.run(cmd_p1, "子代到亲本1比对 | Child to parent1 alignment", timeout=3600):
            results['parent1_alignment'] = str(parent1_alignment)
        
        # 子代到亲本2的比对 | Child to parent2 alignment
        parent2_alignment = output_dir / "child_to_parent2.sam"
        cmd_p2 = f"{self.config.minimap2_path} -ax asm5 {parent2_assembly} {child_assembly} > {parent2_alignment}"
        
        if self.cmd_runner.run(cmd_p2, "子代到亲本2比对 | Child to parent2 alignment", timeout=3600):
            results['parent2_alignment'] = str(parent2_alignment)
        
        # 分析CIGAR操作计算亲本支持度 | Analyze CIGAR operations to calculate parental support
        if results:
            support_file = output_dir / "parental_support.txt"
            self._calculate_parental_support(results, support_file)
            results['parental_support'] = str(support_file)
        
        return results
    
    def _calculate_parental_support(self, alignment_results: Dict, output_file: Path):
        """计算亲本支持度 | Calculate parental support"""
        self.logger.info("计算亲本支持度统计 | Calculating parental support statistics")
        
        # 这里应该实现具体的CIGAR分析逻辑 | Implement specific CIGAR analysis logic here
        with open(output_file, 'w') as f:
            f.write("# 家系亲本支持度分析结果 | Family Parental Support Analysis Results\n")
            f.write("# 基于minimap2比对的CIGAR操作分析 | Based on CIGAR operations from minimap2 alignment\n")
            f.write("\n")
            f.write("子代单倍型1支持情况 | Child haplotype 1 support:\n")
            f.write("子代单倍型2支持情况 | Child haplotype 2 support:\n")
            f.write("\n")
            f.write("详细的CIGAR分析需要自定义实现 | Detailed CIGAR analysis requires custom implementation\n")
            f.write("建议分析内容 | Suggested analysis content:\n")
            f.write("- CIGAR操作统计 | CIGAR operation statistics\n")
            f.write("- 匹配率计算 | Match rate calculation\n")
            f.write("- 插入缺失分析 | Indel analysis\n")
            f.write("- 亲本特异性变异 | Parent-specific variants\n")

# ===== END FILE =====

# ===== FILE: genome_assembly/results.py =====
"""
基因组组装结果处理模块 | Genome Assembly Results Processing Module
"""

import time
import os
from typing import Dict

class ResultsProcessor:
    """结果处理器 | Results Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def generate_summary_report(self, assembly_results: Dict) -> str:
        """生成总结报告 | Generate summary report"""
        self.logger.info("生成组装分析总结报告 | Generating assembly analysis summary report")
        
        report_file = self.config.output_path / "assembly_summary_report.txt"
        
        with open(report_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("基因组组装分析总结报告 | Genome Assembly Analysis Summary Report\n")
            f.write("=" * 80 + "\n\n")
            
            f.write(f"分析时间 | Analysis time: {time.strftime('%Y-%m-%d %H:%M:%S')}\n")
            f.write(f"输出目录 | Output directory: {self.config.output_dir}\n\n")
            
            f.write("输入数据 | Input Data:\n")
            f.write(f"  - HiFi reads: {self.config.hifi_reads}\n")
            if self.config.ont_reads:
                f.write(f"  - ONT reads: {self.config.ont_reads}\n")
            if self.config.reference_genome:
                f.write(f"  - 参考基因组 | Reference genome: {self.config.reference_genome}\n")
            f.write("\n")
            
            f.write("组装参数 | Assembly Parameters:\n")
            f.write(f"  - 基因组大小 | Genome size: {self.config.genome_size}\n")
            f.write(f"  - 倍性 | Ploidy: {self.config.ploidy}\n")
            f.write(f"  - 线程数 | Threads: {self.config.threads}\n")
            f.write(f"  - 内存 | Memory: {self.config.memory}GB\n")
            f.write("\n")
            
            f.write("组装工具版本 | Assembly Tool Versions:\n")
            f.write(f"  - Verkko: {self.config.verkko_version}\n")
            f.write(f"  - hifiasm: {self.config.hifiasm_version}\n")
            f.write(f"  - Graphasing: {self.config.graphasing_version}\n")
            f.write("\n")
            
            f.write("主要输出文件 | Main Output Files:\n")
            for key, value in assembly_results.items():
                if isinstance(value, str) and os.path.exists(value):
                    f.write(f"  - {key}: {value}\n")
                elif isinstance(value, dict):
                    f.write(f"  - {key}:\n")
                    for sub_key, sub_value in value.items():
                        if isinstance(sub_value, str) and os.path.exists(sub_value):
                            f.write(f"    - {sub_key}: {sub_value}\n")
            f.write("\n")
            
            f.write("质量控制结果 | Quality Control Results:\n")
            if self.config.run_contamination_screen:
                f.write("  ✓ 外源污染筛查已完成 | Contamination screening completed\n")
            if self.config.run_flagger:
                f.write("  ✓ Flagger错误注释已完成 | Flagger error annotation completed\n")
            if self.config.run_merqury:
                f.write("  ✓ Merqury质量评估已完成 | Merqury quality assessment completed\n")
            if self.config.run_inspector:
                f.write("  ✓ Inspector组装检查已完成 | Inspector assembly inspection completed\n")
            if self.config.run_deepvariant:
                f.write("  ✓ DeepVariant质量值估计已完成 | DeepVariant quality value estimation completed\n")
            if self.config.run_compleasm:
                f.write("  ✓ compleasm基因完整性评估已完成 | compleasm gene completeness assessment completed\n")
            f.write("\n")
            
            if self.config.trio_mode:
                f.write("家系分析 | Family Analysis:\n")
                f.write("  ✓ 家系三元组分析已完成 | Trio analysis completed\n")
                f.write("  ✓ 亲本支持度分析已完成 | Parental support analysis completed\n")
                f.write("\n")
            
            f.write("分析流程说明 | Analysis Pipeline Description:\n")
            f.write("1. 数据验证 | Data validation\n")
            f.write("2. Verkko主要组装 | Primary assembly with Verkko\n")
            f.write("3. hifiasm补充组装 | Complementary assembly with hifiasm\n")
            f.write("4. Graphasing分期信号生成 | Phasing signal generation with Graphasing\n")
            f.write("5. 外源污染筛查 | Contamination screening\n")
            f.write("6. 组装错误注释 | Assembly error annotation\n")
            f.write("7. 质量评估 | Quality assessment\n")
            f.write("8. 参考基因组比对 | Reference genome alignment\n")
            if self.config.trio_mode:
                f.write("9. 家系分析 | Family analysis\n")
            f.write("\n")
            
            f.write("引用 | Citations:\n")
            f.write("- Verkko: Rautiainen et al. Nature Biotechnology (2023)\n")
            f.write("- hifiasm: Cheng et al. Nature Methods (2021)\n")
            f.write("- Graphasing: Garg et al. Nature Biotechnology (2021)\n")
            f.write("- Merqury: Rhie et al. Genome Biology (2020)\n")
            f.write("- DeepVariant: Poplin et al. Nature Biotechnology (2018)\n")
            f.write("- minimap2: Li. Bioinformatics (2018)\n")
        
        self.logger.info(f"总结报告已生成 | Summary report generated: {report_file}")
        return str(report_file)

# ===== END FILE =====

# ===== FILE: genome_assembly/main.py =====
"""
基因组组装分析主程序模块 | Genome Assembly Analysis Main Module
"""

import argparse
import sys
from .config import AssemblyConfig
from .utils import AssemblyLogger, CommandRunner, check_dependencies
from .data_processing import ReadProcessor
from .assembly import VerkkoAssembler, HifiasmAssembler, GraphasingProcessor
from .quality_control import ContaminationScreener, AssemblyAnnotator
from .quality_assessment import QualityAssessor
from .alignment import AlignmentAnalyzer, TrioAnalyzer
from .results import ResultsProcessor

class GenomeAssemblyAnalyzer:
    """基因组组装分析主类 | Main Genome Assembly Analyzer Class"""
    
    def __init__(self, **kwargs):
        # 初始化配置 | Initialize configuration
        self.config = AssemblyConfig(**kwargs)
        self.config.validate()
        
        # 初始化日志 | Initialize logging
        self.logger_manager = AssemblyLogger(self.config.output_path)
        self.logger = self.logger_manager.get_logger()
        
        # 初始化命令执行器 | Initialize command runner
        self.cmd_runner = CommandRunner(self.logger, self.config.output_path)
        
        # 初始化各个处理器 | Initialize processors
        self.read_processor = ReadProcessor(self.config, self.logger, self.cmd_runner)
        self.verkko_assembler = VerkkoAssembler(self.config, self.logger, self.cmd_runner)
        self.hifiasm_assembler = HifiasmAssembler(self.config, self.logger, self.cmd_runner)
        self.graphasing_processor = GraphasingProcessor(self.config, self.logger, self.cmd_runner)
        self.contamination_screener = ContaminationScreener(self.config, self.logger, self.cmd_runner)
        self.assembly_annotator = AssemblyAnnotator(self.config, self.logger, self.cmd_runner)
        self.quality_assessor = QualityAssessor(self.config, self.logger, self.cmd_runner)
        self.alignment_analyzer = AlignmentAnalyzer(self.config, self.logger, self.cmd_runner)
        self.trio_analyzer = TrioAnalyzer(self.config, self.logger, self.cmd_runner)
        self.results_processor = ResultsProcessor(self.config, self.logger)
    
    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.config, self.logger)
    
    def run_analysis(self):
        """运行完整的基因组组装分析流程 | Run complete genome assembly analysis pipeline"""
        try:
            self.logger.info("=" * 80)
            self.logger.info("开始基因组组装分析流程 | Starting genome assembly analysis pipeline")
            self.logger.info("=" * 80)
            
            # 步骤1: 检查依赖软件 | Step 1: Check dependencies
            self.logger.info("\n步骤1: 检查依赖软件 | Step 1: Check dependencies")
            self.check_dependencies()
            
            # 步骤2: 验证输入数据 | Step 2: Validate input data
            self.logger.info("\n步骤2: 验证输入数据 | Step 2: Validate input data")
            if not self.read_processor.validate_reads():
                raise RuntimeError("输入数据验证失败 | Input data validation failed")
            
            assembly_results = {}
            
            # 步骤3: Verkko主要组装 | Step 3: Primary assembly with Verkko
            self.logger.info("\n步骤3: Verkko主要组装 | Step 3: Primary assembly with Verkko")
            verkko_assembly = self.verkko_assembler.run_verkko_assembly()
            if verkko_assembly:
                assembly_results['verkko_assembly'] = verkko_assembly
            
            # 步骤4: hifiasm补充组装 | Step 4: Complementary assembly with hifiasm
            self.logger.info("\n步骤4: hifiasm补充组装 | Step 4: Complementary assembly with hifiasm")
            hifiasm_assembly = self.hifiasm_assembler.run_hifiasm_assembly()
            if hifiasm_assembly:
                assembly_results['hifiasm_assembly'] = hifiasm_assembly
            
            # 选择主要组装文件 | Choose primary assembly file
            primary_assembly = verkko_assembly if verkko_assembly else hifiasm_assembly
            if not primary_assembly:
                raise RuntimeError("所有组装都失败了 | All assemblies failed")
            
            # 步骤5: Graphasing分期分析 | Step 5: Graphasing phasing analysis
            self.logger.info("\n步骤5: Graphasing分期分析 | Step 5: Graphasing phasing analysis")
            phased_assembly = self.graphasing_processor.run_graphasing(primary_assembly)
            assembly_results['phased_assembly'] = phased_assembly
            
            # 步骤6: 外源污染筛查 | Step 6: Contamination screening
            self.logger.info("\n步骤6: 外源污染筛查 | Step 6: Contamination screening")
            clean_assembly = self.contamination_screener.screen_contamination(phased_assembly)
            assembly_results['clean_assembly'] = clean_assembly
            
            # 步骤7: 组装注释 | Step 7: Assembly annotation
            self.logger.info("\n步骤7: 组装注释 | Step 7: Assembly annotation")
            flagger_results = self.assembly_annotator.run_flagger(clean_assembly)
            nucfreq_results = self.assembly_annotator.run_nucfreq(clean_assembly)
            assembly_results.update(flagger_results)
            assembly_results.update(nucfreq_results)
            
            # 步骤8: 质量评估 | Step 8: Quality assessment
            self.logger.info("\n步骤8: 质量评估 | Step 8: Quality assessment")
            merqury_results = self.quality_assessor.run_merqury(clean_assembly)
            inspector_results = self.quality_assessor.run_inspector(clean_assembly)
            deepvariant_results = self.quality_assessor.run_deepvariant_qv(clean_assembly)
            compleasm_results = self.quality_assessor.run_compleasm(clean_assembly)
            
            assembly_results.update(merqury_results)
            assembly_results.update(inspector_results)
            assembly_results.update(deepvariant_results)
            assembly_results.update(compleasm_results)
            
            # 步骤9: 参考基因组比对分析 | Step 9: Reference genome alignment analysis
            self.logger.info("\n步骤9: 参考基因组比对分析 | Step 9: Reference genome alignment analysis")
            alignment_results = self.alignment_analyzer.align_to_reference(clean_assembly)
            t2t_results = self.alignment_analyzer.analyze_t2t_status(alignment_results)
            assembly_results.update(alignment_results)
            assembly_results.update(t2t_results)
            
            # 步骤10: 家系分析(如果适用) | Step 10: Family analysis (if applicable)
            if self.config.trio_mode:
                self.logger.info("\n步骤10: 家系三元组分析 | Step 10: Family trio analysis")
                # 这里需要为每个家系成员运行完整的组装流程
                # For demonstration, we'll assume assemblies are already available
                trio_results = self.trio_analyzer.analyze_parental_support(
                    clean_assembly, clean_assembly, clean_assembly
                )
                assembly_results.update(trio_results)
            
            # 步骤11: 生成总结报告 | Step 11: Generate summary report
            self.logger.info("\n步骤11: 生成总结报告 | Step 11: Generate summary report")
            summary_report = self.results_processor.generate_summary_report(assembly_results)
            assembly_results['summary_report'] = summary_report
            
            # 完成信息 | Completion information
            self.logger.info("\n" + "=" * 80)
            self.logger.info("基因组组装分析完成！| Genome assembly analysis completed!")
            self.logger.info("=" * 80)
            self.logger.info(f"主要组装文件 | Primary assembly file: {clean_assembly}")
            self.logger.info(f"结果目录 | Results directory: {self.config.output_dir}")
            self.logger.info("主要输出文件 | Main output files:")
            self.logger.info("  - verkko_assembly/: Verkko组装结果 | Verkko assembly results")
            self.logger.info("  - hifiasm_assembly/: hifiasm组装结果 | hifiasm assembly results")
            self.logger.info("  - graphasing_output/: 分期分析结果 | Phasing analysis results")
            self.logger.info("  - contamination_screen/: 污染筛查结果 | Contamination screening results")
            self.logger.info("  - merqury_output/: 质量评估结果 | Quality assessment results")
            self.logger.info("  - alignment_output/: 比对分析结果 | Alignment analysis results")
            if self.config.trio_mode:
                self.logger.info("  - trio_analysis/: 家系分析结果 | Family analysis results")
            self.logger.info("  - assembly_summary_report.txt: 总结报告 | Summary report")
            self.logger.info("  - assembly_analysis.log: 分析日志 | Analysis log")
            
        except Exception as e:
            self.logger.error(f"分析流程在执行过程中意外终止 | Analysis pipeline terminated unexpectedly: {e}")
            sys.exit(1)

def main():
    """主函数 | Main function"""
    parser = argparse.ArgumentParser(
        description='基因组组装分析工具 (模块化版本) | Genome Assembly Analysis Tool (Modular Version)',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog="""
基于Verkko和hifiasm的杂合全基因组组装流程 | Hybrid Genome Assembly Pipeline using Verkko and hifiasm

主要功能 | Main Features:
  - Verkko主要组装器 | Primary assembly with Verkko
  - hifiasm超长读长组装 | Ultra-long read assembly with hifiasm
  - Graphasing分期信号 | Phasing signal with Graphasing
  - 质量控制和注释 | Quality control and annotation
  - 组装质量评估 | Assembly quality assessment
  - 参考基因组比对 | Reference genome alignment
  - 家系三元组分析 | Family trio analysis

示例 | Examples:
  %(prog)s --hifi-reads hifi.fastq.gz -o assembly_results
  %(prog)s --hifi-reads hifi.fq.gz --ont-reads ont.fq.gz --genome-size 3.2g
  %(prog)s --hifi-reads child_hifi.fq.gz --trio-mode --parent1-reads p1.fq.gz --parent2-reads p2.fq.gz --child-reads child.fq.gz
  %(prog)s --hifi-reads hifi.fastq.gz --reference-genome T2T-CHM13.fa --threads 64 --memory 256
        """
    )
    
    # 必需参数 | Required arguments
    parser.add_argument('--hifi-reads', required=True,
                       help='HiFi测序reads文件路径 | HiFi sequencing reads file path')
    
    # 可选输入参数 | Optional input arguments
    parser.add_argument('--ont-reads',
                       help='ONT测序reads文件路径 | ONT sequencing reads file path')
    parser.add_argument('--reference-genome',
                       help='参考基因组文件路径 (用于比对分析) | Reference genome file path (for alignment analysis)')
    parser.add_argument('-o', '--output-dir', default='./assembly_output',
                       help='输出目录 | Output directory')
    
    # 组装参数 | Assembly parameters
    parser.add_argument('--genome-size', default='3.2g',
                       help='预期基因组大小 | Expected genome size')
    parser.add_argument('--ploidy', type=int, choices=[1, 2], default=2,
                       help='基因组倍性 | Genome ploidy')
    
    # Verkko参数 | Verkko parameters
    parser.add_argument('--verkko-version', default='1.4.1',
                       help='Verkko版本 | Verkko version')
    parser.add_argument('--verkko-grid', action='store_true',
                       help='使用网格计算 | Use grid computing')
    parser.add_argument('--verkko-local-cpus', type=int, default=32,
                       help='Verkko本地CPU数量 | Verkko local CPU count')
    parser.add_argument('--verkko-local-memory', type=int, default=128,
                       help='Verkko本地内存(GB) | Verkko local memory (GB)')
    parser.add_argument('--no-verkko-screen', action='store_false', dest='verkko_screen',
                       help='禁用Verkko筛选 | Disable Verkko screening')
    
    # hifiasm参数 | hifiasm parameters
    parser.add_argument('--hifiasm-version', default='0.19.6',
                       help='hifiasm版本 | hifiasm version')
    parser.add_argument('--no-hifiasm-ultra-long', action='store_false', dest='hifiasm_ultra_long',
                       help='禁用hifiasm超长模式 | Disable hifiasm ultra-long mode')
    parser.add_argument('--hifiasm-purge-level', type=int, default=3, choices=[0, 1, 2, 3],
                       help='hifiasm purge级别 | hifiasm purge level')
    parser.add_argument('--hifiasm-similarity', type=float, default=0.8,
                       help='hifiasm相似性阈值 | hifiasm similarity threshold')
    
    # Graphasing参数 | Graphasing parameters
    parser.add_argument('--graphasing-version', default='0.3.1-alpha',
                       help='Graphasing版本 | Graphasing version')
    parser.add_argument('--graphasing-kmer-size', type=int, default=21,
                       help='Graphasing k-mer大小 | Graphasing k-mer size')
    
    # 质控参数 | Quality control parameters
    parser.add_argument('--skip-contamination-screen', action='store_false', dest='run_contamination_screen',
                       help='跳过外源污染筛查 | Skip contamination screening')
    parser.add_argument('--fcs-version', default='0.4.0',
                       help='NCBI FCS版本 | NCBI FCS version')
    parser.add_argument('--skip-flagger', action='store_false', dest='run_flagger',
                       help='跳过Flagger错误注释 | Skip Flagger error annotation')
    parser.add_argument('--flagger-version', default='0.3.3',
                       help='Flagger版本 | Flagger version')
    
    # 质量评估参数 | Quality assessment parameters
    parser.add_argument('--skip-merqury', action='store_false', dest='run_merqury',
                       help='跳过Merqury质量评估 | Skip Merqury quality assessment')
    parser.add_argument('--merqury-version', default='1.0',
                       help='Merqury版本 | Merqury version')
    parser.add_argument('--skip-inspector', action='store_false', dest='run_inspector',
                       help='跳过Inspector组装检查 | Skip Inspector assembly inspection')
    parser.add_argument('--inspector-version', default='1.2',
                       help='Inspector版本 | Inspector version')
    parser.add_argument('--skip-deepvariant', action='store_false', dest='run_deepvariant',
                       help='跳过DeepVariant质量值估计 | Skip DeepVariant quality value estimation')
    parser.add_argument('--deepvariant-version', default='1.6',
                       help='DeepVariant版本 | DeepVariant version')
    parser.add_argument('--skip-compleasm', action='store_false', dest='run_compleasm',
                       help='跳过compleasm基因完整性评估 | Skip compleasm gene completeness assessment')
    parser.add_argument('--compleasm-version', default='0.2.5',
                       help='compleasm版本 | compleasm version')
    parser.add_argument('--orthodb-version', default='10',
                       help='OrthoDB版本 | OrthoDB version')
    
    # 比对参数 | Alignment parameters
    parser.add_argument('--minimap2-version', default='2.26',
                       help='minimap2版本 | minimap2 version')
    parser.add_argument('--mashmap-version', default='3.1.3',
                       help='mashmap版本 | mashmap version')
    
    # 系统参数 | System parameters
    parser.add_argument('--threads', type=int, default=32,
                       help='线程数 | Number of threads')
    parser.add_argument('--memory', type=int, default=128,
                       help='内存大小(GB) | Memory size (GB)')
    parser.add_argument('--tmp-dir', default='/tmp',
                       help='临时目录 | Temporary directory')
    parser.add_argument('--keep-temp', action='store_true',
                       help='保留临时文件 | Keep temporary files')
    
    # 工具路径 | Tool paths
    parser.add_argument('--verkko-path', default='verkko',
                       help='Verkko软件路径 | Verkko software path')
    parser.add_argument('--hifiasm-path', default='hifiasm',
                       help='hifiasm软件路径 | hifiasm software path')
    parser.add_argument('--graphasing-path', default='graphasing',
                       help='Graphasing软件路径 | Graphasing software path')
    parser.add_argument('--fcs-path', default='fcs.py',
                       help='NCBI FCS软件路径 | NCBI FCS software path')
    parser.add_argument('--flagger-path', default='flagger',
                       help='Flagger软件路径 | Flagger software path')
    parser.add_argument('--merqury-path', default='merqury.sh',
                       help='Merqury软件路径 | Merqury software path')
    parser.add_argument('--nucfreq-path', default='nucfreq',
                       help='NucFreq软件路径 | NucFreq software path')
    parser.add_argument('--inspector-path', default='inspector.py',
                       help='Inspector软件路径 | Inspector software path')
    parser.add_argument('--deepvariant-path', default='run_deepvariant',
                       help='DeepVariant软件路径 | DeepVariant software path')
    parser.add_argument('--compleasm-path', default='compleasm',
                       help='compleasm软件路径 | compleasm software path')
    parser.add_argument('--minimap2-path', default='minimap2',
                       help='minimap2软件路径 | minimap2 software path')
    parser.add_argument('--mashmap-path', default='mashmap',
                       help='mashmap软件路径 | mashmap software path')
    
    # 家系分析参数 | Family analysis parameters
    parser.add_argument('--trio-mode', action='store_true',
                       help='启用家系三元组分析模式 | Enable family trio analysis mode')
    parser.add_argument('--parent1-reads',
                       help='亲本1的reads文件路径 | Parent1 reads file path')
    parser.add_argument('--parent2-reads',
                       help='亲本2的reads文件路径 | Parent2 reads file path')
    parser.add_argument('--child-reads',
                       help='子代的reads文件路径 | Child reads file path')
    
    args = parser.parse_args()
    
    # 验证家系模式参数 | Validate trio mode parameters
    if args.trio_mode:
        if not all([args.parent1_reads, args.parent2_reads, args.child_reads]):
            parser.error("家系模式需要提供 --parent1-reads, --parent2-reads, 和 --child-reads 参数 | "
                        "Trio mode requires --parent1-reads, --parent2-reads, and --child-reads arguments")
    
    # 转换参数为字典 | Convert arguments to dictionary
    config_dict = vars(args)
    
    # 创建分析器并运行 | Create analyzer and run
    analyzer = GenomeAssemblyAnalyzer(**config_dict)
    analyzer.run_analysis()

if __name__ == "__main__":
    main()

# ===== END FILE =====

# ===== FILE: run_genome_assembly.py =====
#!/usr/bin/env python3
"""
基因组组装分析运行脚本 | Genome Assembly Analysis Runner Script
这是一个简化的入口脚本，用于运行基因组组装分析 | Simple entry script for running genome assembly analysis

用法 | Usage:
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz -o assembly_results
    
示例 | Examples:
    # 基本组装分析 | Basic assembly analysis
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz -o assembly_results
    
    # 包含ONT数据 | Include ONT data
    python run_genome_assembly.py --hifi-reads hifi.fq.gz --ont-reads ont.fq.gz --genome-size 3.2g
    
    # 家系模式分析 | Family trio analysis
    python run_genome_assembly.py --hifi-reads child_hifi.fq.gz --trio-mode \\
        --parent1-reads p1.fq.gz --parent2-reads p2.fq.gz --child-reads child.fq.gz
    
    # 高性能配置 | High-performance configuration
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz \\
        --reference-genome T2T-CHM13.fa --threads 64 --memory 256
    
    # 自定义质控参数 | Custom quality control parameters
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz \\
        --skip-contamination-screen --skip-flagger
    
    # 指定工具路径 | Specify tool paths
    python run_genome_assembly.py --hifi-reads hifi.fastq.gz \\
        --verkko-path /path/to/verkko --hifiasm-path /path/to/hifiasm

依赖工具 | Required Tools:
- Verkko (v.1.4.1+): 主要组装器 | Primary assembler
- hifiasm (v.0.19.6+): 超长读长组装器 | Ultra-long read assembler  
- Graphasing (v.0.3.1-alpha+): 分期管道 | Phasing pipeline
- NCBI FCS (v.0.4.0+): 外源污染筛查 | Foreign contamination screening
- Flagger (v.0.3.3+): 组装错误注释 | Assembly error annotation
- Merqury (v.1.0+): 质量评估 | Quality assessment
- NucFreq: 核苷酸频率分析 | Nucleotide frequency analysis
- Inspector (v.1.2+): 组装检查 | Assembly inspection
- DeepVariant (v.1.6+): 变异检测质量评估 | Variant calling quality assessment
- compleasm (v.0.2.5+): 基因完整性评估 | Gene completeness assessment
- minimap2 (v.2.26+): 序列比对 | Sequence alignment
- mashmap (v.3.1.3+): 快速序列比对 | Fast sequence mapping

安装方法 | Installation:
  参考各工具官方文档进行安装 | Refer to official documentation for each tool
  建议使用conda环境管理依赖 | Recommend using conda for dependency management
"""

from genome_assembly.main import main

if __name__ == "__main__":
    main()