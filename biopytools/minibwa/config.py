"""
Minibwa配置管理模块|Minibwa Configuration Management Module
"""

import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

from ..common.paths import expand_path, get_tool_path


# 支持的运行模式|Supported run modes
MODE_STANDARD = 'standard'   # 标准短读长PE/SE|Standard short-read PE/SE
MODE_HIC = 'hic'             # Hi-C短读长|Hi-C short reads
MODE_METH = 'meth'           # 方向性BS-seq|Directional BS-seq
MODE_LONG = 'long'           # 准确长读长|Accurate long reads

VALID_MODES = (MODE_STANDARD, MODE_HIC, MODE_METH, MODE_LONG)


@dataclass
class MinibwaConfig:
    """Minibwa配置类|Minibwa Configuration Class"""

    # 必需参数|Required parameters
    genome: str                 # 参考基因组FASTA|Reference genome FASTA
    input_dir: str              # FASTQ输入目录|FASTQ input directory
    pattern: str = '_1.fq.gz'   # R1文件匹配模式|R1 file pattern

    # 输出参数|Output parameters
    output_dir: str = './minibwa_output'

    # 比对模式|Alignment mode (standard/hic/meth/long)
    mode: str = MODE_STANDARD

    # 性能参数|Performance parameters
    threads: int = 12

    # Minibwa map 共用参数|Minibwa map common parameters
    preset: str = 'adap'        # -x 预设: sr/lr/adap|preset
    min_seed: int = 19          # -k 最小种子长度|min seed length
    max_occ: int = 250          # -c 最大种子出现次数|max seed occurrences
    max_gap: int = 100          # -g 最大gap|max gap size
    bandwidth: int = 100        # -w 带宽|bandwidth
    bandwidth_long: int = 20000 # -W 长读带宽|long bandwidth
    min_chain_score: int = 25   # -m 最小链接分数|min chaining score
    sec_ratio: float = 0.5      # -p 次要/主要得分比|min secondary-to-primary ratio
    max_sec: int = 50           # -N 保留次要比对数|retain N secondary alignments
    min_dp_score: int = 30      # -s 最小DP得分*s*{-A}|min DP score (multiplied by -A)

    # Minibwa map 比对打分参数|Minibwa map scoring parameters
    match_score: int = 2        # -A 匹配得分|matching score
    mismatch_penalty: int = 8   # -B 错配罚分|mismatching penalty
    gap_open: str = '12,23'     # -O gap开放罚分|gap open penalty
    gap_ext: str = '2,1'        # -E gap延伸罚分|gap extension penalty

    # 输入输出选项|I/O options
    read_group: Optional[str] = None  # -R SAM read group line
    no_unmap: bool = False      # -u 不输出未比对read|don't output unmapped reads
    out_n: int = 0              # --outn 输出次要比对上限|max secondary to output
    copy_comment: bool = False  # -y 复制FASTA/Q注释|copy comments
    soft_clip_supp: bool = False  # -Y 软剪切补充比对|soft clipping for supp
    batch_size: str = '100m,1g' # -K 批处理大小|batch size

    # 后处理参数|Post-processing parameters
    markdup: bool = False       # 标记重复|Mark duplicates
    remove_dup: bool = False    # 移除重复（需开启markdup）|Remove duplicates

    # 覆盖度参数|Coverage parameters
    min_base_quality: int = 0   # samtools depth -q
    min_mapping_quality: int = 0  # samtools depth -Q
    max_depth: int = 0          # samtools depth -d (0=无限|unlimited)

    # 滑窗覆盖度参数|Windowed coverage parameters
    window_size: int = 1000000  # 窗口大小|Window size in bp
    step_size: int = 100000     # 步长|Step size in bp

    # 运行控制|Run control
    resume: bool = False        # 断点续传|Resume (skip completed samples)
    skip_coverage: bool = False  # 跳过覆盖度分析|Skip coverage analysis

    # 工具路径|Tool paths
    minibwa_path: str = field(
        default_factory=lambda: get_tool_path(
            'minibwa',
            '~/software/minibwa/minibwa',
            'MINIBWA_PATH'
        )
    )
    samtools_path: str = field(
        default_factory=lambda: get_tool_path(
            'samtools',
            '~/.local/bin/samtools',
            'SAMTOOLS_PATH'
        )
    )

    def __post_init__(self):
        """初始化后处理|Post-initialization processing"""
        # 关键：展开所有含~的路径|CRITICAL: Expand all ~ paths
        self.genome = expand_path(self.genome)
        self.input_dir = expand_path(self.input_dir)
        self.output_dir = expand_path(self.output_dir)
        self.minibwa_path = expand_path(self.minibwa_path)
        self.samtools_path = expand_path(self.samtools_path)

        # 输出根目录|Output root
        self.output_path = Path(self.output_dir)
        self.output_path.mkdir(parents=True, exist_ok=True)

        # 流程元数据目录|Pipeline metadata directory
        self.pipeline_info_dir = self.output_path / '00_pipeline_info'

        # 索引目录（流程构建一次，多样本共享）|Index directory (built once, shared)
        self.index_dir = self.output_path / '01_index'
        self.index_prefix = self.index_dir / 'ref'

        # 全局日志目录|Global log directory
        self.log_dir = self.output_path / '99_logs'

        for d in [self.pipeline_info_dir, self.index_dir, self.log_dir]:
            d.mkdir(parents=True, exist_ok=True)

        # 衍生路径|Derived paths
        self.genome_name = Path(self.genome).name

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        # 必需文件存在性|Required file existence
        if not os.path.exists(self.genome):
            errors.append(f"基因组文件不存在|Genome file not found: {self.genome}")
        if not os.path.isdir(self.input_dir):
            errors.append(f"输入目录不存在|Input directory not found: {self.input_dir}")

        # 工具可执行性|Tool executability
        if not os.path.exists(self.minibwa_path):
            errors.append(f"minibwa不存在|minibwa not found: {self.minibwa_path}")
        if not os.path.exists(self.samtools_path):
            errors.append(f"samtools不存在|samtools not found: {self.samtools_path}")

        # 参数范围|Parameter ranges
        if self.threads <= 0:
            errors.append(f"线程数必须为正|Threads must be positive: {self.threads}")
        if self.mode not in VALID_MODES:
            errors.append(
                f"模式无效|Invalid mode: {self.mode}; "
                f"有效值|valid: {', '.join(VALID_MODES)}"
            )
        if self.window_size <= 0 or self.step_size <= 0:
            errors.append("窗口大小和步长必须为正|Window size and step must be positive")
        if self.step_size > self.window_size:
            errors.append("步长不能大于窗口大小|Step cannot be larger than window size")
        if self.remove_dup and not self.markdup:
            errors.append("移除重复需要先开启markdup|--remove-dup requires --markdup")

        if errors:
            raise ValueError("\n".join(errors))

        return True

    def get_map_options(self) -> list:
        """生成minibwa map命令参数列表|Generate minibwa map argument list"""
        args = []

        # 线程|Threads
        args.extend(['-t', str(self.threads)])

        # 模式相关选项|Mode-specific options
        if self.mode == MODE_HIC:
            args.append('--hic')
        elif self.mode == MODE_METH:
            args.append('--meth')
        elif self.mode == MODE_LONG:
            args.append('-f')   # PAF输出|PAF output

        # 预设|Preset
        args.extend(['-x', self.preset])

        # 共用比对参数|Common mapping parameters
        args.extend(['-k', str(self.min_seed)])
        args.extend(['-c', str(self.max_occ)])
        args.extend(['-g', str(self.max_gap)])
        args.extend(['-w', str(self.bandwidth)])
        args.extend(['-W', str(self.bandwidth_long)])
        args.extend(['-m', str(self.min_chain_score)])
        args.extend(['-p', str(self.sec_ratio)])
        args.extend(['-N', str(self.max_sec)])
        args.extend(['-s', str(self.min_dp_score)])

        # 打分参数|Scoring parameters
        args.extend(['-A', str(self.match_score)])
        args.extend(['-B', str(self.mismatch_penalty)])
        args.extend(['-O', self.gap_open])
        args.extend(['-E', self.gap_ext])

        # IO参数|I/O options
        if self.no_unmap:
            args.append('-u')
        if self.out_n > 0:
            args.extend(['--outn', str(self.out_n)])
        if self.copy_comment:
            args.append('-y')
        if self.soft_clip_supp:
            args.append('-Y')
        args.extend(['-K', self.batch_size])

        # Read group|Read group
        if self.read_group:
            args.extend(['-R', self.read_group])

        return args

    def is_meth_mode(self) -> bool:
        """是否BS-seq模式|Is BS-seq mode"""
        return self.mode == MODE_METH

    def is_long_mode(self) -> bool:
        """是否长读模式|Is long-read mode"""
        return self.mode == MODE_LONG

    def to_param_dict(self) -> dict:
        """导出参数字典（用于软件版本记录）|Export param dict for version log"""
        return {
            'genome': self.genome,
            'input_dir': self.input_dir,
            'pattern': self.pattern,
            'mode': self.mode,
            'threads': self.threads,
            'preset': self.preset,
            'min_seed': self.min_seed,
            'max_occ': self.max_occ,
            'max_gap': self.max_gap,
            'bandwidth': self.bandwidth,
            'bandwidth_long': self.bandwidth_long,
            'min_chain_score': self.min_chain_score,
            'sec_ratio': self.sec_ratio,
            'max_sec': self.max_sec,
            'min_dp_score': self.min_dp_score,
            'match_score': self.match_score,
            'mismatch_penalty': self.mismatch_penalty,
            'gap_open': self.gap_open,
            'gap_ext': self.gap_ext,
            'markdup': self.markdup,
            'remove_dup': self.remove_dup,
            'window_size': self.window_size,
            'step_size': self.step_size,
        }
