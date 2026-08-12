"""
CentIER着丝粒鉴定核心分析模块|CentIER Centromere Identification Core Analysis Module
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Optional
from .config import CentIERConfig
from .utils import run_command, format_number


class CentIERAnalyzer:
    """CentIER着丝粒鉴定分析器|CentIER Centromere Identification Analyzer"""

    def __init__(self, config: CentIERConfig, logger):
        """
        初始化分析器|Initialize analyzer

        Args:
            config: CentIER配置对象|CentIERConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def _check_chrnaming(self) -> bool:
        """检查染色体命名是否符合 CentIER 要求的 ChrN 格式|Check chr naming matches CentIER's ChrN format

        CentIER 实际接受"以 ChrN 结尾"的 id(如 So120_Chr1 实测可跑),Chr_1/scaffold_1 会崩。|
        CentIER accepts ids ending in ChrN (e.g. So120_Chr1 works); Chr_1/scaffold_1 crash.
        默认 WARNING + 继续;strict_chrname=True 时中止。|Warns by default; aborts if strict_chrname.

        Returns:
            bool: 是否可继续(False 仅在 strict 模式命名不符时)|Whether to proceed
        """
        import re
        # 允许 So120_Chr1 / Chr1 等以 Chr+字母数字结尾的 id;拒绝 scaffold_1 / Chr_1|
        # allow ids ending in Chr+alnum (So120_Chr1, Chr1); reject scaffold_1, Chr_1
        pattern = re.compile(r'Chr[a-zA-Z0-9]+$')
        bad = []
        with open(self.config.genome_fasta) as f:
            for line in f:
                if line.startswith('>'):
                    seqid = line[1:].strip().split()[0]
                    if not pattern.search(seqid):
                        bad.append(seqid)
        if bad:
            shown = bad[:5] + (['...'] if len(bad) > 5 else [])
            msg = (f"染色体命名不符合 ChrN 格式,CentIER 可能报错|Chr naming not ChrN format "
                   f"({len(bad)} 条): {shown}; 建议用 chr_rename/fasta_id_renamer 规范化|"
                   f"consider chr_rename/fasta_id_renamer")
            if self.config.strict_chrname:
                self.logger.error(msg)
                return False
            self.logger.warning(msg)
        else:
            self.logger.info("染色体命名检查通过|Chr naming check passed (ChrN format)")
        return True

    def _resolve_hicpro_matrices(self) -> bool:
        """从 HiC-Pro 输出定位 100kb/20kb 矩阵+bed,回填 config|Locate matrices from HiC-Pro output, fill config

        HiC-Pro 输出结构: hic_results/matrix/{sample}/{raw,iced}/{bin}/{sample}_{bin}[_iced].matrix
        bed 始终在 raw/ 下(与归一化无关)。|bed always under raw/ (independent of normalization).

        Returns:
            bool: 是否成功定位全部 4 个文件|Whether all 4 files located
        """
        sample = self.config.sample_name
        base = (self.config.output_path / '01_hic_mapping' / 'hicpro_output'
                / 'hic_results' / 'matrix' / sample)
        mtype = self.config.hic_matrix_type

        def matrix_path(bin_size: str):
            if mtype == 'iced':
                return base / 'iced' / bin_size / f'{sample}_{bin_size}_iced.matrix'
            return base / 'raw' / bin_size / f'{sample}_{bin_size}.matrix'

        def bed_path(bin_size: str):
            # bed 始终在 raw/ 下|bed always under raw/
            return base / 'raw' / bin_size / f'{sample}_{bin_size}_abs.bed'

        targets = {
            'matrix1': matrix_path('100000'),
            'bed1': bed_path('100000'),
            'matrix2': matrix_path('20000'),
            'bed2': bed_path('20000'),
        }
        missing = [k for k, p in targets.items() if not p.exists()]
        if missing:
            self.logger.error(f"HiC-Pro 输出缺少矩阵文件|Missing HiC-Pro matrix files: {missing}")
            self.logger.error(f"HiC-Pro 输出目录|HiC-Pro output base: {base}")
            return False
        self.config.matrix1 = str(targets['matrix1'])
        self.config.matrix2 = str(targets['matrix2'])
        self.config.bed1 = str(targets['bed1'])
        self.config.bed2 = str(targets['bed2'])
        for k, p in targets.items():
            self.logger.info(f"定位 Hi-C 矩阵|Located {k}: {p}")
        return True

    def _run_hicpro(self) -> bool:
        """跑 HiC-Pro 产出矩阵(复用 hic_heatmap 的 HiCProPipeline)|Run HiC-Pro via hic_heatmap pipeline

        输出落到 output_path/01_hic_mapping/。|Output goes to output_path/01_hic_mapping/.

        Returns:
            bool: 是否成功|Whether successful
        """
        from ..hic_heatmap.config import HiCProConfig
        from ..hic_heatmap.hicpro_pipeline import HiCProPipeline

        hicpro_out = self.config.output_path / '01_hic_mapping'
        hicpro_out.mkdir(parents=True, exist_ok=True)

        hc_config = HiCProConfig(
            genome=self.config.genome_fasta,
            output_dir=str(hicpro_out),
            genome_id=self.config.genome_id,
            fastq_r1=self.config.fastq_r1,
            fastq_r2=self.config.fastq_r2,
            threads=self.config.threads,
            max_memory_gb=self.config.max_memory_gb,
            restriction_enzyme=self.config.restriction_enzyme,
            bowtie2_idx=self.config.bowtie2_idx,
            bin_sizes=self.config.bin_sizes,
            skip_existing=not self.config.force_hicpro,
            force_hicpro=self.config.force_hicpro,
            require_plothic=False,  # centier 不画热图|centier does not plot
        )
        try:
            hc_config.validate()
        except ValueError as e:
            self.logger.error(f"HiC-Pro 配置校验失败|HiC-Pro config validation failed: {e}")
            return False

        self.logger.info("启动 HiC-Pro 产生 Hi-C 矩阵|Starting HiC-Pro to generate Hi-C matrices")
        pipeline = HiCProPipeline(hc_config, self.logger)
        return pipeline.run_hicpro_only()

    def run(self) -> Dict:
        """
        运行完整的着丝粒鉴定流程|Run complete centromere identification pipeline

        Returns:
            dict: 分析结果统计|Analysis result statistics
        """
        self.logger.info("=" * 60)
        self.logger.info("开始CentIER着丝粒鉴定|Starting CentIER centromere identification")
        self.logger.info("=" * 60)
        self.logger.info(f"输入基因组|Input genome: {self.config.genome_fasta}")
        self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
        self.logger.info(f"线程数|Threads: {self.config.threads}")
        self.logger.info(f"K-mer大小|K-mer size: {self.config.kmer_size}")
        self.logger.info(f"中心容差|Center tolerance: {self.config.center_tolerance}")
        self.logger.info(f"步长|Step length: {self.config.step_len}")

        # Hi-C 自动模式:ChrN 预检 → HiC-Pro → 矩阵定位|Hi-C auto mode preprocessing
        if self.config.fastq_r1:
            self.logger.info("检测到 Hi-C FASTQ 输入,启用自动模式|Hi-C FASTQ detected, auto mode enabled")
            if not self._check_chrnaming():
                return {'success': False, 'error': 'Chr naming check failed (strict mode)'}
            if not self._run_hicpro():
                return {'success': False, 'error': 'HiC-Pro failed'}
            if not self._resolve_hicpro_matrices():
                return {'success': False, 'error': 'HiC-Pro matrix resolution failed'}
            self.logger.info("Hi-C 矩阵已就绪,继续 CentIER|Hi-C matrices ready, proceeding to CentIER")

        try:
            # 临时修改CentIER脚本以使用自定义线程数|Temporarily modify CentIER script to use custom thread count
            centier_script = self.config.get_centier_script_path()
            modified_script = self._patch_centier_script(centier_script)

            # 构建命令|Build command
            cmd, env = self._build_command(modified_script)

            self.logger.info(f"运行CentIER|Running CentIER...")
            self.logger.info(f"命令|Command: {' '.join(cmd)}")
            self.logger.debug(f"工作目录|Working directory: {os.getcwd()}")

            # 记录开始时间用于版本信息|Record start time for version info
            import time
            start_time = time.time()

            # 运行命令（不设置cwd，使用当前目录）|Run command (without cwd, use current directory)
            try:
                result = subprocess.run(
                    cmd,
                    capture_output=True,
                    text=True,
                    check=False,
                    env=env
                )
            finally:
                # 清理临时脚本|Clean up temporary script
                if modified_script != centier_script and os.path.exists(modified_script):
                    os.remove(modified_script)

            # 输出标准输出|Output stdout
            if result.stdout:
                for line in result.stdout.split('\n'):
                    if line.strip():
                        self.logger.info(line)

            # 输出标准错误|Output stderr
            if result.stderr:
                for line in result.stderr.split('\n'):
                    if line.strip():
                        self.logger.warning(line)

            # 检查是否有关键输出文件生成，即使返回码非0也可能成功
            # CentIER可能在某些警告后仍生成结果文件
            output_files = self._check_outputs()

            # 清理CentIER在当前目录生成的临时文件|Clean up CentIER temp files in CWD
            self._cleanup_temp_files()

            # 如果有输出文件生成，视为成功；否则检查返回码
            if result.returncode != 0 and not output_files:
                self.logger.error(f"CentIER运行失败|CentIER run failed with return code: {result.returncode}")
                self.logger.error(f"未生成任何输出文件|No output files generated")
                return {
                    'success': False,
                    'error': f"Return code: {result.returncode}",
                    'stderr': result.stderr
                }

            if result.returncode != 0 and output_files:
                self.logger.warning(f"CentIER返回非零退出码({result.returncode})，但已生成输出文件|"
                                  f"CentIER returned non-zero exit code ({result.returncode}), but output files were generated")

            # 生成软件版本信息|Generate software versions
            self._generate_software_versions()

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("CentIER着丝粒鉴定完成|CentIER centromere identification completed")
            self.logger.info("=" * 60)

            return {
                'success': True,
                'output_files': output_files,
                'genome': self.config.genome_fasta,
                'output_dir': str(self.config.get_centier_output_dir())
            }

        except Exception as e:
            self.logger.error(f"着丝粒鉴定失败|Centromere identification failed: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return {
                'success': False,
                'error': str(e)
            }

    def _patch_centier_script(self, original_script: str) -> str:
        """
        创建临时修改的CentIER脚本，使用自定义线程数|Create temporary modified CentIER script with custom thread count

        Args:
            original_script: 原始脚本路径|Original script path

        Returns:
            str: 修改后的脚本路径|Modified script path
        """
        import tempfile
        import shutil

        # 读取原始脚本|Read original script
        with open(original_script, 'r') as f:
            script_content = f.read()

        # 获取CentIER目录|Get CentIER directory
        centier_dir = os.path.dirname(original_script)

        # 替换硬编码的60线程|Replace hardcoded 60 threads
        script_content = script_content.replace(
            "'-threads','60'",
            f"'-threads','{self.config.threads}'"
        )

        # 修改script_path变量，使其指向原始CentIER目录
        # 这样即使临时脚本在/tmp/，工具路径仍然正确
        script_content = script_content.replace(
            "script_path = os.path.dirname(os.path.abspath(__file__))",
            f"script_path = '{centier_dir}'"
        )

        # 创建临时文件|Create temporary file
        tmp_root = os.path.join(str(self.config.get_centier_output_dir()), 'tmp')
        os.makedirs(tmp_root, exist_ok=True)
        fd, temp_path = tempfile.mkstemp(suffix='_centier.py', text=True, dir=tmp_root)
        with os.fdopen(fd, 'w') as f:
            f.write(script_content)

        # 设置可执行权限|Set executable permission
        os.chmod(temp_path, 0o755)

        self.logger.debug(f"创建临时脚本|Created temporary script: {temp_path}")
        self.logger.debug(f"CentIER目录|CentIER directory: {centier_dir}")
        return temp_path

    def _build_command(self, script_path: str = None):
        """
        构建CentIER命令|Build CentIER command

        Args:
            script_path: 脚本路径（可选）|Script path (optional)

        Returns:
            tuple: (命令列表|Command list, 环境变量字典|Environment variables dict)
        """
        import os
        if script_path is None:
            script_path = self.config.get_centier_script_path()

        # 确保使用绝对路径|Ensure absolute paths
        genome_abs = os.path.abspath(self.config.genome_fasta)
        output_abs = os.path.abspath(str(self.config.get_centier_output_dir()))

        cmd = [sys.executable, script_path, genome_abs]

        # 添加输出目录|Add output directory
        cmd.extend(['-o', output_abs])

        # 添加GFF注释（如果有）|Add GFF annotation (if provided)
        if self.config.gff_annotation:
            cmd.extend(['--gff', os.path.abspath(self.config.gff_annotation)])

        # 添加参数|Add parameters
        cmd.extend(['-k', str(self.config.kmer_size)])
        cmd.extend(['-c', str(self.config.center_tolerance)])
        cmd.extend(['--step_len', str(self.config.step_len)])

        # 添加多着丝粒选项|Add multiple centromeres option
        if self.config.mul_cents:
            cmd.append('--mul_cents')

        # 添加Hi-C数据（如果有）|Add Hi-C data (if provided)
        if self.config.matrix1 and self.config.matrix2 and self.config.bed1 and self.config.bed2:
            cmd.extend(['--matrix1', os.path.abspath(self.config.matrix1)])
            cmd.extend(['--matrix2', os.path.abspath(self.config.matrix2)])
            cmd.extend(['--bed1', os.path.abspath(self.config.bed1)])
            cmd.extend(['--bed2', os.path.abspath(self.config.bed2)])
            cmd.extend(['--MINGAP', str(self.config.mingap)])
            cmd.extend(['--SIGNAL_THRESHOLD', str(self.config.signal_threshold)])

        # 设置环境变量|Set environment variables
        import os as os_module
        env = os_module.environ.copy()

        # 添加CentIER目录到PYTHONPATH，以便找到translate_seq等模块
        # Add CentIER directory to PYTHONPATH to find translate_seq and other modules
        centier_dir = os.path.dirname(self.config.get_centier_script_path())
        if 'PYTHONPATH' in env:
            env['PYTHONPATH'] = f"{centier_dir}:{env['PYTHONPATH']}"
        else:
            env['PYTHONPATH'] = centier_dir

        return cmd, env

    def _check_outputs(self) -> Dict[str, str]:
        """
        检查输出文件|Check output files

        Returns:
            dict: 输出文件路径字典|Output file path dictionary
        """
        prefix = os.path.basename(self.config.genome_fasta)
        output_files = {}

        # 主要输出文件|Main output files
        expected_files = [
            ('centromere_range', f'{prefix}_centromere_range.txt'),
            ('centromere_seq', f'{prefix}_all_centromere_seq.txt'),
            ('monomer_seq', f'{prefix}_monomer_seq.txt'),
            ('monomer_in_centromere', f'{prefix}_monomer_in_centromere.txt'),
            ('ltr_position', f'{prefix}_ltr_position.txt'),
            ('ltr_statistics', f'{prefix}_LTR_statistics.txt'),
            ('visualization', f'{prefix}_draw_cen.svg')
        ]

        for file_type, file_name in expected_files:
            file_path = self.config.get_centier_output_dir() / file_name
            if file_path.exists():
                output_files[file_type] = str(file_path)
                self.logger.info(f"找到输出文件|Output file found: {file_type} -> {file_path}")
            else:
                self.logger.warning(f"输出文件未找到|Output file not found: {file_type} -> {file_path}")

        return output_files

    def _cleanup_temp_files(self):
        """
        清理CentIER在当前工作目录生成的临时文件|Clean up CentIER temp files in CWD

        CentIER上游脚本在CWD生成大量临时文件(以基因组basename为前缀)|
        CentIER upstream script generates many temp files in CWD (prefixed with genome basename)
        """
        import glob
        import shutil

        prefix = os.path.basename(self.config.genome_fasta)
        patterns = [
            f"{prefix}.*",                  # TRF输出|TRF output
            f"{prefix}_LTR_database*",      # LTR数据库|LTR database
            f"{prefix}_l1.txt",             # ltr_finder输出|LTR finder output
            f"{prefix}_centromere_seq.txt", # 中间序列|Intermediate sequence
        ]

        tmp_dir = os.path.join(str(self.config.get_centier_output_dir()), 'tmp')
        moved_count = 0

        for pattern in patterns:
            for file_path in glob.glob(pattern):
                if os.path.isfile(file_path):
                    try:
                        os.makedirs(tmp_dir, exist_ok=True)
                        dest = os.path.join(tmp_dir, os.path.basename(file_path))
                        shutil.move(file_path, dest)
                        moved_count += 1
                        self.logger.debug(f"移动临时文件|Moved temp file: {file_path} -> {dest}")
                    except OSError as e:
                        self.logger.warning(f"无法移动临时文件|Cannot move temp file {file_path}: {e}")
                elif os.path.isdir(file_path) and 'LTR_database' in file_path:
                    try:
                        os.makedirs(tmp_dir, exist_ok=True)
                        dest = os.path.join(tmp_dir, os.path.basename(file_path))
                        if os.path.exists(dest):
                            shutil.rmtree(dest)
                        shutil.move(file_path, dest)
                        moved_count += 1
                        self.logger.debug(f"移动临时目录|Moved temp dir: {file_path} -> {dest}")
                    except OSError as e:
                        self.logger.warning(f"无法移动临时目录|Cannot move temp dir {file_path}: {e}")

        if moved_count > 0:
            self.logger.info(f"已移动{format_number(moved_count)}个临时文件到tmp/|"
                           f"Moved {format_number(moved_count)} temp files to tmp/")

    def _generate_software_versions(self):
        """
        生成软件版本信息文件|Generate software versions info file

        记录CentIER版本、参数和运行环境|Record CentIER version, parameters, and runtime environment
        """
        import time
        import importlib

        pipeline_info_dir = self.config.output_path / '00_pipeline_info'
        pipeline_info_dir.mkdir(parents=True, exist_ok=True)

        versions = {
            'pipeline': {
                'name': 'biopytools centier',
                'version': '1.0.0'
            },
            'centier': {
                'version': '3.0.1',
                'path': self.config.get_centier_script_path(),
            },
            'python_packages': {}
        }

        # 记录Python包版本|Record Python package versions
        for pkg in ['pyfastx', 'numpy', 'pandas', 'scipy']:
            try:
                mod = importlib.import_module(pkg)
                if hasattr(mod, '__version__'):
                    pkg_version = mod.__version__
                elif hasattr(mod, 'version'):
                    pkg_version = mod.version
                else:
                    pkg_version = 'unknown'
                versions['python_packages'][pkg] = pkg_version
            except ImportError:
                versions['python_packages'][pkg] = 'not found'

        # 记录运行参数|Record run parameters
        versions['parameters'] = {
            'genome_fasta': os.path.basename(self.config.genome_fasta),
            'output_dir': str(self.config.output_path),
            'threads': self.config.threads,
            'kmer_size': self.config.kmer_size,
            'center_tolerance': self.config.center_tolerance,
            'step_len': self.config.step_len,
            'mul_cents': self.config.mul_cents,
        }

        if self.config.gff_annotation:
            versions['parameters']['gff_annotation'] = os.path.basename(self.config.gff_annotation)

        # 写入YAML文件|Write YAML file
        import yaml
        version_file = pipeline_info_dir / 'software_versions.yml'
        try:
            with open(version_file, 'w') as f:
                yaml.dump(versions, f, default_flow_style=False, allow_unicode=True)
            self.logger.info(f"软件版本信息已保存|Software versions saved to: {version_file}")
        except Exception as e:
            self.logger.warning(f"无法写入版本文件|Cannot write versions file: {e}")

    def run_step(self, step: int) -> bool:
        """
        运行单个步骤|Run single step

        注意: CentIER原版不支持单步运行，此函数仅用于兼容性
        Note: Original CentIER does not support single-step execution,
              this function is for compatibility only

        Args:
            step: 步骤编号|Step number (1-6)

        Returns:
            bool: 是否成功|Success
        """
        self.logger.warning(f"CentIER不支持单步运行，将运行完整流程|"
                           "CentIER does not support single-step execution, running full pipeline")
        result = self.run()
        return result.get('success', False)

    def get_summary(self) -> Dict:
        """
        获取分析结果摘要|Get analysis result summary

        Returns:
            dict: 结果摘要|Result summary
        """
        prefix = os.path.basename(self.config.genome_fasta)

        summary = {
            'genome_fasta': self.config.genome_fasta,
            'output_dir': str(self.config.output_path),
            'parameters': {
                'kmer_size': self.config.kmer_size,
                'center_tolerance': self.config.center_tolerance,
                'step_len': self.config.step_len,
                'mul_cents': self.config.mul_cents
            }
        }

        # 读取着丝粒范围文件|Read centromere range file
        range_file = self.config.get_centier_output_dir() / f'{prefix}_centromere_range.txt'
        if range_file.exists():
            centromeres = []
            with open(range_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 3:
                            centromeres.append({
                                'chromosome': parts[0],
                                'start': int(parts[1]),
                                'end': int(parts[2])
                            })
            summary['centromeres'] = centromeres
            summary['centromere_count'] = len(centromeres)

        return summary
