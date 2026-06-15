"""
CentIER着丝粒鉴定核心分析模块|CentIER Centromere Identification Core Analysis Module
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Dict, Optional
from .config import CentIERConfig
from .utils import run_command


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

        try:
            # 临时修改CentIER脚本以使用自定义线程数|Temporarily modify CentIER script to use custom thread count
            centier_script = self.config.get_centier_script_path()
            modified_script = self._patch_centier_script(centier_script)

            # 构建命令|Build command
            cmd, env = self._build_command(modified_script)

            self.logger.info(f"运行CentIER|Running CentIER...")
            self.logger.debug(f"工作目录|Working directory: {os.getcwd()}")

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

            # 完成|Complete
            self.logger.info("=" * 60)
            self.logger.info("CentIER着丝粒鉴定完成|CentIER centromere identification completed")
            self.logger.info("=" * 60)

            return {
                'success': True,
                'output_files': output_files,
                'genome': self.config.genome_fasta,
                'output_dir': str(self.config.output_path)
            }

        except Exception as e:
            self.logger.error(f"着丝粒鉴定失败|Centromere identification failed: {e}")
            import traceback
            traceback.print_exc()
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
        fd, temp_path = tempfile.mkstemp(suffix='_centier.py', text=True)
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
        output_abs = os.path.abspath(str(self.config.output_path))

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
            file_path = self.config.output_path / file_name
            if file_path.exists():
                output_files[file_type] = str(file_path)
                self.logger.info(f"找到输出文件|Output file found: {file_type} -> {file_path}")
            else:
                self.logger.warning(f"输出文件未找到|Output file not found: {file_type} -> {file_path}")

        return output_files

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
        range_file = self.config.output_path / f'{prefix}_centromere_range.txt'
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
