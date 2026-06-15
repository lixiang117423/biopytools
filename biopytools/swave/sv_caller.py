"""
Swave结构变异检测核心模块|Swave Structural Variant Detection Core Module
"""

import os
import sys
import subprocess
import re
from typing import Optional, List, Dict
from .config import SwaveConfig
from .utils import SwaveLogger, CommandRunner


class SwaveSVCaller:
    """Swave结构变异检测器|Swave Structural Variant Caller"""

    def __init__(self, config: SwaveConfig, logger, cmd_runner: CommandRunner):
        """
        初始化SV检测器|Initialize SV caller

        Args:
            config: Swave配置对象|Swave config object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    @staticmethod
    def build_call_args(config) -> List[str]:
        """
        构建Swave call子命令的参数列表|Build Swave call subcommand argument list

        Args:
            config: SwaveConfig配置对象|SwaveConfig configuration object

        Returns:
            Swave.py的参数列表（不含python和脚本路径）|Argument list for Swave.py
        """
        args = ['call']

        # 必需参数|Required parameters
        args.extend(['--input_path', config.assemblies_tsv])
        args.extend(['--ref_path', config.ref_fasta])
        args.extend(['--gfa_path', config.gfa_file])
        args.extend(['--gfa_source', config.gfa_source])
        args.extend(['--output_path', config.output_dir])

        # 可选参数|Optional parameters
        if config.decomposed_vcf:
            args.extend(['--decomposed_vcf', config.decomposed_vcf])

        if config.output_mode != 'auto':
            args.extend(['--output_mode', config.output_mode])

        if config.spec_samples and config.spec_samples != ["all"]:
            args.extend(['--spec_samples'] + config.spec_samples)

        if config.spec_snarl:
            args.extend(['--spec_snarl', config.spec_snarl])

        if config.spec_path:
            args.extend(['--spec_path', config.spec_path])

        # SV检测参数|SV detection parameters
        args.extend(['--min_sv_size', str(config.min_sv_size)])
        args.extend(['--max_sv_size', str(config.max_sv_size)])
        args.extend(['--max_sv_comps', str(config.max_sv_comps)])

        # 处理选项|Processing options
        if config.dup_to_ins:
            args.append('--dup_to_ins')

        if config.remove_small:
            args.append('--remove_small')

        if config.force_reverse:
            args.append('--force_reverse')

        # 性能参数|Performance parameters
        args.extend(['--thread_num', str(config.threads)])

        # 外部工具路径|External tool paths
        if config.minigraph_path != 'minigraph':
            args.extend(['--minigraph', config.minigraph_path])

        if config.gfatools_path != 'gfatools':
            args.extend(['--gfatools', config.gfatools_path])

        return args

    @staticmethod
    def _get_swave_env() -> tuple:
        """获取swave conda环境的Python路径和环境变量|Get swave conda env Python path and env vars

        Returns:
            (python路径, 环境变量dict)|(python path, env vars dict)
        """
        conda_base = os.environ.get('CONDA_EXE', '')
        if conda_base:
            env_bin = os.path.join(
                os.path.dirname(os.path.dirname(conda_base)), 'envs', 'swave_v.1.2', 'bin')
            python_path = os.path.join(env_bin, 'python')
            if os.path.exists(python_path):
                env = os.environ.copy()
                # 将conda env的bin目录添加到PATH首位，确保Swave内部调用minigraph/gfatools可找到
                env['PATH'] = env_bin + ':' + env.get('PATH', '')
                return python_path, env
        return 'python', None

    def _build_swave_command(self) -> List[str]:
        """
        构建Swave命令|Build Swave command

        Returns:
            命令列表|Command list
        """
        swave_script = os.path.join(self.config.swave_path, 'Swave.py')
        swave_args = self.build_call_args(self.config)
        python_path, _ = self._get_swave_env()
        return [python_path, swave_script] + swave_args

    def run_call(self) -> bool:
        """
        运行SV检测|Run SV detection

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始Swave结构变异检测|Starting Swave structural variant detection")

        # 构建命令|Build command
        cmd = self._build_swave_command()
        _, env = self._get_swave_env()

        # 执行命令|Execute command
        success, stdout, stderr = self.cmd_runner.run(
            cmd,
            description="Swave SV检测|Swave SV detection",
            cwd=self.config.swave_path,
            env=env
        )

        if not success:
            self.logger.error(f"Swave执行失败|Swave execution failed")
            if stderr:
                self.logger.error(f"错误信息|Error message: {stderr}")
            return False

        # 检查输出文件|Check output files
        output_files = self._check_output_files()

        if output_files:
            self.logger.info("Swave执行成功|Swave execution completed successfully")
            self.logger.info(f"输出文件|Output files:")
            for file_type, file_path in output_files.items():
                self.logger.info(f"  {file_type}: {file_path}")
            return True
        else:
            self.logger.warning("Swave执行完成但未找到预期输出文件|"
                              "Swave execution completed but expected output files not found")
            return False

    def _check_output_files(self) -> dict:
        """
        检查输出文件|Check output files

        Returns:
            输出文件字典|Output files dictionary
        """
        output_files = {}

        # 检查主要输出文件|Check main output files
        expected_files = {
            'log': 'swave.log',
            'vcf_multi': 'swave.sample_level.vcf',
            'vcf_bi': 'swave.sample_level.split.vcf'
        }

        for file_type, file_name in expected_files.items():
            file_path = os.path.join(self.config.output_dir, file_name)
            if os.path.exists(file_path):
                output_files[file_type] = file_path
                self.logger.debug(f"找到输出文件|Found output file: {file_path}")
            else:
                self.logger.debug(f"未找到输出文件|Output file not found: {file_path}")

        return output_files

    def get_output_summary(self) -> dict:
        """
        获取输出摘要|Get output summary

        Returns:
            摘要信息字典|Summary information dictionary
        """
        summary = {
            'output_dir': self.config.output_dir,
            'gfa_source': self.config.gfa_source,
            'output_files': {}
        }

        # 检查VCF文件并统计|Check VCF files and count
        vcf_files = {
            'vcf_multi': 'swave.sample_level.vcf',
            'vcf_bi': 'swave.sample_level.split.vcf'
        }

        for file_type, file_name in vcf_files.items():
            file_path = os.path.join(self.config.output_dir, file_name)
            if os.path.exists(file_path):
                # 统计VCF记录数|Count VCF records
                try:
                    with open(file_path, 'r') as f:
                        record_count = sum(1 for line in f if not line.startswith('#'))
                    summary['output_files'][file_type] = {
                        'path': file_path,
                        'record_count': record_count
                    }
                except Exception as e:
                    self.logger.warning(f"无法统计VCF记录数|Cannot count VCF records: {e}")
                    summary['output_files'][file_type] = {
                        'path': file_path,
                        'record_count': 'unknown'
                    }

        return summary

    def _auto_rename_ref_fasta(self, vcf_path: str) -> str:
        """自动检测VCF contig与FASTA序列名的映射关系，必要时创建重命名的FASTA

        swave的convert_seq使用pysam按contig名从FASTA取序列，
        当VCF contig名（如T2T#0#Chr01）与FASTA序列名（如Chr01）不一致时会报错。
        此方法自动建立映射并创建重命名的FASTA。

        Args:
            vcf_path: VCF文件路径（用于读取contig名）

        Returns:
            可直接用于convert_seq的FASTA路径（原FASTA或重命名后的FASTA）
        """
        # 读取VCF contig名
        vcf_contigs = []
        try:
            import pysam
            with pysam.VariantFile(vcf_path) as vcf:
                vcf_contigs = list(vcf.header.contigs.keys())
        except Exception as e:
            self.logger.warning(f"无法读取VCF contig|Cannot read VCF contigs: {e}")
            return self.config.ref_fasta

        if not vcf_contigs:
            return self.config.ref_fasta

        # 读取FASTA序列名
        import pysam
        try:
            ref_file = pysam.FastaFile(self.config.ref_fasta)
            ref_names = list(ref_file.references)
            ref_file.close()
        except Exception as e:
            self.logger.warning(f"无法读取参考FASTA|Cannot read ref FASTA: {e}")
            return self.config.ref_fasta

        # 检查是否需要重命名：任一ref_name不在vcf_contigs中就需要
        if all(name in vcf_contigs for name in ref_names):
            self.logger.info("FASTA序列名与VCF contig名一致，无需重命名|"
                             "FASTA names match VCF contigs, no rename needed")
            return self.config.ref_fasta

        # 建立映射：对每个ref_name，在vcf_contig的'#'分割片段中查找匹配
        name_mapping: Dict[str, str] = {}
        for ref_name in ref_names:
            matched_contigs = []
            for contig in vcf_contigs:
                segments = contig.split('#')
                if ref_name in segments:
                    matched_contigs.append(contig)
            if len(matched_contigs) == 1:
                name_mapping[ref_name] = matched_contigs[0]
            elif len(matched_contigs) == 0:
                self.logger.warning(f"FASTA序列'{ref_name}'未在VCF contig中找到匹配|"
                                    f"FASTA seq '{ref_name}' has no match in VCF contigs")
            else:
                self.logger.warning(f"FASTA序列'{ref_name}'匹配到多个VCF contig: {matched_contigs}|"
                                    f"FASTA seq '{ref_name}' matched multiple VCF contigs")

        if not name_mapping:
            self.logger.warning("未建立任何名称映射，使用原始FASTA|No name mapping built, using original FASTA")
            return self.config.ref_fasta

        # 创建重命名的FASTA
        renamed_fasta = os.path.join(self.config.output_dir, 'ref.renamed.fa')
        self.logger.info(f"创建重命名FASTA: {renamed_fasta}|Creating renamed FASTA")

        try:
            ref_file = pysam.FastaFile(self.config.ref_fasta)
            with open(renamed_fasta, 'w') as fout:
                for ref_name in ref_names:
                    new_name = name_mapping.get(ref_name, ref_name)
                    seq = ref_file.fetch(ref_name)
                    fout.write(f">{new_name}\n{seq}\n")
            ref_file.close()
            self.logger.info(f"名称映射示例|Name mapping example: "
                             f"{list(name_mapping.items())[0][0]} -> {list(name_mapping.items())[0][1]}")
            return renamed_fasta
        except Exception as e:
            self.logger.error(f"创建重命名FASTA失败|Failed to create renamed FASTA: {e}")
            return self.config.ref_fasta

    def run_convert_seq(self) -> bool:
        """在call成功后自动将VCF中的图路径编号转换为实际碱基序列

        对swave.sample_level.vcf和swave.sample_level.split.vcf分别调用convert_seq，
        自动处理contig名不匹配问题，输出*.with_seq.vcf文件。

        Returns:
            是否全部成功|Whether all conversions succeeded
        """
        import pysam

        self.logger.info("开始自动序列转换|Starting automatic sequence conversion")

        # 待转换的VCF文件列表
        vcf_files = [
            'swave.sample_level.vcf',
            'swave.sample_level.split.vcf',
        ]

        # 检查哪些VCF文件存在
        existing_vcfs = []
        for vcf_name in vcf_files:
            vcf_path = os.path.join(self.config.output_dir, vcf_name)
            if os.path.exists(vcf_path):
                existing_vcfs.append((vcf_name, vcf_path))

        if not existing_vcfs:
            self.logger.warning("未找到可转换的VCF文件|No VCF files found for conversion")
            return False

        # 用第一个VCF文件来检测是否需要重命名FASTA
        _, first_vcf = existing_vcfs[0]
        actual_ref = self._auto_rename_ref_fasta(first_vcf)

        # 创建临时输出目录（swave convert_seq的output_path是目录）
        tmp_output = os.path.join(self.config.output_dir, 'tmp_convert_seq')
        os.makedirs(tmp_output, exist_ok=True)

        swave_script = os.path.join(self.config.swave_path, 'Swave.py')
        python_path, env = self._get_swave_env()

        all_success = True
        need_cleanup_renamed = (actual_ref != self.config.ref_fasta)

        for vcf_name, vcf_path in existing_vcfs:
            self.logger.info(f"转换序列: {vcf_name}|Converting sequences: {vcf_name}")

            args = [python_path, swave_script, 'convert_seq',
                    '--vcf_path', vcf_path,
                    '--gfa_path', self.config.gfa_file,
                    '--ref_path', actual_ref,
                    '--output_path', tmp_output]

            try:
                result = subprocess.run(
                    args, shell=False, check=False,
                    cwd=self.config.swave_path, env=env
                )

                # swave输出文件名: {原名}.converted.vcf
                converted_name = vcf_name.replace('.vcf', '.converted.vcf')
                converted_path = os.path.join(tmp_output, converted_name)
                # 目标文件名: {原名}.with_seq.vcf
                target_name = vcf_name.replace('.vcf', '.with_seq.vcf')
                target_path = os.path.join(self.config.output_dir, target_name)

                if os.path.exists(converted_path):
                    # 重命名为.with_seq.vcf
                    os.rename(converted_path, target_path)
                    self.logger.info(f"序列转换完成: {target_name}|Sequence conversion done: {target_name}")
                else:
                    self.logger.warning(f"转换输出文件未生成: {converted_name}|"
                                        f"Conversion output not found: {converted_name}")
                    all_success = False
            except Exception as e:
                self.logger.error(f"序列转换失败 {vcf_name}: {e}|Sequence conversion failed: {e}")
                all_success = False

        # 清理临时文件
        try:
            import shutil
            if os.path.exists(tmp_output):
                shutil.rmtree(tmp_output)
            if need_cleanup_renamed and os.path.exists(actual_ref):
                os.remove(actual_ref)
        except Exception as e:
            self.logger.warning(f"清理临时文件失败|Failed to cleanup temp files: {e}")

        if all_success:
            self.logger.info("所有VCF序列转换完成|All VCF sequence conversions completed")
        else:
            self.logger.warning("部分VCF序列转换失败|Some VCF sequence conversions failed")

        return all_success


class SwaveConverter:
    """Swave格式转换器|Swave Format Converter"""

    def __init__(self, config: SwaveConfig, logger, cmd_runner: CommandRunner):
        """
        初始化转换器|Initialize converter

        Args:
            config: Swave配置对象|Swave config object
            logger: 日志器|Logger
            cmd_runner: 命令执行器|Command runner
        """
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def convert_seq(self, vcf_path: str, output_path: Optional[str] = None,
                    force_pangenie: bool = False) -> bool:
        """
        转换图路径为VCF的REF和ALT序列|Convert graph paths to REF and ALT sequences in VCF

        Args:
            vcf_path: VCF文件路径|VCF file path
            output_path: 输出路径|Output path
            force_pangenie: 强制输出满足pangenie要求的序列|Force output sequences to meet pangenie requirements

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("开始序列转换|Starting sequence conversion")

        swave_script = os.path.join(self.config.swave_path, 'Swave.py')

        python_path, env = self._get_swave_env()
        args = [python_path, swave_script, 'convert_seq',
                '--vcf_path', vcf_path,
                '--gfa_path', self.config.gfa_file,
                '--ref_path', self.config.ref_fasta]

        if output_path:
            args.extend(['--output_path', output_path])

        if force_pangenie:
            args.append('--force_pangenie')

        success, stdout, stderr = self.cmd_runner.run(
            args,
            description="序列转换|Sequence conversion",
            cwd=self.config.swave_path,
            env=env
        )

        return success

    def extract_sample(self, vcf_path: str, spec_samples: List[str],
                       output_path: Optional[str] = None) -> bool:
        """
        从VCF提取特定样本的SV|Extract SVs for specific samples from VCF

        Args:
            vcf_path: VCF文件路径|VCF file path
            spec_samples: 指定样本列表|Specific sample list
            output_path: 输出路径|Output path

        Returns:
            是否成功|Whether successful
        """
        self.logger.info(f"提取样本SV|Extracting SVs for samples: {spec_samples}")

        swave_script = os.path.join(self.config.swave_path, 'Swave.py')

        python_path, env = self._get_swave_env()
        args = [python_path, swave_script, 'extract_sample',
                '--vcf_path', vcf_path,
                '--spec_samples'] + spec_samples

        if output_path:
            args.extend(['--output_path', output_path])

        success, stdout, stderr = self.cmd_runner.run(
            args,
            description="提取样本SV|Extract sample SVs",
            cwd=self.config.swave_path,
            env=env
        )

        return success
