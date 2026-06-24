"""
rMVP核心分析类|rMVP Core Analyzer Class
"""

import os
import subprocess
import time
from pathlib import Path
from typing import List, Dict, Optional

from .config import RMVPConfig
from .utils import RMVPLogger, detect_r_environment, check_r_mvp, validate_vcf_file, validate_phenotype_file, prepare_vcf_for_rmvp, CommandRunner, build_conda_command
from .rscript_generator import RMVPRScriptGenerator


class RMVPAnalyzer:
    """rMVP GWAS分析器|rMVP GWAS Analyzer"""

    def __init__(self, config: RMVPConfig):
        """
        初始化分析器|Initialize analyzer

        Args:
            config: RMVPConfig配置对象|RMVPConfig object
        """
        self.config = config

        # 初始化日志|Initialize logger
        log_file = self.config.output_dir / f"{self.config.output_prefix}.log"
        logger_manager = RMVPLogger(log_file, self.config.log_level)
        self.logger = logger_manager.get_logger()

        # R脚本生成器|R script generator
        self.script_generator = RMVPRScriptGenerator(config)

        # 临时文件列表|Temporary files list
        self.temp_files: List[Path] = []

    def check_dependencies(self) -> bool:
        """
        检查依赖软件|Check dependencies

        Returns:
            是否所有依赖都满足|Whether all dependencies are satisfied
        """
        self.logger.info(" 检查依赖软件|Checking dependencies")

        # 检查R环境|Check R environment
        self.logger.info("  检测R环境|Detecting R environment")
        r_success, r_env_name, r_type = detect_r_environment(
            self.config.r_env,
            self.config.r_path,
            self.logger
        )

        if not r_success:
            self.logger.error("  R环境检测失败|R environment detection failed")
            return False

        self.logger.info(f"  R环境检测成功|R environment detected: {r_env_name} (type: {r_type})")

        # 保存R环境名称和类型|Save R env name and type
        self.config.r_env_name = r_env_name
        self.config.r_env_type = r_type

        # 检查rMVP包|Check rMVP package
        self.logger.info("  检查rMVP包|Checking rMVP package")
        r_mvp_installed = check_r_mvp(r_env_name, r_type, self.logger)

        if not r_mvp_installed:
            self.logger.error("  rMVP包未安装，请安装后重试|rMVP package not installed, please install and retry")
            self.logger.error("  安装命令|Install command: conda run -n {} R -e \"install.packages('rMVP', repos='https://cloud.r-project.org')\"".format(r_env_name))
            return False

        self.logger.info("  rMVP包已安装|rMVP package is installed")

        # 检查PLINK（LD去连锁开启时）|Check PLINK (when LD pruning enabled)
        if self.config.ld_pruning:
            self.logger.info("  检查PLINK|Checking PLINK")
            if not Path(self.config.plink_path).exists():
                self.logger.error(f"  PLINK不存在|PLINK not found: {self.config.plink_path}")
                self.logger.error("  请设置 --plink-path 或环境变量 PLINK_PATH|Set --plink-path or env PLINK_PATH")
                return False
            self.logger.info(f"  PLINK可用|PLINK available: {self.config.plink_path}")

        return True

    def validate_input_files(self) -> bool:
        """
        验证输入文件|Validate input files

        Returns:
            是否所有文件都有效|Whether all files are valid
        """
        self.logger.info(" 验证输入文件|Validating input files")

        # 验证VCF文件|Validate VCF file
        if not validate_vcf_file(Path(self.config.vcf_file), self.logger):
            return False

        # 验证表型文件|Validate phenotype file
        valid, n_cols = validate_phenotype_file(Path(self.config.pheno_file), self.logger)

        if not valid:
            return False

        # 保存表型数量|Save number of traits
        self.n_traits = n_cols - 1

        return True

    def _is_step_completed(self, output_file: Path) -> bool:
        """
        检查步骤是否已完成（通过输出文件存在性判断）|Check if step is done

        Args:
            output_file: 该步骤的关键输出文件|Key output file for this step

        Returns:
            是否已完成|Whether completed
        """
        return Path(output_file).exists()

    def run_data_conversion(self) -> bool:
        """
        运行数据转换（VCF → rMVP格式）|Run data conversion (VCF → rMVP format)

        Returns:
            是否成功|Whether successful
        """
        # 检查关键输出文件|Check key output file
        output_dir_abs = Path(self.config.output_dir).resolve()
        key_output_file = output_dir_abs / f"{self.config.output_prefix}.geno.desc"

        if self._is_step_completed(key_output_file):
            self.logger.info(" 开始数据转换|Starting data conversion")
            self.logger.info(f"   跳过已完成步骤|Skipping completed step: 数据转换|Data conversion")
            self.logger.info(f"   输出文件已存在|Output file exists: {key_output_file}")
            return True

        self.logger.info(" 开始数据转换|Starting data conversion")

        try:
            # 预处理VCF文件（处理bgzip格式不兼容）|Pre-process VCF file (handle bgzip incompatibility)
            original_vcf = self.config.vcf_file
            prepared_vcf = prepare_vcf_for_rmvp(
                self.config.vcf_file, str(self.config.output_dir), self.logger
            )
            if prepared_vcf != original_vcf:
                self.config.vcf_file = prepared_vcf

            # 生成R脚本|Generate R script
            script_content = self.script_generator.generate_data_conversion_script()

            # 确保输出目录存在|Ensure output directory exists
            output_dir_abs.mkdir(parents=True, exist_ok=True)

            # 使用绝对路径保存脚本|Save script using absolute path
            script_file = output_dir_abs / f"{self.config.output_prefix}_convert.R"

            with open(script_file, 'w', encoding='utf-8') as f:
                f.write(script_content)

            self.logger.info(f"   R脚本已生成|R script generated: {script_file}")

            # 执行R脚本|Execute R script
            return self._execute_r_script(script_file, "数据转换|Data conversion")

        except Exception as e:
            self.logger.error(f" 数据转换失败|Data conversion failed: {e}")
            return False

    def run_ld_pruning(self) -> bool:
        """
        LD去连锁：用PLINK修剪连锁SNP，产出去连锁VCF供rMVP计算K/PCA
        LD pruning: prune linked SNPs with PLINK, produce pruned VCF for rMVP K/PCA

        kinship/PCA在去连锁SNP上计算，GWAS仍用全部SNP
        Kinship/PCA computed on pruned SNPs, GWAS still uses all SNPs

        Returns:
            是否成功|Whether successful
        """
        output_dir_abs = Path(self.config.output_dir).resolve()
        prefix = self.config.output_prefix
        plink_prefix = str(output_dir_abs / f"{prefix}_plink")
        ld_prefix = str(output_dir_abs / f"{prefix}_ld")
        pruned_prefix = str(output_dir_abs / f"{prefix}_pruned")
        pruned_vcf = output_dir_abs / f"{prefix}_pruned.vcf"
        self.pruned_vcf = str(pruned_vcf)

        # 断点续传：去连锁VCF已存在则跳过|Checkpoint: skip if pruned VCF exists
        if pruned_vcf.exists():
            self.logger.info(f"   跳过已完成步骤|Skipping completed step: LD去连锁|LD pruning")
            self.logger.info(f"   去连锁VCF已存在|Pruned VCF exists: {pruned_vcf}")
            return True

        self.logger.info(f"   LD去连锁参数|LD pruning params: --indep-pairwise {self.config.ld_window} {self.config.ld_step} {self.config.ld_r2}")
        self.logger.info(f"   线程数|Threads: {self.config.ncpus}")

        # 公共PLINK参数|Common PLINK args
        common = [
            '--keep-allele-order', '--real-ref-alleles',
            '--allow-extra-chr', '--double-id',
            '--no-sex', '--no-parents', '--no-fid',
            '--threads', str(self.config.ncpus),
        ]
        cmd_runner = CommandRunner(self.logger, output_dir_abs)

        # 第1步：VCF → PLINK bed|Step 1: VCF to PLINK bed
        if not cmd_runner.run(
            build_conda_command(self.config.plink_path,
                ['--vcf', self.config.vcf_file, '--make-bed', '--out', plink_prefix] + common),
            "LD去连锁步骤1/4：VCF转PLINK|LD pruning step 1/4: VCF to PLINK"):
            return False

        # 第2步：LD修剪|Step 2: LD pruning (indep-pairwise)
        if not cmd_runner.run(
            build_conda_command(self.config.plink_path,
                ['--bfile', plink_prefix,
                 '--indep-pairwise', self.config.ld_window, str(self.config.ld_step), str(self.config.ld_r2),
                 '--out', ld_prefix] + common),
            "LD去连锁步骤2/4：LD修剪|LD pruning step 2/4: indep-pairwise"):
            return False

        # 第3步：提取非连锁SNP|Step 3: extract pruned SNPs
        if not cmd_runner.run(
            build_conda_command(self.config.plink_path,
                ['--bfile', plink_prefix, '--extract', f"{ld_prefix}.prune.in",
                 '--make-bed', '--out', pruned_prefix] + common),
            "LD去连锁步骤3/4：提取非连锁SNP|LD pruning step 3/4: extract pruned SNPs"):
            return False

        # 第4步：去连锁bed → VCF（rMVP输入）|Step 4: pruned bed to VCF (rMVP input)
        if not cmd_runner.run(
            build_conda_command(self.config.plink_path,
                ['--bfile', pruned_prefix, '--recode', 'vcf-iid', '--out', pruned_prefix] + common),
            "LD去连锁步骤4/4：去连锁bed转VCF|LD pruning step 4/4: pruned bed to VCF"):
            return False

        if not pruned_vcf.exists():
            self.logger.error(f"   去连锁VCF未生成|Pruned VCF not generated: {pruned_vcf}")
            return False

        self.logger.info(f"   去连锁VCF已生成|Pruned VCF generated: {pruned_vcf}")
        return True

    def run_batch_analysis(self, trait_names: List[str]) -> bool:
        """
        运行批量GWAS分析|Run batch GWAS analysis

        Args:
            trait_names: 表型名称列表|List of trait names

        Returns:
            是否成功|Whether successful
        """
        self.logger.info(" 开始批量GWAS分析|Starting batch GWAS analysis")
        self.logger.info(f"   表型数量|Number of traits: {len(trait_names)}")
        self.logger.info(f"   分析模型|Models: {', '.join(self.config.models)}")

        # 检查所有输出文件是否都已存在|Check if all output files exist
        output_dir_abs = Path(self.config.output_dir).resolve()
        all_completed = True
        missing_files = []

        for trait_name in trait_names:
            for model in self.config.models:
                # 输出文件格式：{prefix}_{trait}.{model}.pmap
                # Output file format: {prefix}_{trait}.{model}.pmap
                output_file = output_dir_abs / f"{self.config.output_prefix}_{trait_name}.{model.lower()}.pmap"
                if not self._is_step_completed(output_file):
                    all_completed = False
                    missing_files.append(str(output_file.name))

        if all_completed:
            self.logger.info(f"   跳过已完成步骤|Skipping completed step: 批量GWAS分析|Batch GWAS analysis")
            self.logger.info(f"   所有输出文件已存在|All output files exist")
            return True

        if missing_files:
            self.logger.info(f"   需要生成|Need to generate: {len(missing_files)} 个文件|files")

        try:
            # 生成批量分析R脚本|Generate batch analysis R script
            script_content = self.script_generator.generate_batch_script(trait_names)

            # 确保输出目录存在|Ensure output directory exists
            output_dir_abs.mkdir(parents=True, exist_ok=True)

            # 使用绝对路径保存脚本|Save script using absolute path
            script_file = output_dir_abs / f"{self.config.output_prefix}_batch.R"

            with open(script_file, 'w', encoding='utf-8') as f:
                f.write(script_content)

            self.logger.info(f"   R脚本已生成|R script generated: {script_file}")

            # 执行R脚本|Execute R script
            success = self._execute_r_script(script_file, "批量GWAS分析|Batch GWAS analysis")

            if success:
                self.logger.info(" 批量分析完成|Batch analysis completed")
            else:
                self.logger.error(" 批量分析失败|Batch analysis failed")

            return success

        except Exception as e:
            self.logger.error(f" 批量分析失败|Batch analysis failed: {e}")
            return False

    def run_analysis(self) -> bool:
        """
        运行完整的GWAS分析流程|Run complete GWAS analysis pipeline

        Returns:
            是否成功|Whether successful
        """
        start_time = time.time()

        self.logger.info("=" * 70)
        self.logger.info(" rMVP GWAS分析流程|rMVP GWAS Analysis Pipeline")
        self.logger.info("=" * 70)

        try:
            # 步骤总数（LD去连锁开启时为6步，否则5步）|Total steps (6 with LD pruning, else 5)
            total = 6 if self.config.ld_pruning else 5
            cur = 0

            # 1. 检查依赖|Check dependencies
            cur += 1
            self.logger.info(f"[{cur}/{total}] 检查依赖|Checking dependencies")
            if not self.check_dependencies():
                return False

            # 2. 验证输入文件|Validate input files
            cur += 1
            self.logger.info(f"[{cur}/{total}] 验证输入文件|Validating input files")
            if not self.validate_input_files():
                return False

            # 读取表型名称|Read trait names
            trait_names = self._get_trait_names()
            self.logger.info(f"   检测到表型|Detected traits: {trait_names}")

            # 3. LD去连锁（开启时）|LD pruning (when enabled)
            if self.config.ld_pruning:
                cur += 1
                self.logger.info(f"[{cur}/{total}] LD去连锁（PLINK）|LD pruning (PLINK)")
                if not self.run_ld_pruning():
                    return False

            # 4. 数据转换|Data conversion
            cur += 1
            self.logger.info(f"[{cur}/{total}] 数据转换|Data conversion")
            if not self.run_data_conversion():
                return False

            # 5. GWAS分析|GWAS analysis
            cur += 1
            self.logger.info(f"[{cur}/{total}] GWAS分析|GWAS analysis")
            if not self.run_batch_analysis(trait_names):
                return False

            # 6. 整合结果|Integrate results
            cur += 1
            self.logger.info(f"[{cur}/{total}] 整合结果|Integrating results")
            # 这个功能在result_parser.py中实现|This will be implemented in result_parser.py

            # 完成|Completed
            elapsed_time = time.time() - start_time
            self.logger.info("=" * 70)
            self.logger.info(f" 分析完成！|Analysis completed!")
            self.logger.info(f" 总运行时间|Total runtime: {elapsed_time:.2f}秒|seconds")
            self.logger.info("=" * 70)

            return True

        except KeyboardInterrupt:
            self.logger.info("\n用户中断操作|User interrupted operation")
            return False
        except Exception as e:
            self.logger.error(f" 分析过程中发生错误|Error during analysis: {e}")
            return False

    def _get_trait_names(self) -> List[str]:
        """
        获取表型名称|Get trait names

        Returns:
            表型名称列表|List of trait names
        """
        try:
            with open(self.config.pheno_file, 'r') as f:
                header = f.readline().strip()
                columns = header.split('\t')
                # 第1列是样本ID，第2列开始是表型|Column 1 is sample ID, columns 2+ are traits
                return columns[1:]
        except Exception as e:
            self.logger.error(f" 读取表型名称失败|Failed to read trait names: {e}")
            return []

    def _execute_r_script(self, script_file: Path, description: str) -> bool:
        """
        执行R脚本|Execute R script

        Args:
            script_file: R脚本文件路径|R script file path
            description: 任务描述|Task description

        Returns:
            是否成功|Whether successful
        """
        try:
            # 构建命令|Build command
            if self.config.r_env_type == "conda":
                # 使用conda环境|Use conda environment
                # 注意：使用 Rscript 而不是 R，因为 Rscript 用于运行脚本文件|Note: use Rscript not R, because Rscript is for running script files
                cmd = f"conda run -n {self.config.r_env_name} Rscript {script_file}"
            else:
                # 直接使用Rscript|Use Rscript directly
                cmd = f"{self.config.r_env_name} Rscript {script_file}"

            self.logger.info(f"   命令|Command: {cmd}")

            # 使用绝对路径作为工作目录|Use absolute path as working directory
            cwd = str(Path(self.config.output_dir).resolve())
            self.logger.info(f"   工作目录|Working directory: {cwd}")

            # 输出文件|Output files
            stdout_file = Path(cwd) / f"{Path(script_file).stem}.out"
            # stderr合并到stdout，不再单独创建.err文件|stderr merged to stdout, no separate .err file

            # 执行命令，重定向输出到文件|Execute command, redirect output to files
            start_time = time.time()

            with open(stdout_file, 'w') as stdout_f:
                result = subprocess.run(
                    cmd,
                    shell=True,
                    cwd=cwd,
                    stdout=stdout_f,
                    stderr=subprocess.STDOUT,  # 合并stderr到stdout|Merge stderr to stdout
                    text=True,
                    timeout=86400  # 24小时超时|24 hours timeout
                )

            elapsed_time = time.time() - start_time

            # 记录返回码|Log return code
            self.logger.info(f"   命令返回码|Command return code: {result.returncode}")

            # 检查执行结果|Check result
            if result.returncode == 0:
                self.logger.info(f"   {description}成功|{description} successful")
                self.logger.info(f"   执行时间|Execution time: {elapsed_time:.2f}秒|seconds")
                return True
            else:
                self.logger.error(f"   {description}失败|{description} failed")
                return False

        except subprocess.TimeoutExpired:
            self.logger.error(f"   {description}超时|{description} timeout")
            return False
        except Exception as e:
            self.logger.error(f"   {description}异常|{description} error: {e}")
            return False
