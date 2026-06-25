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

    def _get_batch_completion_files(self, trait_names: List[str]) -> List[Path]:
        """
        返回批量GWAS的完成标志文件列表|Return batch GWAS completion marker files

        rMVP 1.4.6 实际产出 {trait}.{model}.{trait}.csv（如 phe.GLM.phe.csv），
        这是用户真正需要的结果文件。以此作为续传完成标志，而非 .RData——.RData 可能在
        结果文件缺失时仍存在（如异常中断），会导致误跳过GWAS。
        |rMVP 1.4.6 emits {trait}.{model}.{trait}.csv; use it as the resume marker
        rather than .RData, which can persist when result files are absent.
        """
        output_dir_abs = Path(self.config.output_dir).resolve()
        files = []
        for t in trait_names:
            for m in self.config.models:
                files.append(output_dir_abs / f"{t}.{m}.{t}.csv")
        return files

    def _is_conversion_completed(self) -> bool:
        """
        数据转换是否完成|Is data conversion complete?

        ld_pruning=True 时，转换步产出两部分：完整VCF的基因型(.geno.*) 和去连锁VCF的
        K/PCA(_pruned.kin.*/_pruned.pc.*)。仅检查 .geno.desc 会漏掉后者，导致后续 batch
        读取 _pruned.kin.desc 失败。|When ld_pruning=True, conversion produces both the
        full-VCF genotype (.geno.*) and pruned-VCF K/PCA (_pruned.kin.*/_pruned.pc.*).
        Checking only .geno.desc misses the latter, causing batch to fail reading
        _pruned.kin.desc.
        """
        output_dir_abs = Path(self.config.output_dir).resolve()
        prefix = self.config.output_prefix
        required = [output_dir_abs / f"{prefix}.geno.desc"]
        if self.config.ld_pruning:
            required.append(output_dir_abs / f"{prefix}_pruned.kin.desc")
            required.append(output_dir_abs / f"{prefix}_pruned.pc.desc")
        return all(self._is_step_completed(f) for f in required)

    def run_data_conversion(self) -> bool:
        """
        运行数据转换（VCF → rMVP格式）|Run data conversion (VCF → rMVP format)

        Returns:
            是否成功|Whether successful
        """
        # 检查关键输出文件|Check key output files
        output_dir_abs = Path(self.config.output_dir).resolve()

        if self._is_conversion_completed():
            self.logger.info(" 开始数据转换|Starting data conversion")
            self.logger.info(f"   跳过已完成步骤|Skipping completed step: 数据转换|Data conversion")
            self.logger.info(f"   输出文件已存在|Output files exist: geno + pruned K/PCA")
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

        output_dir_abs = Path(self.config.output_dir).resolve()

        # 检查所有表型是否都已完成|Check if all traits are completed
        # 完成标志：{trait}.{model}.{trait}.csv（rMVP真实结果文件，见 _get_batch_completion_files）
        # Marker: {trait}.{model}.{trait}.csv (real rMVP result; see _get_batch_completion_files)
        completion_files = self._get_batch_completion_files(trait_names)
        missing_files = [f for f in completion_files if not self._is_step_completed(f)]

        if not missing_files:
            self.logger.info(f"   跳过已完成步骤|Skipping completed step: 批量GWAS分析|Batch GWAS analysis")
            self.logger.info(f"   所有表型结果已存在|All trait results exist")
            return True

        self.logger.info(f"   需要生成|Need to generate: {len(missing_files)} 个表型结果|trait results")

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

    def _build_r_command(self, script_file: Path) -> List[str]:
        """
        构建R脚本执行命令（列表形式，配合shell=False）|Build R script command (list form, for shell=False)

        conda环境：conda run -n <env> --no-capture-output Rscript <script>
          --no-capture-output 避免 conda 缓冲输出（§13.2.0），让 R 的 flush.console()
          输出能实时透传到本进程|R batch run via conda with --no-capture-output (§13.2.0) so
          flush.console() output streams through in real time.
        direct：直接调用 Rscript|direct: invoke Rscript directly

        Args:
            script_file: R脚本文件路径|R script file path

        Returns:
            命令列表（配合 subprocess shell=False）|Command list (for subprocess shell=False)
        """
        script_file = str(script_file)
        if self.config.r_env_type == "conda":
            return ['conda', 'run', '-n', self.config.r_env_name,
                    '--no-capture-output', 'Rscript', script_file]
        else:
            # direct 模式：r_env_name 为 R 可执行文件路径|direct: r_env_name is the R executable path
            return [self.config.r_env_name, 'Rscript', script_file]

    def _execute_r_script(self, script_file: Path, description: str) -> bool:
        """
        执行R脚本，流式转发输出到日志|Execute R script, stream output to logger

        R 的 stdout 逐行读取并转发到 logger（实时可见于调度系统 .out 和模块 .log），
        同时 tee 到 {stem}.out 保留原始 R 输出。配合生成脚本中的 flush.console()，
        解决 R 块缓冲导致的进度黑洞。|R stdout is read line-by-line and forwarded to the
        logger (visible in real time in the scheduler .out and module .log), and teed to
        {stem}.out to keep raw R output. Combined with flush.console() in the generated
        script, this fixes the progress black hole caused by R block buffering.

        Args:
            script_file: R脚本文件路径|R script file path
            description: 任务描述|Task description

        Returns:
            是否成功|Whether successful
        """
        try:
            cmd = self._build_r_command(script_file)
            self.logger.info(f"   命令|Command: {' '.join(cmd)}")

            # 使用绝对路径作为工作目录|Use absolute path as working directory
            cwd = str(Path(self.config.output_dir).resolve())
            self.logger.info(f"   工作目录|Working directory: {cwd}")

            # 原始R输出文件（保留兼容）|Raw R output file (kept for compatibility)
            stdout_file = Path(cwd) / f"{Path(script_file).stem}.out"

            start_time = time.time()

            # 列表命令 + shell=False，逐行流式读取|List command + shell=False, line-streamed
            proc = subprocess.Popen(
                cmd,
                shell=False,
                cwd=cwd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,  # 合并stderr到stdout|Merge stderr to stdout
                text=True,
                bufsize=1,  # 行缓冲|Line buffered
            )

            # 兜底超时（仅在两次输出之间检查；主要靠调度器墙钟限制）|Backstop timeout
            # (checked only between output lines; the scheduler enforces wall-time primarily)
            timeout_seconds = 86400  # 24小时|24 hours
            deadline = time.time() + timeout_seconds
            timed_out = False

            try:
                with open(stdout_file, 'w', encoding='utf-8') as out_f:
                    while True:
                        line = proc.stdout.readline()
                        if not line:
                            break
                        # tee：原始输出写文件|Tee: raw output to file
                        out_f.write(line)
                        out_f.flush()
                        # 转发到日志：过滤 \r 进度条减少噪声|Forward to log, filtering \r progress bars
                        clean = line.rstrip('\n')
                        if clean and '\r' not in clean:
                            self.logger.info(f"   {clean}")
                        if time.time() > deadline and proc.poll() is None:
                            timed_out = True
                            break
            finally:
                if timed_out:
                    proc.kill()
                proc.stdout.close()
                proc.wait()
                returncode = proc.returncode

            elapsed_time = time.time() - start_time

            if timed_out:
                self.logger.error(f"   {description}超时|{description} timeout (>{timeout_seconds}s)")
                return False

            # 记录返回码|Log return code
            self.logger.info(f"   命令返回码|Command return code: {returncode}")

            # 检查执行结果|Check result
            if returncode == 0:
                self.logger.info(f"   {description}成功|{description} successful")
                self.logger.info(f"   执行时间|Execution time: {elapsed_time:.2f}秒|seconds")
                return True
            else:
                self.logger.error(f"   {description}失败|{description} failed (return code {returncode})")
                return False

        except Exception as e:
            self.logger.error(f"   {description}异常|{description} error: {e}")
            return False
