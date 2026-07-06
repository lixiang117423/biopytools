"""
LAI指数评估模块|LAI Index Evaluation Module

使用EDTA流程计算LAI指数
Calculate LAI index using EDTA pipeline
"""

import os
import subprocess
import re
import glob
from pathlib import Path
from typing import Dict, Any, Optional, List
from .utils import get_conda_env_from_path, build_conda_command


class LAIEvaluator:
    """LAI指数评估器|LAI Index Evaluator"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
        self.working_dir = Path(self.config.lai_output_dir)
        self.working_dir.mkdir(parents=True, exist_ok=True)

        # 从conda环境路径提取环境名|Extract env names from conda env paths
        self.edta_env_name = get_conda_env_from_path(self.config.conda_env_edta)
        self.ltr_harvest_env_name = get_conda_env_from_path(self.config.conda_env_ltr_harvest)
        self.ltr_retriever_env_name = get_conda_env_from_path(self.config.conda_env_ltr_retriever)

        # EDTA实际消费的基因组路径(默认原基因组,需改名时替换为改名副本)
        # Genome path EDTA actually consumes (defaults to original; replaced by renamed copy when needed)
        self.edta_genome = self.config.genome

    # EDTA对序列名有13字符硬限制|EDTA has a hard 13-char limit on seq IDs
    # RepeatMasker/Census不支持超过13字符的序列名,EDTA内置两轮改名(截断/提数字)对HiC_scaffold_XX类名字会撞名失败
    # RepeatMasker/Census reject >13-char seq IDs; EDTA's 2 built-in rename attempts (truncate/extract-number) collide on names like HiC_scaffold_XX
    EDTA_SEQ_ID_MAX_LEN = 13

    def _prepare_genome_for_edta(self) -> str:
        """
        准备EDTA兼容的基因组|Prepare an EDTA-compatible genome

        EDTA要求序列名<=13字符且唯一|EDTA requires seq IDs <=13 chars and unique
        检测到超长或重复序列名时,生成唯一短名副本(如s0001)供EDTA使用,并保留旧名→新名映射表
        When long or duplicate IDs are detected, generates a uniquely-renamed copy (e.g. s0001) for EDTA, keeping an old->new map

        Returns:
            EDTA实际消费的基因组路径|Genome path EDTA will actually consume
        """
        genome = self.config.genome
        seq_ids = self._read_seq_ids(genome)

        if not seq_ids:
            self.logger.warning(f"未读取到序列名,使用原基因组|No seq IDs read, using original genome: {genome}")
            return genome

        max_len = max(len(sid) for sid in seq_ids)
        all_unique = len(set(seq_ids)) == len(seq_ids)
        if max_len <= self.EDTA_SEQ_ID_MAX_LEN and all_unique:
            self.logger.info(
                f"序列名最长{max_len}字符且唯一,EDTA兼容,无需改名|"
                f"Max seq ID length {max_len} <= {self.EDTA_SEQ_ID_MAX_LEN} and unique, no rename needed"
            )
            return genome

        # 需要改名: 生成唯一短名副本|Need rename: create uniquely-renamed copy
        genome_stem = Path(genome).stem
        renamed_path = self.working_dir / f"{genome_stem}.edta_renamed.fa"
        map_path = self.working_dir / f"{genome_stem}.edta_renamed.map.tsv"

        # 断点续传: 已存在则直接复用|Checkpoint resume: reuse if already exists
        if renamed_path.exists() and map_path.exists():
            self.logger.info(
                f"序列名最长{max_len}字符超限,复用已有改名副本|"
                f"Seq IDs up to {max_len} chars exceed limit, reusing renamed copy: {renamed_path}"
            )
            return str(renamed_path)

        self.logger.warning(
            f"序列名最长{max_len}字符超过EDTA {self.EDTA_SEQ_ID_MAX_LEN}字符限制,生成改名副本|"
            f"Seq IDs up to {max_len} chars exceed EDTA {self.EDTA_SEQ_ID_MAX_LEN}-char limit, creating renamed copy"
        )
        self.logger.info(f"改名映射表(旧名→新名)|Rename map (old->new): {map_path}")

        new_ids = self._generate_short_ids(len(seq_ids))
        self._write_renamed_fasta(genome, renamed_path, new_ids, map_path)

        self.logger.info(
            f"改名副本已生成|Renamed genome created: {renamed_path} "
            f"(共{len(seq_ids)}条序列|{len(seq_ids)} sequences)"
        )
        return str(renamed_path)

    def _read_seq_ids(self, genome_path: str) -> List[str]:
        """
        读取FASTA所有序列名(取>后第一个空白分隔token,与EDTA一致)|Read all seq IDs (first token after >, matching EDTA)

        Args:
            genome_path: 基因组FASTA路径|Genome FASTA path

        Returns:
            序列名列表(按文件顺序)|List of seq IDs (in file order)
        """
        seq_ids: List[str] = []
        with open(genome_path, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    # 取>后第一个token,与EDTA的(split)[0]一致|First token after >, matching EDTA's (split)[0]
                    header = line[1:].strip()
                    seq_ids.append(header.split()[0] if header else '')
        return seq_ids

    def _generate_short_ids(self, count: int) -> List[str]:
        """
        生成唯一短序列名(如s0001),保证<=13字符|Generate unique short seq IDs (e.g. s0001), guaranteed <=13 chars

        Args:
            count: 序列数量|Number of sequences

        Returns:
            短序列名列表|List of short seq IDs
        """
        # 宽度按序列数决定,前缀's',总长=1+width,远小于13|Width by count, prefix 's', total = 1+width << 13
        width = max(1, len(str(count)))
        return [f"s{str(i + 1).zfill(width)}" for i in range(count)]

    def _write_renamed_fasta(self, src: str, dst: Path, new_ids: List[str], map_path: Path) -> None:
        """
        写入改名后的FASTA及映射表|Write renamed FASTA and mapping table

        Args:
            src: 原基因组路径|Source genome path
            dst: 改名副本路径|Destination renamed FASTA path
            new_ids: 新序列名列表(顺序与原文件头一致)|New seq IDs (same order as source headers)
            map_path: 映射表路径(旧名\\t新名)|Map path (old_id\\tnew_id per line)

        Raises:
            IndexError: 序列数与预期不符|Sequence count mismatch
        """
        seq_idx = 0
        with open(src, 'r') as fin, open(dst, 'w') as fout, open(map_path, 'w') as fmap:
            for line in fin:
                if line.startswith('>'):
                    old_header = line[1:].strip()
                    old_id = old_header.split()[0] if old_header else ''
                    if seq_idx >= len(new_ids):
                        raise IndexError(
                            f"序列数({seq_idx + 1})超过预期({len(new_ids)}),请检查FASTA头是否变动|"
                            f"Sequence count ({seq_idx + 1}) exceeds expected ({len(new_ids)}), check FASTA headers"
                        )
                    new_id = new_ids[seq_idx]
                    fmap.write(f"{old_id}\t{new_id}\n")
                    fout.write(f">{new_id}\n")
                    seq_idx += 1
                else:
                    fout.write(line)

    def evaluate(self) -> Optional[Dict[str, Any]]:
        """运行LAI评估|Run LAI evaluation"""
        if self.config.skip_lai:
            self.logger.info("跳过LAI评估|Skipping LAI evaluation")
            return None

        self.logger.info("开始LAI指数评估|Starting LAI index evaluation")

        # 检查是否已完成|Check if already completed
        # LAI文件可能在多个位置，需要智能查找
        # LAI file may be in multiple locations, need intelligent search
        if self.config.resume:
            lai_file = self._find_existing_lai_file()
            if lai_file:
                self.logger.info(f"LAI评估已完成，跳过|LAI evaluation already completed, skipping: {lai_file}")
                return self._parse_lai_result(lai_file)

        # 运行LAI流程|Run LAI pipeline
        try:
            if self.config.lai_mode == 'edta':
                return self._run_edta_pipeline()
            else:
                return self._run_legacy_pipeline()

        except Exception as e:
            self.logger.error(f"LAI评估异常|LAI evaluation error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _run_edta_pipeline(self) -> Optional[Dict[str, Any]]:
        """
        使用EDTA流程运行LAI计算|Run LAI calculation using EDTA pipeline

        EDTA是一个完整的TE注释流程，集成了LTR_FINDER、LTRharvest和LTR_retriever
        EDTA is a complete TE annotation pipeline integrating LTR_FINDER, LTRharvest, and LTR_retriever
        """
        self.logger.info("使用EDTA流程计算LAI|Calculating LAI using EDTA pipeline")

        # 准备EDTA兼容的基因组(序列名<=13字符,必要时生成改名副本)|Prepare EDTA-compatible genome (seq IDs <=13 chars; rename if needed)
        self.edta_genome = self._prepare_genome_for_edta()

        # 获取EDTA.pl脚本路径|Get EDTA.pl script path
        # EDTA.pl通常在conda环境的bin目录中|EDTA.pl is usually in the bin directory of conda env
        edta_pl_path = os.path.join(self.config.conda_env_edta, "bin", "EDTA.pl")

        # 如果不存在，尝试查找其他可能的位置|If not found, try other possible locations
        if not os.path.exists(edta_pl_path):
            # 尝试share目录|Try share directory
            edta_pl_path = os.path.join(self.config.conda_env_edta, "share", "edta", "EDTA.pl")

        if not os.path.exists(edta_pl_path):
            # 尝试直接使用EDTA.pl命令（如果在PATH中）|Try using EDTA.pl command directly (if in PATH)
            edta_pl_path = "EDTA.pl"
            self.logger.warning(f"未找到EDTA.pl在标准路径，尝试使用系统PATH|EDTA.pl not found in standard path, trying system PATH")
        else:
            self.logger.info(f"使用EDTA.pl|Using EDTA.pl: {edta_pl_path}")

        # 构建EDTA命令|Build EDTA command
        # 参考|Reference: https://github.com/oushujun/EDTA
        #
        # EDTA参数说明|EDTA parameters:
        # --genome: 基因组FASTA文件|Genome FASTA file
        # --species: 物种类型 (Rice|Maize|others)|Species type
        # --step: 运行步骤 (all|filter|final|anno)|Which steps to run
        # --overwrite: 覆盖已有结果|Overwrite previous results
        # --threads: 线程数|Number of threads
        # --anno: 执行全基因组TE注释|Perform whole-genome TE annotation
        # --evaluate: 评估注释一致性|Evaluate annotation consistency

        genome_basename = Path(self.edta_genome).stem

        # EDTA线程: 内存随线程近线性增长, 默认封顶避免OOM(见config.EDTA_DEFAULT_THREAD_CAP)
        # EDTA threads: memory grows ~linearly with threads; capped by default to avoid OOM (see config.EDTA_DEFAULT_THREAD_CAP)
        edta_threads = self.config.edta_threads
        if self.config.lai_threads > edta_threads:
            self.logger.warning(
                f"EDTA线程受限|EDTA threads capped: {self.config.lai_threads} -> {edta_threads} "
                f"(EDTA内存随线程近线性增长, 300Mb基因组88线程实测峰值1.2TB易OOM|"
                f"EDTA memory scales ~linearly with threads, observed 1.2TB peak at 88 threads on 300Mb genome)"
            )
            self.logger.warning(
                f"如需更多线程且内存充足, 用 --edta-threads N 覆盖|"
                f"Override with --edta-threads N if memory is ample"
            )

        edta_args = [
            edta_pl_path,
            "--genome", self.edta_genome,
            "--species", "others",  # 非水稻/玉米物种|Non-rice/maize species
            "--step", "all",  # 运行完整流程|Run entire pipeline
            "--overwrite", "1",  # 覆盖已有结果|Overwrite previous results
            "--threads", str(edta_threads),
            "--anno", "1",  # 需要执行注释才能生成LAI文件|Need annotation to generate LAI file
            "--evaluate", "0",  # 跳过耗时的评估步骤|Skip time-consuming evaluation step
        ]

        # 使用conda run运行EDTA|Use conda run to execute EDTA
        cmd = build_conda_command(self.edta_env_name, "perl", edta_args)

        self.logger.info(f"运行EDTA流程|Running EDTA pipeline:")
        self.logger.info(f"  工作目录|Working directory: {self.working_dir}")
        self.logger.info(f"  命令|Command: {' '.join(cmd)}")

        # 运行EDTA|Run EDTA
        result = subprocess.run(
            cmd,
            cwd=str(self.working_dir),
            capture_output=True,
            text=True,
            check=False
        )

        # 记录输出|Log output
        if result.stdout:
            self.logger.debug(f"EDTA STDOUT:\n{result.stdout[:1000]}")
        if result.stderr:
            self.logger.warning(f"EDTA STDERR:\n{result.stderr[:1000]}")

        if result.returncode != 0:
            self.logger.error(f"EDTA运行失败|EDTA failed, return code: {result.returncode}")
            self.logger.error(f"错误信息|Error message: {result.stderr[:500]}")

            # 检查是否有LTR结果（即使EDTA失败）|Check if LTR results exist (even if EDTA failed)
            self.logger.info("检查EDTA是否生成了LTR结果|Checking if EDTA generated LTR results")
            passed_list = None
            for candidate in self.working_dir.glob("**/*.pass.list"):
                if "nmtf" not in candidate.name:
                    passed_list = candidate
                    self.logger.info(f"找到passed.list|Found passed.list: {passed_list}")
                    break

            if passed_list:
                self.logger.info("EDTA失败但找到了LTR结果，尝试手动完成注释|EDTA failed but LTR results found, attempting manual annotation")
                return self._complete_edta_manually(passed_list)

            return None

        self.logger.info("EDTA流程完成|EDTA pipeline completed")

        # 检查输出文件|Check output files
        # EDTA会在工作目录下创建一个以基因组名开头的目录
        # 例如: genome.fa.mod.EDTA.conda/
        edta_output_dir = None
        for item in self.working_dir.iterdir():
            if item.is_dir() and "EDTA" in item.name:
                edta_output_dir = item
                break

        if not edta_output_dir:
            self.logger.warning(f"未找到EDTA输出目录|EDTA output directory not found")
            # 尝试直接在工作目录中查找LAI文件
            # Try to find LAI file directly in working directory
            possible_lai_files = list(self.working_dir.glob("*.LAI"))
            if possible_lai_files:
                lai_file = possible_lai_files[0]
                self.logger.info(f"找到LAI文件|Found LAI file: {lai_file}")
                return self._parse_lai_result(lai_file)
            return None

        self.logger.info(f"找到EDTA输出目录|Found EDTA output directory: {edta_output_dir}")

        # 查找LAI文件|Find LAI file
        # EDTA会在EDTA_raw/LTR/目录下生成LAI文件
        # EDTA generates LAI file under EDTA_raw/LTR/
        lai_file = edta_output_dir / "EDTA_raw" / "LTR" / f"{genome_basename}.out.LAI"

        if not lai_file.exists():
            # 尝试其他可能的LAI文件位置
            # Try other possible LAI file locations
            possible_lai_files = list(edta_output_dir.glob("**/*.LAI"))
            if possible_lai_files:
                lai_file = possible_lai_files[0]
                self.logger.info(f"找到LAI文件|Found LAI file: {lai_file}")
                return self._parse_lai_result(lai_file)
            else:
                self.logger.warning(f"未找到LAI文件，尝试手动计算|LAI file not found, attempting manual calculation")
                # 尝试手动运行LAI计算
                # Try to run LAI calculation manually
                return self._run_lai_after_edta(edta_output_dir)

        # 解析LAI结果|Parse LAI results
        return self._parse_lai_result(lai_file)

    def _complete_edta_manually(self, passed_list: Path) -> Optional[Dict[str, Any]]:
        """
        EDTA失败后手动完成注释（运行RepeatMasker）|Manually complete annotation after EDTA failure (run RepeatMasker)

        Args:
            passed_list: passed.list文件路径|Path to passed.list file

        Returns:
            LAI结果字典|LAI results dictionary
        """
        self.logger.info("=" * 60)
        self.logger.info("手动完成EDTA注释|Manually completing EDTA annotation")
        self.logger.info("=" * 60)

        genome_mod = self.working_dir / f"{Path(self.edta_genome).name}.mod"
        if not genome_mod.exists():
            self.logger.error(f"未找到genome.fa.mod文件|Cannot find genome.fa.mod file: {genome_mod}")
            return None

        # 运行RepeatMasker生成.out文件|Run RepeatMasker to generate .out file
        self.logger.info("运行RepeatMasker生成注释文件|Running RepeatMasker to generate annotation file")

        # 不指定-species参数，让RepeatMasker使用默认库（适用于非模式生物）|Don't specify -species to use default library (for non-model organisms)
        # 如果需要特定物种库，用户应先配置RepeatMasker库|For specific species, users should pre-configure RepeatMasker library
        # -pa与EDTA --threads同理: 内存随并行度近线性增长, 用edta_threads(已封顶)而非lai_threads
        # -pa same as EDTA --threads: memory scales ~linearly with parallelism, use edta_threads (capped) not lai_threads
        cmd_str = (
            f"conda run -n {self.edta_env_name} --no-capture-output "
            f"RepeatMasker -pa {self.config.edta_threads} -gff -dir {self.working_dir} {genome_mod}"
        )

        self.logger.info(f"命令|Command: {cmd_str}")

        result = subprocess.run(
            ["bash", "-c", cmd_str],
            cwd=str(self.working_dir),
            capture_output=True,
            text=True,
            check=False,
            env=os.environ.copy()
        )

        if result.returncode != 0:
            self.logger.error(f"RepeatMasker运行失败|RepeatMasker failed: {result.stderr[:500]}")
            return None

        self.logger.info("RepeatMasker完成|RepeatMasker completed")

        # 查找生成的.out文件|Find generated .out file
        out_file = self.working_dir / f"{genome_mod.name}.out"
        if not out_file.exists():
            self.logger.error(f"未找到RepeatMasker输出文件|Cannot find RepeatMasker output: {out_file}")
            return None

        self.logger.info(f"找到.out文件|Found .out file: {out_file}")

        # 继续LAI计算|Continue with LAI calculation
        return self._run_lai_calculation(passed_list, out_file, genome_mod)

    def _run_lai_after_edta(self, edta_output_dir: Path) -> Optional[Dict[str, Any]]:
        """
        EDTA运行后手动计算LAI|Manually calculate LAI after EDTA run

        Args:
            edta_output_dir: EDTA输出目录|EDTA output directory

        Returns:
            LAI结果字典|LAI results dictionary
        """
        self.logger.info("手动计算LAI|Manually calculating LAI")

        # 在整个工作目录中查找必要的文件|Search for necessary files in entire working directory
        # 因为EDTA输出可能分散在多个目录中|Because EDTA outputs may be scattered in multiple directories

        # 1. 查找passed.list - 完整LTR列表|Find intact LTR list
        self.logger.info("查找passed.list文件|Searching for passed.list file")
        passed_list = None
        for candidate in self.working_dir.glob("**/*.pass.list"):
            # 排除nmtf.pass.list，优先使用普通的pass.list
            # Exclude nmtf.pass.list, prefer regular pass.list
            if "nmtf" not in candidate.name:
                passed_list = candidate
                self.logger.info(f"找到passed.list|Found passed.list: {passed_list}")
                break

        if not passed_list:
            # 如果没有找到普通pass.list，使用nmtf.pass.list
            # If no regular pass.list found, use nmtf.pass.list
            for candidate in self.working_dir.glob("**/*nmtf.pass.list"):
                passed_list = candidate
                self.logger.info(f"找到nmtf.pass.list|Found nmtf.pass.list: {passed_list}")
                break

        if not passed_list:
            self.logger.error("未找到passed.list文件|Cannot find passed.list file")
            return None

        # 2. 查找.out文件 - 所有LTR候选|Find all LTR candidates
        self.logger.info("查找.out文件|Searching for .out file")
        out_file = None
        for candidate in self.working_dir.glob("**/*.out"):
            # 优先选择EDTA.anno目录中的.out文件
            # Prefer .out file in EDTA.anno directory
            if "EDTA.anno" in str(candidate):
                out_file = candidate
                self.logger.info(f"找到.out文件|Found .out file: {out_file}")
                break

        # 如果没有在EDTA.anno中找到，查找任何.out文件
        # If not found in EDTA.anno, search for any .out file
        if not out_file:
            for candidate in self.working_dir.glob("*.mod.out"):
                out_file = candidate
                self.logger.info(f"找到.out文件|Found .out file: {out_file}")
                break

        if not out_file:
            self.logger.error("未找到.out文件|Cannot find .out file")
            return None

        # 运行LAI计算|Run LAI calculation
        genome_mod = self.working_dir / f"{Path(self.edta_genome).name}.mod"
        return self._run_lai_calculation(passed_list, out_file, genome_mod)

    def _run_lai_calculation(self, passed_list: Path, out_file: Path, genome_mod: Path) -> Optional[Dict[str, Any]]:
        """
        运行LAI计算|Run LAI calculation

        Args:
            passed_list: passed.list文件路径|Path to passed.list file
            out_file: .out文件路径|Path to .out file
            genome_mod: genome.fa.mod文件路径|Path to genome.fa.mod file

        Returns:
            LAI结果字典|LAI results dictionary
        """
        self.logger.info("=" * 60)
        self.logger.info("运行LAI计算|Running LAI calculation")
        self.logger.info("=" * 60)

        # 验证输入文件|Validate input files
        if not genome_mod.exists():
            self.logger.error(f"未找到genome.fa.mod文件|Cannot find genome.fa.mod file: {genome_mod}")
            return None

        # LAI输出文件名基于输入的.out文件|LAI output filename based on input .out file
        # 输入: genome.fa.mod.EDTA.TEanno.out -> 输出: genome.fa.mod.EDTA.TEanno.out.LAI
        output_lai = Path(str(out_file) + ".LAI")

        # 根据配置决定是否使用-qq快速模式|Decide whether to use -qq quick mode based on config
        # -qq: 跳过blastn计算（快速，5秒）|Skip blastn calculation (fast, 5 seconds)
        # 无-qq: 运行完整blastn（数小时，用于种间比较）|Run full blastn (hours, for interspecies comparison)
        quick_mode_flag = "-qq" if self.config.lai_quick_mode else ""

        # 记录模式选择|Log mode selection
        if self.config.lai_quick_mode:
            self.logger.info("使用LAI快速模式(-qq)，适用于种内比较|Using LAI quick mode (-qq) for within-species comparison")
        else:
            self.logger.warning("使用LAI完整模式（无-qq），将运行blastn计算，可能需要数小时|Using LAI full mode (no -qq), will run blastn calculation, may take hours")

        self.logger.info(f"运行LAI计算|Running LAI calculation")

        # 使用conda run调用LAI（从EDTA环境）|Use conda run to invoke LAI (from EDTA environment)
        # 注意：LAI工具位于EDTA环境，不在LTR_retriever环境中
        # Note: LAI tool is in EDTA environment, not in LTR_retriever environment
        cmd_str = (
            f"conda run -n {self.edta_env_name} --no-capture-output "
            f"LAI -genome {genome_mod} -intact {passed_list} -all {out_file} "
            f"-t {self.config.lai_threads} {quick_mode_flag}"
        )

        self.logger.debug(f"命令|Command: {cmd_str}")

        result = subprocess.run(
            ["bash", "-c", cmd_str],
            cwd=str(self.working_dir),
            capture_output=True,
            text=True,
            check=False,
            env=os.environ.copy()
        )

        if result.stdout:
            self.logger.debug(f"LAI计算STDOUT|LAI calculation STDOUT:\n{result.stdout[:500]}")
        if result.stderr:
            self.logger.warning(f"LAI计算STDERR|LAI calculation STDERR:\n{result.stderr[:500]}")

        # 检查LAI是否不适用于当前基因组|Check if LAI is not applicable
        if "LAI is not applicable" in result.stdout or "Total LTR sequence content" in result.stdout:
            self.logger.warning("LAI不适用于当前基因组（LTR含量低于5%）|LAI is not applicable for this genome (LTR content < 5%)")
            return {
                'lai_score': None,
                'raw_lai': None,
                'error': 'LAI not applicable: LTR content too low (< 5%)',
                'note': 'This genome has insufficient LTR content for LAI calculation'
            }

        if result.returncode != 0:
            self.logger.error(f"LAI计算失败|LAI calculation failed")
            if result.stderr:
                self.logger.error(f"错误信息|Error: {result.stderr}")
            if result.stdout:
                self.logger.error(f"输出信息|Output: {result.stdout}")
            return None

        # LAI脚本输出文件可能在working_dir根目录或anno子目录
        # LAI script may output to working_dir root or anno subdirectory
        # 尝试多个可能的位置|Try multiple possible locations
        possible_locations = [
            output_lai,  # anno子目录|anno subdirectory
            self.working_dir / f"{Path(out_file).name}.LAI",  # working_dir根目录|working_dir root
        ]

        for lai_file in possible_locations:
            if lai_file.exists():
                self.logger.info(f"找到LAI文件|Found LAI file: {lai_file}")
                return self._parse_lai_result(lai_file)

        # 如果上述位置都找不到，尝试glob递归查找
        # If not found in above locations, try glob recursive search
        lai_files = glob.glob(str(self.working_dir / "*.LAI"), recursive=False)
        if lai_files:
            self.logger.info(f"通过glob找到LAI文件|Found LAI file via glob: {lai_files[0]}")
            return self._parse_lai_result(Path(lai_files[0]))

        self.logger.error(f"LAI输出文件未生成|LAI output file not generated, searched in: {possible_locations}")
        return None

    def _run_legacy_pipeline(self) -> Optional[Dict[str, Any]]:
        """
        运行传统LTR_retriever流程（已弃用）|Run legacy LTR_retriever pipeline (deprecated)

        保留此方法以保持向后兼容性
        Keep this method for backward compatibility
        """
        self.logger.warning("使用传统LTR_retriever流程（已弃用，建议使用EDTA模式）|Using legacy LTR_retriever pipeline (deprecated, recommend using EDTA mode)")

        # 1. LTRharvest
        if self.config.lai_mode in ['full', 'harvest']:
            harvest_file = self._run_ltrharvest()
            if not harvest_file:
                return None

        # 2. LTR_retriever
        if self.config.lai_mode in ['full', 'retrieve']:
            retriever_file = self._run_ltr_retriever_pipeline()
            if not retriever_file:
                return None

        # 3. LAI计算
        if self.config.lai_mode in ['full', 'calculate']:
            return self._run_lai_calculation()

        return None

    def _run_ltrharvest(self) -> Optional[Path]:
        """运行LTRharvest识别LTR|Run LTRharvest to identify LTRs"""
        self.logger.info("运行LTRharvest|Running LTRharvest")

        # 建立索引|Build index
        genome_link = self.working_dir / "genome.fa"
        if not genome_link.exists():
            genome_link.symlink_to(self.config.genome)

        index_name = str(genome_link)
        suffixerator_args = [
            "suffixerator",
            "-db", str(genome_link),
            "-indexname", index_name,
            "-tis", "-suf", "-lcp", "-des", "-ssp", "-sds", "-dna",
            "-memlimit", "3500MB"
        ]
        cmd_index = build_conda_command(self.ltr_harvest_env_name, "gt", suffixerator_args)

        result = subprocess.run(
            cmd_index,
            cwd=str(self.working_dir),
            capture_output=True,
            text=True,
            check=False
        )
        if result.returncode != 0:
            self.logger.error(f"索引建立失败|Index building failed: {result.stderr}")
            return None

        # 运行harvest|Run harvest
        output_file = self.working_dir / f"{Path(self.config.genome).name}.harvest.scn"
        harvest_args = [
            "ltrharvest",
            "-index", index_name,
            "-minlenltr", "100",
            "-maxlenltr", "7000",
            "-mintsd", "4",
            "-maxtsd", "6",
            "-motif", "TGCA",
            "-motifmis", "1",
            "-similar", "85",
            "-vic", "10",
            "-seed", "20",
            "-seqids", "yes"
        ]
        cmd_harvest = build_conda_command(self.ltr_harvest_env_name, "gt", harvest_args)

        with open(output_file, 'w') as f:
            result = subprocess.run(
                cmd_harvest,
                cwd=str(self.working_dir),
                stdout=f,
                stderr=subprocess.PIPE,
                text=True
            )

        if result.returncode != 0:
            self.logger.error(f"LTRharvest运行失败|LTRharvest failed: {result.stderr}")
            return None

        self.logger.info(f"LTRharvest完成|LTRharvest completed: {output_file}")
        return output_file

    def _run_ltr_retriever_pipeline(self) -> Optional[Path]:
        """运行LTR_retriever流程|Run LTR_retriever pipeline"""
        self.logger.info("运行LTR_retriever|Running LTR_retriever")

        harvest_file = self.working_dir / f"{Path(self.config.genome).name}.harvest.scn"
        if not harvest_file.exists():
            self.logger.error(f"Harvest文件不存在|Harvest file not found: {harvest_file}")
            return None

        retriever_args = [
            "-adapter", "vt",
            "-genome", self.config.genome,
            "-harvest", str(harvest_file),
            "-outdir", str(self.working_dir),
            "-threads", str(self.config.lai_threads)
        ]
        cmd = build_conda_command(self.ltr_retriever_env_name, "LTR_retriever", retriever_args)

        result = subprocess.run(
            cmd,
            cwd=str(self.working_dir),
            capture_output=True,
            text=True,
            check=False
        )

        if result.stderr:
            self.logger.warning(f"LTR_retriever STDERR:\n{result.stderr}")

        if result.returncode != 0:
            self.logger.error(f"LTR_retriever运行失败|LTR_retriever failed: {result.stderr}")
            return None

        passed_list = self.working_dir / f"{Path(self.config.genome).name}.passed.list"
        if not passed_list.exists():
            self.logger.warning(f"LTR_retriever输出文件不存在|LTR_retriever output not found: {passed_list}")
            return None

        return passed_list

    def _parse_lai_result(self, lai_file: Path) -> Optional[Dict[str, Any]]:
        """解析LAI结果文件|Parse LAI result file

        LAI文件格式|LAI file format:
        Chr    From    To    Intact    Total    raw_LAI    LAI
        whole_genome    1    181288361    0.0517    0.2273    22.77    NA
        """
        try:
            self.logger.info(f"解析LAI结果文件|Parsing LAI result file: {lai_file}")

            with open(lai_file, 'r') as f:
                content = f.read()
                self.logger.debug(f"LAI文件内容|LAI file content:\n{content[:500]}")

            with open(lai_file, 'r') as f:
                for line in f:
                    # 查找包含whole_genome的行|Look for lines containing whole_genome
                    if "whole_genome" in line:
                        parts = line.split()
                        self.logger.debug(f"解析LAI行|Parsing LAI line: {line.strip()}")

                        # LAI文件格式|LAI file format:
                        # 0      1      2  3       4      5         6
                        # Chr    From   To  Intact  Total  raw_LAI   LAI
                        # 索引5是raw_LAI，索引6是LAI|Index 5 is raw_LAI, index 6 is LAI
                        if len(parts) >= 6:
                            try:
                                # 解析raw_LAI（索引5）|Parse raw_LAI (index 5)
                                raw_lai = float(parts[5])

                                # 尝试解析LAI（索引6，如果存在且不是NA）|Try to parse LAI (index 6, if exists and not NA)
                                lai_score = None
                                if len(parts) >= 7:
                                    try:
                                        if parts[6].upper() != 'NA':
                                            lai_score = float(parts[6])
                                    except (ValueError, IndexError):
                                        pass

                                result = {
                                    'raw_lai': raw_lai,
                                    'lai_score': lai_score if lai_score else raw_lai,
                                }

                                if lai_score:
                                    self.logger.info(f"全基因组LAI|Whole-genome LAI: {lai_score}, raw_LAI: {raw_lai}")
                                else:
                                    self.logger.info(f"全基因组raw_LAI|Whole-genome raw_LAI: {raw_lai} (快速模式|-qq mode)")

                                return result

                            except (ValueError, IndexError) as e:
                                self.logger.warning(f"无法解析LAI值|Cannot parse LAI value from: {line.strip()}, error: {e}")
                                continue

            # 如果没有找到whole_genome，尝试其他格式|If whole_genome not found, try other formats
            with open(lai_file, 'r') as f:
                for line in f:
                    if line.startswith("LAI"):
                        parts = line.split()
                        if len(parts) >= 2:
                            try:
                                lai_score = float(parts[1])
                                self.logger.info(f"从LAI行解析|Parsed from LAI line: {lai_score}")
                                return {
                                    'raw_lai': lai_score,
                                    'lai_score': lai_score,
                                }
                            except (ValueError, IndexError):
                                pass

            self.logger.error("无法解析LAI值|Cannot parse LAI value")
            return None

        except Exception as e:
            self.logger.error(f"解析LAI结果失败|Failed to parse LAI results: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _find_existing_lai_file(self) -> Optional[Path]:
        """
        查找已存在的LAI文件|Find existing LAI file

        Returns:
            找到的LAI文件路径，未找到返回None|Found LAI file path, None if not found
        """
        # 可能的LAI文件位置|Possible LAI file locations
        # 递归查找working_dir下所有.LAI文件(兼容改名前后的EDTA输出目录名)
        # Recursive search for all .LAI files (works with EDTA output dirs before/after renaming)
        lai_files = glob.glob(str(self.working_dir / "**" / "*.LAI"), recursive=True)

        if lai_files:
            # 如果找到多个，选择最新的|If multiple found, choose the latest
            lai_files.sort(key=lambda x: os.path.getmtime(x), reverse=True)
            return Path(lai_files[0])

        return None
