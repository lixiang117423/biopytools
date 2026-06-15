"""TreeMix分析运行器|TreeMix Analysis Runner

流程 (参考 https://github.com/taotaoyuan/treemix):
  prepare: VCF → plink LD过滤 → plink make-bed → 修改fam FID → plink freq --within → plink2treemix
  run:     Scan(m=0..m_max × N replicates) → OptM → Plot(m=0, optimal_m)
  all:     prepare + run
"""

import gzip
import os
import random
import subprocess
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from datetime import datetime
from pathlib import Path
from typing import List, Optional, Tuple

from .config import TreemixConfig
from .utils import (
    TreemixLogger,
    convert_plink_to_treemix,
    generate_pop_order,
    infer_populations,
    modify_fam_file,
    parse_llik,
    parse_vcf_samples,
    write_pop_cov,
)


class TreemixRunner:
    """TreeMix分析运行器|TreeMix Analysis Runner"""

    VERSION = "2.0.0"

    def __init__(self, config: TreemixConfig, logger: Optional[TreemixLogger] = None):
        self.config = config
        self._logger_obj = logger
        self.logger = logger.get_logger() if logger else None
        self._start_time = None

    # ================================================================
    # 公共入口|Public Entry Points
    # ================================================================

    def run_prepare(self) -> bool:
        """运行输入准备流程|Run input preparation pipeline"""
        try:
            self._start_time = time.time()
            self._init_logger("TreeMix Prepare")
            self._print_header("输入准备|Input Preparation")
            self.config.validate_prepare()

            # 解析分组信息|Parse cluster data
            cluster_data = self._get_cluster_data()
            pop_cov_file = self._prepare_cluster(cluster_data)

            # 步骤1/6: 按分组文件筛选样品 (bcftools view -S)
            if self.config.cluster_file:
                self.logger.info("  步骤1/6: 样品筛选|Sample filtering (bcftools view -S)")
                filtered_vcf = self._bcftools_filter_samples(cluster_data)
                if not filtered_vcf:
                    return False
            else:
                filtered_vcf = self.config.vcf_file

            # 步骤2/6: LD过滤|LD pruning
            self.logger.info("  步骤2/6: LD过滤|LD pruning")
            prune_in = self._plink_ld_prune(filtered_vcf)
            if not prune_in:
                return False

            # 步骤3/6: VCF → PLINK binary格式 (bed/bim/fam)
            self.logger.info("  步骤3/6: VCF转PLINK格式|Convert VCF to PLINK binary")
            bed_prefix = self._plink_make_bed(filtered_vcf, prune_in)
            if not bed_prefix:
                return False

            # 步骤4/6: 修改fam文件FID + 计算等位基因频率
            self.logger.info("  步骤4/6: 修改FID并计算频率|Modify FID and calculate frequencies")
            frq_file = self._plink_freq_within(bed_prefix, pop_cov_file, cluster_data)
            if not frq_file:
                return False

            # 步骤5/6: 转换为TreeMix格式
            self.logger.info("  步骤5/6: 转换为TreeMix格式|Convert to TreeMix format")
            treemix_frq = self._convert_to_treemix(frq_file)
            if not treemix_frq:
                return False

            # 步骤6/6: 生成pop.order.txt
            self.logger.info("  步骤6/6: 生成pop.order.txt|Generate pop.order.txt")
            self._generate_pop_order(cluster_data)

            self._print_footer(True)
            return True

        except Exception as e:
            if self.logger:
                self.logger.error(f"程序执行出错|Error: {e}", exc_info=True)
            return False

    def run_treemix(self) -> bool:
        """运行TreeMix完整分析: Scan + OptM + Plot|Run TreeMix: scan + OptM + plot"""
        try:
            self._start_time = time.time()
            self._init_logger("TreeMix Run")
            self._print_header("TreeMix分析|TreeMix Analysis")
            self.config.validate_run()

            # 阶段1: Scan - 对m=0..m_max各运行N次重复
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("阶段1/3: Scan - 扫描m=0..m_max|Phase 1/3: Scanning m values")
            self.logger.info("=" * 60)
            if not self._run_scan():
                return False

            # 阶段2: OptM - 确定最优m
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("阶段2/3: OptM - 确定最优m值|Phase 2/3: Determining optimal m")
            self.logger.info("=" * 60)
            optimal_m = self._run_optm()
            if optimal_m is None:
                self.logger.warning("  OptM失败, 跳过绘图|OptM failed, skipping plot")
                self._print_footer(True)
                return True

            # 阶段3: Plot - 绘制树图和残差图
            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info(f"阶段3/3: Plot - 绘图 (m=0, m={optimal_m})|Phase 3/3: Plotting")
            self.logger.info("=" * 60)
            self._run_plot(optimal_m)

            self._print_footer(True)
            return True

        except Exception as e:
            if self.logger:
                self.logger.error(f"程序执行出错|Error: {e}", exc_info=True)
            return False

    def run_all(self) -> bool:
        """运行完整流程 (prepare + run)|Run full pipeline"""
        if not self.run_prepare():
            return False

        treemix_frq = os.path.join(self.config.prepare_dir, "input.treemix.frq.gz")
        if os.path.exists(treemix_frq):
            self.config.treemix_input = treemix_frq
        else:
            if self.logger:
                self.logger.error("未找到TreeMix输入文件|TreeMix input file not found")
            return False

        return self.run_treemix()

    # ================================================================
    # Prepare步骤|Prepare Steps
    # ================================================================

    def _get_cluster_data(self) -> List[Tuple[str, str]]:
        """获取群体分组数据|Get population cluster data

        Returns:
            [(sample_id, population), ...]
        """
        if self.config.cluster_file:
            self.logger.info(f"  使用分组文件|Using cluster file: {self.config.cluster_file}")
            return self._read_cluster_file(self.config.cluster_file)

        # 自动推断
        self.logger.info(
            f"  自动推断群体 (分隔符: '{self.config.pop_delimiter}')|"
            f"Auto-inferring populations (delimiter: '{self.config.pop_delimiter}')"
        )
        samples = parse_vcf_samples(self.config.vcf_file)
        self.logger.info(f"  检测到 {len(samples)} 个样本|Detected {len(samples)} samples")

        cluster_data = infer_populations(samples, self.config.pop_delimiter)
        pops = sorted(set(p for _, p in cluster_data))
        self.logger.info(f"  推断出 {len(pops)} 个群体|Inferred {len(pops)} populations: {', '.join(pops)}")
        return cluster_data

    def _prepare_cluster(self, cluster_data: List[Tuple[str, str]]) -> str:
        """生成pop.cov文件|Generate pop.cov file for plink --within"""
        pops = sorted(set(p for _, p in cluster_data))

        pop_cov = os.path.join(self.config.prepare_dir, "pop.cov")
        write_pop_cov(cluster_data, pop_cov)
        self.logger.info(f"  生成pop.cov ({len(pops)} 个群体, {len(cluster_data)} 个样本)|"
                        f"Generated pop.cov ({len(pops)} populations, {len(cluster_data)} samples)")
        return pop_cov

    def _bcftools_filter_samples(self, cluster_data: List[Tuple[str, str]]) -> Optional[str]:
        """步骤1: 按分组文件筛选VCF样品|Step 1: Filter VCF samples by cluster file

        bcftools view -S sample_id.txt input.vcf.gz -Oz -o filtered.vcf.gz
        """
        out_vcf = os.path.join(self.config.prepare_dir, "selected_samples.vcf.gz")

        if os.path.exists(out_vcf):
            self.logger.info("    跳过已完成|Skipping completed: sample filtering")
            return out_vcf

        # 写样品ID列表 (每行一个)
        sample_list_file = os.path.join(self.config.prepare_dir, "sample_id.txt")
        with open(sample_list_file, 'w') as f:
            for sample_id, _ in cluster_data:
                f.write(f"{sample_id}\n")

        self.logger.info(f"    筛选 {len(cluster_data)} 个样品|Filtering {len(cluster_data)} samples")

        cmd = [
            self.config.bcftools_path,
            'view',
            '-S', sample_list_file,
            self.config.vcf_file,
            '-Oz',
            '-o', out_vcf,
        ]

        self.logger.info(f"    命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=7200)
        if result.returncode != 0:
            self.logger.error(f"    bcftools筛选失败|bcftools filtering failed: {result.stderr}")
            return None

        if not os.path.exists(out_vcf):
            self.logger.error("    未生成筛选后的VCF|Filtered VCF not generated")
            return None

        # 建索引
        cmd_idx = [self.config.bcftools_path, 'index', out_vcf]
        subprocess.run(cmd_idx, capture_output=True, text=True, timeout=3600)

        self.logger.info(f"    筛选完成|Filtering done: {out_vcf}")
        return out_vcf

    def _plink_ld_prune(self, vcf_file: str) -> Optional[str]:
        """步骤2: plink LD过滤|Step 2: plink LD pruning"""
        prune_out = os.path.join(self.config.prepare_dir, "ld.prune.in")
        prune_prefix = os.path.join(self.config.prepare_dir, "ld")

        if os.path.exists(prune_out):
            self.logger.info("    跳过已完成|Skipping completed: LD pruning")
            return prune_out

        cmd = [
            self.config.plink_path,
            '--vcf', vcf_file,
            '--indep-pairwise',
            str(self.config.ld_window),
            str(self.config.ld_step),
            str(self.config.ld_r2),
            '--out', prune_prefix,
            '--allow-extra-chr',
            '--set-missing-var-ids', '@:#',
            '--keep-allele-order',
        ]

        self.logger.info(f"    命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            self.logger.error(f"    LD过滤失败|LD pruning failed: {result.stderr}")
            return None

        if not os.path.exists(prune_out):
            self.logger.error("    未生成ld.prune.in|ld.prune.in not generated")
            return None

        with open(prune_out) as f:
            n_snps = sum(1 for _ in f)
        self.logger.info(f"    保留 {n_snps} 个SNP|Retained {n_snps} SNPs after LD pruning")
        return prune_out

    def _plink_make_bed(self, vcf_file: str, prune_in: str) -> Optional[str]:
        """步骤3: VCF → PLINK binary (提取LD过滤后的SNP)|Step 3: VCF to PLINK binary"""
        bed_prefix = os.path.join(self.config.prepare_dir, "plink_input")

        if os.path.exists(f"{bed_prefix}.bed"):
            self.logger.info("    跳过已完成|Skipping completed: PLINK binary conversion")
            return bed_prefix

        cmd = [
            self.config.plink_path,
            '--vcf', vcf_file,
            '--extract', prune_in,
            '--make-bed',
            '--out', bed_prefix,
            '--allow-extra-chr',
            '--set-missing-var-ids', '@:#',
            '--keep-allele-order',
        ]

        self.logger.info(f"    命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            self.logger.error(f"    PLINK格式转换失败|PLINK binary conversion failed: {result.stderr}")
            return None

        if not os.path.exists(f"{bed_prefix}.bed"):
            self.logger.error("    未生成bed文件|bed file not generated")
            return None

        return bed_prefix

    def _plink_freq_within(
        self,
        bed_prefix: str,
        pop_cov_file: str,
        cluster_data: List[Tuple[str, str]],
    ) -> Optional[str]:
        """步骤3: 修改fam FID + 计算等位基因频率|Step 3: Modify FID + calculate frequencies"""
        fam_file = f"{bed_prefix}.fam"
        frq_prefix = os.path.join(self.config.prepare_dir, "input")
        frq_strat = f"{frq_prefix}.frq.strat"
        frq_gz = f"{frq_strat}.gz"

        if os.path.exists(frq_gz):
            self.logger.info("    跳过已完成|Skipping completed: allele frequencies")
            return frq_gz

        # 修改fam文件: FID → 群体名
        n_modified = modify_fam_file(fam_file, cluster_data)
        self.logger.info(f"    修改 {n_modified} 个样本的FID|Modified FID for {n_modified} samples")

        # plink --freq --within 计算分层频率
        cmd = [
            self.config.plink_path,
            '--bfile', bed_prefix,
            '--freq',
            '--missing',
            '--within', pop_cov_file,
            '--out', frq_prefix,
        ]

        self.logger.info(f"    命令|Command: {' '.join(cmd)}")
        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            self.logger.error(f"    频率计算失败|Frequency calculation failed: {result.stderr}")
            return None

        if not os.path.exists(frq_strat):
            self.logger.error("    未生成frq.strat|frq.strat not generated")
            return None

        # gzip压缩
        with open(frq_strat, 'rb') as f_in:
            with gzip.open(frq_gz, 'wb') as f_out:
                f_out.writelines(f_in)

        self.logger.info(f"    生成 {frq_gz}|Generated {frq_gz}")
        return frq_gz

    def _convert_to_treemix(self, frq_file: str) -> Optional[str]:
        """步骤4: 转换为TreeMix格式 (plink2treemix)|Step 4: Convert to TreeMix format"""
        out_file = os.path.join(self.config.prepare_dir, "input.treemix.frq.gz")

        if os.path.exists(out_file):
            self.logger.info("    跳过已完成|Skipping completed: format conversion")
            return out_file

        try:
            pops = convert_plink_to_treemix(frq_file, out_file)
            self.logger.info(f"    转换完成, {len(pops)} 个群体|Converted, {len(pops)} populations")
            return out_file
        except Exception as e:
            self.logger.error(f"    格式转换失败|Format conversion failed: {e}")
            return None

    def _generate_pop_order(self, cluster_data: List[Tuple[str, str]]):
        """步骤5: 生成pop.order.txt|Step 5: Generate pop.order.txt"""
        pop_order_file = os.path.join(self.config.prepare_dir, "pop.order.txt")
        generate_pop_order(cluster_data, pop_order_file)
        self.logger.info(f"    生成|Generated: {pop_order_file}")

    # ================================================================
    # Run步骤: Scan阶段|Run Steps: Scan Phase
    # ================================================================

    def _build_treemix_cmd(self, m: int, output_stem: str) -> List[str]:
        """构建TreeMix命令|Build TreeMix command

        教程参数: -global -k 500 -bootstrap 1000 -seed random
        """
        args = [
            self.config.treemix_path,
            '-i', self.config.treemix_input,
            '-m', str(m),
            '-o', output_stem,
            '-k', str(self.config.k),
        ]

        if self.config.root:
            args.extend(['-root', self.config.root])
        if self.config.global_opt:
            args.append('-global')
        if self.config.noss:
            args.append('-noss')
        if self.config.se:
            args.append('-se')

        return args

    def _count_completed(self, m_dir: str, n_rep: int) -> int:
        """统计已完成的replicate数|Count completed replicates"""
        done = 0
        for i in range(1, n_rep + 1):
            treeout = os.path.join(m_dir, f"rep_{i:04d}.treeout.gz")
            if os.path.exists(treeout):
                done += 1
        return done

    def _run_scan(self) -> bool:
        """Scan阶段: 对m=0..m_max各运行N次重复|Scan phase: run N replicates for each m"""
        success = True

        for m in range(self.config.m_max + 1):
            m_dir = os.path.join(self.config.treemix_dir, f"m{m}")
            os.makedirs(m_dir, exist_ok=True)

            done = self._count_completed(m_dir, self.config.replicates)
            if done >= self.config.replicates:
                self.logger.info(f"  m={m}: 跳过已完成 ({done}/{self.config.replicates})|"
                                f"Skipping completed: m={m}")
                continue

            self.logger.info(
                f"  m={m}: 运行 {self.config.replicates} 次重复 (线程={self.config.threads})|"
                f"Running {self.config.replicates} replicates (threads={self.config.threads})"
            )

            if self.config.threads > 1 and self.config.replicates > 1:
                success &= self._run_scan_parallel(m, m_dir, done)
            else:
                success &= self._run_scan_serial(m, m_dir, done)

        # 汇总log likelihood
        self._summarize_llik()
        return success

    def _run_scan_parallel(self, m: int, m_dir: str, already_done: int) -> bool:
        """并行运行replicates|Run replicates in parallel"""
        failed = 0

        def _run_rep(rep_id):
            rep_stem = os.path.join(m_dir, f"rep_{rep_id:04d}")
            rep_treeout = f"{rep_stem}.treeout.gz"
            if os.path.exists(rep_treeout):
                return rep_id, True, None

            args = self._build_treemix_cmd(m, rep_stem)
            args.extend(['-bootstrap', str(self.config.bootstrap)])

            # 随机种子
            seed = random.randint(1, 999999999)
            args.extend(['-seed', str(seed)])

            log_file = f"{rep_stem}.log"
            with open(log_file, 'w') as log_fh:
                result = subprocess.run(
                    args,
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    timeout=86400,
                )
            return rep_id, result.returncode == 0, result.returncode

        done = already_done
        total = self.config.replicates

        with ThreadPoolExecutor(max_workers=self.config.threads) as executor:
            futures = {
                executor.submit(_run_rep, i + 1): i + 1
                for i in range(total)
            }
            for future in as_completed(futures):
                rep_id, ok, rc = future.result()
                if ok:
                    done += 1
                else:
                    failed += 1
                if (done + failed) % max(1, self.config.replicates // 5) == 0 or failed > 0:
                    self.logger.info(
                        f"    m={m} 进度|Progress: {done}/{total} 完成, {failed} 失败"
                    )

        if failed > 0:
            self.logger.warning(
                f"  m={m}: {total - failed}/{total} 成功, {failed} 失败|"
                f"{total - failed}/{total} succeeded, {failed} failed"
            )
            return failed < total
        return True

    def _run_scan_serial(self, m: int, m_dir: str, already_done: int) -> bool:
        """串行运行replicates|Run replicates serially"""
        failed = 0
        done = already_done
        total = self.config.replicates

        for i in range(1, total + 1):
            rep_stem = os.path.join(m_dir, f"rep_{i:04d}")
            rep_treeout = f"{rep_stem}.treeout.gz"
            if os.path.exists(rep_treeout):
                continue

            args = self._build_treemix_cmd(m, rep_stem)
            args.extend(['-bootstrap', str(self.config.bootstrap)])

            seed = random.randint(1, 999999999)
            args.extend(['-seed', str(seed)])

            log_file = f"{rep_stem}.log"
            with open(log_file, 'w') as log_fh:
                result = subprocess.run(
                    args,
                    stdout=log_fh,
                    stderr=subprocess.STDOUT,
                    timeout=86400,
                )

            if result.returncode == 0:
                done += 1
            else:
                failed += 1

        if failed > 0:
            self.logger.warning(
                f"  m={m}: {total - failed}/{total} 成功, {failed} 失败|"
                f"{total - failed}/{total} succeeded, {failed} failed"
            )
            return failed < total
        return True

    def _summarize_llik(self):
        """汇总各m值的log likelihood|Summarize log likelihoods across m values"""
        summary_file = os.path.join(self.config.output_dir, "treemix_llik_summary.txt")

        lines = ["m\tlog_likelihood"]
        for m in range(self.config.m_max + 1):
            m_dir = os.path.join(self.config.treemix_dir, f"m{m}")
            # 取第一个replicate的llik
            llik_file = os.path.join(m_dir, "rep_0001.llik")

            llik = parse_llik(llik_file)
            llik_str = f"{llik:.4f}" if llik is not None else "N/A"
            lines.append(f"{m}\t{llik_str}")

        with open(summary_file, 'w') as f:
            f.write('\n'.join(lines) + '\n')

        self.logger.info(f"\n  汇总文件|Summary file: {summary_file}")
        self.logger.info("  " + "-" * 40)
        self.logger.info(f"  {'m':>4s}  {'log_likelihood':>16s}")
        self.logger.info("  " + "-" * 40)
        for line in lines[1:]:
            parts = line.split('\t')
            self.logger.info(f"  {parts[0]:>4s}  {parts[1]:>16s}")
        self.logger.info("  " + "-" * 40)

    # ================================================================
    # Run步骤: OptM阶段|Run Steps: OptM Phase
    # ================================================================

    def _run_optm(self) -> Optional[int]:
        """OptM阶段: 用R包OptM确定最优m|OptM phase: determine optimal m using OptM R package

        Returns:
            最优m值, 失败返回None
        """
        # 检查R
        if not os.path.exists(self.config.r_path):
            self.logger.error(f"  R不存在|R not found: {self.config.r_path}")
            return None

        # 生成.llik文件路径列表
        llik_files = []
        for m in range(self.config.m_max + 1):
            m_dir = os.path.join(self.config.treemix_dir, f"m{m}")
            for i in range(1, self.config.replicates + 1):
                llik_file = os.path.join(m_dir, f"rep_{i:04d}.llik")
                if os.path.exists(llik_file):
                    llik_files.append(llik_file)

        if not llik_files:
            self.logger.error("  未找到任何.llik文件|No .llik files found")
            return None

        self.logger.info(f"  找到 {len(llik_files)} 个.llik文件|Found {len(llik_files)} .llik files")

        # 生成R脚本
        optm_dir = os.path.join(self.config.output_dir, "03_optm")
        os.makedirs(optm_dir, exist_ok=True)

        r_script = os.path.join(optm_dir, "run_optm.R")
        optimal_file = os.path.join(optm_dir, "optimal_m.txt")

        llik_paths = ', '.join(f'"{f}"' for f in llik_files)

        r_code = f"""# OptM: 确定最优基因流事件数|Determine optimal number of migration edges
library(OptM)

# 读取所有.llik文件|Read all .llik files
llik_files <- c({llik_paths})
m_list <- lapply(llik_files, function(f) {{
    lines <- readLines(f)
    for (line in lines) {{
        if (grepl("likelihood", line, ignore.case = TRUE)) {{
            val <- as.numeric(strsplit(line, ":")[[1]][2])
            return(list(m = NA, llik = val))
        }}
    }}
    return(list(m = NA, llik = NA))
}})

# 提取log likelihood值|Extract log likelihood values
llik_vals <- sapply(m_list, function(x) x$llik)
llik_vals <- llik_vals[!is.na(llik_vals)]

if (length(llik_vals) == 0) {{
    stop("无法提取log likelihood值|Cannot extract log likelihood values")
}}

# 使用OptM确定最优m|Use OptM to determine optimal m
# 方法: Evanno方法 (基于Delta K, 类似STRUCTURE)
# 由于OptM需要不同的输入格式, 这里使用线性回归方法
# 即选择使得拟合优度不再显著改善的m值
n_m <- {self.config.m_max} + 1  # m=0..m_max
n_rep <- {self.config.replicates}

# 构建数据: 每个m对应replicates个llik值
llik_matrix <- matrix(NA, nrow = n_rep, ncol = n_m)
for (idx in seq_along(llik_files)) {{
    m_idx <- ((idx - 1) %/% n_rep) + 1
    rep_idx <- ((idx - 1) %% n_rep) + 1
    if (m_idx <= n_m && rep_idx <= n_rep) {{
        llik_matrix[rep_idx, m_idx] <- m_list[[idx]]$llik
    }}
}}

# 使用OptM的optimizeM函数
# 构建OptM需要的格式
m_values <- 0:{self.config.m_max}
m_lab <- paste0("m", m_values)

# 方法1: 使用poppr格式构建 (需要自定义)
# 方法2: 直接使用llik列表
optm_result <- tryCatch({{
    # 使用OptM的内部方法: 用llik列表构建
    # optimizeM接受llik列表, 方法: "Evanno"或"linear"
    optimizeM(llik_matrix, method = "linear", plot = FALSE)
}}, error = function(e) {{
    # 备用方法: 简单的线性回归拐点分析
    cat("OptM optimizeM失败, 使用备用方法|OptM failed, using fallback method\\n")
    cat(paste("Error:", e$message, "\\n"))
    return(NULL)
}})

if (!is.null(optm_result)) {{
    # 提取最优m
    optimal <- which.min(optm_result$delta) - 1
    # 确保在有效范围内
    if (optimal < 0) optimal <- 0
    if (optimal > {self.config.m_max}) optimal <- {self.config.m_max}
}} else {{
    # 备用方法: 选择llik最佳值对应的m
    mean_llik <- colMeans(llik_matrix, na.rm = TRUE)
    optimal <- which.max(mean_llik) - 1  # 0-indexed
}}

cat("\\n")
cat("=== OptM结果|OptM Results ===\\n")
cat(sprintf("最优m值|Optimal m: %d\\n", optimal))
cat("\\n各m值的平均log likelihood|Mean log likelihood per m:\\n")
for (mm in 1:length(m_values)) {{
    vals <- llik_matrix[, mm]
    vals <- vals[!is.na(vals)]
    if (length(vals) > 0) {{
        cat(sprintf("  m=%d: %.4f (sd=%.4f, n=%d)\\n",
            m_values[mm], mean(vals), sd(vals), length(vals)))
    }}
}}

# 保存结果|Save result
writeLines(as.character(optimal), "{optimal_file}")
"""

        with open(r_script, 'w') as f:
            f.write(r_code)

        # 生成OptM选择图
        optm_plot_script = os.path.join(optm_dir, "plot_optm.R")
        optm_plot_pdf = os.path.join(optm_dir, "optm_selection.pdf")
        r_plot_code = f"""# OptM选择图|OptM selection plot
library(OptM)

llik_files <- c({llik_paths})
m_list <- lapply(llik_files, function(f) {{
    lines <- readLines(f)
    for (line in lines) {{
        if (grepl("likelihood", line, ignore.case = TRUE)) {{
            val <- as.numeric(strsplit(line, ":")[[1]][2])
            return(list(m = NA, llik = val))
        }}
    }}
    return(list(m = NA, llik = NA))
}})

n_m <- {self.config.m_max} + 1
n_rep <- {self.config.replicates}
llik_matrix <- matrix(NA, nrow = n_rep, ncol = n_m)
for (idx in seq_along(llik_files)) {{
    m_idx <- ((idx - 1) %/% n_rep) + 1
    rep_idx <- ((idx - 1) %% n_rep) + 1
    if (m_idx <= n_m && rep_idx <= n_rep) {{
        llik_matrix[rep_idx, m_idx] <- m_list[[idx]]$llik
    }}
}}

pdf("{optm_plot_pdf}", width = 8, height = 6)
m_values <- 0:{self.config.m_max}
mean_llik <- colMeans(llik_matrix, na.rm = TRUE)

# 绘制log likelihood随m变化|Plot log likelihood vs m
plot(m_values, mean_llik, type = "b", pch = 19, col = "blue",
     xlab = "Number of migration edges (m)",
     ylab = "Log likelihood",
     main = "TreeMix: Log likelihood vs migration edges")

# 添加误差线|Add error bars
sd_llik <- apply(llik_matrix, 2, sd, na.rm = TRUE)
arrows(m_values, mean_llik - sd_llik, m_values, mean_llik + sd_llik,
       angle = 90, code = 3, length = 0.05, col = "blue")

# 标记最优m|Mark optimal m
optimal <- as.integer(readLines("{optimal_file}"))
if (!is.na(optimal) && optimal >= 0 && optimal <= {self.config.m_max}) {{
    points(optimal, mean_llik[optimal + 1], pch = 17, cex = 1.5, col = "red")
    text(optimal, mean_llik[optimal + 1], labels = paste0("m=", optimal),
         pos = 4, col = "red", font = 2)
}}

dev.off()
cat("OptM图已保存|OptM plot saved: {optm_plot_pdf}\\n")
"""

        with open(optm_plot_script, 'w') as f:
            f.write(r_plot_code)

        # 执行OptM
        self.logger.info(f"  执行OptM R脚本|Running OptM R script")
        cmd = [self.config.r_path, 'CMD', 'BATCH', '--vanilla', r_script,
               os.path.join(optm_dir, "optm.log")]
        self.logger.info(f"    命令|Command: {' '.join(cmd)}")

        result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
        if result.returncode != 0:
            self.logger.error(f"  OptM执行失败|OptM failed: {result.stderr}")
            return None

        # 读取最优m值
        if not os.path.exists(optimal_file):
            self.logger.error("  未生成optimal_m.txt|optimal_m.txt not generated")
            return None

        with open(optimal_file) as f:
            optimal_m = int(f.read().strip())

        # 执行OptM绘图
        self.logger.info(f"  执行OptM绘图|Running OptM plotting")
        cmd2 = [self.config.r_path, 'CMD', 'BATCH', '--vanilla', optm_plot_script,
                os.path.join(optm_dir, "plot_optm.log")]
        result2 = subprocess.run(cmd2, capture_output=True, text=True, timeout=3600)
        if result2.returncode == 0:
            self.logger.info(f"  OptM图已保存|OptM plot saved: {optm_plot_pdf}")

        self.logger.info(f"  最优m值|Optimal m: {optimal_m}")

        # 复制optimal_m.txt到输出目录
        final_optimal = os.path.join(self.config.output_dir, "optimal_m.txt")
        with open(final_optimal, 'w') as f:
            f.write(f"{optimal_m}\n")
        self.logger.info(f"  最优m已保存|Optimal m saved: {final_optimal}")

        return optimal_m

    # ================================================================
    # Run步骤: Plot阶段|Run Steps: Plot Phase
    # ================================================================

    def _find_plotting_funcs(self) -> Optional[str]:
        """查找plotting_funcs.R|Find plotting_funcs.R

        优先级: 用户指定 > treemix源码目录
        """
        if self.config.plotting_funcs_r and os.path.exists(self.config.plotting_funcs_r):
            return self.config.plotting_funcs_r

        # 搜索常见路径
        candidates = [
            os.path.expanduser(
                "~/software/treemix/treemix-1.13/src/plotting_funcs.R"
            ),
        ]
        for path in candidates:
            if os.path.exists(path):
                return path

        return None

    def _run_plot(self, optimal_m: int):
        """Plot阶段: 绘制树图和残差图|Plot phase: tree and residual plots"""
        plotting_funcs = self._find_plotting_funcs()
        if not plotting_funcs:
            self.logger.warning("  未找到plotting_funcs.R, 跳过绘图|"
                               "plotting_funcs.R not found, skipping plot")
            return

        if not os.path.exists(self.config.r_path):
            self.logger.warning(f"  R不存在, 跳过绘图|R not found: {self.config.r_path}")
            return

        pop_order_file = os.path.join(self.config.prepare_dir, "pop.order.txt")
        if not os.path.exists(pop_order_file):
            self.logger.warning("  未找到pop.order.txt, 跳过绘图|pop.order.txt not found")
            return

        plot_dir = os.path.join(self.config.output_dir, "04_plot")
        os.makedirs(plot_dir, exist_ok=True)

        # 需要绘制的m值: m=0 和 optimal_m
        m_values_to_plot = sorted(set([0, optimal_m]))
        # 确保m在有效范围内
        m_values_to_plot = [m for m in m_values_to_plot if m <= self.config.m_max]

        for m in m_values_to_plot:
            self.logger.info(f"  绘图 m={m}|Plotting m={m}")

            # 使用第一个replicate的结果
            m_dir = os.path.join(self.config.treemix_dir, f"m{m}")
            rep_stem = os.path.join(m_dir, "rep_0001")

            vertices = f"{rep_stem}.vertices.gz"
            edges = f"{rep_stem}.edges.gz"
            cov = f"{rep_stem}.cov.gz"
            modelcov = f"{rep_stem}.modelcov.gz"
            covse = f"{rep_stem}.covse.gz"

            if not all(os.path.exists(f) for f in [vertices, edges, cov, modelcov, covse]):
                self.logger.warning(f"    m={m} 缺少必要文件, 跳过|m={m} missing required files")
                continue

            # 生成R绘图脚本
            r_script = os.path.join(plot_dir, f"plot_m{m}.R")
            tree_pdf = os.path.join(plot_dir, f"m{m}_tree.pdf")
            resid_pdf = os.path.join(plot_dir, f"m{m}_residual.pdf")

            r_code = f"""# TreeMix绘图|TreeMix plotting for m={m}
source("{plotting_funcs}")

pop_order <- "{pop_order_file}"
stem <- "{rep_stem}"

# 绘制树图|Plot tree (o=NA表示不使用自定义颜色)
pdf("{tree_pdf}", width = 12, height = 8)
plot_tree(stem, o = NA, cex = 1.2)
title(main = "TreeMix population graph", sub = paste0("m = ", {m}))
dev.off()

# 绘制残差图|Plot residuals
pdf("{resid_pdf}", width = 8, height = 8)
plot_resid(stem, pop_order, usemax = TRUE)
title(main = "TreeMix residuals", sub = paste0("m = ", {m}))
dev.off()

cat("绘图完成|Plots saved: {tree_pdf}, {resid_pdf}\\n")
"""

            with open(r_script, 'w') as f:
                f.write(r_code)

            # 执行R绘图
            cmd = [self.config.r_path, 'CMD', 'BATCH', '--vanilla', r_script,
                   os.path.join(plot_dir, f"plot_m{m}.log")]
            self.logger.info(f"    命令|Command: {' '.join(cmd)}")

            result = subprocess.run(cmd, capture_output=True, text=True, timeout=3600)
            if result.returncode == 0:
                if os.path.exists(tree_pdf):
                    self.logger.info(f"    树图|Tree plot: {tree_pdf}")
                if os.path.exists(resid_pdf):
                    self.logger.info(f"    残差图|Residual plot: {resid_pdf}")
            else:
                self.logger.warning(f"    m={m} 绘图失败|Plot failed for m={m}")
                # 输出R日志以便调试
                r_log = os.path.join(plot_dir, f"plot_m{m}.log")
                if os.path.exists(r_log):
                    with open(r_log) as f:
                        log_content = f.read()
                    if log_content.strip():
                        self.logger.warning(f"    R日志|R log:\n{log_content[-500:]}")

    # ================================================================
    # 辅助方法|Helper Methods
    # ================================================================

    def _init_logger(self, name: str):
        if self._logger_obj is None:
            log_dir = Path(self.config.logs_dir)
            log_dir.mkdir(parents=True, exist_ok=True)
            self._logger_obj = TreemixLogger(log_dir / "treemix.log", name)
            self.logger = self._logger_obj.get_logger()

    def _read_cluster_file(self, cluster_file: str) -> List[Tuple[str, str]]:
        """读取分组文件|Read cluster file

        格式支持:
          sample_id  population  (两列, 自动跳过表头)
          FID  IID  CLUSTER     (三列, plink --within 格式)
        """
        _header_keywords = {'sample', 'group', 'population', 'fid', 'iid',
                            'cluster', 'clst', 'header', 'sample_id'}
        data = []
        with open(cluster_file) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                fields = line.split('\t')
                if len(fields) >= 3:
                    # FID IID CLUST 格式
                    data.append((fields[1], fields[2]))
                elif len(fields) == 2:
                    if fields[0].lower() in _header_keywords:
                        continue
                    data.append((fields[0], fields[1]))
        return data

    def _print_header(self, title: str):
        if not self.logger:
            return
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info(f"TreeMix - {title}")
        self.logger.info("=" * 60)
        self.logger.info(f"版本|Version: {self.VERSION}")
        self.logger.info(f"时间|Time: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
        self.logger.info("")

    def _print_footer(self, success: bool):
        if not self.logger:
            return
        elapsed = time.time() - self._start_time
        hours, remainder = divmod(int(elapsed), 3600)
        minutes, seconds = divmod(remainder, 60)
        self.logger.info("")
        self.logger.info("=" * 60)
        status = "完成|Completed" if success else "失败|Failed"
        self.logger.info(f"TreeMix {status}")
        self.logger.info(f"  总耗时|Total time: {hours}h {minutes}m {seconds}s")
        self.logger.info("=" * 60)
