"""
 ADMIXTURE分析核心模块|ADMIXTURE Analysis Core Module

运行各K值的ADMIXTURE/ADAMIXTURE、按CV误差/LogLikelihood选最优K、处理Q/P矩阵。
所有命令走 build_conda_command 传完整工具路径(§13.6);执行前记录完整命令(§2.2.1);
per-K支持断点续传(§10.2);dry_run只记录命令不执行。
"""

import os
import subprocess
import pandas as pd
from .utils import CommandRunner, build_conda_command


class AdmixtureAnalyzer:
    """ADMIXTURE分析器|ADMIXTURE Analyzer"""

    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def run_admixture_analysis(self, plink_prefix: str):
        """ 运行ADMIXTURE分析|Run ADMIXTURE analysis

        Returns:
            最优K值;dry_run或无法判定时返回None|Best K value; None on dry_run or if undetermined
        """
        self.logger.info("开始ADMIXTURE分析|Starting ADMIXTURE analysis")

        bed_file = f"{plink_prefix}.bed"
        if not os.path.exists(bed_file):
            raise FileNotFoundError(f" PLINK文件不存在|PLINK file not found: {bed_file}")

        #  为每个K值运行ADMIXTURE|Run ADMIXTURE for each K value
        for k in range(self.config.min_k, self.config.max_k + 1):
            self._run_single_k(bed_file, k)

        # dry_run时跳过最优K判定|Skip best-K determination on dry_run
        if self.config.dry_run:
            self.logger.info("模拟运行,跳过最优K判定|Dry run, skipping best-K determination")
            return None

        #  计算交叉验证误差/LogLikelihood选最优K|Determine best K
        best_k = self._find_best_k()

        self.logger.info(f" ADMIXTURE分析完成|ADMIXTURE analysis completed")
        self.logger.info(f" 最优K值|Best K value: {best_k}")

        return best_k

    def _run_single_k(self, bed_file: str, k: int):
        """ 为单个K值运行ADMIXTURE或ADAMIXTURE|Run ADMIXTURE or ADAMIXTURE for single K value"""
        bed_basename = os.path.basename(bed_file)
        base_prefix = os.path.splitext(bed_basename)[0]

        # 断点续传:该K的Q已存在则跳过(§10.2)|Checkpoint: skip if this K's Q already exists
        q_out = os.path.join(self.config.admixture_dir, f"{base_prefix}.{k}.Q")
        if (not self.config.force) and os.path.exists(q_out) and os.path.getsize(q_out) > 0:
            self.logger.info(f"跳过已完成K|Skipping completed K={k} ({q_out})")
            return

        if self.config.method == "adamixture":
            self._run_adamixture_k(bed_file, k, base_prefix)
        else:
            self._run_admixture_k(bed_file, k, base_prefix)

    def _run_admixture_k(self, bed_file: str, k: int, base_prefix: str):
        """ 运行ADMIXTURE|Run ADMIXTURE (传完整路径,§13.6)"""
        self.logger.info(f"运行ADMIXTURE K={k}|Running ADMIXTURE K={k}")

        log_file = os.path.join(self.config.admixture_dir, f"log_{k}.out")

        # 使用conda wrapper构建命令(传完整admixture路径)|Build cmd with full tool path
        args = [
            '--cv', str(self.config.cv_folds),
            '-j' + str(self.config.threads),
            bed_file,   # 绝对路径;ADMIXTURE按basename命名输出,写cwd|absolute path
            str(k),
        ]
        cmd = build_conda_command(self.config.admixture_path, args)

        # 记录完整命令(§2.2.1)|Log full command
        self.logger.info(f"执行|Executing: ADMIXTURE K={k}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        # dry_run只记录不执行|Dry run: log only
        if self.config.dry_run:
            self.logger.info(f"模拟运行,跳过执行|Dry run, skipping execution: ADMIXTURE K={k}")
            return

        try:
            with open(log_file, 'w') as log_f:
                subprocess.run(
                    cmd,
                    shell=False,
                    check=True,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    cwd=self.config.admixture_dir,
                )

            expected_q_file = os.path.join(self.config.admixture_dir, f"{base_prefix}.{k}.Q")
            if os.path.exists(expected_q_file):
                self.logger.info(f"K={k} 分析成功完成|K={k} analysis completed successfully")
            else:
                self.logger.warning(f"K={k} 分析可能未完全成功，Q文件不存在|K={k} analysis may be incomplete, Q file missing")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"K={k}分析失败|K={k} analysis failed: {e}")
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    self.logger.error(f"ADMIXTURE日志|ADMIXTURE log:\n{f.read()}")
            raise

    def _run_adamixture_k(self, bed_file: str, k: int, base_prefix: str):
        """ 运行ADAMIXTURE|Run ADAMIXTURE (传完整路径,§13.6)"""
        self.logger.info(f"运行ADAMIXTURE K={k}|Running ADAMIXTURE K={k}")

        log_file = os.path.join(self.config.admixture_dir, f"adamixture_{k}.log")

        # 使用conda wrapper构建命令(传完整adamixture路径,禁止basename提取)|Full path, no basename
        args = [
            '--k', str(k),
            '--data_path', bed_file,
            '--save_dir', self.config.admixture_dir,
            '--name', base_prefix,
            '--threads', str(self.config.threads),
            '--seed', str(self.config.adamixture_seed),
            '--lr', str(self.config.adamixture_lr),
            '--beta1', str(self.config.adamixture_beta1),
            '--beta2', str(self.config.adamixture_beta2),
            '--max_iter', str(self.config.adamixture_max_iter),
        ]
        cmd = build_conda_command(self.config.adamixture_path, args)

        self.logger.info(f"执行|Executing: ADAMIXTURE K={k}")
        self.logger.info(f"命令|Command: {' '.join(cmd)}")

        if self.config.dry_run:
            self.logger.info(f"模拟运行,跳过执行|Dry run, skipping execution: ADAMIXTURE K={k}")
            return

        try:
            with open(log_file, 'w') as log_f:
                subprocess.run(
                    cmd,
                    shell=False,
                    check=True,
                    stdout=log_f,
                    stderr=subprocess.STDOUT,
                    cwd=self.config.admixture_dir,
                )

            # ADAMIXTURE可能输出两种格式|ADAMIXTURE may output two name formats
            q_f1 = os.path.join(self.config.admixture_dir, f"{base_prefix}.{k}.Q")
            q_f2 = os.path.join(self.config.admixture_dir, f"{base_prefix}.{k}.{k}.Q")
            if os.path.exists(q_f1) or os.path.exists(q_f2):
                self.logger.info(f"K={k} 分析成功完成|K={k} analysis completed successfully")
            else:
                self.logger.warning(f"K={k} 分析可能未完全成功，Q文件不存在|K={k} analysis may be incomplete, Q file missing")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"K={k}分析失败|K={k} analysis failed: {e}")
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    self.logger.error(f"ADAMIXTURE日志|ADAMIXTURE log:\n{f.read()}")
            raise

    def _find_best_k(self):
        """ 找到最优K值|Find best K value"""
        if self.config.method == "adamixture":
            return self._find_best_k_adamixture()
        return self._find_best_k_admixture()

    def _find_best_k_admixture(self):
        """ 使用CV误差找到最优K值|Find best K using CV error (ADMIXTURE)"""
        cv_results = []

        for k in range(self.config.min_k, self.config.max_k + 1):
            log_file = os.path.join(self.config.admixture_dir, f"log_{k}.out")
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    for line in f:
                        if 'CV error' in line:
                            cv_error = float(line.strip().split()[-1])
                            cv_results.append({'K': k, 'CV_error': cv_error})
                            break

        if not cv_results:
            raise ValueError(" 未找到CV误差信息|No CV error information found")

        cv_df = pd.DataFrame(cv_results)
        best_k = cv_df.loc[cv_df['CV_error'].idxmin(), 'K']

        cv_file = os.path.join(self.config.admixture_dir, "cv_results.csv")
        cv_df.to_csv(cv_file, index=False)

        return int(best_k)

    def _find_best_k_adamixture(self):
        """ 使用log-likelihood找到最优K值|Find best K using log-likelihood (ADAMIXTURE)"""
        ll_results = []

        for k in range(self.config.min_k, self.config.max_k + 1):
            log_file = os.path.join(self.config.admixture_dir, f"adamixture_{k}.log")
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    for line in f:
                        if 'Log-likelihood' in line or 'log-likelihood' in line:
                            try:
                                ll_value = float(line.split(':')[-1].strip())
                                ll_results.append({'K': k, 'LogLikelihood': ll_value})
                                self.logger.info(f"K={k} Log-likelihood: {ll_value}")
                                break
                            except (ValueError, IndexError):
                                continue

        if not ll_results:
            self.logger.warning(" 未找到Log-likelihood信息，返回最大K值|No log-likelihood found, returning max K")
            return self.config.max_k

        ll_df = pd.DataFrame(ll_results)
        best_k = ll_df.loc[ll_df['LogLikelihood'].idxmax(), 'K']

        ll_file = os.path.join(self.config.admixture_dir, "loglikelihood_results.csv")
        ll_df.to_csv(ll_file, index=False)

        self.logger.info(f"最佳K值|Best K value: {best_k} (Log-likelihood: {ll_df.loc[ll_df['LogLikelihood'].idxmax(), 'LogLikelihood']})")

        return int(best_k)


class ResultsProcessor:
    """ 结果处理器|Results Processor"""

    def __init__(self, config, logger):
        self.config = config
        self.logger = logger

    def process_results(self, best_k: int):
        """ 处理分析结果|Process analysis results"""
        self.logger.info("处理分析结果|Processing analysis results")

        q_file, p_file = self._locate_q_p_files(best_k)

        #  读取Q文件（个体祖先成分）|Read Q file (individual ancestry proportions)
        q_data = self._read_q_file(q_file, best_k)

        #  读取P文件（等位基因频率）|Read P file (allele frequencies)
        p_data = self._read_p_file(p_file, best_k)

        #  计算统计信息|Calculate statistics
        stats = self._calculate_statistics(q_data, best_k)

        #  保存处理后的结果|Save processed results
        self._save_processed_results(q_data, p_data, stats, best_k)

        return q_data, p_data, stats

    def _locate_q_p_files(self, best_k: int):
        """ 定位Q/P文件(03_admixture/,兼容ADAMIXTURE双格式)|Locate Q/P files (in 03_admixture/)"""
        base = self.config.base_name
        q_f1 = os.path.join(self.config.admixture_dir, f"{base}.{best_k}.Q")
        q_f2 = os.path.join(self.config.admixture_dir, f"{base}.{best_k}.{best_k}.Q")
        p_f1 = os.path.join(self.config.admixture_dir, f"{base}.{best_k}.P")
        p_f2 = os.path.join(self.config.admixture_dir, f"{base}.{best_k}.{best_k}.P")

        q_file = q_f1 if os.path.exists(q_f1) else q_f2
        p_file = p_f1 if os.path.exists(p_f1) else p_f2

        if not os.path.exists(q_file):
            raise FileNotFoundError(
                f"Q文件不存在|Q file not found: tried {q_f1} and {q_f2}")
        return q_file, p_file

    def _read_q_file(self, q_file: str, k: int):
        """ 读取Q文件|Read Q file"""
        if not os.path.exists(q_file):
            raise FileNotFoundError(f" Q文件不存在|Q file not found: {q_file}")

        q_data = pd.read_csv(q_file, sep=r'\s+', header=None)
        q_data.columns = [f"Pop{i+1}" for i in range(k)]

        #  添加个体信息(从02_plink/的fam)|Add individual info from fam (in 02_plink/)
        fam_file = os.path.join(self.config.plink_dir, f"{self.config.base_name}.fam")
        if os.path.exists(fam_file):
            fam_data = pd.read_csv(fam_file, sep=r'\s+', header=None)
            q_data['FID'] = fam_data.iloc[:, 0].values
            q_data['IID'] = fam_data.iloc[:, 1].values

        return q_data

    def _read_p_file(self, p_file: str, k: int):
        """ 读取P文件|Read P file"""
        if not os.path.exists(p_file):
            self.logger.warning(f"P文件不存在|P file not found: {p_file}")
            return None

        p_data = pd.read_csv(p_file, sep=r'\s+', header=None)
        p_data.columns = [f"Pop{i+1}" for i in range(k)]
        return p_data

    def _calculate_statistics(self, q_data: pd.DataFrame, k: int):
        """ 计算统计信息|Calculate statistics"""
        pop_cols = [f"Pop{i+1}" for i in range(k)]
        max_ancestry = q_data[pop_cols].max(axis=1)
        admixture_level = 1 - max_ancestry

        stats = {
            'total_individuals': len(q_data),
            'highly_admixed': int((admixture_level > 0.3).sum()),
            'pure_individuals': int((max_ancestry > 0.9).sum()),
            'mean_admixture_level': float(admixture_level.mean()),
            'mean_max_ancestry': float(max_ancestry.mean()),
        }
        return stats

    def _save_processed_results(self, q_data: pd.DataFrame, p_data: pd.DataFrame,
                               stats: dict, best_k: int):
        """ 保存处理后的结果到04_results/|Save processed results to 04_results/"""
        q_output = os.path.join(self.config.results_dir, "admixture_proportions.csv")
        q_data.to_csv(q_output, index=False)

        stats_output = os.path.join(self.config.results_dir, "admixture_statistics.txt")
        with open(stats_output, 'w') as f:
            f.write(" ADMIXTURE分析统计信息|ADMIXTURE Analysis Statistics\n")
            f.write("=" * 80 + "\n\n")
            f.write(f" 最优K值|Best K value: {best_k}\n")
            f.write(f"总个体数|Total individuals: {stats['total_individuals']}\n")
            f.write(f" 高度混合个体数|Highly admixed individuals: {stats['highly_admixed']}\n")
            f.write(f"纯合个体数|Pure individuals: {stats['pure_individuals']}\n")
            f.write(f"平均混合程度|Mean admixture level: {stats['mean_admixture_level']:.3f}\n")
            f.write(f" 平均最大祖先成分|Mean max ancestry: {stats['mean_max_ancestry']:.3f}\n")

        self.logger.info(f" 结果已保存|Results saved:")
        self.logger.info(f"  - 个体祖先成分|Individual ancestry proportions: {q_output}")
        self.logger.info(f"  - 统计信息|Statistics: {stats_output}")
