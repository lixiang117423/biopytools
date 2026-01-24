"""
 ADMIXTURE分析核心模块|ADMIXTURE Analysis Core Module 
"""

import os
import glob
import pandas as pd
from pathlib import Path
from .utils import CommandRunner

class AdmixtureAnalyzer:
    """ADMIXTURE分析器|ADMIXTURE Analyzer"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_admixture_analysis(self, plink_prefix: str):
        """ 运行ADMIXTURE分析|Run ADMIXTURE analysis"""
        self.logger.info("  开始ADMIXTURE分析|Starting ADMIXTURE analysis")
        
        bed_file = f"{plink_prefix}.bed"
        if not os.path.exists(bed_file):
            raise FileNotFoundError(f" PLINK文件不存在|PLINK file not found: {bed_file}")
        
        #  为每个K值运行ADMIXTURE|Run ADMIXTURE for each K value
        for k in range(self.config.min_k, self.config.max_k + 1):
            self._run_single_k(bed_file, k)
        
        #  计算交叉验证误差|Calculate cross-validation error
        best_k = self._find_best_k()
        
        self.logger.info(f" ADMIXTURE分析完成|ADMIXTURE analysis completed")
        self.logger.info(f" 最优K值|Best K value: {best_k}")
        
        return best_k
    
    def _run_single_k(self, bed_file: str, k: int):
        """ 为单个K值运行ADMIXTURE|Run ADMIXTURE for single K value"""
        self.logger.info(f" 运行ADMIXTURE K={k}|Running ADMIXTURE K={k}")
        
        #  设置输出文件名|Set output file names
        log_file = os.path.join(self.config.output_dir, f"log_{k}.out")
        
        #  构建命令，使用相对路径避免长路径问题
        bed_basename = os.path.basename(bed_file)
        cmd = (
            f"cd {self.config.output_dir} && "
            f"admixture --cv={self.config.cv_folds} -j{self.config.threads} "
            f"{bed_basename} {k} > log_{k}.out 2>&1"
        )
        
        #  运行命令|Run command
        try:
            result = self.cmd_runner.run(cmd, f"ADMIXTURE分析 K={k}|ADMIXTURE analysis K={k}")
            
            #  检查输出文件是否生成成功
            expected_q_file = os.path.join(self.config.output_dir, f"{os.path.splitext(bed_basename)[0]}.{k}.Q")
            if os.path.exists(expected_q_file):
                self.logger.info(f"  K={k} 分析成功完成|K={k} analysis completed successfully")
            else:
                self.logger.warning(f" K={k} 分析可能未完全成功，Q文件不存在|K={k} analysis may not be fully successful, Q file missing")
                
        except Exception as e:
            self.logger.error(f" K={k}分析失败|K={k} analysis failed: {e}")
            # 尝试读取日志文件了解错误原因
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    log_content = f.read()
                    self.logger.error(f" ADMIXTURE日志|ADMIXTURE log:\n{log_content}")
            raise
    
    def _find_best_k(self):
        """ 找到最优K值|Find best K value"""
        cv_results = []
        
        # 从log文件中提取CV误差|Extract CV error from log files
        for k in range(self.config.min_k, self.config.max_k + 1):
            log_file = os.path.join(self.config.output_dir, f"log_{k}.out")
            
            if os.path.exists(log_file):
                with open(log_file, 'r') as f:
                    for line in f:
                        if 'CV error' in line:
                            cv_error = float(line.strip().split()[-1])
                            cv_results.append({'K': k, 'CV_error': cv_error})
                            break
        
        if not cv_results:
            raise ValueError(" 未找到CV误差信息|No CV error information found")
        
        # 找到最小CV误差对应的K值|Find K with minimum CV error
        cv_df = pd.DataFrame(cv_results)
        best_k = cv_df.loc[cv_df['CV_error'].idxmin(), 'K']
        
        #  保存CV结果|Save CV results
        cv_file = os.path.join(self.config.output_dir, "cv_results.csv")
        cv_df.to_csv(cv_file, index=False)
        
        return int(best_k)

class ResultsProcessor:
    """ 结果处理器|Results Processor"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def process_results(self, best_k: int):
        """ 处理分析结果|Process analysis results"""
        self.logger.info(" 处理分析结果|Processing analysis results")
        
        #  读取Q文件（个体祖先成分）| Read Q file (individual ancestry proportions)
        q_file = os.path.join(self.config.output_dir, f"{self.config.base_name}.{best_k}.Q")
        q_data = self._read_q_file(q_file, best_k)
        
        #  读取P文件（等位基因频率）| Read P file (allele frequencies)
        p_file = os.path.join(self.config.output_dir, f"{self.config.base_name}.{best_k}.P")
        p_data = self._read_p_file(p_file, best_k)
        
        #  计算统计信息|Calculate statistics
        stats = self._calculate_statistics(q_data, best_k)
        
        #  保存处理后的结果|Save processed results
        self._save_processed_results(q_data, p_data, stats, best_k)
        
        return q_data, p_data, stats
    
    def _read_q_file(self, q_file: str, k: int):
        """ 读取Q文件|Read Q file"""
        if not os.path.exists(q_file):
            raise FileNotFoundError(f" Q文件不存在|Q file not found: {q_file}")
        
        q_data = pd.read_csv(q_file, sep=r'\s+', header=None)
        q_data.columns = [f"Pop{i+1}" for i in range(k)]
        
        #  添加个体信息|Add individual information
        fam_file = os.path.join(self.config.output_dir, f"{self.config.base_name}.fam")
        if os.path.exists(fam_file):
            fam_data = pd.read_csv(fam_file, sep=r'\s+', header=None)
            q_data['FID'] = fam_data.iloc[:, 0]
            q_data['IID'] = fam_data.iloc[:, 1]
        
        return q_data
    
    def _read_p_file(self, p_file: str, k: int):
        """ 读取P文件|Read P file"""
        if not os.path.exists(p_file):
            self.logger.warning(f" P文件不存在|P file not found: {p_file}")
            return None
        
        p_data = pd.read_csv(p_file, sep=r'\s+', header=None)
        p_data.columns = [f"Pop{i+1}" for i in range(k)]
        
        return p_data
    
    def _calculate_statistics(self, q_data: pd.DataFrame, k: int):
        """ 计算统计信息|Calculate statistics"""
        pop_cols = [f"Pop{i+1}" for i in range(k)]
        
        #  计算每个个体的最大祖先成分|Calculate max ancestry for each individual
        max_ancestry = q_data[pop_cols].max(axis=1)
        
        # 计算混合程度|Calculate admixture level
        admixture_level = 1 - max_ancestry
        
        #  统计信息|Statistics
        stats = {
            'total_individuals': len(q_data),
            'highly_admixed': sum(admixture_level > 0.3),
            'pure_individuals': sum(max_ancestry > 0.9),
            'mean_admixture_level': admixture_level.mean(),
            'mean_max_ancestry': max_ancestry.mean()
        }
        
        return stats
    
    def _save_processed_results(self, q_data: pd.DataFrame, p_data: pd.DataFrame, 
                               stats: dict, best_k: int):
        """ 保存处理后的结果|Save processed results"""
        # 保存详细Q矩阵|Save detailed Q matrix
        q_output = os.path.join(self.config.output_dir, "admixture_proportions.csv")
        q_data.to_csv(q_output, index=False)
        
        #  保存统计信息|Save statistics
        stats_output = os.path.join(self.config.output_dir, "admixture_statistics.txt")
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
        self.logger.info(f"  -  个体祖先成分|Individual ancestry proportions: {q_output}")
        self.logger.info(f"  -  统计信息|Statistics: {stats_output}")