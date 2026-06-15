"""
LAI模块主程序|LAI Module Main Program
"""

import argparse
import sys
from pathlib import Path
from .config import LAIConfig, LTRHarvestConfig, LTRFinderConfig, LTRRetrieverConfig, LAICalculateConfig
from .utils import LAILogger, check_dependencies
from .ltrharvest import LTRHarvestRunner
from .ltrfinder import LTRFinderRunner
from .ltr_retriever import LTRRetrieverRunner
from .lai import LAICalculatorRunner


class LAICalculator:
    """LAI计算器主类|LAI Calculator Main Class"""

    def __init__(self, **kwargs):
        # 初始化配置|Initialize configuration
        self.config = LAIConfig(**kwargs)
        self.config.validate()

        # 初始化子配置|Initialize sub-configurations
        self.harvest_config = LTRHarvestConfig()
        self.finder_config = LTRFinderConfig()
        self.retriever_config = LTRRetrieverConfig()
        self.lai_config = LAICalculateConfig(threads=kwargs.get('threads', 64))

        # 初始化日志|Initialize logging
        log_file = self.config.output_path / "lai.log"
        self.logger_manager = LAILogger(log_file)
        self.logger = self.logger_manager.get_logger()

    def run(self):
        """运行LAI计算流程|Run LAI calculation pipeline"""
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始LAI分析流程|Starting LAI analysis pipeline")
            self.logger.info("=" * 60)

            # 检查依赖|Check dependencies
            self.logger.info("检查依赖软件|Checking dependencies")
            if not check_dependencies(self.config, self.logger):
                sys.exit(1)

            # 根据模式运行|Run based on mode
            if self.config.mode in ['full', 'harvest']:
                # 步骤1-2: 识别LTR候选|Steps 1-2: Identify LTR candidates
                harvest_result = self._run_harvest()
                if not harvest_result:
                    raise RuntimeError("LTR候选识别失败|LTR candidate identification failed")

            if self.config.mode in ['full', 'retrieve']:
                # 步骤3: 筛选LTR|Step 3: Filter LTRs
                retriever_result = self._run_retriever()
                if not retriever_result:
                    raise RuntimeError("LTR筛选失败|LTR filtering failed")

            if self.config.mode in ['full', 'calculate']:
                # 步骤4: 计算LAI|Step 4: Calculate LAI
                lai_result = self._run_lai()
                if not lai_result:
                    raise RuntimeError("LAI计算失败|LAI calculation failed")

            self.logger.info("=" * 60)
            self.logger.info("LAI分析流程完成|LAI analysis pipeline completed")
            self.logger.info("=" * 60)
            self.logger.info(f"输出目录|Output directory: {self.config.output_dir}")
            self.logger.info(f"日志文件|Log file: {self.logger_manager.log_file}")

            sys.exit(0)

        except Exception as e:
            self.logger.error(f"流程执行失败|Pipeline execution failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            sys.exit(1)

    def _run_harvest(self) -> dict:
        """
        运行LTR候选识别|Run LTR candidate identification

        Returns:
            dict: 候选文件路径|Candidate file paths
        """
        try:
            # 检查输出文件是否已存在|Check if output files already exist
            prefix = self.config.genome_path.name
            harvest_file = self.config.output_path / f"{prefix}.harvest.scn"
            finder_file = self.config.output_path / f"{prefix}.finder.combine.scn"
            raw_ltr_file = self.config.output_path / f"{prefix}.rawLTR.scn"

            # 运行LTRharvest|Run LTRharvest
            harvest_runner = LTRHarvestRunner(
                self.config,
                self.logger,
                self.harvest_config
            )

            # 检查harvest结果是否已存在|Check if harvest result already exists
            if self.config.skip_completed and harvest_file.exists():
                self.logger.info("=" * 60)
                self.logger.info("检测到已完成的LTRharvest结果，跳过此步骤|Detected completed LTRharvest results, skipping this step")
                self.logger.info("=" * 60)
                self.logger.info(f"使用现有文件|Using existing file: {harvest_file}")
                harvest_result = harvest_file
            else:
                harvest_result = harvest_runner.run()
                if not harvest_result:
                    return None

            # 运行LTR_FINDER|Run LTR_FINDER
            finder_runner = LTRFinderRunner(
                self.config,
                self.logger,
                self.finder_config
            )

            # 检查finder结果是否已存在|Check if finder result already exists
            if self.config.skip_completed and finder_file.exists():
                self.logger.info("=" * 60)
                self.logger.info("检测到已完成的LTR_FINDER结果，跳过此步骤|Detected completed LTR_FINDER results, skipping this step")
                self.logger.info("=" * 60)
                self.logger.info(f"使用现有文件|Using existing file: {finder_file}")
                finder_result = finder_file
            else:
                finder_result = finder_runner.run()
                if not finder_result:
                    return None

            # 合并结果|Merge results
            # 检查合并文件是否已存在|Check if merged file already exists
            if self.config.skip_completed and raw_ltr_file.exists():
                self.logger.info("=" * 60)
                self.logger.info("检测到已完成的合并结果，跳过此步骤|Detected completed merge results, skipping this step")
                self.logger.info("=" * 60)
                self.logger.info(f"使用现有文件|Using existing file: {raw_ltr_file}")
                merged_result = raw_ltr_file
            else:
                merged_result = self._merge_candidates(harvest_result, finder_result)
                if not merged_result:
                    return None

            return {
                'harvest': harvest_result,
                'finder': finder_result,
                'merged': merged_result
            }

        except Exception as e:
            self.logger.error(f"LTR候选识别失败|LTR candidate identification failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return None

    def _run_retriever(self) -> dict:
        """
        运行LTR筛选|Run LTR filtering

        Returns:
            dict: 输出文件路径|Output file paths
        """
        try:
            # 查找合并的候选文件|Find merged candidate file
            prefix = self.config.genome_path.name
            raw_ltr_file = self.config.output_path / f"{prefix}.rawLTR.scn"

            if not raw_ltr_file.exists():
                self.logger.error(f"未找到合并的候选文件|Merged candidate file not found: {raw_ltr_file}")
                return None

            # 检查输出文件是否已存在|Check if output files already exist
            pass_list = self.config.output_path / f"{prefix}.pass.list"
            out_file = self.config.output_path / f"{prefix}.out"

            # 如果输出文件存在且启用了跳过选项|If output files exist and skip is enabled
            if self.config.skip_completed:
                if pass_list.exists() and out_file.exists():
                    self.logger.info("=" * 60)
                    self.logger.info("检测到已完成的LTR筛选结果，跳过此步骤|Detected completed LTR filtering results, skipping this step")
                    self.logger.info("=" * 60)
                    self.logger.info(f"使用现有文件|Using existing files:")
                    self.logger.info(f"  - pass.list: {pass_list}")
                    self.logger.info(f"  - .out: {out_file}")
                    return {
                        'pass_list': pass_list,
                        'out': out_file
                    }

            # 运行LTR_retriever|Run LTR_retriever
            retriever_runner = LTRRetrieverRunner(
                self.config,
                self.logger,
                self.retriever_config
            )
            result = retriever_runner.run(raw_ltr_file)

            return result

        except Exception as e:
            self.logger.error(f"LTR筛选失败|LTR filtering failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return None

    def _run_lai(self) -> Path:
        """
        运行LAI计算|Run LAI calculation

        Returns:
            Path: LAI输出文件路径|LAI output file path
        """
        try:
            # 查找必需的文件|Find required files
            prefix = self.config.genome_path.name
            pass_list = self.config.output_path / f"{prefix}.pass.list"
            out_file = self.config.output_path / f"{prefix}.out"
            lai_file = self.config.output_path / f"{prefix}.LAI"

            if not pass_list.exists():
                self.logger.error(f"未找到pass.list文件|pass.list file not found: {pass_list}")
                return None

            if not out_file.exists():
                self.logger.error(f"未找到.out文件|.out file not found: {out_file}")
                return None

            # 检查输出文件是否已存在|Check if output file already exists
            if self.config.skip_completed and lai_file.exists():
                self.logger.info("=" * 60)
                self.logger.info("检测到已完成的LAI计算结果，跳过此步骤|Detected completed LAI calculation results, skipping this step")
                self.logger.info("=" * 60)
                self.logger.info(f"使用现有文件|Using existing file: {lai_file}")
                return lai_file

            # 运行LAI计算|Run LAI calculation
            lai_runner = LAICalculatorRunner(
                self.config,
                self.logger,
                self.lai_config
            )
            result = lai_runner.run(pass_list, out_file)

            return result

        except Exception as e:
            self.logger.error(f"LAI计算失败|LAI calculation failed: {e}")
            import traceback
            self.logger.debug(traceback.format_exc())
            return None

    def _merge_candidates(self, harvest_file: Path, finder_file: Path) -> Path:
        """
        合并LTR候选文件|Merge LTR candidate files

        Args:
            harvest_file: LTRharvest输出文件|LTRharvest output file
            finder_file: LTR_FINDER输出文件|LTR_FINDER output file

        Returns:
            Path: 合并后的文件路径|Merged file path
        """
        try:
            prefix = self.config.genome_path.name
            merged_file = self.config.output_path / f"{prefix}.rawLTR.scn"

            self.logger.info("合并LTR候选文件|Merging LTR candidate files")

            # 使用cat命令合并文件|Use cat command to merge files
            with open(merged_file, 'w') as outfile:
                # 写入harvest结果|Write harvest results
                with open(harvest_file, 'r') as infile:
                    outfile.write(infile.read())

                # 写入finder结果|Write finder results
                with open(finder_file, 'r') as infile:
                    outfile.write(infile.read())

            self.logger.info(f"合并完成|Merge completed: {merged_file}")

            return merged_file

        except Exception as e:
            self.logger.error(f"合并文件失败|Failed to merge files: {e}")
            return None


def main():
    """主函数|Main function"""
    parser = argparse.ArgumentParser(
        description='LAI计算工具|LTR Assembly Index Calculator',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
示例|Examples:
  %(prog)s -g genome.fa -o output_dir
  %(prog)s -g genome.fa -o output_dir -m harvest
  %(prog)s -g genome.fa -o output_dir -m calculate
        """
    )

    # 必需参数|Required arguments
    parser.add_argument('-g', '--genome', required=True,
                       help='基因组FASTA文件|Genome FASTA file')

    parser.add_argument('-o', '--output', required=True,
                       help='输出目录|Output directory')

    # 运行参数|Runtime parameters
    parser.add_argument('-t', '--threads', type=int, default=64,
                       help='线程数|Number of threads')

    parser.add_argument('-m', '--mode', default='full',
                       choices=['full', 'harvest', 'retrieve', 'calculate'],
                       help='运行模式|Run mode (default: full)')

    # Conda环境路径|Conda environment path
    parser.add_argument('--conda-env',
                       default='~/miniforge3/envs/ltr_retriever_v.3.0.1',
                       help='Conda环境路径|Conda environment path')

    args = parser.parse_args()

    # 创建计算器并运行|Create calculator and run
    calculator = LAICalculator(
        genome=args.genome,
        output_dir=args.output,
        threads=args.threads,
        mode=args.mode,
        conda_env=args.conda_env
    )

    calculator.run()


if __name__ == "__main__":
    main()
