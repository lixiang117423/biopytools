"""
EDTA转座子注释运行器|EDTA TE Annotation Runner
"""

import os
import sys
import subprocess
from pathlib import Path
from typing import Optional

from .config import EDTAConfig
from .utils import EDTALogger, CommandRunner, check_edta_dependencies


class EDTARunner:
    """EDTA转座子注释运行器|EDTA TE Annotation Runner"""

    def __init__(self, config: EDTAConfig):
        """初始化EDTA运行器|Initialize EDTA Runner

        Args:
            config: EDTA配置对象|EDTA configuration object
        """
        self.config = config
        self.logger_manager = EDTALogger(config.output_path)
        self.logger = self.logger_manager.get_logger()
        self.cmd_runner = CommandRunner(self.logger, config.output_path)

        # 检查EDTA依赖|Check EDTA dependencies
        try:
            check_edta_dependencies(config, self.logger)
        except RuntimeError as e:
            self.logger.error(f"依赖检查失败|Dependency check failed: {e}")
            raise

    def build_command(self) -> list:
        """构建EDTA命令|Build EDTA command

        Returns:
            EDTA命令列表|EDTA command list
        """
        cmd = []

        # 查找EDTA.pl|Find EDTA.pl
        edta_pl = self._find_edta_pl()

        if not edta_pl or not os.path.exists(edta_pl):
            raise RuntimeError(f"未找到EDTA.pl|EDTA.pl not found at: {edta_pl}")

        # 使用perl执行EDTA.pl|Use perl to execute EDTA.pl
        cmd = ["perl", edta_pl]

        # 必需参数|Required parameters
        cmd.extend(["--genome", self.config.genome])

        # 可选参数|Optional parameters
        cmd.extend(["--species", self.config.species])
        cmd.extend(["--step", self.config.step])
        cmd.extend(["--overwrite", str(self.config.overwrite)])
        cmd.extend(["--maxdiv", str(self.config.maxdiv)])
        cmd.extend(["--u", str(self.config.u)])
        cmd.extend(["--threads", str(self.config.threads)])
        cmd.extend(["--debug", str(self.config.debug)])

        # 布尔参数|Boolean parameters
        if self.config.sensitive:
            cmd.extend(["--sensitive", "1"])
        if self.config.anno:
            cmd.extend(["--anno", "1"])
        if self.config.evaluate:
            cmd.extend(["--evaluate", "1"])
        if self.config.force:
            cmd.extend(["--force", "1"])

        # 文件参数|File parameters
        if self.config.cds:
            cmd.extend(["--cds", self.config.cds])
        if self.config.curatedlib:
            cmd.extend(["--curatedlib", self.config.curatedlib])
        if self.config.rmlib:
            cmd.extend(["--rmlib", self.config.rmlib])
        if self.config.rmout:
            cmd.extend(["--rmout", self.config.rmout])
        if self.config.exclude:
            cmd.extend(["--exclude", self.config.exclude])

        # 路径参数|Path parameters
        if self.config.repeatmasker_path:
            cmd.extend(["--repeatmasker", self.config.repeatmasker_path])
        if self.config.repeatmodeler_path:
            cmd.extend(["--repeatmodeler", self.config.repeatmodeler_path])
        if self.config.ltrretriever_path:
            cmd.extend(["--ltrretriever", self.config.ltrretriever_path])
        if self.config.annosine_path:
            cmd.extend(["--annosine", self.config.annosine_path])

        return cmd

    def _find_edta_pl(self) -> Optional[str]:
        """查找EDTA.pl脚本|Find EDTA.pl script

        Returns:
            EDTA.pl路径或None|EDTA.pl path or None
        """
        # 1. 检查配置中的edta_path|Check edta_path in config
        if self.config.edta_path:
            edta_pl = os.path.join(self.config.edta_path, "EDTA.pl")
            if os.path.exists(edta_pl):
                return edta_pl

        # 2. 检查conda环境|Check conda environment
        conda_prefix = os.environ.get('CONDA_PREFIX', '')
        if conda_prefix:
            edta_pl = os.path.join(conda_prefix, 'share', 'EDTA', 'EDTA.pl')
            if os.path.exists(edta_pl):
                return edta_pl

        # 3. 检查系统PATH|Check system PATH
        try:
            result = subprocess.run(
                ["which", "EDTA.pl"],
                capture_output=True,
                text=True,
                timeout=5
            )
            if result.returncode == 0 and result.stdout.strip():
                return result.stdout.strip()
        except Exception:
            pass

        # 4. 使用默认路径|Use default path
        default_path = '~/miniforge3/envs/EDTA_v.2.2.2/share/EDTA/EDTA.pl'
        if os.path.exists(default_path):
            return default_path

        return None

    def run(self) -> bool:
        """运行EDTA分析|Run EDTA analysis

        Returns:
            是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始EDTA转座子注释分析|Starting EDTA TE annotation analysis")
        self.logger.info("=" * 60)

        try:
            # 构建命令|Build command
            cmd = self.build_command()

            # 打印参数信息|Print parameter information
            self._log_parameters()

            # 执行EDTA|Execute EDTA
            self.logger.info("开始执行EDTA.pl|Starting EDTA.pl execution")
            try:
                result = self.cmd_runner.run(
                    cmd,
                    description="EDTA转座子注释|EDTA TE annotation",
                    check=True,
                    show_progress=True
                )
            except subprocess.CalledProcessError as e:
                # 检查是否是TE缺失相关的错误（SINE, TIR, LINE, Helitron, LTR）
                error_output = e.stdout if hasattr(e, 'stdout') else str(e)

                te_types = ['SINE', 'TIR', 'LINE', 'Helitron', 'LTR']
                missing_te = None

                for te_type in te_types:
                    if te_type in error_output and 'not found' in error_output:
                        missing_te = te_type
                        break

                if missing_te:
                    self.logger.warning(f"检测到{missing_te}查找失败|{missing_te} detection failed")
                    self.logger.info("尝试创建空文件并使用--force 1继续|Attempting to create empty file and continue with --force 1")

                    # 对于SINE和LINE，可以创建空文件；对于TIR和Helitron，需要使用--force
                    if missing_te in ['SINE', 'LINE']:
                        if not self._create_empty_te_file(missing_te):
                            self.logger.error(f"无法创建空{missing_te}文件|Failed to create empty {missing_te} file")
                            return False

                    # 使用filter步骤和force参数继续运行
                    return self._run_filter_step()
                else:
                    # 其他错误，直接返回失败
                    raise

            # 检查结果|Check results
            if result.returncode == 0:
                self.logger.info("EDTA分析完成|EDTA analysis completed")
                self._check_outputs()
                return True
            else:
                self.logger.error(f"EDTA分析失败|EDTA analysis failed with return code: {result.returncode}")
                return False

        except Exception as e:
            self.logger.error(f"EDTA运行出错|EDTA runtime error: {str(e)}")
            return False

    def _create_empty_te_file(self, te_type: str) -> bool:
        """创建空的TE文件|Create empty TE file

        Args:
            te_type: TE类型 (SINE, LINE等)

        Returns:
            是否成功|Whether successful
        """
        try:
            genome_name = Path(self.config.genome).stem
            raw_dir = Path(self.config.genome).parent / f"{genome_name}.mod.EDTA.raw"
            te_file = raw_dir / f"{genome_name}.mod.{te_type}.raw.fa"

            # 创建目录
            raw_dir.mkdir(parents=True, exist_ok=True)

            # 创建空的FASTA文件
            with open(te_file, 'w') as f:
                f.write("")

            self.logger.info(f"已创建空{te_type}文件|Created empty {te_type} file: {te_file}")
            return True
        except Exception as e:
            self.logger.error(f"创建空{te_type}文件失败|Failed to create empty {te_type} file: {e}")
            return False

    def _run_filter_step(self) -> bool:
        """运行filter步骤继续EDTA分析，使用--force参数|Run filter step with --force parameter"""
        try:
            # 修改配置为filter步骤，并启用force
            original_step = self.config.step
            original_force = self.config.force
            self.config.step = 'filter'
            self.config.force = 1  # 启用force参数以自动填充缺失的TE文件

            # 重新构建命令
            cmd = self.build_command()

            self.logger.info("继续执行EDTA filter步骤（使用--force 1）|Continuing with EDTA filter step (--force 1)")
            self.logger.info("--force参数将自动填充缺失的TIR/LINE/Helitron文件|--force will auto-fill missing TIR/LINE/Helitron files")

            result = self.cmd_runner.run(
                cmd,
                description="EDTA filter步骤|EDTA filter step",
                check=True,
                show_progress=True
            )

            # 恢复原始配置
            self.config.step = original_step
            self.config.force = original_force

            # 检查结果
            if result.returncode == 0:
                self.logger.info("EDTA filter步骤完成|EDTA filter step completed")
                self._check_outputs()
                return True
            else:
                self.logger.error(f"EDTA filter步骤失败|EDTA filter step failed with return code: {result.returncode}")
                return False

        except Exception as e:
            self.logger.error(f"EDTA filter步骤运行出错|EDTA filter step runtime error: {str(e)}")
            return False

    def _log_parameters(self):
        """记录参数信息|Log parameter information"""
        self.logger.info("EDTA运行参数|EDTA running parameters:")
        self.logger.info(f"  基因组文件|Genome: {self.config.genome}")
        self.logger.info(f"  物种|Species: {self.config.species}")
        self.logger.info(f"  步骤|Step: {self.config.step}")
        self.logger.info(f"  覆盖|Overwrite: {self.config.overwrite}")
        self.logger.info(f"  线程数|Threads: {self.config.threads}")
        self.logger.info(f"  敏感模式|Sensitive: {self.config.sensitive}")
        self.logger.info(f"  注释|Annotation: {self.config.anno}")
        self.logger.info(f"  评估|Evaluate: {self.config.evaluate}")
        self.logger.info(f"  最大分歧度|Max divergence: {self.config.maxdiv}")
        self.logger.info(f"  输出目录|Output directory: {self.config.output_dir}")

        if self.config.cds:
            self.logger.info(f"  CDS文件|CDS: {self.config.cds}")
        if self.config.curatedlib:
            self.logger.info(f"  筛选库|Curated library: {self.config.curatedlib}")
        if self.config.exclude:
            self.logger.info(f"  排除区域|Exclude regions: {self.config.exclude}")

    def _check_outputs(self):
        """检查输出文件|Check output files"""
        self.logger.info("检查EDTA输出文件|Checking EDTA output files")

        genome_name = Path(self.config.genome).stem
        expected_files = [
            f"{genome_name}.mod.EDTA.TElib.fa",
            f"{genome_name}.mod.EDTA.TElib.fa.mod.tsv"
        ]

        if self.config.anno:
            expected_files.extend([
                f"{genome_name}.mod.EDTA.TEanno.gff3",
                f"{genome_name}.mod.EDTA.TEanno.sum"
            ])

        found_files = []
        missing_files = []

        for filename in expected_files:
            filepath = Path(self.config.genome).parent / filename
            if filepath.exists():
                found_files.append(filename)
                self.logger.info(f"  找到|Found: {filename}")
            else:
                missing_files.append(filename)
                self.logger.warning(f"  未找到|Not found: {filename}")

        if found_files:
            self.logger.info(f"共找到{len(found_files)}个输出文件|Found {len(found_files)} output files")

        if missing_files:
            self.logger.warning(f"缺失{len(missing_files)}个输出文件|Missing {len(missing_files)} output files")
