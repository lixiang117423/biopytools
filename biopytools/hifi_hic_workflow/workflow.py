"""
HiFi+Hi-C工作流主流程模块|HiFi+Hi-C Workflow Main Module

整合4个子模块，实现完整的基因组组装流程
Integrate 4 sub-modules for complete genome assembly workflow
"""

import os
import sys
import logging
from datetime import datetime
from typing import Optional
from pathlib import Path


class HifiHicWorkflow:
    """HiFi+Hi-C工作流主类|HiFi+Hi-C Workflow Main Class"""

    def __init__(self, config, logger: logging.Logger):
        """
        初始化工作流|Initialize workflow

        Args:
            config: HifiHicWorkflowConfig配置对象|HifiHicWorkflowConfig object
            logger: 日志记录器|Logger
        """
        self.config = config
        self.logger = logger

        # 导入数据传递管理器|Import data transfer manager
        from .data_transfer import DataTransferManager
        self.data_transfer = DataTransferManager(config, logger)

    def run_workflow(self) -> bool:
        """
        运行完整工作流|Run complete workflow

        Returns:
            bool: 是否成功|Whether successful
        """
        start_time = datetime.now()

        try:
            self.logger.info("=" * 80)
            self.logger.info("开始HiFi+Hi-C基因组组装与挂载流程|Starting HiFi+Hi-C Genome Assembly and Scaffolding Workflow")
            self.logger.info("=" * 80)

            # 显示配置信息|Display configuration
            self._print_config()

            # Step 1: HiFi组装|Step 1: HiFi Assembly
            hifi_output = None
            if not self.config.skip_hifi_hic:
                hifi_output = self._run_hifi_hic()
                if not hifi_output:
                    self.logger.error("HiFi组装失败|HiFi assembly failed")
                    return False
            else:
                self.logger.info("跳过HiFi组装|Skipping HiFi assembly")
                hifi_output = self.data_transfer.get_hifi_hic_output()
                if not hifi_output:
                    self.logger.error("未找到HiFi组装输出|HiFi assembly output not found")
                    return False

            # Step 2: HapHiC挂载|Step 2: HapHiC Scaffolding
            haphic_output = None
            if not self.config.skip_haphic:
                haphic_output = self._run_haphic(hifi_output)
                if not haphic_output:
                    self.logger.error("HapHiC挂载失败|HapHiC scaffolding failed")
                    return False
            else:
                self.logger.info("跳过HapHiC挂载|Skipping HapHiC scaffolding")
                haphic_output = self.data_transfer.get_haphic_output()
                if not haphic_output:
                    self.logger.error("未找到HapHiC输出|HapHiC output not found")
                    return False

            # Step 3: 染色体重命名|Step 3: Chromosome Rename
            renamed_output = None
            if not self.config.skip_rename:
                renamed_output = self._run_rename_chromosomes(haphic_output)
                if not renamed_output:
                    self.logger.error("染色体重命名失败|Chromosome rename failed")
                    return False
            else:
                self.logger.info("跳过染色体重命名|Skipping chromosome rename")
                renamed_output = self.data_transfer.get_rename_output()
                if not renamed_output or not os.path.exists(renamed_output):
                    self.logger.error("未找到重命名输出|Renamed output not found")
                    return False

            # Step 4: Hi-C热图|Step 4: Hi-C Heatmap
            if not self.config.skip_heatmap:
                heatmap_output = self._run_hic_heatmap(renamed_output)
                if not heatmap_output:
                    self.logger.error("Hi-C热图生成失败|Hi-C heatmap generation failed")
                    return False
            else:
                self.logger.info("跳过Hi-C热图生成|Skipping Hi-C heatmap generation")

            # 生成流程报告|Generate workflow report
            end_time = datetime.now()
            self._generate_workflow_report(start_time, end_time)

            self.logger.info("=" * 80)
            self.logger.info("HiFi+Hi-C工作流完成|HiFi+Hi-C workflow completed successfully")
            self.logger.info("=" * 80)

            return True

        except Exception as e:
            self.logger.error(f"工作流执行异常|Workflow execution error: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _run_hifi_hic(self) -> Optional[str]:
        """
        运行HiFi组装步骤|Run HiFi assembly step

        Returns:
            str: 输出FASTA文件路径|Output FASTA file path, or None if failed
        """
        try:
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤1: HiFi组装|Step 1: HiFi Assembly")
            self.logger.info("=" * 60)

            # 检查是否已完成|Check if already completed
            if self.config.resume and self.data_transfer.check_step_completion("hifi_hic"):
                self.logger.info("HiFi组装已完成，跳过|HiFi assembly already completed, skipping")
                return self.data_transfer.get_hifi_hic_output()

            # 获取参数|Get parameters
            params = self.data_transfer.get_hifi_hic_params()

            # 导入并运行hifi_hic模块|Import and run hifi_hic module
            from biopytools.hifi_hic.main import GenomeAssembler

            assembler = GenomeAssembler(**params)
            assembler.run_assembly()

            # 检查输出|Check output
            output_file = self.data_transfer.get_hifi_hic_output()
            if not output_file:
                self.logger.error("未找到HiFi组装输出|HiFi assembly output not found")
                return None

            # 保存步骤信息|Save step info
            self.data_transfer.create_step_info_file("hifi_hic", {
                "output_file": output_file,
                "use_ngs_polish": self.config.use_ngs_polish
            })

            self.logger.info(f"步骤1完成|Step 1 completed: {output_file}")
            return output_file

        except Exception as e:
            self.logger.error(f"HiFi组装失败|HiFi assembly failed: {e}")
            return None

    def _run_haphic(self, asm_file: str) -> Optional[str]:
        """
        运行HapHiC挂载步骤|Run HapHiC scaffolding step

        Args:
            asm_file: 组装文件路径|Assembly file path

        Returns:
            str: 输出FASTA文件路径|Output FASTA file path, or None if failed
        """
        try:
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤2: HapHiC挂载|Step 2: HapHiC Scaffolding")
            self.logger.info("=" * 60)

            # 检查是否已完成|Check if already completed
            if self.config.resume and self.data_transfer.check_step_completion("haphic"):
                self.logger.info("HapHiC挂载已完成，跳过|HapHiC scaffolding already completed, skipping")
                return self.data_transfer.get_haphic_output()

            # 获取参数|Get parameters
            params = self.data_transfer.get_haphic_params(asm_file)

            # 导入并运行haphic模块|Import and run haphic module
            from biopytools.haphic.main import HapHiCProcessor

            assembler = HapHiCProcessor(**params)
            success = assembler.run_pipeline()

            if not success:
                self.logger.error("HapHiC流程执行失败|HapHiC pipeline failed")
                return None

            # 检查输出|Check output
            output_file = self.data_transfer.get_haphic_output()
            if not output_file:
                self.logger.error("未找到HapHiC输出|HapHiC output not found")
                return None

            # 保存步骤信息|Save step info
            self.data_transfer.create_step_info_file("haphic", {
                "input_file": asm_file,
                "output_file": output_file,
                "nchrs": self.config.nchrs
            })

            self.logger.info(f"步骤2完成|Step 2 completed: {output_file}")
            return output_file

        except Exception as e:
            self.logger.error(f"HapHiC挂载失败|HapHiC scaffolding failed: {e}")
            return None

    def _run_rename_chromosomes(self, input_fa: str) -> Optional[str]:
        """
        运行染色体重命名步骤|Run chromosome rename step

        Args:
            input_fa: 输入FASTA文件|Input FASTA file

        Returns:
            str: 输出FASTA文件路径|Output FASTA file path, or None if failed
        """
        try:
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤3: 染色体重命名|Step 3: Chromosome Rename")
            self.logger.info("=" * 60)

            # 检查是否已完成|Check if already completed
            if self.config.resume and self.data_transfer.check_step_completion("rename"):
                self.logger.info("染色体重命名已完成，跳过|Chromosome rename already completed, skipping")
                return self.data_transfer.get_rename_output()

            # 获取参数|Get parameters
            params = self.data_transfer.get_rename_params(input_fa)
            output_fa = params["output_file"]

            # 使用参考基因组引导命名|Use reference genome guided naming
            from .reference_namer import ReferenceGenomeNamer

            namer = ReferenceGenomeNamer(self.config, self.logger)
            success = namer.run_rename(input_fa, output_fa)

            if not success:
                self.logger.error("染色体重命名失败|Chromosome rename failed")
                return None

            # 检查输出|Check output
            if not os.path.exists(output_fa):
                self.logger.error(f"未找到重命名输出|Renamed output not found: {output_fa}")
                return None

            # 保存步骤信息|Save step info
            self.data_transfer.create_step_info_file("rename", {
                "input_file": input_fa,
                "output_file": output_fa,
                "reference_genome": self.config.reference_genome,
                "nchrs": self.config.nchrs
            })

            self.logger.info(f"步骤3完成|Step 3 completed: {output_fa}")
            return output_fa

        except Exception as e:
            self.logger.error(f"染色体重命名失败|Chromosome rename failed: {e}")
            return None

    def _run_hic_heatmap(self, genome_fa: str) -> Optional[str]:
        """
        运行Hi-C热图生成步骤|Run Hi-C heatmap generation step

        Args:
            genome_fa: 基因组FASTA文件|Genome FASTA file

        Returns:
            str: 热图文件路径|Heatmap file path, or None if failed
        """
        try:
            self.logger.info("\n" + "=" * 60)
            self.logger.info("步骤4: Hi-C热图生成|Step 4: Hi-C Heatmap Generation")
            self.logger.info("=" * 60)

            # 检查是否已完成|Check if already completed
            if self.config.resume and self.data_transfer.check_step_completion("heatmap"):
                self.logger.info("Hi-C热图已生成，跳过|Hi-C heatmap already generated, skipping")
                return self.data_transfer.get_heatmap_output()

            # 获取参数|Get parameters
            params = self.data_transfer.get_heatmap_params(genome_fa)

            # 导入并运行hic_heatmap模块|Import and run hic_heatmap module
            from biopytools.hic_heatmap.main import main as hic_heatmap_main

            # 临时修改sys.argv以传递参数|Temporarily modify sys.argv to pass parameters
            original_argv = sys.argv

            try:
                # 构建参数列表|Build argument list
                sys.argv = [
                    'hic_heatmap',
                    '-i', genome_fa,
                    '-g', self.config.prefix,
                    '-1', self.config.hic_r1,
                    '-2', self.config.hic_r2,
                    '-o', self.config.heatmap_output_dir,
                    '-t', str(self.config.threads),
                    '--max-memory', str(self.config.hicpro_max_memory_gb),
                    '--restriction-enzyme', self.config.hicpro_restriction_enzyme,
                    '--bin-sizes', self.config.hicpro_bin_sizes,
                    '--resolution', str(self.config.heatmap_resolution),
                    '--color-map', self.config.heatmap_color_map,
                    '--dpi', str(self.config.heatmap_dpi),
                    '--format', self.config.heatmap_format,
                    '--bar-max', str(self.config.heatmap_bar_max),
                ]

                if self.config.hicpro_path:
                    sys.argv.extend(['--hicpro-sif', self.config.hicpro_path])

                if self.config.plothic_path:
                    sys.argv.extend(['--plothic-path', self.config.plothic_path])

                if self.config.force_rerun:
                    sys.argv.append('--force')

                # 运行hic_heatmap|Run hic_heatmap
                exit_code = hic_heatmap_main()

                if exit_code != 0:
                    self.logger.error("Hi-C热图生成失败|Hi-C heatmap generation failed")
                    return None

            finally:
                # 恢复原始argv|Restore original argv
                sys.argv = original_argv

            # 检查输出|Check output
            output_file = self.data_transfer.get_heatmap_output()
            if not output_file or not os.path.exists(output_file):
                self.logger.warning(f"未找到热图输出|Heatmap output not found: {output_file}")
                # 尝试查找可能的输出文件|Try to find possible output files
                import glob
                pattern = os.path.join(self.config.heatmap_output_dir, "plot", f"*.{self.config.heatmap_format}")
                matching_files = glob.glob(pattern)
                if matching_files:
                    output_file = matching_files[0]
                    self.logger.info(f"找到热图文件|Found heatmap file: {output_file}")
                else:
                    output_file = None

            # 保存步骤信息|Save step info
            if output_file:
                self.data_transfer.create_step_info_file("heatmap", {
                    "genome_file": genome_fa,
                    "output_file": output_file,
                    "hic_r1": self.config.hic_r1,
                    "hic_r2": self.config.hic_r2
                })

                self.logger.info(f"步骤4完成|Step 4 completed: {output_file}")
                return output_file
            else:
                return None

        except Exception as e:
            self.logger.error(f"Hi-C热图生成失败|Hi-C heatmap generation failed: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None

    def _print_config(self):
        """显示配置信息|Display configuration"""
        self.logger.info("\n工作流配置|Workflow Configuration:")
        self.logger.info("-" * 60)

        summary = self.config.get_summary()

        self.logger.info(f"样本前缀|Sample prefix: {summary['prefix']}")
        self.logger.info(f"线程数|Threads: {summary['threads']}")
        self.logger.info(f"染色体数|Chromosomes: {summary['nchrs']}")
        self.logger.info(f"工作目录|Work directory: {summary['work_dir']}")

        self.logger.info(f"\n输入文件|Input files:")
        self.logger.info(f"  HiFi reads|HiFi reads: {summary['hifi_reads']}")
        self.logger.info(f"  Hi-C R1|Hi-C R1: {summary['hic_r1']}")
        self.logger.info(f"  Hi-C R2|Hi-C R2: {summary['hic_r2']}")
        self.logger.info(f"  参考基因组|Reference genome: {summary['reference_genome']}")
        if summary['ngs_data']:
            self.logger.info(f"  NGS数据|NGS data: {summary['ngs_data']}")

        self.logger.info(f"\n流程控制|Workflow control:")
        self.logger.info(f"  跳过HiFi组装|Skip HiFi assembly: {summary['skip_hifi_hic']}")
        self.logger.info(f"  跳过HapHiC|Skip HapHiC: {summary['skip_haphic']}")
        self.logger.info(f"  跳过重命名|Skip rename: {summary['skip_rename']}")
        self.logger.info(f"  跳过热图|Skip heatmap: {summary['skip_heatmap']}")
        self.logger.info(f"  使用NGS polish|Use NGS polish: {summary['use_ngs_polish']}")
        self.logger.info(f"  断点续传|Resume: {summary['resume']}")

        self.logger.info("-" * 60)

    def _generate_workflow_report(self, start_time: datetime, end_time: datetime):
        """
        生成工作流报告|Generate workflow report

        Args:
            start_time: 开始时间|Start time
            end_time: 结束时间|End time
        """
        try:
            report_file = os.path.join(self.config.work_dir, "workflow_report.txt")

            duration = (end_time - start_time).total_seconds()

            with open(report_file, 'w') as f:
                f.write("=" * 60 + "\n")
                f.write("HiFi+Hi-C工作流报告|HiFi+Hi-C Workflow Report\n")
                f.write("=" * 60 + "\n\n")

                f.write(f"开始时间|Start time: {start_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"结束时间|End time: {end_time.strftime('%Y-%m-%d %H:%M:%S')}\n")
                f.write(f"总耗时|Total duration: {duration:.1f}秒 ({duration/60:.1f}分钟)\n\n")

                f.write(f"工作目录|Work directory: {self.config.work_dir}\n\n")

                f.write("步骤完成状态|Step completion status:\n")
                f.write("-" * 60 + "\n")

                steps = [
                    ("hifi_hic", "HiFi组装|HiFi assembly"),
                    ("haphic", "HapHiC挂载|HapHiC scaffolding"),
                    ("rename", "染色体重命名|Chromosome rename"),
                    ("heatmap", "Hi-C热图|Hi-C heatmap")
                ]

                for step_key, step_name in steps:
                    completed = self.data_transfer.check_step_completion(step_key)
                    status = " 完成|completed" if completed else " 未完成|not completed"
                    f.write(f"  [{status}] {step_name}\n")

                f.write("\n")

                # 输出文件列表|Output file list
                f.write("主要输出文件|Main output files:\n")
                f.write("-" * 60 + "\n")

                if not self.config.skip_hifi_hic:
                    hifi_output = self.data_transfer.get_hifi_hic_output()
                    if hifi_output:
                        f.write(f"  HiFi组装|HiFi assembly: {hifi_output}\n")

                if not self.config.skip_haphic:
                    haphic_output = self.data_transfer.get_haphic_output()
                    if haphic_output:
                        f.write(f"  HapHiC挂载|HapHiC scaffolds: {haphic_output}\n")

                if not self.config.skip_rename:
                    rename_output = self.data_transfer.get_rename_output()
                    if rename_output and os.path.exists(rename_output):
                        f.write(f"  重命名基因组|Renamed genome: {rename_output}\n")

                if not self.config.skip_heatmap:
                    heatmap_output = self.data_transfer.get_heatmap_output()
                    if heatmap_output and os.path.exists(heatmap_output):
                        f.write(f"  Hi-C热图|Hi-C heatmap: {heatmap_output}\n")

            self.logger.info(f"工作流报告已保存|Workflow report saved: {report_file}")

        except Exception as e:
            self.logger.warning(f"生成工作流报告失败|Failed to generate workflow report: {e}")
