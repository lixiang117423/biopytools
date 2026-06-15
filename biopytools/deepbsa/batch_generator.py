"""
DeepBSA批量命令生成器|DeepBSA Batch Command Generator
生成批量处理的任务脚本|Generate batch processing task scripts
"""

import logging
from pathlib import Path
from typing import List


class BatchGenerator:
    """DeepBSA批量命令生成器|DeepBSA Batch Command Generator"""

    def __init__(self, config, logger: logging.Logger):
        """初始化生成器|Initialize generator

        Args:
            config: BatchConfig配置对象|BatchConfig configuration object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger

    def generate(self) -> bool:
        """生成批量任务脚本|Generate batch task script

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            self.logger.info("=" * 60)
            self.logger.info("开始生成批量任务脚本|Starting batch script generation")
            self.logger.info("=" * 60)
            self.logger.info(f"输入文件|Input file: {self.config.input_path}")
            self.logger.info(f"输出目录|Output directory: {self.config.output_path}")
            self.logger.info(f"脚本文件|Script file: {self.config.script_name}")
            self.logger.info(f"运行方法|Methods: {', '.join(self.config.methods_list)}")
            self.logger.info(f"线程数|Threads: {self.config.threads}")
            self.logger.info("=" * 60)
            self.logger.info("")

            # 生成命令|Generate commands
            self.logger.info("生成运行脚本|Generate run scripts")
            commands = self._generate_commands()
            self._write_run_script(commands)

            self.logger.info("")
            self.logger.info("=" * 60)
            self.logger.info("批量任务脚本生成完成|Batch script generation completed")
            self.logger.info("=" * 60)
            self._show_usage_guide()

            return True

        except Exception as e:
            self.logger.error(f"生成失败|Generation failed: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return False

    def _generate_commands(self) -> List[str]:
        """生成任务命令|Generate task commands

        Returns:
            list: 命令列表|Command list
        """
        commands = []

        for method in self.config.methods_list:
            # 每个方法一个独立的输出目录|Each method has independent output directory
            output_dir = self.config.output_path / method

            # 构建run命令|Build run command
            cmd = (
                f"biopytools deepbsa run "
                f"-i {self.config.input_path} "
                f"-o {output_dir} "
                f"--methods {method} "
                f"--threads {self.config.threads} "
                f"--smooth-func {self.config.smooth_func} "
                f"--no-parallel "  # 强制串行模式，因为只运行一个方法|Force serial mode, only running one method
                f"--skip-merge"  # 跳过合并，最后统一合并|Skip merge, merge later
            )

            commands.append(cmd)
            self.logger.debug(f"生成命令|Generated command: {method}")

        return commands

    def _write_run_script(self, commands: List[str]):
        """写入运行脚本|Write run script

        Args:
            commands: 命令列表|Command list
        """
        self.logger.info(f"写入运行脚本|Writing run script: {self.config.script_path}")
        with open(self.config.script_path, 'w') as f:
            # 添加shebang|Add shebang
            f.write("#!/bin/bash\n")
            f.write("# DeepBSA批量任务脚本|DeepBSA Batch Task Script\n")
            f.write(f"# 输入文件|Input file: {self.config.input_path}\n")
            f.write(f"# 方法数量|Method count: {len(commands)}\n")
            f.write(f"# 线程数|Threads: {self.config.threads}\n")
            f.write("\n")

            # 写入每个命令|Write each command
            for i, cmd in enumerate(commands, 1):
                f.write(f"# 方法|Method {i}: {self.config.methods_list[i-1]}\n")
                f.write(cmd + "\n")

        # 设置可执行权限|Set executable permission
        self.config.script_path.chmod(0o755)
        self.logger.info(f"运行脚本已保存|Run script saved: {self.config.script_path}")

    def _show_usage_guide(self):
        """显示使用指南|Show usage guide"""
        self.logger.info("")
        self.logger.info("=" * 60)
        self.logger.info("执行方式|Execution methods:")
        self.logger.info("=" * 60)
        self.logger.info("步骤1|Step 1: 运行各方法|Run methods")
        self.logger.info(f"  使用您的批量投递脚本逐行执行 {self.config.script_name}")
        self.logger.info(f"  Use your batch submission script to execute {self.config.script_name} line by line")
        self.logger.info("")
        self.logger.info("步骤2|Step 2: 合并结果|Merge results (所有方法完成后执行|After all methods complete)")
        self.logger.info(f"  biopytools deepbsa merge -i {self.config.output_path} -o {self.config.output_path}/merged_results")
        self.logger.info("")
        self.logger.info("输出目录结构|Output directory structure:")
        self.logger.info(f"  {self.config.output_path}/")
        for method in self.config.methods_list:
            self.logger.info(f"  ├── {method}/          # 方法{method}的结果|Results for method {method}")
        self.logger.info(f"  └── merged_results/   # 合并后的结果|Merged results")
        self.logger.info("=" * 60)
