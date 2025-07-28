# biopytools/core/runner.py

import subprocess
from pathlib import Path

class CommandRunner:
    def __init__(self, logger, working_dir: Path, param_db: dict = None):
        self.logger = logger
        self.working_dir = Path(working_dir).resolve()
        self.param_db = param_db if param_db is not None else {}

    def _explain_command(self, cmd: list):
        if not cmd or cmd[0] not in self.param_db: return

        command_name = cmd[0]
        explanations = self.param_db[command_name]
        
        self.logger.icon_info("param", "=" * 60)
        self.logger.icon_info("param", f"命令参数详解 | Command Parameter Explanation for '{command_name}'")
        self.logger.icon_info("param", "=" * 60)

        i = 0
        while i < len(cmd):
            param = cmd[i]
            if param in explanations:
                desc = explanations[param]
                value_str = ""
                # 检查参数后面是否有值
                if i + 1 < len(cmd) and not str(cmd[i+1]).startswith('-'):
                    value_str = f" {cmd[i+1]}"
                    i += 1
                self.logger.icon_info("info", f"  {param}{value_str}")
                self.logger.icon_info("info", f"    💡 {desc}")
            i += 1
        self.logger.icon_info("param", "=" * 60)

    def run(self, cmd, description: str, check: bool = True):
        self.logger.start(description)
        
        cmd_str = ' '.join(str(c) for c in cmd)
        self.logger.icon_info("run", f"命令 | Command: {cmd_str}")
        self.logger.icon_info("data", f"工作目录 | Working directory: {self.working_dir}")

        if self.param_db: self._explain_command(cmd)

        try:
            result = subprocess.run(cmd, capture_output=True, text=True, check=check,
                                    cwd=str(self.working_dir), encoding='utf-8')
            self.logger.success(f"完成: {description}")
            if result.stdout: self.logger.debug(f"标准输出 | Stdout:\n{result.stdout}")
            if result.stderr: self.logger.debug(f"标准错误 | Stderr:\n{result.stderr}")
            return result
        except FileNotFoundError:
            self.logger.error(f"命令未找到: '{cmd[0]}'. 请确保它已安装并配置在系统PATH中。")
            if check: raise
            return None
        except subprocess.CalledProcessError as e:
            self.logger.error(f"命令执行失败: {description}")
            self.logger.error(f"错误信息 | Error output:\n{e.stderr}")
            if check: raise
            return e