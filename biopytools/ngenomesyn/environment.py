"""
🔍 环境检查模块 | Environment Check Module
"""

import subprocess
import logging
from typing import Dict

class EnvironmentChecker:
    """🔍 环境检查器 | Environment Checker"""
    
    def __init__(self, logger):
        self.logger = logger
        self.required_tools = ["NGenomeSyn"]
        self.optional_tools = {
            "minimap2": "🧬 Minimap2比对器 | Minimap2 aligner",
            "nucmer": "🧮 MUMmer nucmer比对器 | MUMmer nucmer aligner",
            "show-coords": "📊 MUMmer show-coords工具 | MUMmer show-coords tool",
            "syri": "🔬 SyRI结构变异检测 | SyRI structural variation detector", 
            "convert": "🎨 ImageMagick转换工具 | ImageMagick convert tool"
        }
        self.command_log = []  # 存储执行的命令
    
    def log_command(self, command: str, description: str = ""):
        """记录执行的命令到日志 | Log executed command"""
        self.logger.info("🔍" + "=" * 60)
        self.logger.info(f"🔍 环境检查命令 | Environment check command: {description}")
        self.logger.info(f"💻 命令内容 | Command: {command}")
        self.logger.info("🔍" + "=" * 60)
        
        # 保存到命令历史
        self.command_log.append({
            "description": description,
            "command": command
        })
    
    def save_environment_check_script(self, output_dir: str = "./"):
        """保存环境检查脚本 | Save environment check script"""
        import os
        script_path = os.path.join(output_dir, "check_environment.sh")
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 🔍 环境检查脚本 | Environment check script\n")
            f.write(f"# ⏰ 生成时间 | Generated at: $(date)\n\n")
            
            f.write("echo \"🔍 开始环境检查 | Starting environment check\"\n\n")
            
            # 检查必需工具
            f.write("# 🎯 检查必需工具 | Check required tools\n")
            f.write("echo \"🎯 检查必需工具 | Checking required tools\"\n")
            f.write("MISSING_TOOLS=()\n\n")
            
            for tool in self.required_tools:
                f.write(f"if command -v {tool} >/dev/null 2>&1; then\n")
                f.write(f"    echo \"✅ {tool} 已安装\"\n")
                f.write("else\n")
                f.write(f"    echo \"❌ {tool} 未安装\"\n")
                f.write(f"    MISSING_TOOLS+=(\"{tool}\")\n")
                f.write("fi\n\n")
            
            # 检查可选工具
            f.write("# 🛠️ 检查可选工具 | Check optional tools\n")
            f.write("echo \"🛠️ 检查可选工具 | Checking optional tools\"\n")
            for tool, description in self.optional_tools.items():
                f.write(f"if command -v {tool} >/dev/null 2>&1; then\n")
                f.write(f"    echo \"✅ {tool} 已安装 - {description}\"\n")
                f.write("else\n")
                f.write(f"    echo \"⚠️ {tool} 未安装 - {description}\"\n")
                f.write("fi\n\n")
            
            # 检查结果
            f.write("# 📊 检查结果 | Check results\n")
            f.write("if [ ${#MISSING_TOOLS[@]} -eq 0 ]; then\n")
            f.write("    echo \"✅ 所有必需工具都已安装，可以开始分析\"\n")
            f.write("    exit 0\n")
            f.write("else\n")
            f.write("    echo \"❌ 缺少以下必需工具: ${MISSING_TOOLS[*]}\"\n")
            f.write("    echo \"📋 请安装缺少的工具后重试\"\n")
            f.write("    exit 1\n")
            f.write("fi\n")
        
        # 设置可执行权限
        os.chmod(script_path, 0o755)
        
        self.logger.info(f"📝 环境检查脚本已保存 | Environment check script saved: {script_path}")
        self.logger.info("🚀 您可以使用以下命令检查环境:")
        self.logger.info(f"bash {script_path}")
        
        return script_path
    
    def check_tool_installation(self, tool_name: str) -> bool:
        """🔍 检查工具是否安装 | Check if tool is installed"""
        try:
            command = f"which {tool_name}"
            self.log_command(command, f"检查 {tool_name} 是否安装")

            result = subprocess.run(
                ["which", tool_name],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                tool_path = result.stdout.strip()
                self.logger.info(f"✅ {tool_name} 已安装 | installed: {tool_path}")
                return True
            else:
                self.logger.warning(f"⚠️ {tool_name} 未安装 | not installed")
                return False
        except (subprocess.TimeoutExpired, FileNotFoundError):
            self.logger.warning(f"❌ {tool_name} 检查失败 | check failed")
            return False

    def check_syri(self) -> bool:
        """🔬 检查SyRI安装 | Check SyRI installation

        SyRI可以以以下方式安装:
        1. 作为'python3 syri'命令（从PATH）
        2. 作为'syri'命令（可执行文件）
        """
        # 方法1: 检查'syri'是否在PATH中
        try:
            command = "which syri"
            self.log_command(command, "检查syri是否在PATH中")

            result = subprocess.run(
                ["which", "syri"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                tool_path = result.stdout.strip()
                self.logger.info(f"✅ SyRI 已安装 | installed: {tool_path}")
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

        # 方法2: 尝试运行'python3 syri --help'
        try:
            command = "python3 syri --help"
            self.log_command(command, "测试python3 syri命令")

            result = subprocess.run(
                ["python3", "syri", "--help"],
                capture_output=True,
                text=True,
                timeout=10
            )
            # SyRI的--help会返回0，即使找不到也会返回错误
            # 只要能调用syri模块就认为安装了
            if "syri" in result.stdout.lower() or "usage" in result.stdout.lower():
                self.logger.info("✅ SyRI 已安装 | python3 syri available")
                return True
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass

        # 方法3: 检查常见安装路径
        import os
        possible_paths = [
            os.path.expanduser("~/bin/syri"),
            "/usr/local/bin/syri",
            os.path.expanduser("~/syri/syri"),
        ]

        for path in possible_paths:
            if os.path.exists(path):
                self.logger.info(f"✅ SyRI 已安装 | found at: {path}")
                return True

        self.logger.warning("⚠️ SyRI 未安装 | SyRI not installed")
        return False
    
    def check_ngenomesyn(self) -> bool:
        """🎯 检查NGenomeSyn安装 | Check NGenomeSyn installation"""
        # 🔍 首先用which命令检查是否在PATH中
        try:
            command = "which NGenomeSyn"
            self.log_command(command, "检查NGenomeSyn是否在PATH中")
            
            result = subprocess.run(
                ["which", "NGenomeSyn"],
                capture_output=True,
                text=True,
                timeout=10
            )
            if result.returncode == 0:
                tool_path = result.stdout.strip()
                self.logger.info(f"✅ NGenomeSyn 已安装 | installed: {tool_path}")
                
                # 进一步验证可以运行（不使用-h参数，因为可能不支持）
                try:
                    test_command = "NGenomeSyn"
                    self.log_command(test_command, "测试NGenomeSyn是否可执行")
                    
                    test_result = subprocess.run(
                        [tool_path],
                        capture_output=True,
                        text=True,
                        timeout=5
                    )
                    # NGenomeSyn运行后无论返回码如何，只要能执行就说明安装正确
                    self.logger.info(f"✅ NGenomeSyn 验证通过 | NGenomeSyn verification passed")
                    return True
                except (subprocess.TimeoutExpired, FileNotFoundError):
                    self.logger.warning(f"⚠️ NGenomeSyn 可执行但验证失败 | NGenomeSyn executable but verification failed")
                    return True  # 即使验证失败，只要which找到了就认为安装了
        except (subprocess.TimeoutExpired, FileNotFoundError):
            pass
        
        # 🔍 检查多个可能的相对路径
        possible_paths = [
            "./bin/NGenomeSyn", 
            "../bin/NGenomeSyn",
            "../../bin/NGenomeSyn"
        ]
        
        for path in possible_paths:
            try:
                # 直接运行不带参数，避免-h参数的问题
                command = f"{path}"
                self.log_command(command, f"测试相对路径中的NGenomeSyn: {path}")
                
                result = subprocess.run(
                    [path],
                    capture_output=True,
                    text=True,
                    timeout=5
                )
                # 只要能找到可执行文件并运行（无论返回码），就认为安装了
                self.logger.info(f"✅ NGenomeSyn 已安装 | installed: {path}")
                return True
            except (subprocess.TimeoutExpired, FileNotFoundError):
                continue
        
        self.logger.error("❌ NGenomeSyn 未安装或不在PATH中 | not installed or not in PATH")
        return False
    
    def check_environment(self, required_aligner: str = "minimap2") -> Dict[str, bool]:
        """🔧 检查完整环境 | Check complete environment"""
        self.logger.info("🔍 开始环境检查 | Starting environment check")
        self.logger.info("=" * 50)

        results = {}

        # 🎯 检查NGenomeSyn (必需)
        results["NGenomeSyn"] = self.check_ngenomesyn()

        # 🔗 检查必需的比对器
        results[required_aligner] = self.check_tool_installation(required_aligner)

        # 🛠️ 检查可选工具
        for tool, description in self.optional_tools.items():
            if tool != required_aligner:  # 避免重复检查
                # 使用专门的SyRI检查方法
                if tool == "syri":
                    results[tool] = self.check_syri()
                else:
                    results[tool] = self.check_tool_installation(tool)

        self.logger.info("=" * 50)

        # ⚠️ 检查关键工具
        critical_missing = []
        if not results["NGenomeSyn"]:
            critical_missing.append("NGenomeSyn")
        if not results.get(required_aligner, False):
            critical_missing.append(required_aligner)

        if critical_missing:
            self.logger.error(f"❌ 关键工具缺失 | Critical tools missing: {', '.join(critical_missing)}")
        else:
            self.logger.info("✅ 环境检查通过 | Environment check passed")

        # 生成环境检查脚本
        try:
            self.save_environment_check_script()
        except Exception as e:
            self.logger.warning(f"⚠️ 保存环境检查脚本失败 | Failed to save environment check script: {e}")

        return results