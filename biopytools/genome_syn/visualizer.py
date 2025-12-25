"""
可视化模块 | Visualizer Module
"""

import os
import subprocess
import shutil
from typing import Dict, List, Optional
from .data_processing import GenomeProcessor

class NGenomeSynVisualizer:
    """NGenomeSyn可视化器 | NGenomeSyn Visualizer"""
    
    def __init__(self, logger):
        self.logger = logger
        self.processor = GenomeProcessor(logger)
        self.command_log = []  # 存储执行的命令
    
    def log_command(self, command: str, working_dir: str = None, description: str = ""):
        """记录执行的命令到日志 | Log executed command"""
        if working_dir is None:
            working_dir = os.getcwd()
        
        self.logger.info("🎨" + "=" * 80)
        self.logger.info(f"🎨 执行命令 | Executing command: {description}")
        self.logger.info(f"📁 工作目录 | Working directory: {working_dir}")
        self.logger.info(f"💻 命令内容 | Command: {command}")
        self.logger.info("🎨" + "=" * 80)
        
        # 保存到命令历史
        self.command_log.append({
            "description": description,
            "command": command,
            "working_dir": working_dir
        })
    
    def save_command_script(self, output_dir: str, script_name: str = "visualization_commands.sh"):
        """保存所有执行过的命令到脚本文件 | Save all executed commands to script file"""
        script_path = os.path.join(output_dir, script_name)
        
        with open(script_path, 'w') as f:
            f.write("#!/bin/bash\n")
            f.write("# 🎨 可视化命令脚本 | Visualization command script\n")
            f.write(f"# ⏰ 生成时间 | Generated at: $(date)\n\n")
            
            for i, cmd_info in enumerate(self.command_log, 1):
                f.write(f"# 📋 步骤 {i}: {cmd_info['description']}\n")
                f.write(f"echo \"🔄 执行可视化步骤 {i}: {cmd_info['description']}\"\n")
                f.write(f"cd \"{cmd_info['working_dir']}\"\n")
                f.write(f"{cmd_info['command']}\n")
                f.write("if [ $? -ne 0 ]; then\n")
                f.write(f"    echo \"❌ 错误: 可视化步骤 {i} 执行失败\"\n")
                f.write("    exit 1\n")
                f.write("fi\n\n")
        
        # 设置可执行权限
        os.chmod(script_path, 0o755)
        
        self.logger.info(f"📝 可视化命令脚本已保存 | Visualization command script saved: {script_path}")
        self.logger.info("🚀 您可以使用以下命令手动执行可视化 | You can manually execute visualization with:")
        self.logger.info(f"bash {script_path}")
        
        return script_path
    
    def log_config_generation(self, config_file: str, output_dir: str):
        """记录配置文件生成步骤 | Log configuration file generation step"""
        config_content_display = f"cat \"{os.path.basename(config_file)}\""
        self.log_command(
            command=config_content_display,
            working_dir=output_dir,
            description="📋 查看生成的NGenomeSyn配置文件内容"
        )
    
    def generate_ngenomesyn_config(self, config_data: Dict, output_dir: str, selected_chromosomes: Optional[List[int]] = None) -> str:
        """生成NGenomeSyn配置文件 | Generate NGenomeSyn configuration file"""
        self.logger.info("📝 生成NGenomeSyn配置文件 | Generating NGenomeSyn configuration file")
        
        if selected_chromosomes:
            self.logger.info(f"🧬 使用染色体过滤 | Using chromosome filter: {selected_chromosomes}")
        
        config_lines = ["SetParaFor = global\n\n"]
        
        # 生成.len文件并添加基因组信息
        for genome in config_data["genomes"]:
            order = genome["order"]
            
            # 根据是否有染色体过滤决定文件名
            if selected_chromosomes:
                chr_suffix = "_chr" + "_".join(map(str, selected_chromosomes))
                len_file = os.path.join(output_dir, f"{genome['name']}{chr_suffix}.len")
            else:
                len_file = os.path.join(output_dir, f"{genome['name']}.len")
            
            # 生成.len文件，传递染色体过滤参数
            self.processor.generate_len_file(
                genome["file"], 
                len_file, 
                selected_chromosomes=selected_chromosomes
            )
            
            # 记录len文件生成命令（模拟）
            len_cmd = f"# 📏 生成基因组长度文件: python -c \"from data_processing import GenomeProcessor; GenomeProcessor().generate_len_file('{genome['file']}', '{len_file}', selected_chromosomes={selected_chromosomes})\""
            self.log_command(
                command=len_cmd,
                working_dir=output_dir,
                description=f"📏 生成基因组 {genome['name']} 的长度文件"
            )
            
            # 检查生成的len文件
            if os.path.exists(len_file) and os.path.getsize(len_file) > 0:
                config_lines.append(f"GenomeInfoFile{order}={os.path.basename(len_file)}\n")
                self.logger.info(f"✅ 生成基因组信息文件 | Generated genome info file: {len_file}")
            else:
                self.logger.error(f"❌ 基因组信息文件生成失败 | Failed to generate genome info file: {len_file}")
                raise RuntimeError(f"❌ 基因组信息文件生成失败 | Failed to generate genome info file")
        
        config_lines.append("\n")
        
        # 添加链接文件信息并验证文件存在
        missing_files = []
        chr_suffix = "_chr" + "_".join(map(str, selected_chromosomes)) if selected_chromosomes else ""
        
        for link in config_data["links"]:
            from_order = link["from"]
            to_order = link["to"]
            
            # 尝试查找link文件或paf文件
            link_file = f"{link['from_name']}_vs_{link['to_name']}{chr_suffix}.link"
            paf_file = f"{link['from_name']}_vs_{link['to_name']}{chr_suffix}.paf"
            
            link_path = os.path.join(output_dir, link_file)
            paf_path = os.path.join(output_dir, paf_file)
            
            if os.path.exists(link_path) and os.path.getsize(link_path) > 0:
                config_lines.append(f"LinkFileRef{from_order}VsRef{to_order}={link_file}\n")
                self.logger.info(f"🔗 使用链接文件 | Using link file: {link_file}")
            elif os.path.exists(paf_path) and os.path.getsize(paf_path) > 0:
                # 如果只有PAF文件，NGenomeSyn可能不支持，需要转换
                self.logger.warning(f"⚠️ 只找到PAF文件，NGenomeSyn可能需要LINK格式 | Only found PAF file, NGenomeSyn may need LINK format: {paf_file}")
                config_lines.append(f"LinkFileRef{from_order}VsRef{to_order}={paf_file}\n")
            else:
                missing_files.append(f"{link['from_name']}_vs_{link['to_name']}{chr_suffix}")
                self.logger.error(f"❌ 缺失比对文件 | Missing alignment file: {link_file} or {paf_file}")
        
        if missing_files:
            raise RuntimeError(f"❌ 缺失以下比对文件 | Missing alignment files: {', '.join(missing_files)}")
        
        config_lines.append("\n")
        
        # 添加可视化设置
        viz_config = config_data["visualization"]
        if viz_config.get("canvas_width"):
            config_lines.append(f"body={viz_config['canvas_width']}\n")
        if viz_config.get("canvas_height"):
            ratio = viz_config["canvas_height"] / viz_config.get("canvas_width", 1200)
            config_lines.append(f"CanvasHeightRatio={ratio:.2f}\n")
        
        config_lines.append("\n")
        
        # 🔧 修复：为每个基因组单独设置染色体名称显示参数
        for genome in config_data["genomes"]:
            order = genome["order"]
            config_lines.append(f"SetParaFor = Genome{order}\n")
            config_lines.append("ChrNameShow=1\n")           # 显示染色体名称
            config_lines.append("ChrNameColor=black\n")      # 染色体名称颜色
            config_lines.append("ChrNameShiftX=0\n")         # X轴偏移
            config_lines.append("ChrNameShiftY=15\n")        # Y轴偏移（稍微向下）
            config_lines.append("\n")
        
        # 如果使用了染色体过滤，添加注释
        if selected_chromosomes:
            config_lines.append(f"# 🧬 分析染色体 | Analyzed chromosomes: {', '.join(map(str, selected_chromosomes))}\n")
        
        # 保存配置文件
        config_filename = f"ngenomesyn{chr_suffix}.conf" if selected_chromosomes else "ngenomesyn.conf"
        ngenomesyn_config = os.path.join(output_dir, config_filename)
        
        with open(ngenomesyn_config, 'w') as f:
            f.writelines(config_lines)
        
        # 记录配置文件内容到日志
        self.log_config_generation(ngenomesyn_config, output_dir)
        
        # 记录生成的配置文件内容用于调试
        self.logger.info(f"📋 NGenomeSyn配置文件内容 | NGenomeSyn config file content:")
        with open(ngenomesyn_config, 'r') as f:
            content = f.read()
            self.logger.info(f"📄 配置内容 | Config content:\n{content}")
        
        self.logger.info(f"✅ NGenomeSyn配置文件已生成 | NGenomeSyn config generated: {ngenomesyn_config}")
        return ngenomesyn_config
    
    def run_ngenomesyn(self, config_file: str, output_dir: str, output_formats: List[str]) -> List[str]:
        """运行NGenomeSyn绘图 | Run NGenomeSyn visualization"""
        self.logger.info("🎨 开始NGenomeSyn可视化 | Starting NGenomeSyn visualization")
        
        output_files = []
        
        # 从配置文件名推断输出文件前缀
        config_basename = os.path.splitext(os.path.basename(config_file))[0]
        if config_basename != "ngenomesyn":
            output_prefix = f"genome_synteny_{config_basename.replace('ngenomesyn', '')}"
        else:
            output_prefix = "genome_synteny"
        
        base_output = os.path.join(output_dir, output_prefix)
        
        # 在运行前检查输入文件
        self.logger.info("🔍 检查输入文件状态 | Checking input file status")
        self.logger.info(f"📁 工作目录 | Working directory: {output_dir}")
        self.logger.info(f"📋 配置文件 | Config file: {config_file}")
        
        # 记录文件检查命令
        check_cmd = f"ls -la \"{output_dir}\""
        self.log_command(
            command=check_cmd,
            working_dir=output_dir,
            description="🔍 检查输出目录中的文件"
        )
        
        # 列出输出目录中的所有文件用于调试
        if os.path.exists(output_dir):
            files = os.listdir(output_dir)
            self.logger.info(f"📂 输出目录文件列表 | Output directory file list: {files}")
        
        # 检查配置文件是否存在
        if not os.path.exists(config_file):
            self.save_command_script(output_dir, "visualization_commands_FAILED.sh")
            raise RuntimeError(f"❌ 配置文件不存在 | Config file does not exist: {config_file}")
        
        try:
            # 构建并运行NGenomeSyn命令
            cmd = f"NGenomeSyn -InConf \"{os.path.basename(config_file)}\" -OutPut \"{os.path.basename(base_output)}\""
            
            self.logger.info(f"🚀 执行命令 | Executing command: {cmd}")
            self.logger.info(f"📁 工作目录 | Working directory: {output_dir}")
            
            # 记录NGenomeSyn测试命令
            test_cmd = "NGenomeSyn"
            self.log_command(
                command=test_cmd,
                working_dir=output_dir,
                description="🔧 测试NGenomeSyn是否可用"
            )
            
            # 先测试NGenomeSyn是否可用
            test_result = subprocess.run(
                ["NGenomeSyn"],
                capture_output=True,
                text=True,
                timeout=10
            )
            self.logger.info(f"🔧 NGenomeSyn测试结果 | NGenomeSyn test result: returncode={test_result.returncode}")
            if test_result.stdout:
                self.logger.info(f"📤 NGenomeSyn输出 | NGenomeSyn output: {test_result.stdout[:500]}")
            if test_result.stderr:
                self.logger.info(f"❌ NGenomeSyn错误 | NGenomeSyn stderr: {test_result.stderr[:500]}")
            
            # 记录实际的NGenomeSyn执行命令
            self.log_command(
                command=cmd,
                working_dir=output_dir,
                description=f"🎨 执行NGenomeSyn生成共线性图: {output_prefix}"
            )
            
            # 运行实际命令
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                cwd=output_dir,
                timeout=300  # 5分钟超时
            )
            
            # 详细记录执行结果
            self.logger.info(f"📊 命令执行结果 | Command execution result:")
            self.logger.info(f"🔢 返回码 | Return code: {result.returncode}")
            
            if result.stdout:
                self.logger.info(f"📤 标准输出 | STDOUT: {result.stdout}")
            else:
                self.logger.info("📤 标准输出为空 | STDOUT is empty")
            
            if result.stderr:
                self.logger.error(f"❌ 标准错误 | STDERR: {result.stderr}")
            else:
                self.logger.info("✅ 标准错误为空 | STDERR is empty")
            
            # 检查生成的文件
            svg_file = f"{base_output}.svg"
            
            # 记录检查输出文件命令
            check_output_cmd = f"ls -la \"{os.path.basename(base_output)}\"*"
            self.log_command(
                command=check_output_cmd,
                working_dir=output_dir,
                description="🔍 检查生成的输出文件"
            )
            
            # 修改错误检查逻辑：优先检查文件是否生成，而不是返回码
            if os.path.exists(svg_file) and os.path.getsize(svg_file) > 0:
                output_files.append(svg_file)
                self.logger.info(f"🎨 SVG文件已生成 | SVG file generated: {svg_file}")
                
                # 如果返回码不为0但文件已生成，只警告而不报错
                if result.returncode != 0:
                    self.logger.warning(f"⚠️ NGenomeSyn返回非零退出码但文件已成功生成 | NGenomeSyn returned non-zero exit code but files were generated successfully")
                
                # 如果需要PNG格式，使用convert转换
                if "png" in output_formats:
                    png_file = f"{base_output}.png"
                    cmd_convert = f"convert \"{os.path.basename(svg_file)}\" \"{os.path.basename(png_file)}\""
                    
                    # 记录PNG转换命令
                    self.log_command(
                        command=cmd_convert,
                        working_dir=output_dir,
                        description=f"🖼️ 将SVG转换为PNG格式: {os.path.basename(svg_file)} -> {os.path.basename(png_file)}"
                    )
                    
                    try:
                        convert_result = subprocess.run(
                            cmd_convert, 
                            shell=True, 
                            check=True, 
                            capture_output=True, 
                            text=True,
                            cwd=output_dir
                        )
                        if os.path.exists(png_file) and os.path.getsize(png_file) > 0:
                            output_files.append(png_file)
                            self.logger.info(f"🖼️ PNG文件已生成 | PNG file generated: {png_file}")
                    except subprocess.CalledProcessError as e:
                        self.logger.warning(f"⚠️ PNG转换失败 | PNG conversion failed: {e}")
                        self.logger.warning("🔧 请检查ImageMagick是否安装 | Please check if ImageMagick is installed")
                
                # 保存成功的命令脚本
                self.save_command_script(output_dir, "visualization_commands.sh")
                
            else:
                # 只有在文件确实没有生成时才报错
                error_msg = f"❌ NGenomeSyn执行失败 | NGenomeSyn failed: return code {result.returncode}"
                if result.stderr:
                    error_msg += f", stderr: {result.stderr}"
                if result.stdout:
                    error_msg += f", stdout: {result.stdout}"
                
                # 列出执行后的文件状态
                if os.path.exists(output_dir):
                    files_after = os.listdir(output_dir)
                    self.logger.info(f"📂 执行后文件列表 | File list after execution: {files_after}")
                
                self.logger.error(f"❌ SVG文件未生成或为空 | SVG file not generated or empty: {svg_file}")
                self.logger.error(error_msg)
                
                # 保存失败的命令脚本
                self.save_command_script(output_dir, "visualization_commands_FAILED.sh")
                raise RuntimeError(error_msg)
            
            return output_files
            
        except subprocess.TimeoutExpired:
            self.logger.error("⏰ NGenomeSyn执行超时 | NGenomeSyn execution timeout")
            self.save_command_script(output_dir, "visualization_commands_TIMEOUT.sh")
            raise RuntimeError("⏰ NGenomeSyn执行超时 | NGenomeSyn execution timeout")
        except Exception as e:
            self.logger.error(f"💥 NGenomeSyn可视化失败 | NGenomeSyn visualization failed: {e}")
            self.save_command_script(output_dir, "visualization_commands_FAILED.sh")
            raise