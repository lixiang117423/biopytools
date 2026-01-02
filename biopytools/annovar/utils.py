"""
ANNOVAR注释工具函数模块 | ANNOVAR Annotation Utility Functions Module
"""

import logging
import os
import subprocess
import sys
from pathlib import Path

class ANNOVARLogger:
    """ANNOVAR注释日志管理器 | ANNOVAR Annotation Logger Manager"""
    
    def __init__(self, output_dir: Path, log_name: str = "annovar_annotation.log"):
        self.output_dir = output_dir
        self.log_file = output_dir / log_name
        self.setup_logging()
    
    def setup_logging(self):
        """设置日志 | Setup logging"""
        if self.log_file.exists():
            self.log_file.unlink()
        
        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s - %(levelname)s - %(message)s',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)
    
    def get_logger(self):
        """获取日志器 | Get logger"""
        return self.logger

class CommandRunner:
    """命令执行器 | Command Runner"""
    
    def __init__(self, logger, working_dir: Path):
        self.logger = logger
        self.working_dir = working_dir
    
    def run(self, cmd: str, description: str = "") -> bool:
        """执行命令 | Execute command"""
        if description:
            self.logger.info(f"⚡ 执行步骤 | Executing step: {description}")
        
        self.logger.info(f"💻 命令 | Command: {cmd}")
        
        try:
            result = subprocess.run(
                cmd, 
                shell=True, 
                capture_output=True, 
                text=True, 
                check=True,
                cwd=self.working_dir
            )
            
            self.logger.info(f"✅ 命令执行成功 | Command executed successfully: {description}")
            
            if result.stdout:
                self.logger.debug(f"📊 标准输出 | Stdout: {result.stdout}")
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 命令执行失败 | Command execution failed: {description}")
            self.logger.error(f"🔢 错误代码 | Error code: {e.returncode}")
            self.logger.error(f"💬 错误信息 | Error message: {e.stderr}")
            return False

class GFF3Validator:
    """GFF3文件验证和修复器 | GFF3 File Validator and Fixer"""
    
    def __init__(self, logger):
        self.logger = logger
    
    def clean_and_fix_gff3(self, gff3_file: str):
        """清理和修复GFF3文件的格式问题 | Clean and fix GFF3 file format issues"""
        self.logger.info(
            "🧹 开始清理GFF3文件格式 | Starting GFF3 file format cleaning"
        )
        
        base_name = os.path.splitext(gff3_file)[0]
        clean_file = f"{base_name}.clean.gff3"
        final_file = f"{base_name}.final_fixed.gff3"
        
        # 步骤1: 清理attributes列，只保留ID、Name、Parent | Step 1: Clean attributes column, keep only ID, Name, Parent
        self.logger.info("🔧 步骤1: 清理attributes列 | Step 1: Cleaning attributes column")
        
        clean_cmd = f"""awk 'BEGIN{{FS=OFS="\\t"}} /^#/ {{print; next}} {{
    n = split($9, attrs, ";");
    new_attrs = "";
    for (i = 1; i <= n; i++) {{
        if (attrs[i] ~ /^ID=/ || attrs[i] ~ /^Name=/ || attrs[i] ~ /^Parent=/) {{
            if (new_attrs == "") {{
                new_attrs = attrs[i];
            }} else {{
                new_attrs = new_attrs ";" attrs[i];
            }}
        }}
    }}
    $9 = new_attrs;
    print $0;
}}' "{gff3_file}" > "{clean_file}" """
        
        try:
            result = subprocess.run(clean_cmd, shell=True, capture_output=True, text=True, check=True)
            self.logger.info(f"✅ 第一步清理完成 | First step cleaning completed: {clean_file}")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 第一步清理失败 | First step cleaning failed: {e.stderr}")
            return False
        
        # 步骤2: 修复坐标问题(start > end) | Step 2: Fix coordinate issues (start > end)
        self.logger.info("🔧 步骤2: 修复坐标问题 | Step 2: Fixing coordinate issues")
        
        fix_cmd = f"""awk 'BEGIN{{FS=OFS="\\t"}} /^#/ {{print; next}} {{
    if ($4 > $5) {{
        tmp = $4;
        $4 = $5;
        $5 = tmp;
    }}
    print $0;
}}' "{clean_file}" > "{final_file}" """
        
        try:
            result = subprocess.run(fix_cmd, shell=True, capture_output=True, text=True, check=True)
            self.logger.info(f"✅ 第二步修复完成 | Second step fixing completed: {final_file}")
            
            # 备份原文件 | Backup original file
            backup_file = f"{gff3_file}.original_backup"
            if not os.path.exists(backup_file):
                os.rename(gff3_file, backup_file)
                self.logger.info(f"💾 原文件已备份 | Original file backed up: {backup_file}")
            
            # 用修复后的文件替换原文件 | Replace original file with fixed file
            os.rename(final_file, gff3_file)
            self.logger.info(f"🔄 已用修复后的文件替换原文件 | Replaced original file with fixed file")
            
            # 清理临时文件 | Clean up temporary file
            if os.path.exists(clean_file):
                os.remove(clean_file)
            
            return True
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ 第二步修复失败 | Second step fixing failed: {e.stderr}")
            return False

    def check_gff3_header(self, gff3_file: str):
        """检查并确保GFF3文件有正确的头部 | Check and ensure GFF3 file has correct header"""
        with open(gff3_file, 'r') as f:
            first_line = f.readline().strip()
        
        if not first_line.startswith('##gff-version'):
            self.logger.warning(
                "⚠️ GFF3文件缺少版本头部，将添加 ##gff-version 3 | "
                "GFF3 file missing version header, will add ##gff-version 3"
            )
            temp_file = gff3_file + '.tmp'
            with open(gff3_file, 'r') as original, open(temp_file, 'w') as temp:
                temp.write('##gff-version 3\n')
                temp.write(original.read())
            
            os.replace(temp_file, gff3_file)
            self.logger.info("✅ 已添加GFF3版本头部 | Added GFF3 version header")
    
    def fix_gff3_cds_phase(self, gff3_file: str):
        """修复GFF3文件中CDS特征缺少phase的问题 | Fix missing phase in CDS features"""
        self.logger.info(
            "🔧 检查并修复GFF3文件中的CDS phase问题 | "
            "Checking and fixing CDS phase issues in GFF3 file"
        )
        
        temp_file = gff3_file + '.phase_fixed'
        lines_fixed = 0
        
        with open(gff3_file, 'r') as input_file, open(temp_file, 'w') as output_file:
            for line_num, line in enumerate(input_file, 1):
                line = line.rstrip('\n\r')
                
                # 跳过注释行和空行 | Skip comment lines and empty lines
                if line.startswith('#') or not line.strip():
                    output_file.write(line + '\n')
                    continue
                
                # 分割GFF行（应该有9列） | Split GFF line (should have 9 columns)
                fields = line.split('\t')
                
                if len(fields) >= 3 and fields[2].upper() == 'CDS':
                    # 这是一个CDS行，检查phase | This is a CDS line, check phase
                    if len(fields) < 8:
                        # 列数不够，补充到8列 | Not enough columns, pad to 8 columns
                        while len(fields) < 8:
                            fields.append('.')
                        fields[7] = '0'  # 设置phase为0 | Set phase to 0
                        lines_fixed += 1
                        self.logger.debug(f"🔧 行 | Line {line_num}: 添加缺失的列并设置phase为0 | Added missing columns and set phase to 0")
                    elif len(fields) == 8:
                        # 正好8列，但缺少phase | Exactly 8 columns, but missing phase
                        fields.append('0')  # 添加phase | Add phase
                        lines_fixed += 1
                        self.logger.debug(f"🔧 行 | Line {line_num}: 添加缺失的phase为0 | Added missing phase as 0")
                    elif len(fields) >= 9:
                        # 已有9列或更多，检查phase是否有效 | 9+ columns, check if phase is valid
                        phase = fields[7]
                        if phase not in ['0', '1', '2']:
                            fields[7] = '0'  # 修复无效的phase | Fix invalid phase
                            lines_fixed += 1
                            self.logger.debug(f"🔧 行 | Line {line_num}: 修复无效phase '{phase}' 为 '0' | Fixed invalid phase '{phase}' to '0'")
                    
                    # 确保至少有9列（包括attributes） | Ensure at least 9 columns (including attributes)
                    if len(fields) == 8:
                        fields.append('.')  # 添加空的attributes列 | Add empty attributes column
                    
                    output_file.write('\t'.join(fields) + '\n')
                else:
                    # 非CDS行，直接输出 | Non-CDS line, output directly
                    output_file.write(line + '\n')
        
        if lines_fixed > 0:
            # 备份原文件 | Backup original file
            backup_file = gff3_file + '.backup'
            os.rename(gff3_file, backup_file)
            os.rename(temp_file, gff3_file)
            self.logger.info(f"✅ 修复了 {lines_fixed} 行CDS的phase问题 | Fixed {lines_fixed} CDS phase issues")
            self.logger.info(f"💾 原文件已备份为 | Original file backed up as: {backup_file}")
        else:
            # 删除临时文件 | Remove temporary file
            os.remove(temp_file)
            self.logger.info("✅ GFF3文件中的CDS phase都正常，无需修复 | All CDS phases in GFF3 file are normal, no fix needed")