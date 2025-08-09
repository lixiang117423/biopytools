"""
KMC3接口封装模块 | KMC3 Interface Wrapper Module 🤖🔧
"""

import os
import subprocess
import tempfile
import shutil
from pathlib import Path
from typing import List, Dict, Optional, Tuple, Union
import logging

class KMCInterface:
    """KMC3接口封装类 🤖"""
    
    def __init__(self, config):
        self.config = config
        self.logger = logging.getLogger(__name__)
        
        # 检查KMC是否可用
        self._check_kmc_availability()
    
    # def _check_kmc_availability(self):
    #     """检查KMC是否安装和可用"""
    #     try:
    #         result = subprocess.run(['kmc'], capture_output=True, text=True)
    #         if "K-Mer Counter (KMC)" not in result.stderr:
    #             raise RuntimeError("KMC not found or not working properly")
    #         self.logger.info("KMC3 is available")
    #     except FileNotFoundError:
    #         raise RuntimeError("KMC not found in PATH. Please install KMC3")

    def _check_kmc_availability(self):
        """检查KMC是否安装和可用 🔍"""
        try:
            # 尝试运行 kmc -h 来获取帮助信息
            result = subprocess.run(['kmc', '-h'], capture_output=True, text=True)
            
            # KMC的帮助信息可能在stdout或stderr中
            output_text = result.stdout + result.stderr
            
            if "K-Mer Counter (KMC)" not in output_text:
                # 如果-h不工作，尝试无参数运行
                result2 = subprocess.run(['kmc'], capture_output=True, text=True)
                output_text2 = result2.stdout + result2.stderr
                
                if "K-Mer Counter (KMC)" not in output_text2:
                    raise RuntimeError("KMC not found or not working properly ❌")
            
            # 提取版本信息（如果可能）
            version_info = "unknown version"
            for line in output_text.split('\n'):
                if "K-Mer Counter (KMC)" in line and "ver." in line:
                    version_info = line.strip()
                    break
            
            self.logger.info(f"KMC3 is available: {version_info} ✅")
            
        except FileNotFoundError:
            raise RuntimeError("KMC not found in PATH. Please install KMC3 ❌")
        except Exception as e:
            raise RuntimeError(f"Error checking KMC availability: {e} ❌")
    
    def count_kmers(self, input_files, output_prefix: str, 
               file_format: str = "fq") -> str:
        """使用KMC计数k-mer，支持单文件或文件列表"""

        # 添加这行调试
        self.logger.info(f"🔧 count_kmers called with {input_files}, {output_prefix}, {file_format}")
        
        # 处理输入参数
        if isinstance(input_files, str):
            file_list = [input_files]
        else:
            file_list = list(input_files)
        
        self.logger.info(f"🔧 KMC input files: {len(file_list)} files")
        for i, f in enumerate(file_list):
            self.logger.info(f"🔧   {i+1}. {os.path.basename(f)}")
        
        # 创建临时文件列表
        temp_list_file = None
        if len(file_list) > 1:
            temp_list_file = self._create_file_list(file_list)
            input_arg = f"@{temp_list_file}"
            self.logger.info(f"🔧 Created file list: {temp_list_file}")
        else:
            input_arg = file_list[0]
        
        # 构建KMC命令
        cmd = self._build_kmc_command(input_arg, output_prefix, file_format)
        
        try:
            # 显示完整的KMC命令
            self.logger.info(f"🚀 Executing KMC command:")
            self.logger.info(f"🚀   {' '.join(cmd)}")
            
            # 执行KMC命令
            result = subprocess.run(cmd, capture_output=True, text=True, check=True)
            
            # 显示KMC输出
            if result.stdout:
                self.logger.info(f"📤 KMC stdout: {result.stdout}")
            if result.stderr:
                self.logger.info(f"📤 KMC stderr: {result.stderr}")
            
            if self.config.verbose:
                self.logger.info(f"KMC return code: {result.returncode}")
            
            self.logger.info("✅ KMC counting completed successfully")
            return output_prefix
            
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ KMC failed with return code {e.returncode}")
            self.logger.error(f"❌ KMC command: {' '.join(cmd)}")
            self.logger.error(f"❌ KMC stdout: {e.stdout}")
            self.logger.error(f"❌ KMC stderr: {e.stderr}")
            raise RuntimeError(f"KMC counting failed: {e.stderr}")
        
        finally:
            # 清理临时文件
            if temp_list_file and os.path.exists(temp_list_file):
                os.unlink(temp_list_file)
                self.logger.info(f"🗑️ Cleaned up temp file: {temp_list_file}")
    
    def _create_file_list(self, files: List[str]) -> str:
        """创建KMC输入文件列表 📝"""
        temp_fd, temp_path = tempfile.mkstemp(suffix='.txt', prefix='kmc_files_')
        
        try:
            with os.fdopen(temp_fd, 'w') as f:
                for file_path in files:
                    f.write(f"{file_path}\n")
            return temp_path
        except:
            os.close(temp_fd)
            if os.path.exists(temp_path):
                os.unlink(temp_path)
            raise
    
    def _build_kmc_command(self, input_arg: str, output_prefix: str, 
                          file_format: str) -> List[str]:
        """构建KMC命令行 🛠️"""
        cmd = ['kmc']
        
        # KMC参数
        kmc_params = self.config.get_kmc_params()
        
        # 添加参数
        cmd.extend([f'-k{kmc_params["k"]}'])
        cmd.extend([f'-m{kmc_params["m"]}'])
        cmd.extend([f'-ci{kmc_params["ci"]}'])
        cmd.extend([f'-cx{kmc_params["cx"]}'])
        cmd.extend([f'-p{kmc_params["p"]}'])
        cmd.extend([f'-t{kmc_params["t"]}'])
        
        # 可选参数
        if kmc_params.get('b', False):
            cmd.append('-b')
        if kmc_params.get('r', False):
            cmd.append('-r')
        if kmc_params.get('sm', False):
            cmd.append('-sm')
        if kmc_params.get('v', False):
            cmd.append('-v')
        
        # 文件格式
        cmd.append(f'-f{file_format}')
        
        # 输入、输出、工作目录
        cmd.extend([input_arg, output_prefix, self.config.temp_dir])
        
        return cmd
    
    def query_kmers(self, kmc_database: str, query_kmers: List[str]) -> Dict[str, int]:
        """
        查询k-mer在数据库中的丰度 🔎
        
        Args:
            kmc_database: KMC数据库文件前缀 🗃️
            query_kmers: 要查询的k-mer列表 🧬
        
        Returns:
            k-mer -> 丰度的字典 📊
        """
        # 这里需要使用KMC的Python API或者kmc_tools
        # 暂时使用kmc_tools dump然后解析
        
        # 创建查询文件
        query_file = tempfile.mktemp(suffix='.txt')
        result_file = tempfile.mktemp(suffix='.txt')
        
        try:
            # 写入查询k-mer ✍️
            with open(query_file, 'w') as f:
                for kmer in query_kmers:
                    f.write(f"{kmer}\n")
            
            # 使用kmc_tools查询 ⚙️
            # cmd = [
            #     'kmc_tools', 'simple',
            #     kmc_database,
            #     '-ci1', '-cx1000000',
            #     'dump', result_file
            # ]
            # 使用kmc_dump查询
            cmd = ['kmc_dump', kmc_database, result_file]
            
            subprocess.run(cmd, check=True, capture_output=True)
            
            # 解析结果 📊
            results = {}
            if os.path.exists(result_file):
                with open(result_file, 'r') as f:
                    for line in f:
                        if line.strip():
                            parts = line.strip().split('\t')
                            if len(parts) >= 2:
                                kmer, count = parts[0], int(parts[1])
                                if kmer in query_kmers:
                                    results[kmer] = count
            
            # 添加未找到的k-mer（丰度为0）👻
            for kmer in query_kmers:
                if kmer not in results:
                    results[kmer] = 0
            
            return results
            
        finally:
            # 清理临时文件 🧹
            for temp_file in [query_file, result_file]:
                if os.path.exists(temp_file):
                    os.unlink(temp_file)
    
    def get_kmer_list(self, kmc_database: str, min_count: int = 1) -> List[Tuple[str, int]]:
        """
        获取KMC数据库中的所有k-mer列表 📋
        
        Args:
            kmc_database: KMC数据库文件前缀 🗃️
            min_count: 最小计数阈值 📉
        
        Returns:
            [(kmer, count), ...] 列表 📊
        """
        result_file = tempfile.mktemp(suffix='.txt')
        
        try:
            # 使用kmc_dump导出k-mer 📤
            cmd = ['kmc_dump', kmc_database, result_file]
            subprocess.run(cmd, check=True, capture_output=True)
            
            # 解析结果 📊
            kmers = []
            with open(result_file, 'r') as f:
                for line in f:
                    if line.strip():
                        parts = line.strip().split('\t')
                        if len(parts) >= 2:
                            kmer, count = parts[0], int(parts[1])
                            if count >= min_count:
                                kmers.append((kmer, count))
            
            self.logger.info(f"Retrieved {len(kmers)} k-mers from database ✅")
            return kmers
            
        finally:
            if os.path.exists(result_file):
                os.unlink(result_file)
    
    def intersect_databases(self, db1: str, db2: str, output: str) -> str:
        """计算两个KMC数据库的交集 🔗"""
        cmd = [
            'kmc_tools', 'simple',
            db1, db2, 'intersect', output
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            return output
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"KMC intersect failed: {e} ❌")
    
    def union_databases(self, db1: str, db2: str, output: str) -> str:
        """计算两个KMC数据库的并集 ➕"""
        cmd = [
            'kmc_tools', 'simple',
            db1, db2, 'union', output
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            return output
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"KMC union failed: {e} ❌")
    
    def subtract_databases(self, db1: str, db2: str, output: str) -> str:
        """计算两个KMC数据库的差集 (db1 - db2) ➖"""
        cmd = [
            'kmc_tools', 'simple',
            db1, db2, 'kmers_subtract', output
        ]
        
        try:
            subprocess.run(cmd, check=True, capture_output=True)
            return output
        except subprocess.CalledProcessError as e:
            raise RuntimeError(f"KMC subtract failed: {e} ❌")