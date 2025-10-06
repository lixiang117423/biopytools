"""
系统发育树构建模块 | Phylogenetic Tree Construction Module
使用scikit-bio构建NJ树的简化实现 | Simplified implementation using scikit-bio for NJ tree construction
"""

import os
from pathlib import Path
from typing import Tuple, List
import numpy as np

class TreeBuilder:
    """系统发育树构建器 | Phylogenetic Tree Builder"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    def build_nj_tree(self) -> bool:
        """构建NJ系统发育树 | Build NJ phylogenetic tree"""
        self.logger.info("开始构建NJ系统发育树 | Starting NJ phylogenetic tree construction")
        
        try:
            # 导入scikit-bio | Import scikit-bio
            try:
                from skbio import DistanceMatrix
                from skbio.tree import nj
                self.logger.info("✓ 成功导入 scikit-bio | Successfully imported scikit-bio")
            except ImportError as e:
                self.logger.error("缺少scikit-bio库 | Missing scikit-bio library")
                self.logger.error("请运行: pip install scikit-bio | Please run: pip install scikit-bio")
                self.logger.error(f"或者: conda install -c conda-forge scikit-bio | Or: conda install -c conda-forge scikit-bio")
                return False
            
            # 读取距离矩阵 | Read distance matrix
            self.logger.info(f"读取距离矩阵文件 | Reading distance matrix file: {self.config.distance_matrix}")
            distance_data, sample_names = self._read_distance_matrix()
            
            self.logger.info(f"成功读取 {len(sample_names)} 个样本的距离矩阵 | Successfully read distance matrix for {len(sample_names)} samples")
            self.logger.info(f"样本列表 | Sample list: {sample_names[:5]}{'...' if len(sample_names) > 5 else ''}")
            
            # 创建scikit-bio DistanceMatrix对象 | Create scikit-bio DistanceMatrix object
            self.logger.info("创建 scikit-bio DistanceMatrix 对象 | Creating scikit-bio DistanceMatrix object")
            dm = DistanceMatrix(distance_data, ids=sample_names)
            
            # 使用NJ算法构建树 | Build tree using NJ algorithm
            self.logger.info("使用NJ算法构建系统发育树 | Building phylogenetic tree using NJ algorithm")
            tree = nj(dm)
            
            # 获取Newick格式字符串 | Get Newick format string
            newick_string = str(tree)
            self.logger.info(f"成功构建NJ树，树长度: {len(newick_string)} 字符 | Successfully built NJ tree, tree length: {len(newick_string)} characters")
            
            # 保存树文件 | Save tree file
            self._save_tree(newick_string)
            
            # 输出树的ASCII表示（如果样本数不太多）| Output ASCII representation of tree (if not too many samples)
            if len(sample_names) <= 20:
                try:
                    self.logger.info("系统发育树 ASCII 表示 | Phylogenetic tree ASCII representation:")
                    ascii_art = tree.ascii_art()
                    for line in ascii_art.split('\n')[:20]:  # 限制输出行数
                        self.logger.info(line)
                except Exception as e:
                    self.logger.warning(f"无法生成ASCII树表示 | Cannot generate ASCII tree representation: {e}")
            
            self.logger.info(f"NJ系统发育树已保存 | NJ phylogenetic tree saved: {self.config.tree_output}")
            self._log_tree_stats(sample_names)
            return True
            
        except Exception as e:
            self.logger.error(f"NJ树构建失败 | NJ tree construction failed: {e}")
            import traceback
            self.logger.error(f"详细错误信息 | Detailed error: {traceback.format_exc()}")
            return False
    
    def _read_distance_matrix(self) -> Tuple[np.ndarray, List[str]]:
        """读取距离矩阵文件 | Read distance matrix file"""
        try:
            with open(self.config.distance_matrix, 'r') as f:
                lines = [line.strip() for line in f.readlines() if line.strip()]
            
            if not lines:
                raise ValueError("距离矩阵文件为空 | Distance matrix file is empty")
            
            # 检查第一行是否是数字（矩阵维度）| Check if first line is a number (matrix dimension)
            first_line = lines[0]
            data_start_idx = 0
            n_samples = None
            
            try:
                n_samples = int(first_line)
                data_start_idx = 1
                self.logger.info(f"检测到矩阵维度信息: {n_samples} | Detected matrix dimension: {n_samples}")
            except ValueError:
                # 第一行不是数字，可能直接是数据 | First line is not a number, might be data directly
                n_samples = len(lines)
                data_start_idx = 0
                self.logger.info(f"未检测到维度信息，推断样本数: {n_samples} | No dimension detected, inferred sample count: {n_samples}")
            
            # 解析数据行 | Parse data lines
            sample_names = []
            distance_matrix = []
            
            data_lines = lines[data_start_idx:data_start_idx + n_samples]
            
            for i, line in enumerate(data_lines):
                parts = line.split()
                if len(parts) < 2:
                    raise ValueError(f"数据行格式错误 | Invalid data line format: line {i+1}")
                
                # 第一列是样本名 | First column is sample name
                sample_name = parts[0]
                sample_names.append(sample_name)
                
                # 其余列是距离值 | Remaining columns are distance values
                distances = []
                for j, val_str in enumerate(parts[1:n_samples+1]):  # 确保只读取n_samples个值
                    try:
                        val = float(val_str)
                        distances.append(val)
                    except (ValueError, IndexError):
                        # 如果缺少值，用0填充 | If missing values, fill with 0
                        distances.append(0.0)
                        self.logger.warning(f"样本 {sample_name} 缺少距离值，用0填充 | Sample {sample_name} missing distance value, filled with 0")
                
                # 确保距离向量长度正确 | Ensure distance vector has correct length
                while len(distances) < n_samples:
                    distances.append(0.0)
                distances = distances[:n_samples]  # 截断多余的值
                
                distance_matrix.append(distances)
            
            # 转换为numpy数组 | Convert to numpy array
            distance_array = np.array(distance_matrix, dtype=float)
            
            # 验证和修正矩阵 | Validate and fix matrix
            distance_array = self._validate_and_fix_matrix(distance_array, sample_names)
            
            return distance_array, sample_names
            
        except Exception as e:
            self.logger.error(f"读取距离矩阵失败 | Failed to read distance matrix: {e}")
            raise
    
    def _validate_and_fix_matrix(self, distance_array: np.ndarray, sample_names: List[str]) -> np.ndarray:
        """验证和修正距离矩阵 | Validate and fix distance matrix"""
        n_samples = len(sample_names)
        
        self.logger.info(f"验证距离矩阵 | Validating distance matrix: {distance_array.shape}")
        
        # 检查矩阵形状 | Check matrix shape
        if distance_array.shape[0] != n_samples:
            raise ValueError(f"距离矩阵行数({distance_array.shape[0]})与样本数({n_samples})不匹配 | "
                           f"Distance matrix rows ({distance_array.shape[0]}) don't match sample count ({n_samples})")
        
        if distance_array.shape[1] != n_samples:
            raise ValueError(f"距离矩阵列数({distance_array.shape[1]})与样本数({n_samples})不匹配 | "
                           f"Distance matrix columns ({distance_array.shape[1]}) don't match sample count ({n_samples})")
        
        # 使矩阵对称 | Make matrix symmetric
        if not np.allclose(distance_array, distance_array.T, rtol=1e-10):
            self.logger.warning("距离矩阵不对称，正在对称化 | Distance matrix is not symmetric, symmetrizing")
            distance_array = (distance_array + distance_array.T) / 2
        
        # 确保对角线为0 | Ensure diagonal is 0
        np.fill_diagonal(distance_array, 0)
        
        # 检查负值 | Check for negative values
        if np.any(distance_array < 0):
            self.logger.warning("距离矩阵包含负值，将转换为绝对值 | Distance matrix contains negative values, converting to absolute values")
            distance_array = np.abs(distance_array)
        
        # 检查是否有无效值 | Check for invalid values
        if np.any(np.isnan(distance_array)) or np.any(np.isinf(distance_array)):
            self.logger.warning("距离矩阵包含NaN或无穷值，正在修复 | Distance matrix contains NaN or infinite values, fixing")
            distance_array = np.nan_to_num(distance_array, nan=0.0, posinf=1e6, neginf=0.0)
        
        self.logger.info(f"距离矩阵验证完成 | Distance matrix validation completed")
        self.logger.info(f"矩阵统计 | Matrix statistics: min={np.min(distance_array):.6f}, max={np.max(distance_array):.6f}, mean={np.mean(distance_array):.6f}")
        
        return distance_array
    
    def _save_tree(self, newick_string: str):
        """保存系统发育树 | Save phylogenetic tree"""
        # 确保输出目录存在 | Ensure output directory exists
        output_dir = Path(self.config.tree_output).parent
        output_dir.mkdir(parents=True, exist_ok=True)
        
        # 保存Newick格式的树 | Save tree in Newick format
        with open(self.config.tree_output, 'w') as f:
            f.write(newick_string)
            if not newick_string.endswith('\n'):
                f.write('\n')
        
        self.logger.info(f"系统发育树已保存到 | Phylogenetic tree saved to: {self.config.tree_output}")
    
    def _log_tree_stats(self, sample_names: List[str]):
        """记录系统发育树统计信息 | Log phylogenetic tree statistics"""
        self.logger.info(f"系统发育树统计 | Phylogenetic tree statistics:")
        self.logger.info(f"  - 样本数 | Sample count: {len(sample_names)}")
        self.logger.info(f"  - 树格式 | Tree format: Newick")
        self.logger.info(f"  - 构建方法 | Construction method: Neighbor-Joining (NJ)")
        self.logger.info(f"  - 使用库 | Library used: scikit-bio")
        
        if os.path.exists(self.config.tree_output):
            size = os.path.getsize(self.config.tree_output)
            self.logger.info(f"  - 文件大小 | File size: {size} bytes")
            
            # 读取并验证树文件 | Read and validate tree file
            try:
                with open(self.config.tree_output, 'r') as f:
                    tree_content = f.read().strip()
                    if tree_content.endswith(';'):
                        self.logger.info(f"  - 树格式验证 | Tree format validation: ✓ 有效的Newick格式 | Valid Newick format")
                    else:
                        self.logger.warning(f"  - 树格式验证 | Tree format validation: ⚠ 可能格式不完整 | Possibly incomplete format")
            except Exception as e:
                self.logger.warning(f"  - 无法验证树文件 | Cannot validate tree file: {e}")
