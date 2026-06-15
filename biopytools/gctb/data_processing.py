"""
GCTB数据处理模块|GCTB Data Processing Module

包含VCF到PLINK转换、质量控制和LD矩阵构建功能
Includes VCF to PLINK conversion, quality control, and LD matrix construction
"""

import os
import pandas as pd
import logging
from pathlib import Path
from typing import Optional, Tuple, List
from .config import GCTBConfig
from .utils import CommandRunner


class DataProcessor:
    """数据处理器|Data Processor"""

    def __init__(self, config: GCTBConfig, logger: logging.Logger,
                 cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 定义文件路径|Define file paths
        self.prefix = self.config.convert_dir / "input"
        self.plink_prefix = str(self.prefix)
        self.qc_prefix = self.config.qc_dir / "input_qc"
        self.plink_qc_prefix = str(self.qc_prefix)

    def vcf_to_plink(self) -> bool:
        """
        VCF文件转换为PLINK格式|Convert VCF to PLINK format

        Returns:
            bool: 转换是否成功|Whether conversion succeeded
        """
        self.logger.info("开始VCF到PLINK格式转换|Starting VCF to PLINK conversion")

        # 构建PLINK命令|Build PLINK command
        args = [
            '--vcf', self.config.vcf_file,
            '--make-bed',
            '--set-missing-var-ids', '@:#',  # 为缺失的SNP ID自动生成唯一ID|Auto-generate unique IDs for missing SNP IDs
            '--allow-extra-chr',  # 允许非常规染色体名|Allow non-standard chromosome names
            '--out', self.plink_prefix,
            '--threads', str(self.config.threads)
        ]

        cmd = self._build_plink_command(args)

        success = self.cmd_runner.run(
            cmd,
            "VCF到PLINK转换|VCF to PLINK conversion"
        )

        if success:
            self.logger.info(f"PLINK文件已生成|PLINK files generated: {self.plink_prefix}")

        return success

    def quality_control(self) -> bool:
        """
        执行质量控制过滤|Perform quality control filtering

        Returns:
            bool: 质控是否成功|Whether QC succeeded
        """
        self.logger.info("开始质量控制|Starting quality control")

        # 构建质控命令|Build QC command
        args = [
            '--bfile', self.plink_prefix,
            '--maf', str(self.config.maf_threshold),
            '--geno', str(self.config.miss_threshold),
            '--hwe', str(self.config.hwe_p),
            '--make-bed',
            '--out', self.plink_qc_prefix,
            '--threads', str(self.config.threads)
        ]

        cmd = self._build_plink_command(args)

        success = self.cmd_runner.run(
            cmd,
            "质量控制|Quality control"
        )

        if success:
            self.logger.info(f"质控后文件已生成|QC filtered files generated: {self.plink_qc_prefix}")

        return success

    def calculate_freq(self) -> bool:
        """
        计算等位基因频率|Calculate allele frequencies

        Returns:
            bool: 计算是否成功|Whether calculation succeeded
        """
        self.logger.info("开始计算等位基因频率|Starting allele frequency calculation")

        freq_file = f"{self.plink_qc_prefix}.frq"

        # 检查是否已存在|Check if already exists
        if os.path.exists(freq_file):
            self.logger.info(f"频率文件已存在，跳过|Frequency file already exists, skipping: {freq_file}")
            return True

        args = [
            '--bfile', self.plink_qc_prefix,
            '--freq',
            '--out', str(self.qc_prefix),
            '--threads', str(self.config.threads)
        ]

        cmd = self._build_plink_command(args)

        success = self.cmd_runner.run(
            cmd,
            "计算等位基因频率|Calculate allele frequencies"
        )

        return success

    def convert_and_split_phenotype(self) -> List[str]:
        """
        转换表型文件格式并拆分为多个单表型文件|Convert phenotype file format and split into multiple single-phenotype files

        支持的输入格式|Supported input formats:
        1. 标准格式（FID, IID, PHENO）|Standard format (FID, IID, PHENO)
        2. 样品名称格式（sample, trait1, trait2, ...）|Sample name format (sample, trait1, trait2, ...)

        Returns:
            List[str]: 转换后的表型文件路径列表|List of converted phenotype file paths
        """
        self.logger.info("开始转换表型文件|Starting phenotype file conversion")

        # 读取表型文件|Read phenotype file
        try:
            df = pd.read_csv(self.config.pheno_file, sep='\t')
        except Exception as e:
            self.logger.error(f"读取表型文件失败|Failed to read phenotype file: {e}")
            return []

        # 检测文件格式|Detect file format
        if len(df.columns) < 2:
            self.logger.error("表型文件列数不足|Phenotype file has insufficient columns")
            return []

        # 判断是否需要转换|Check if conversion is needed
        first_col = df.columns[0].lower()
        if first_col in ['fid', 'family']:
            # 已经是标准格式，检查是否有多个表型列|Already in standard format, check for multiple phenotype columns
            if 'PHENO' in df.columns or 'pheno' in df.columns:
                self.logger.info("表型文件已是标准格式|Phenotype file is already in standard format")
                # 如果有多个表型列，需要拆分|If multiple phenotype columns, need to split
                trait_cols = [col for col in df.columns if col.lower() not in ['fid', 'iid', 'family']]
                if len(trait_cols) > 1:
                    return self._split_standard_phenotype(df, trait_cols)
                return [self.config.pheno_file]
            else:
                # 找到表型列|Find phenotype columns
                trait_cols = [col for col in df.columns if col.lower() not in ['fid', 'iid', 'family']]
                if len(trait_cols) == 0:
                    self.logger.error("未找到表型列|No phenotype column found")
                    return []
                return self._split_standard_phenotype(df, trait_cols)
        else:
            # 样品名称格式，需要转换|Sample name format, need conversion
            return self._convert_sample_format(df)

    def _split_standard_phenotype(self, df: pd.DataFrame, trait_cols: List[str]) -> List[str]:
        """
        拆分标准格式的多表型文件|Split standard format multi-phenotype file

        Args:
            df: 数据框|Dataframe
            trait_cols: 表型列名列表|List of phenotype column names

        Returns:
            List[str]: 拆分后的表型文件路径列表|List of split phenotype file paths
        """
        self.logger.info(f"检测到 {len(trait_cols)} 个表型，将分别处理|Detected {len(trait_cols)} phenotypes, will process separately")

        converted_files = []
        output_dir = self.config.output_path / "converted_phenotypes"
        output_dir.mkdir(parents=True, exist_ok=True)

        for trait_col in trait_cols:
            # 创建单表型数据框|Create single-phenotype dataframe
            single_trait_df = df[['FID' if 'FID' in df.columns else df.columns[0],
                                   'IID' if 'IID' in df.columns else df.columns[1] if len(df.columns) > 1 else df.columns[0],
                                   trait_col]].copy()
            single_trait_df.columns = ['FID', 'IID', 'PHENO']

            # 保存|Save
            output_file = output_dir / f"{trait_col}.pheno"
            single_trait_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
            converted_files.append(str(output_file))
            self.logger.info(f"生成表型文件|Generated phenotype file: {output_file}")

        return converted_files

    def _convert_sample_format(self, df: pd.DataFrame) -> List[str]:
        """
        转换样品名称格式的表型文件|Convert sample name format phenotype file

        Args:
            df: 数据框|Dataframe

        Returns:
            List[str]: 转换后的表型文件路径列表|List of converted phenotype file paths
        """
        self.logger.info("检测到样品名称格式，正在转换|Detected sample name format, converting")

        # 第一列是样品名称，其余是表型|First column is sample name, rest are phenotypes
        sample_col = df.columns[0]
        trait_cols = [col for col in df.columns[1:]]

        self.logger.info(f"检测到 {len(trait_cols)} 个表型|Detected {len(trait_cols)} phenotypes")

        converted_files = []
        output_dir = self.config.output_path / "converted_phenotypes"
        output_dir.mkdir(parents=True, exist_ok=True)

        for trait_col in trait_cols:
            # 创建标准格式|Create standard format
            # FID和IID都使用样品名称，与PLINK的fam文件保持一致
            # Both FID and IID use sample name to match PLINK fam file
            single_trait_df = pd.DataFrame({
                'FID': df[sample_col],
                'IID': df[sample_col],
                'PHENO': df[trait_col]
            })

            # 移除缺失值|Remove missing values
            single_trait_df = single_trait_df.dropna(subset=['PHENO'])

            # 保存|Save
            output_file = output_dir / f"{trait_col}.pheno"
            single_trait_df.to_csv(output_file, sep='\t', index=False, float_format='%.6f')
            converted_files.append(str(output_file))
            self.logger.info(f"生成表型文件|Generated phenotype file: {output_file} ({len(single_trait_df)} samples)")

        return converted_files

    def build_sparse_ld_matrix(self) -> bool:
        """
        构建稀疏LD矩阵|Build sparse LD matrix

        Returns:
            bool: 构建是否成功|Whether construction succeeded
        """
        self.logger.info("开始构建稀疏LD矩阵|Starting sparse LD matrix construction")

        ld_prefix = str(self.config.ld_dir / "ld_matrix")

        # 使用GCTB构建稀疏LD矩阵|Use GCTB to build sparse LD matrix
        args = [
            '--bfile', self.plink_qc_prefix,
            '--make-sparse-ldm',
            '--out', ld_prefix
        ]

        cmd = self._build_gctb_command(args)

        success = self.cmd_runner.run(
            cmd,
            "构建稀疏LD矩阵|Build sparse LD matrix"
        )

        if success:
            self.logger.info(f"稀疏LD矩阵已生成|Sparse LD matrix generated: {ld_prefix}.ldm.sparse")

        return success

    def build_block_ld_matrix(self) -> bool:
        """
        构建分块LD矩阵|Build block LD matrix

        Returns:
            bool: 构建是否成功|Whether construction succeeded
        """
        self.logger.info("开始构建分块LD矩阵|Starting block LD matrix construction")

        ld_prefix = str(self.config.ld_dir / "ld_matrix")

        # 使用GCTB构建分块LD矩阵|Use GCTB to build block LD matrix
        args = [
            '--bfile', self.plink_qc_prefix,
            '--make-block-ldm',
            '--out', ld_prefix
        ]

        cmd = self._build_gctb_command(args)

        success = self.cmd_runner.run(
            cmd,
            "构建分块LD矩阵|Build block LD matrix"
        )

        if success:
            self.logger.info(f"分块LD矩阵已生成|Block LD matrix generated: {ld_prefix}.ldm.block")

        return success

    def eigen_decompose_ld_matrix(self) -> bool:
        """
        特征值分解LD矩阵|Eigen-decompose LD matrix

        Returns:
            bool: 分解是否成功|Whether decomposition succeeded
        """
        self.logger.info("开始LD矩阵特征值分解|Starting LD matrix eigen-decomposition")

        # 首先需要构建分块LD矩阵|First need to build block LD matrix
        ld_prefix = str(self.config.ld_dir / "ld_matrix")

        if not os.path.exists(f"{ld_prefix}.ldm.block.info"):
            self.logger.info("需要先构建分块LD矩阵|Need to build block LD matrix first")
            if not self.build_block_ld_matrix():
                return False

        # 然后进行特征值分解|Then perform eigen-decomposition
        args = [
            '--ldm', f"{ld_prefix}.ldm.block",
            '--make-ldm-eigen',
            '--out', ld_prefix
        ]

        cmd = self._build_gctb_command(args)

        success = self.cmd_runner.run(
            cmd,
            "LD矩阵特征值分解|LD matrix eigen-decomposition"
        )

        if success:
            self.logger.info(f"特征值分解LD矩阵已生成|Eigen-decomposed LD matrix generated: {ld_prefix}")

        return success

    def _build_plink_command(self, args: list) -> list:
        """构建PLINK命令|Build PLINK command"""
        from .utils import build_conda_command
        return build_conda_command(self.config.plink_path, args)

    def _build_gctb_command(self, args: list) -> list:
        """构建GCTB命令|Build GCTB command"""
        from .utils import build_conda_command
        return build_conda_command(self.config.gctb_path, args)

    def get_processed_files(self) -> dict:
        """
        获取处理后的文件列表|Get list of processed files

        Returns:
            dict: 文件类型到路径的映射|Mapping from file type to path
        """
        files = {
            'plink_prefix': self.plink_qc_prefix,
            'bed': f"{self.plink_qc_prefix}.bed",
            'bim': f"{self.plink_qc_prefix}.bim",
            'fam': f"{self.plink_qc_prefix}.fam",
            'freq': f"{self.plink_qc_prefix}.frq"
        }

        # 添加LD矩阵文件|Add LD matrix files
        ld_prefix = str(self.config.ld_dir / "ld_matrix")

        if self.config.ld_matrix_type == "sparse":
            files['ld_matrix'] = f"{ld_prefix}.ldm.sparse"
        elif self.config.ld_matrix_type == "block":
            files['ld_matrix'] = f"{ld_prefix}.ldm.block"
        elif self.config.ld_matrix_type == "eigen":
            files['ld_matrix'] = ld_prefix

        return files
