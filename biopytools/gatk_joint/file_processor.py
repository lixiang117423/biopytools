"""
文件处理模块 | File Processing Module
"""

import os
import glob
from pathlib import Path

class FileTypeDetector:
    """文件类型检测器 | File Type Detector"""
    
    def __init__(self, config, logger):
        self.config = config
        self.logger = logger
    
    # def detect_file_type(self):
    #     """检测输入文件类型 | Detect input file type"""
    #     self.logger.info("🔎 检测输入文件类型 | Detecting input file type")
        
    #     # 获取所有vcf相关文件 | Get all vcf-related files
    #     vcf_files = []
    #     for ext in ['*.vcf', '*.vcf.gz', '*.g.vcf', '*.g.vcf.gz', '*.gvcf', '*.gvcf.gz']:
    #         vcf_files.extend(glob.glob(os.path.join(self.config.input_dir, ext)))
        
    #     if not vcf_files:
    #         raise ValueError(f"❌ 未找到VCF/GVCF文件 | No VCF/GVCF files found in: {self.config.input_dir}")
        
    #     self.logger.info(f"📊 找到 {len(vcf_files)} 个文件 | Found {len(vcf_files)} files")
        
    #     # 检测文件类型 | Detect file type
    #     gvcf_count = sum(1 for f in vcf_files if '.g.vcf' in f or '.gvcf' in f)
    #     vcf_count = len(vcf_files) - gvcf_count
        
    #     if gvcf_count > 0 and vcf_count == 0:
    #         self.config.file_type = 'gvcf'
    #         self.logger.info(f"✅ 检测到GVCF文件 | Detected GVCF files: {gvcf_count} files")
    #     elif vcf_count > 0 and gvcf_count == 0:
    #         self.config.file_type = 'vcf'
    #         self.logger.info(f"✅ 检测到VCF文件 | Detected VCF files: {vcf_count} files")
    #     else:
    #         self.config.file_type = 'mixed'
    #         self.logger.warning(f"⚠️ 检测到混合文件类型 | Detected mixed file types: {gvcf_count} GVCF + {vcf_count} VCF")
    #         self.logger.warning("将按GVCF流程处理 | Will process as GVCF")
        
    #     return vcf_files

    def detect_file_type(self):
        """检测输入文件类型 | Detect input file type"""
        self.logger.info("🔎 检测输入文件类型 | Detecting input file type")
        
        # 🔥 优先检测GVCF文件（包含所有位点）| Check GVCF files first (contains all sites)
        gvcf_files = []
        for ext in ['*.g.vcf.gz', '*.g.vcf', '*.gvcf.gz', '*.gvcf']:
            gvcf_files.extend(glob.glob(os.path.join(self.config.input_dir, ext)))
        
        # 过滤掉索引文件 | Filter out index files
        gvcf_files = [f for f in gvcf_files if not f.endswith(('.tbi', '.idx', '.csi'))]
        
        if gvcf_files:
            # 找到GVCF文件，使用GVCF流程 | Found GVCF files, use GVCF pipeline
            vcf_files = gvcf_files
            self.config.file_type = 'gvcf'
            self.logger.info(f"✅ 检测到GVCF文件 | Detected GVCF files: {len(vcf_files)} files")
        else:
            # 没有GVCF，检测普通VCF文件（仅变异位点）| No GVCF, check regular VCF files (variant sites only)
            vcf_files = []
            for ext in ['*.vcf.gz', '*.vcf']:
                vcf_files.extend(glob.glob(os.path.join(self.config.input_dir, ext)))
            
            # 过滤掉索引文件 | Filter out index files
            vcf_files = [f for f in vcf_files if not f.endswith(('.tbi', '.idx', '.csi'))]
            
            if vcf_files:
                self.config.file_type = 'vcf'
                self.logger.info(f"✅ 检测到VCF文件 | Detected VCF files: {len(vcf_files)} files")
            else:
                raise ValueError(f"❌ 未找到VCF/GVCF文件 | No VCF/GVCF files found in: {self.config.input_dir}")
        
        self.logger.info(f"📊 共找到 {len(vcf_files)} 个文件 | Total {len(vcf_files)} files found")
        
        # 排序保证结果可重复 | Sort for reproducibility
        vcf_files.sort()
        
        return vcf_files
    
    def create_sample_map(self, vcf_files):
        """创建样本映射文件 | Create sample map file"""
        sample_map_file = self.config.output_path / "sample_map.txt"
        
        self.logger.info(f"📝 创建样本映射文件 | Creating sample map file")
        
        with open(sample_map_file, 'w') as f:
            for vcf_file in vcf_files:
                sample_name = Path(vcf_file).stem
                # 移除.g或.gvcf后缀 | Remove .g or .gvcf suffix
                sample_name = sample_name.replace('.g.vcf', '').replace('.g', '').replace('.vcf', '')
                f.write(f"{sample_name}\t{vcf_file}\n")
        
        self.logger.info(f"✅ 样本映射文件已创建 | Sample map created: {sample_map_file}")
        return str(sample_map_file)
