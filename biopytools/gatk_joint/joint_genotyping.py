"""
Joint Genotyping 主流程模块 | Joint Genotyping Main Pipeline Module
"""

import os
from pathlib import Path
from .utils import CommandRunner

class JointGenotyper:
    """Joint Genotyping 执行器 | Joint Genotyping Executor"""
    
    def __init__(self, config, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
    
    def run_genotype_gvcfs(self):
        """运行GenotypeGVCFs | Run GenotypeGVCFs"""
        self.logger.info("🧬 开始联合分型 | Starting joint genotyping")
        
        output_vcf = self.config.output_path / f"{self.config.base_name}_raw.vcf.gz"
        
        # 构建命令 | Build command
        cmd = f"{self.config.gatk_path} --java-options '-Xmx{self.config.memory}' GenotypeGVCFs"
        cmd += f" -R {self.config.reference}"
        cmd += f" -V gendb://{self.config.genomicsdb_path}"
        cmd += f" -O {output_vcf}"
        
        # 添加区间参数 | Add interval parameter
        if self.config.intervals:
            cmd += f" -L {self.config.intervals}"
        
        success = self.cmd_runner.run(
            cmd,
            "联合分型 | Joint genotyping"
        )
        
        if success:
            self.config.raw_vcf = str(output_vcf)
            self.logger.info(f"✅ 原始VCF已生成 | Raw VCF generated: {output_vcf}")
        
        return success

    def run_genomicsdb_import(self, sample_map_file):
        """运行GenomicsDBImport | Run GenomicsDBImport"""
        self.logger.info("📦 开始GenomicsDB导入 | Starting GenomicsDB import")
        
        genomicsdb_workspace = self.config.output_path / "genomicsdb_workspace"
        
        # 🔥 如果workspace已存在，删除它 | Remove existing workspace if it exists
        if genomicsdb_workspace.exists():
            self.logger.warning(f"⚠️ GenomicsDB workspace已存在，正在删除... | Workspace exists, removing...")
            import shutil
            try:
                shutil.rmtree(genomicsdb_workspace)
                self.logger.info(f"✅ 已删除旧workspace | Old workspace removed")
            except Exception as e:
                self.logger.error(f"❌ 删除workspace失败 | Failed to remove workspace: {e}")
                raise RuntimeError(f"无法删除已存在的workspace | Cannot remove existing workspace: {e}")
        
        # 🔥 如果用户没有指定区间，自动生成全基因组区间 | Auto-generate intervals if not specified
        if not self.config.intervals:
            self.logger.info("⚠️ 未指定区间，自动从参考基因组提取 | No intervals specified, extracting from reference")
            intervals_file = self._generate_intervals_from_reference()
            if not intervals_file:
                raise RuntimeError("❌ 无法生成区间文件 | Failed to generate intervals file")
            self.config.intervals = intervals_file
        
        # 构建命令 | Build command
        cmd = f"{self.config.gatk_path} --java-options '-Xmx{self.config.memory}' GenomicsDBImport"
        cmd += f" --genomicsdb-workspace-path {genomicsdb_workspace}"
        cmd += f" --sample-name-map {sample_map_file}"
        cmd += f" --reader-threads {self.config.threads}"
        
        # 添加区间参数 | Add interval parameter
        cmd += f" -L {self.config.intervals}"
        
        success = self.cmd_runner.run(
            cmd, 
            "GenomicsDB导入 | GenomicsDB import"
        )
        
        if success:
            self.config.genomicsdb_path = str(genomicsdb_workspace)
        
        return success

    def _generate_intervals_from_reference(self):
        """从参考基因组生成区间文件 | Generate intervals from reference genome"""
        import subprocess
        import os
        
        # 🔥 正确构建.dict和.fai文件路径 | Correctly build .dict and .fai paths
        ref_base = self.config.reference
        # 移除所有可能的参考基因组后缀 | Remove all possible reference suffixes
        for suffix in ['.fasta', '.fa', '.fna', '.fasta.gz', '.fa.gz', '.fna.gz']:
            if ref_base.endswith(suffix):
                ref_base = ref_base[:-len(suffix)]
                break
        
        dict_file = f"{ref_base}.dict"
        fai_file = f"{self.config.reference}.fai"
        intervals_file = self.config.output_path / "intervals.list"
        
        self.logger.info(f"📍 参考基因组 | Reference: {self.config.reference}")
        self.logger.info(f"📍 字典文件路径 | Dict path: {dict_file}")
        self.logger.info(f"📍 索引文件路径 | Fai path: {fai_file}")
        
        try:
            # 第1步：检查并创建.dict文件 | Step 1: Check and create .dict file
            if not os.path.exists(dict_file):
                self.logger.info("⚠️ 未找到字典文件，正在创建... | Dict file not found, creating...")
                cmd = f"{self.config.gatk_path} CreateSequenceDictionary -R {self.config.reference}"
                result = subprocess.run(cmd, shell=True, capture_output=True, text=True, cwd=self.config.output_path)
                
                if result.returncode == 0:
                    self.logger.info("✅ 字典文件创建成功 | Dict file created successfully")
                else:
                    self.logger.error(f"❌ 字典文件创建失败 | Dict creation failed: {result.stderr}")
                    # 尝试使用.fai文件 | Try using .fai file
                    if not os.path.exists(fai_file):
                        self.logger.info("⚠️ 尝试创建索引文件... | Trying to create fai file...")
                        fai_cmd = f"samtools faidx {self.config.reference}"
                        fai_result = subprocess.run(fai_cmd, shell=True, capture_output=True, text=True)
                        if fai_result.returncode != 0:
                            self.logger.error(f"❌ 索引文件创建失败 | Fai creation failed: {fai_result.stderr}")
                            return None
            
            # 第2步：从.dict文件提取染色体 | Step 2: Extract chromosomes from .dict
            if os.path.exists(dict_file):
                self.logger.info(f"✅ 使用字典文件 | Using dict file: {dict_file}")
                with open(dict_file, 'r') as f:
                    chromosomes = []
                    for line in f:
                        if line.startswith('@SQ'):
                            # 提取SN字段 | Extract SN field
                            for field in line.split('\t'):
                                if field.startswith('SN:'):
                                    chrom = field.replace('SN:', '').strip()
                                    chromosomes.append(chrom)
                    
                    if chromosomes:
                        with open(intervals_file, 'w') as out:
                            for chrom in chromosomes:
                                out.write(f"{chrom}\n")
                        self.logger.info(f"✅ 生成区间文件 | Generated intervals: {len(chromosomes)} chromosomes")
                        self.logger.info(f"📄 区间文件位置 | Intervals file: {intervals_file}")
                        return str(intervals_file)
            
            # 第3步：如果.dict不可用，尝试.fai | Step 3: Try .fai if .dict unavailable
            if os.path.exists(fai_file):
                self.logger.info(f"✅ 使用索引文件 | Using fai file: {fai_file}")
                with open(fai_file, 'r') as f:
                    chromosomes = []
                    for line in f:
                        if line.strip():
                            chrom = line.split('\t')[0]
                            chromosomes.append(chrom)
                    
                    if chromosomes:
                        with open(intervals_file, 'w') as out:
                            for chrom in chromosomes:
                                out.write(f"{chrom}\n")
                        self.logger.info(f"✅ 生成区间文件 | Generated intervals: {len(chromosomes)} chromosomes")
                        self.logger.info(f"📄 区间文件位置 | Intervals file: {intervals_file}")
                        return str(intervals_file)
            
            # 第4步：都失败了 | Step 4: Everything failed
            self.logger.error("❌ 无法从任何源提取染色体信息 | Cannot extract chromosome info from any source")
            return None
        
        except Exception as e:
            self.logger.error(f"❌ 生成区间文件失败 | Failed to generate intervals: {e}")
            import traceback
            self.logger.error(traceback.format_exc())
            return None