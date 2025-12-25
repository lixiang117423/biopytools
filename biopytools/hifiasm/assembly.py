"""
HiFiasm组装模块 | HiFiasm Assembly Module
"""

import os
import time
import logging
from pathlib import Path
from typing import List, Dict, Optional, Tuple

class HifiasmAssembler:
    """HiFiasm组装器 | HiFiasm Assembler"""
    
    def __init__(self, config, logger: logging.Logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.assembly_dir = Path(config.output_dir) / 'assembly'
        self.assembly_dir.mkdir(exist_ok=True)
    
    def run_assembly(self) -> bool:
        """运行HiFiasm组装 | Run HiFiasm assembly"""
        try:
            self.logger.info("准备HiFiasm组装 | Preparing HiFiasm assembly")
            
            # 设置工作目录
            original_working_dir = self.cmd_runner.working_dir
            self.cmd_runner.working_dir = self.assembly_dir
            
            try:
                # 构建HiFiasm命令
                hifiasm_cmd = self._build_hifiasm_command()
                
                # 记录命令参数解释
                self._log_parameter_explanations(hifiasm_cmd)
                
                # 运行组装
                if not self._execute_hifiasm(hifiasm_cmd):
                    return False
                
                # 检查输出文件
                if not self._verify_assembly_outputs():
                    return False
                
                self.logger.success("HiFiasm组装完成 | HiFiasm assembly completed")
                return True
                
            finally:
                # 恢复工作目录
                self.cmd_runner.working_dir = original_working_dir
                
        except Exception as e:
            self.logger.error(f"HiFiasm组装失败 | HiFiasm assembly failed: {e}")
            return False
    
    def _build_hifiasm_command(self) -> List[str]:
        """构建HiFiasm命令 | Build HiFiasm command"""
        cmd = [
            self.config.hifiasm_path,
            '-o', self.config.prefix,
            '-t', str(self.config.threads)
        ]
        
        # 基因组大小估计
        genome_size = self.config.estimate_genome_size()
        cmd.extend(['--hg-size', genome_size])
        
        # Purge参数
        cmd.extend(['-l', str(self.config.purge_level)])
        cmd.extend(['--purge-max', str(self.config.purge_max)])
        
        # 相似性阈值
        cmd.extend(['-s', str(self.config.similarity_threshold)])
        
        # ONT数据
        if self.config.ont_reads:
            cmd.extend(['--ul', self.config.ont_reads])
        
        # Hi-C数据
        if self.config.hi_c_1 and self.config.hi_c_2:
            cmd.extend(['--h1', self.config.hi_c_1])
            cmd.extend(['--h2', self.config.hi_c_2])
        
        # 根据组装类型调整参数
        if self.config.assembly_type == 'triploid':
            cmd.extend(['--trio'])
        elif self.config.assembly_type == 'polyploid':
            cmd.extend(['--polyploid'])
        
        # 额外参数
        if self.config.extra_hifiasm_args:
            extra_args = self.config.extra_hifiasm_args.split()
            cmd.extend(extra_args)
        
        # 输入文件（最后添加）
        cmd.append(self.config.input_reads)
        
        return cmd
    
    def _log_parameter_explanations(self, cmd: List[str]):
        """记录参数解释 | Log parameter explanations"""
        self.logger.info("="*60)
        self.logger.info("HiFiasm参数详解 | HiFiasm Parameter Explanations")
        self.logger.info("="*60)
        
        param_explanations = {
            '-o': '输出文件前缀 | Output file prefix',
            '-t': '线程数 | Number of threads',
            '--hg-size': '基因组大小估计，用于优化内存使用 | Genome size estimation for memory optimization',
            '-l': 'Purge级别(0-3)，控制重复序列处理 | Purge level (0-3) for repetitive sequence handling',
            '--purge-max': '最大purge覆盖度阈值 | Maximum purge coverage threshold',
            '-s': '相似性阈值，用于单倍型分离 | Similarity threshold for haplotype separation',
            '--ul': 'ONT长读长数据，用于辅助组装 | ONT long-read data for assembly assistance',
            '--h1': 'Hi-C第一端数据，用于染色体级别组装 | Hi-C first-end data for chromosome-level assembly',
            '--h2': 'Hi-C第二端数据，用于染色体级别组装 | Hi-C second-end data for chromosome-level assembly',
            '--trio': '三倍体模式组装 | Triploid mode assembly',
            '--polyploid': '多倍体模式组装 | Polyploid mode assembly'
        }
        
        i = 0
        while i < len(cmd):
            param = cmd[i]
            if param in param_explanations:
                desc = param_explanations[param]
                value_str = ""
                # 检查参数后面是否有值
                if i + 1 < len(cmd) and not str(cmd[i+1]).startswith('-'):
                    value_str = f" {cmd[i+1]}"
                    i += 1
                self.logger.info(f"  {param}{value_str}")
                self.logger.info(f"    💡 {desc}")
            i += 1
        
        self.logger.info("="*60)
    
    def _execute_hifiasm(self, cmd: List[str]) -> bool:
        """执行HiFiasm命令 | Execute HiFiasm command"""
        try:
            self.logger.info("开始HiFiasm组装，这可能需要数小时... | Starting HiFiasm assembly, this may take several hours...")
            
            start_time = time.time()
            
            # 执行命令并监控进度
            result = self.cmd_runner.run_with_progress(
                cmd=cmd,
                description="HiFiasm基因组组装 | HiFiasm genome assembly",
                progress_pattern="[M::"  # HiFiasm的进度输出模式
            )
            
            end_time = time.time()
            assembly_time = end_time - start_time
            
            self.logger.info(f"HiFiasm组装耗时: {assembly_time/3600:.1f}小时 | HiFiasm assembly time: {assembly_time/3600:.1f} hours")
            
            return result.returncode == 0
            
        except Exception as e:
            self.logger.error(f"HiFiasm执行异常 | HiFiasm execution exception: {e}")
            return False
    
    def _verify_assembly_outputs(self) -> bool:
        """验证组装输出文件 | Verify assembly output files"""
        self.logger.info("验证HiFiasm输出文件 | Verifying HiFiasm output files")
        
        expected_files = self._get_expected_output_files()
        missing_files = []
        existing_files = []
        
        for file_type, file_path in expected_files.items():
            full_path = self.assembly_dir / file_path
            if full_path.exists() and full_path.stat().st_size > 0:
                existing_files.append((file_type, full_path))
                file_size = full_path.stat().st_size / (1024*1024)  # MB
                self.logger.info(f"✓ {file_type}: {full_path.name} ({file_size:.1f} MB)")
            else:
                missing_files.append((file_type, full_path))
        
        if missing_files:
            self.logger.warning("以下预期文件未找到 | Following expected files not found:")
            for file_type, file_path in missing_files:
                self.logger.warning(f"  ✗ {file_type}: {file_path.name}")
        
        # 至少需要有primary contigs
        primary_files = [f for f in existing_files if 'primary' in f[0].lower() or 'p_ctg' in str(f[1])]
        if not primary_files:
            self.logger.error("未找到primary assembly文件 | Primary assembly files not found")
            return False
        
        self.logger.success(f"找到 {len(existing_files)} 个有效输出文件 | Found {len(existing_files)} valid output files")
        return True
    
    def _get_expected_output_files(self) -> Dict[str, str]:
        """获取预期输出文件列表 | Get expected output files list"""
        files = {
            'Primary contigs (GFA)': f'{self.config.prefix}.bp.p_ctg.gfa',
            'Alternate contigs (GFA)': f'{self.config.prefix}.bp.a_ctg.gfa',
            'Assembly log': f'{self.config.prefix}.log'
        }
        
        # 根据组装类型添加特定文件
        if self.config.assembly_type == 'triploid':
            files.update({
                'Haplotype 1 (GFA)': f'{self.config.prefix}.bp.hap1.p_ctg.gfa',
                'Haplotype 2 (GFA)': f'{self.config.prefix}.bp.hap2.p_ctg.gfa'
            })
        
        # 如果有Hi-C数据，会生成scaffold文件
        if self.config.hi_c_1 and self.config.hi_c_2:
            files.update({
                'Scaffolds (GFA)': f'{self.config.prefix}.hic.p_ctg.gfa'
            })
        
        return files
    
    def get_assembly_files(self) -> Dict[str, Path]:
        """获取组装文件路径 | Get assembly file paths"""
        files = {}
        
        # 查找所有GFA文件
        for gfa_file in self.assembly_dir.glob(f'{self.config.prefix}*.gfa'):
            file_type = self._classify_gfa_file(gfa_file.name)
            files[file_type] = gfa_file
        
        return files
    
    def _classify_gfa_file(self, filename: str) -> str:
        """分类GFA文件类型 | Classify GFA file type"""
        if '.bp.p_ctg.gfa' in filename:
            return 'primary_contigs'
        elif '.bp.a_ctg.gfa' in filename:
            return 'alternate_contigs'
        elif '.hic.p_ctg.gfa' in filename:
            return 'hi_c_scaffolds'
        elif 'hap1' in filename:
            return 'haplotype1'
        elif 'hap2' in filename:
            return 'haplotype2'
        else:
            return 'unknown'

class GFAConverter:
    """GFA格式转换器 | GFA Format Converter"""
    
    def __init__(self, config, logger: logging.Logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.assembly_dir = Path(config.output_dir) / 'assembly'
        self.fasta_dir = Path(config.output_dir) / 'final_results'
        self.fasta_dir.mkdir(exist_ok=True)
    
    def convert_gfa_to_fasta(self) -> bool:
        """转换GFA文件到FASTA格式 | Convert GFA files to FASTA format"""
        try:
            self.logger.info("转换GFA文件到FASTA格式 | Converting GFA files to FASTA format")
            
            # 获取组装文件
            assembler = HifiasmAssembler(self.config, self.logger, self.cmd_runner)
            assembly_files = assembler.get_assembly_files()
            
            if not assembly_files:
                self.logger.error("未找到GFA文件 | No GFA files found")
                return False
            
            converted_files = {}
            
            for file_type, gfa_path in assembly_files.items():
                fasta_path = self._convert_single_gfa(gfa_path, file_type)
                if fasta_path:
                    converted_files[file_type] = fasta_path
            
            # 生成合并的primary assembly
            if 'primary_contigs' in converted_files:
                self._create_primary_assembly(converted_files)
            
            self.logger.success(f"成功转换 {len(converted_files)} 个GFA文件 | Successfully converted {len(converted_files)} GFA files")
            return True
            
        except Exception as e:
            self.logger.error(f"GFA转换失败 | GFA conversion failed: {e}")
            return False
    
    def _convert_single_gfa(self, gfa_path: Path, file_type: str) -> Optional[Path]:
        """转换单个GFA文件 | Convert single GFA file"""
        try:
            self.logger.info(f"转换 {file_type}: {gfa_path.name}")
            
            # 生成输出文件名
            fasta_name = gfa_path.stem.replace('.bp', '').replace('.hic', '') + '.fasta'
            fasta_path = self.fasta_dir / fasta_name
            
            # 使用awk命令转换GFA到FASTA
            cmd = [
                'awk',
                '/^S/{print ">"$2"\\n"$3}',
                str(gfa_path)
            ]
            
            result = self.cmd_runner.run(
                cmd=cmd,
                description=f"转换GFA到FASTA | Convert GFA to FASTA: {gfa_path.name}"
            )
            
            # 将结果写入文件
            with open(fasta_path, 'w') as f:
                f.write(result.stdout)
            
            # 验证输出文件
            if fasta_path.exists() and fasta_path.stat().st_size > 0:
                file_size = fasta_path.stat().st_size / (1024*1024)
                self.logger.info(f"✓ 生成FASTA文件: {fasta_path.name} ({file_size:.1f} MB)")
                return fasta_path
            else:
                self.logger.error(f"FASTA文件生成失败 | FASTA file generation failed: {fasta_path}")
                return None
                
        except Exception as e:
            self.logger.error(f"转换GFA文件失败 | GFA file conversion failed: {e}")
            return None
    
    def _create_primary_assembly(self, converted_files: Dict[str, Path]):
        """创建主要组装文件 | Create primary assembly file"""
        try:
            primary_path = self.fasta_dir / f"{self.config.prefix}.primary.fasta"
            
            # 如果有primary_contigs文件，直接复制
            if 'primary_contigs' in converted_files:
                import shutil
                shutil.copy2(converted_files['primary_contigs'], primary_path)
                self.logger.info(f"创建主要组装文件 | Created primary assembly file: {primary_path.name}")
            
            # 如果需要，也可以创建单倍型特异的组装
            if self.config.analyze_haplotypes:
                self._create_haplotype_assemblies(converted_files)
                
        except Exception as e:
            self.logger.error(f"创建主要组装文件失败 | Primary assembly file creation failed: {e}")
    
    def _create_haplotype_assemblies(self, converted_files: Dict[str, Path]):
        """创建单倍型组装文件 | Create haplotype assembly files"""
        try:
            for hap_type in ['haplotype1', 'haplotype2']:
                if hap_type in converted_files:
                    hap_path = self.fasta_dir / f"{self.config.prefix}.{hap_type}.fasta"
                    import shutil
                    shutil.copy2(converted_files[hap_type], hap_path)
                    self.logger.info(f"创建单倍型文件 | Created haplotype file: {hap_path.name}")
                    
        except Exception as e:
            self.logger.error(f"创建单倍型文件失败 | Haplotype file creation failed: {e}")
    
    def get_converted_files(self) -> Dict[str, Path]:
        """获取转换后的FASTA文件 | Get converted FASTA files"""
        fasta_files = {}
        
        for fasta_file in self.fasta_dir.glob('*.fasta'):
            if 'primary' in fasta_file.name:
                fasta_files['primary'] = fasta_file
            elif 'haplotype1' in fasta_file.name or 'hap1' in fasta_file.name:
                fasta_files['haplotype1'] = fasta_file
            elif 'haplotype2' in fasta_file.name or 'hap2' in fasta_file.name:
                fasta_files['haplotype2'] = fasta_file
            elif 'alternate' in fasta_file.name:
                fasta_files['alternate'] = fasta_file
        
        return fasta_files