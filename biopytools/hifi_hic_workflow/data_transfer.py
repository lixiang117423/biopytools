"""
数据传递管理模块|Data Transfer Management Module

负责各步骤之间的输入输出文件定位和传递
Responsible for input/output file location and transfer between steps
"""

import os
import glob
import logging
from typing import Optional, Tuple, Dict, Any
from pathlib import Path


class DataTransferManager:
    """数据传递管理器|Data Transfer Manager"""

    def __init__(self, config, logger: logging.Logger):
        """
        初始化数据传递管理器|Initialize data transfer manager

        Args:
            config: HifiHicWorkflowConfig配置对象|HifiHicWorkflowConfig object
            logger: 日志记录器|Logger
        """
        self.config = config
        self.logger = logger

    def get_hifi_hic_params(self) -> Dict[str, Any]:
        """
        获取hifi_hic模块的参数|Get hifi_hic module parameters

        Returns:
            dict: hifi_hic参数字典|hifi_hic parameter dictionary
        """
        params = {
            "hifi_data": self.config.hifi_reads,
            "hic_r1": self.config.hic_r1,
            "hic_r2": self.config.hic_r2,
            "prefix": self.config.prefix,
            "threads": self.config.threads,
            "genome_size": self.config.genome_size,
            "n_hap": self.config.n_hap,
            "base_dir": self.config.hifi_hic_output_dir,
            "resume": self.config.resume,
        }

        # 可选参数|Optional parameters
        if self.config.purge_level is not None:
            params["purge_level"] = self.config.purge_level

        if self.config.hom_cov is not None:
            params["hom_cov"] = self.config.hom_cov

        # NGS polish参数|NGS polish parameters
        if self.config.use_ngs_polish and self.config.ngs_data:
            params["ngs_data"] = self.config.ngs_data
            params["high_cov"] = self.config.ngs_high_cov
            params["ngs_pattern"] = self.config.ngs_pattern

        self.logger.info(f"hifi_hic参数|hifi_hic parameters: base_dir={params['base_dir']}, use_ngs={self.config.use_ngs_polish}")
        return params

    def get_hifi_hic_output(self) -> Optional[str]:
        """
        获取hifi_hic的输出FASTA文件|Get hifi_hic output FASTA file

        Returns:
            str: primary.fa文件路径|Path to primary.fa file, or None if not found
        """
        self.logger.info("定位hifi_hic输出文件|Locating hifi_hic output file")

        # 根据是否使用NGS polish选择不同的路径|Choose different path based on NGS polish
        if self.config.use_ngs_polish and self.config.ngs_data:
            # NGS polish后的输出路径|Output path after NGS polish
            polished_fa = os.path.join(
                self.config.hifi_hic_output_dir,
                self.config.prefix,
                "03.ngs_polish",
                "04.reassembly",
                "02.fasta",
                f"{self.config.prefix}.primary.fa"
            )

            if os.path.exists(polished_fa):
                self.logger.info(f"找到NGS polish后的基因组|Found NGS polished genome: {polished_fa}")
                return polished_fa
            else:
                self.logger.warning(f"NGS polish输出不存在|NGS polish output not found: {polished_fa}")

        # 标准输出路径（无NGS或NGS polish失败）|Standard output path (no NGS or NGS polish failed)
        primary_fa = os.path.join(
            self.config.hifi_hic_output_dir,
            self.config.prefix,
            "02.fasta",
            f"{self.config.prefix}.primary.fa"
        )

        if os.path.exists(primary_fa):
            self.logger.info(f"找到HiFi组装基因组|Found HiFi assembled genome: {primary_fa}")
            return primary_fa
        else:
            self.logger.error(f"未找到HiFi组装输出|HiFi assembly output not found: {primary_fa}")
            return None

    def get_haphic_params(self, asm_file: str) -> Dict[str, Any]:
        """
        获取haphic模块的参数|Get haphic module parameters

        Args:
            asm_file: 组装文件路径|Assembly file path

        Returns:
            dict: haphic参数字典|haphic parameter dictionary
        """
        params = {
            "asm_file": asm_file,
            "hic_file": self.config.hic_r1,
            "hic_file_type": "fastq",  # 使用FASTQ模式|Use FASTQ mode
            "nchrs": self.config.nchrs,
            "prefix": self.config.prefix,
            "output_dir": self.config.haphic_output_dir,
            "threads": self.config.haphic_threads,
            "memory_limit": self.config.memory_limit,
            "force_rerun": self.config.force_rerun,
        }

        # 工具路径|Tool paths
        if self.config.haphic_bin != "haphic":
            params["haphic_bin"] = self.config.haphic_bin

        if self.config.bwa_bin != "bwa":
            params["bwa_bin"] = self.config.bwa_bin

        if self.config.samtools_bin != "samtools":
            params["samtools_bin"] = self.config.samtools_bin

        # Hi-C数据处理参数|Hi-C data processing parameters
        params["mapq_threshold"] = self.config.mapq_threshold
        params["edit_distance"] = self.config.edit_distance
        params["min_re_sites"] = self.config.min_re_sites

        # 聚类参数|Clustering parameters
        params["min_inflation"] = self.config.min_inflation
        params["max_inflation"] = self.config.max_inflation
        params["inflation_step"] = self.config.inflation_step
        params["nx"] = self.config.nx
        params["min_group_len"] = self.config.min_group_len

        # 排序和定向参数|Ordering and orientation parameters
        params["processes"] = self.config.processes
        params["fast_sorting"] = self.config.fast_sorting
        params["allhic_optimization"] = self.config.allhic_optimization

        # 组装校正参数|Assembly correction parameters
        params["correct_nrounds"] = self.config.correct_nrounds
        params["correct_min_coverage"] = self.config.correct_min_coverage

        if self.config.correct_resolution is not None:
            params["correct_resolution"] = self.config.correct_resolution

        # 高级选项|Advanced options
        if self.config.quick_view:
            params["quick_view"] = True

        if self.config.re_sites != "GATC":
            params["re_sites"] = self.config.re_sites

        if self.config.skip_clustering:
            params["skip_clustering"] = True

        # 可视化参数|Visualization parameters
        if self.config.generate_plots:
            params["generate_plots"] = True
            params["bin_size"] = self.config.bin_size
            params["min_len"] = self.config.min_len
            params["separate_plots"] = self.config.separate_plots

        self.logger.info(f"haphic参数|haphic parameters: nchrs={params['nchrs']}, output_dir={params['output_dir']}")
        return params

    def get_haphic_output(self) -> Optional[str]:
        """
        获取haphic的输出FASTA文件|Get haphic output FASTA file

        Returns:
            str: scaffold FASTA文件路径|Path to scaffold FASTA file, or None if not found
        """
        self.logger.info("定位haphic输出文件|Locating haphic output file")

        # HapHiC的最终输出在04.build目录|HapHiC final output is in 04.build directory
        scaffold_fa = os.path.join(
            self.config.haphic_output_dir,
            "04.build",
            f"{self.config.prefix}.fa"
        )

        if os.path.exists(scaffold_fa):
            self.logger.info(f"找到HapHiC挂载结果|Found HapHiC scaffolds: {scaffold_fa}")
            return scaffold_fa
        else:
            self.logger.error(f"未找到HapHiC输出|HapHiC output not found: {scaffold_fa}")
            return None

    def get_rename_params(self, input_fa: str) -> Dict[str, Any]:
        """
        获取rename_chromosomes模块的参数|Get rename_chromosomes module parameters

        Args:
            input_fa: 输入FASTA文件|Input FASTA file

        Returns:
            dict: rename参数字典|rename parameter dictionary
        """
        output_fa = os.path.join(
            self.config.rename_output_dir,
            f"{self.config.prefix}.renamed.fa"
        )

        params = {
            "input_file": input_fa,
            "output_file": output_fa,
            "chromosome_number": self.config.nchrs,
            "keep_all": self.config.rename_keep_all,
        }

        self.logger.info(f"rename_chromosomes参数|rename_chromosomes parameters: nchrs={params['chromosome_number']}, output={output_fa}")
        return params

    def get_rename_output(self) -> str:
        """
        获取rename_chromosomes的输出FASTA文件|Get rename_chromosomes output FASTA file

        Returns:
            str: 重命名后的FASTA文件路径|Path to renamed FASTA file
        """
        renamed_fa = os.path.join(
            self.config.rename_output_dir,
            f"{self.config.prefix}.renamed.fa"
        )

        self.logger.info(f"rename_chromosomes输出|rename_chromosomes output: {renamed_fa}")
        return renamed_fa

    def get_heatmap_params(self, genome_fa: str) -> Dict[str, Any]:
        """
        获取hic_heatmap模块的参数|Get hic_heatmap module parameters

        Args:
            genome_fa: 基因组FASTA文件|Genome FASTA file

        Returns:
            dict: hic_heatmap参数字典|hic_heatmap parameter dictionary
        """
        # 提取genome_id|Extract genome_id
        genome_id = self.config.prefix

        params = {
            "genome": genome_fa,
            "genome_id": genome_id,
            "fastq_r1": self.config.hic_r1,
            "fastq_r2": self.config.hic_r2,
            "output_dir": self.config.heatmap_output_dir,
            "threads": self.config.threads,
            "max_memory_gb": self.config.hicpro_max_memory_gb,
            "restriction_enzyme": self.config.hicpro_restriction_enzyme,
            "bin_sizes": self.config.hicpro_bin_sizes,
            "resolution": self.config.heatmap_resolution,
            "color_map": self.config.heatmap_color_map,
            "dpi": self.config.heatmap_dpi,
            "output_format": self.config.heatmap_format,
            "bar_max": self.config.heatmap_bar_max,
            "skip_existing": self.config.resume,
            "force_hicpro": self.config.force_rerun,
        }

        # 工具路径|Tool paths
        if self.config.hicpro_path:
            params["hicpro_path"] = self.config.hicpro_path

        if self.config.hicpro_sif:
            params["hicpro_sif"] = self.config.hicpro_sif

        if self.config.plothic_path:
            params["plothic_path"] = self.config.plothic_path

        self.logger.info(f" hic_heatmap参数| hic_heatmap parameters: genome_id={genome_id}, resolution={params['resolution']}")
        return params

    def get_heatmap_output(self) -> str:
        """
        获取hic_heatmap的输出文件|Get hic_heatmap output file

        Returns:
            str: 热图文件路径|Path to heatmap file
        """
        heatmap_file = os.path.join(
            self.config.heatmap_output_dir,
            "plot",
            f"{self.config.prefix}_hic_heatmap.{self.config.heatmap_format}"
        )

        self.logger.info(f" hic_heatmap输出| hic_heatmap output: {heatmap_file}")
        return heatmap_file

    def check_step_completion(self, step_name: str) -> bool:
        """
        检查指定步骤是否已完成|Check if specified step is completed

        Args:
            step_name: 步骤名称|Step name (hifi_hic, haphic, rename, heatmap)

        Returns:
            bool: 是否完成|Whether completed
        """
        if step_name == "hifi_hic":
            output_file = self.get_hifi_hic_output()
        elif step_name == "haphic":
            output_file = self.get_haphic_output()
        elif step_name == "rename":
            output_file = self.get_rename_output()
        elif step_name == "heatmap":
            output_file = self.get_heatmap_output()
        else:
            self.logger.warning(f"未知步骤名|Unknown step name: {step_name}")
            return False

        if output_file and os.path.exists(output_file):
            # 检查文件是否为空|Check if file is empty
            if os.path.getsize(output_file) > 0:
                self.logger.info(f"步骤{step_name}已完成|Step {step_name} completed: {output_file}")
                return True
            else:
                self.logger.warning(f"步骤{step_name}输出文件为空|Step {step_name} output file is empty: {output_file}")
                return False
        else:
            self.logger.info(f"步骤{step_name}未完成|Step {step_name} not completed")
            return False

    def create_step_info_file(self, step_name: str, info: Dict[str, Any]):
        """
        创建步骤信息文件|Create step info file

        Args:
            step_name: 步骤名称|Step name
            info: 步骤信息|Step information
        """
        step_info_dir = os.path.join(self.config.work_dir, "00.pipeline_info")
        Path(step_info_dir).mkdir(parents=True, exist_ok=True)

        info_file = os.path.join(step_info_dir, f"{step_name}_info.txt")

        with open(info_file, 'w') as f:
            f.write(f"步骤|Step: {step_name}\n")
            f.write(f"状态|Status: completed\n")
            for key, value in info.items():
                f.write(f"{key}: {value}\n")

        self.logger.info(f"步骤信息已保存|Step info saved: {info_file}")
