"""
PlotSR核心流程模块|PlotSR Core Pipeline Module
"""

import os
import subprocess
from typing import List, Tuple
import logging


class PlotSRPipeline:
    """PlotSR流程执行类|PlotSR Pipeline Executor"""

    def __init__(self, config, logger):
        """
        初始化流程|Initialize pipeline

        Args:
            config: PlotSRConfig配置对象|PlotSRConfig object
            logger: 日志器|Logger
        """
        self.config = config
        self.logger = logger
        self.original_genomes = config.genomes.copy()  # 保存原始基因组路径|Save original genome paths
        self.temp_genomes = []  # 临时基因组文件|Temporary genome files

    def run(self) -> bool:
        """
        运行完整流程|Run complete pipeline

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始PlotSR流程|Starting PlotSR pipeline")
        self.logger.info("=" * 60)

        # 步骤1：检查依赖|Step 1: Check dependencies
        if not self._check_dependencies():
            return False

        # 步骤2：如果指定了染色体，先提取染色体序列|Step 2: If chromosomes specified, extract first
        # 注意：由于PlotSR对提取的FASTA文件兼容性问题，暂时禁用此优化
        # Note: Disabled due to PlotSR compatibility issues with extracted FASTA files
        # if self.config.chromosomes:
        #     self.logger.info("步骤0: 提取指定染色体序列|Step 0: Extracting specified chromosomes")
        #     if not self._extract_chromosomes():
        #         return False
        if self.config.chromosomes:
            self.logger.info("步骤0: 指定染色体|Step 0: Chromosomes specified (will use --chr filter in PlotSR)")
            resolved_chrs = self._resolve_chromosome_names_for_extraction(self.config.chromosomes)
            self.logger.info(f"将分析染色体|Will analyze chromosomes: {', '.join(resolved_chrs)}")

        # 步骤3：提取染色体长度|Step 3: Extract chromosome lengths
        self.logger.info("步骤1: 提取染色体长度|Step 1: Extracting chromosome lengths")
        chr_length_files = self._extract_chromosome_lengths()

        # 步骤4：基因组比对|Step 4: Genome alignment
        self.logger.info("步骤2: 基因组比对|Step 2: Genome alignment")
        bam_files = self._run_alignments()

        # 步骤5：SyRI结构注释|Step 5: SyRI structural annotation
        self.logger.info("步骤3: SyRI结构注释|Step 3: SyRI structural annotation")
        syri_files = self._run_syri(bam_files)

        # 检查SyRI是否成功|Check if SyRI succeeded
        if not syri_files:
            self.logger.error("SyRI分析失败，流程终止|SyRI analysis failed, pipeline terminated")
            self.logger.error("建议：检查基因组染色体方向是否一致，或使用RectChr等工具校正|Suggestion: check chromosome orientation consistency or use RectChr for correction")
            return False

        # 步骤6：生成PlotSR配置|Step 6: Generate PlotSR configuration
        self.logger.info("步骤4: 生成PlotSR配置|Step 4: Generating PlotSR configuration")
        self._generate_plotsr_config(chr_length_files, syri_files)

        # 步骤7：运行PlotSR可视化|Step 7: Run PlotSR visualization
        self.logger.info("步骤5: PlotSR可视化|Step 5: PlotSR visualization")
        success = self._run_plotsr(syri_files)

        if success:
            self.logger.info("=" * 60)
            self.logger.info("PlotSR流程完成|PlotSR pipeline completed")
            self.logger.info(f"结果保存在|Results saved to: {self.config.output_dir}")
            self.logger.info("=" * 60)

        return success

    def _check_dependencies(self) -> bool:
        """
        检查依赖工具|Check dependency tools

        Returns:
            bool: 是否全部可用|Whether all available
        """
        from .utils import check_dependencies

        self.logger.info("检查依赖工具|Checking dependency tools...")
        all_available, missing = check_dependencies()

        if not all_available:
            self.logger.error(
                f"缺少依赖工具|Missing required tools: {', '.join(missing)}"
            )
            return False

        self.logger.info("所有依赖工具可用|All dependency tools available")
        return True

    def _extract_chromosomes(self) -> bool:
        """
        从基因组中提取指定染色体|Extract specified chromosomes from genomes

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            # 先解析染色体名称（将数字转换为名称）|Resolve chromosome names first
            resolved_chrs = self._resolve_chromosome_names_for_extraction(self.config.chromosomes)
            self.logger.info(f"提取染色体|Extracting chromosomes: {', '.join(resolved_chrs)}")

            # 为每个基因组提取指定染色体|Extract specified chromosomes for each genome
            new_genomes = []
            for i, (genome_file, name) in enumerate(zip(self.config.genomes, self.config.names)):
                self.logger.info(f"处理基因组|Processing genome [{i+1}]: {name}")

                # 输出文件|Output file
                extracted_fa = os.path.join(
                    self.config.output_dir,
                    'plotsr',
                    f"{name}_extracted.fa"
                )

                # 检查是否已存在|Check if already exists
                if self.config.skip_existing and os.path.exists(extracted_fa) and os.path.getsize(extracted_fa) > 0:
                    self.logger.info(f"跳过|Skipping: 提取文件已存在|Extracted file exists")
                    new_genomes.append(extracted_fa)
                    self.temp_genomes.append(extracted_fa)
                    continue

                # 使用samtools提取染色体|Use samtools to extract chromosomes
                chr_list = ','.join(resolved_chrs)
                cmd = ['samtools', 'faidx', genome_file, chr_list]

                self.logger.debug(f"运行|Running: {' '.join(cmd)}")

                with open(extracted_fa, 'w') as f_out:
                    result = subprocess.run(
                        cmd,
                        stdout=f_out,
                        stderr=subprocess.PIPE,
                        text=True,
                        check=True
                    )

                # 验证文件|Verify file
                if not os.path.exists(extracted_fa) or os.path.getsize(extracted_fa) == 0:
                    self.logger.error(f"染色体提取失败|Chromosome extraction failed: {name}")
                    return False

                new_genomes.append(extracted_fa)
                self.temp_genomes.append(extracted_fa)
                self.logger.info(f"染色体已提取|Chromosomes extracted: {extracted_fa}")

            # 更新config中的基因组路径|Update genome paths in config
            self.config.genomes = new_genomes

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"染色体提取失败|Chromosome extraction failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"染色体提取异常|Chromosome extraction error: {e}")
            return False

    def _resolve_chromosome_names_for_extraction(self, chr_list: List[str]) -> List[str]:
        """
        解析染色体名称用于提取（支持数字索引）|Resolve chromosome names for extraction (support numeric index)

        Args:
            chr_list: 染色体列表（数字或名称）|List of chromosomes (numbers or names)

        Returns:
            list: 解析后的染色体名称列表|List of resolved chromosome names
        """
        if not chr_list:
            return None

        # 读取第一个基因组的.fai文件|Read first genome .fai file
        first_genome = self.original_genomes[0]
        fai_file = f"{first_genome}.fai"

        # 确保.fai文件存在|Ensure .fai file exists
        if not os.path.exists(fai_file):
            self.logger.warning(f"索引文件不存在，创建索引|Index file not found, creating index: {fai_file}")
            try:
                subprocess.run(['samtools', 'faidx', first_genome], check=True, capture_output=True)
            except subprocess.CalledProcessError as e:
                self.logger.error(f"创建索引失败|Failed to create index: {e}")
                return chr_list

        # 读取所有染色体名称|Read all chromosome names
        chr_names = []
        with open(fai_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 1:
                    chr_names.append(fields[0])

        resolved = []
        for item in chr_list:
            # 检查是否为数字|Check if it's a number
            try:
                chr_idx = int(item) - 1  # 转换为0-based索引|Convert to 0-based index
                if 0 <= chr_idx < len(chr_names):
                    resolved.append(chr_names[chr_idx])
                else:
                    self.logger.warning(f"染色体索引超出范围|Chromosome index out of range: {item} (1-{len(chr_names)})")
                    resolved.append(item)  # 保留原值|Keep original value
            except ValueError:
                # 不是数字，直接使用|Not a number, use directly
                resolved.append(item)

        return resolved

    def _extract_chromosome_lengths(self) -> List[str]:
        """
        提取染色体长度|Extract chromosome lengths

        Returns:
            list: 染色体长度文件列表|List of chromosome length files
        """
        from .utils import extract_chromosome_lengths

        chr_length_files = []

        for i, (genome, name) in enumerate(zip(self.config.genomes, self.config.names)):
            output_file = os.path.join(
                self.config.output_dir,
                'plotsr',
                f"{name}.chrlen"
            )

            # 检查是否已存在|Check if already exists
            if self.config.skip_existing and os.path.exists(output_file) and os.path.getsize(output_file) > 0:
                self.logger.info(f"跳过|Skipping [{i+1}]: {name} (文件已存在|File exists)")
                chr_length_files.append(output_file)
                continue

            self.logger.info(f"处理基因组|Processing genome [{i+1}]: {name}")
            extract_chromosome_lengths(genome, output_file)
            chr_length_files.append(output_file)
            self.logger.info(f"染色体长度已保存|Chromosome lengths saved: {output_file}")

        return chr_length_files

    def _run_alignments(self) -> List[str]:
        """
        运行minimap2比对|Run minimap2 alignments

        Returns:
            list: BAM文件列表|List of BAM files
        """
        bam_files = []

        # 按顺序比对相邻基因组|Align adjacent genomes in order
        for i in range(len(self.config.genomes) - 1):
            ref_genome = self.config.genomes[i]
            query_genome = self.config.genomes[i + 1]
            ref_name = self.config.names[i]
            query_name = self.config.names[i + 1]

            # 输出文件名|Output file name
            bam_file = os.path.join(
                self.config.output_dir,
                'alignment',
                f"{ref_name}_vs_{query_name}.bam"
            )
            bam_index_file = f"{bam_file}.bai"

            # 检查是否已存在|Check if already exists
            if self.config.skip_existing and os.path.exists(bam_file) and os.path.exists(bam_index_file):
                if os.path.getsize(bam_file) > 0 and os.path.getsize(bam_index_file) > 0:
                    self.logger.info(f"跳过|Skipping: {ref_name} vs {query_name} (BAM文件已存在|BAM exists)")
                    bam_files.append(bam_file)
                    continue

            self.logger.info(f"比对|Aligning: {ref_name} vs {query_name}")

            # 运行minimap2|Run minimap2
            if self._run_minimap2(ref_genome, query_genome, bam_file):
                bam_files.append(bam_file)
                self.logger.info(f"比对完成|Alignment completed: {bam_file}")
            else:
                self.logger.error(f"比对失败|Alignment failed: {ref_name} vs {query_name}")
                return []

        return bam_files

    def _run_minimap2(self, ref: str, query: str, output_bam: str) -> bool:
        """
        运行minimap2比对|Run minimap2 alignment

        Args:
            ref: 参考基因组|Reference genome
            query: 查询基因组|Query genome
            output_bam: 输出BAM文件|Output BAM file

        Returns:
            bool: 是否成功|Whether successful
        """
        threads = self.config.threads
        preset = self.config.minimap2_preset

        try:
            # 使用管道方式：minimap2 | samtools sort | samtools index
            # Use pipeline: minimap2 | samtools sort | samtools index
            cmd = f"""
                minimap2 -x {preset} -t {threads} -a --eqx {ref} {query} |
                samtools sort -@ {threads} -O BAM -o {output_bam} -
            """

            self.logger.debug(f"运行|Running: minimap2 -x {preset} {ref} {query}")

            result = subprocess.run(
                cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )

            if result.stderr:
                # 检查是否有警告|Check for warnings
                for line in result.stderr.split('\n'):
                    if line.strip() and not line.strip().startswith('[W::hts_set_opt]'):
                        self.logger.debug(f"minimap2|minimap2: {line}")

            # 索引BAM文件|Index BAM file
            self.logger.debug(f"索引BAM文件|Indexing BAM file: {output_bam}")

            subprocess.run(
                ['samtools', 'index', '-@', str(threads), output_bam],
                check=True,
                capture_output=True
            )

            # 验证BAM文件|Verify BAM file
            if not os.path.exists(output_bam) or os.path.getsize(output_bam) == 0:
                self.logger.error(f"BAM文件生成失败|BAM file generation failed: {output_bam}")
                return False

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"minimap2执行失败|minimap2 execution failed: {e.stderr if hasattr(e, 'stderr') else str(e)}")
            return False
        except Exception as e:
            self.logger.error(f"比对异常|Alignment error: {e}")
            return False

    def _run_syri(self, bam_files: List[str]) -> List[str]:
        """
        运行SyRI结构注释|Run SyRI structural annotation

        Args:
            bam_files: BAM文件列表|List of BAM files

        Returns:
            list: SyRI输出文件列表|List of SyRI output files
        """
        syri_files = []

        for i, bam_file in enumerate(bam_files):
            ref_name = self.config.names[i]
            query_name = self.config.names[i + 1]
            ref_genome = self.config.genomes[i]
            query_genome = self.config.genomes[i + 1]

            # 输出目录和前缀|Output directory and prefix
            syri_dir = os.path.join(self.config.output_dir, 'syri')
            prefix_name = f"{ref_name}_vs_{query_name}"

            # 确定输出文件|Determine output file
            if self.config.syri_filter:
                syri_out = os.path.join(syri_dir, f"{prefix_name}syri.filtered.out")
            else:
                syri_out = os.path.join(syri_dir, f"{prefix_name}syri.out")

            # 检查是否已存在|Check if already exists
            if self.config.skip_existing and os.path.exists(syri_out) and os.path.getsize(syri_out) > 0:
                self.logger.info(f"跳过|Skipping: {ref_name} vs {query_name} (SyRI文件已存在|SyRI exists)")
                syri_files.append(syri_out)
                continue

            self.logger.info(f"SyRI分析|SyRI analysis: {ref_name} vs {query_name}")

            # 运行SyRI|Run SyRI
            if self._run_syri_single(bam_file, ref_genome, query_genome, syri_dir, prefix_name):
                raw_syri_out = os.path.join(syri_dir, f"{prefix_name}syri.out")

                # 可选：过滤|Optional: filtering
                if self.config.syri_filter:
                    filtered_out = os.path.join(syri_dir, f"{prefix_name}syri.filtered.out")
                    self._filter_syri_output(raw_syri_out, filtered_out)
                    syri_files.append(filtered_out)
                else:
                    syri_files.append(raw_syri_out)

                self.logger.info(f"SyRI完成|SyRI completed: {syri_files[-1]}")
            else:
                self.logger.error(f"SyRI失败|SyRI failed: {ref_name} vs {query_name}")
                return []

        return syri_files

    def _run_syri_single(self, bam: str, ref: str, query: str, syri_dir: str, prefix_name: str) -> bool:
        """
        运行单个SyRI分析|Run single SyRI analysis

        Args:
            bam: BAM文件|BAM file
            ref: 参考基因组|Reference genome
            query: 查询基因组|Query genome
            syri_dir: 输出目录|Output directory
            prefix_name: 文件名前缀|Filename prefix

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            cmd = [
                'syri',
                '-c', bam,
                '-r', ref,
                '-q', query,
                '-F', 'B',
                '--dir', syri_dir,
                '--prefix', prefix_name
            ]

            self.logger.debug(f"运行|Running: {' '.join(cmd)}")

            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            # 检查输出文件是否生成|Check if output file was generated
            syri_out = os.path.join(syri_dir, f"{prefix_name}syri.out")
            if not os.path.exists(syri_out):
                self.logger.error(f"SyRI输出文件未生成|SyRI output file not generated: {syri_out}")
                self.logger.error(f"可能原因：染色体链方向不一致或共线性区域太少|Possible reason: chromosome strand mismatch or insufficient syntenic regions")
                return False

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"SyRI执行失败|SyRI execution failed: {e.stderr}")
            return False

    def _filter_syri_output(self, input_file: str, output_file: str):
        """
        过滤SyRI输出|Filter SyRI output

        Args:
            input_file: 输入文件|Input file
            output_file: 输出文件|Output file
        """
        # 读取并过滤|Read and filter
        with open(input_file, 'r') as f_in:
            lines = f_in.readlines()

        # 写入过滤后的结果|Write filtered results
        with open(output_file, 'w') as f_out:
            for line in lines:
                # 移除注释行|Remove comment lines
                if not line.startswith('#'):
                    f_out.write(line)

    def _generate_plotsr_config(self, chr_length_files: List[str], syri_files: List[str]):
        """
        生成PlotSR配置文件|Generate PlotSR configuration files

        Args:
            chr_length_files: 染色体长度文件列表|List of chromosome length files
            syri_files: SyRI输出文件列表|List of SyRI output files
        """
        # 生成genomes.txt|Generate genomes.txt
        genomes_txt = os.path.join(self.config.output_dir, 'plotsr', 'genomes.txt')

        # 使用原始FASTA文件而不是chrlen文件，避免PlotSR读取chrlen文件时的问题
        # Use original FASTA files instead of chrlen files to avoid PlotSR reading issues
        with open(genomes_txt, 'w') as f:
            for genome_file, name in zip(self.config.genomes, self.config.names):
                f.write(f"{genome_file}\t{name}\tlw:1.5\n")

        self.logger.info(f"配置文件已生成|Config file generated: {genomes_txt}")

    def _resolve_chromosome_names(self, chr_list: List[str]) -> List[str]:
        """
        解析染色体名称（支持数字索引）|Resolve chromosome names (support numeric index)

        Args:
            chr_list: 染色体列表（数字或名称）|List of chromosomes (numbers or names)

        Returns:
            list: 解析后的染色体名称列表|List of resolved chromosome names
        """
        if not chr_list:
            return None

        # 读取第一个基因组的染色体信息|Read first genome chromosome info
        first_genome_chrlen = os.path.join(
            self.config.output_dir,
            'plotsr',
            f"{self.config.names[0]}.chrlen"
        )

        if not os.path.exists(first_genome_chrlen):
            self.logger.warning(f"染色体长度文件不存在|Chromosome length file not found: {first_genome_chrlen}")
            return chr_list

        # 读取所有染色体名称|Read all chromosome names
        chr_names = []
        with open(first_genome_chrlen, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) >= 1:
                    chr_names.append(fields[0])

        resolved = []
        for item in chr_list:
            # 检查是否为数字|Check if it's a number
            try:
                chr_idx = int(item) - 1  # 转换为0-based索引|Convert to 0-based index
                if 0 <= chr_idx < len(chr_names):
                    resolved.append(chr_names[chr_idx])
                else:
                    self.logger.warning(f"染色体索引超出范围|Chromosome index out of range: {item} (1-{len(chr_names)})")
                    resolved.append(item)  # 保留原值|Keep original value
            except ValueError:
                # 不是数字，直接使用|Not a number, use directly
                resolved.append(item)

        return resolved

    def _run_plotsr(self, syri_files: List[str]) -> bool:
        """
        运行PlotSR可视化|Run PlotSR visualization

        Args:
            syri_files: SyRI输出文件列表|List of SyRI output files

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            # 检查输出文件是否已存在|Check if output file already exists
            plot_output = os.path.join(self.config.output_dir, f'plot.{self.config.output_format}')

            if self.config.skip_existing and os.path.exists(plot_output) and os.path.getsize(plot_output) > 0:
                self.logger.info(f"跳过|Skipping: PlotSR可视化已完成 (文件已存在|File exists: {plot_output})")
                return True

            # 查找所有SyRI输出文件|Find all SyRI output files
            import glob

            syri_pattern = os.path.join(self.config.output_dir, 'syri', '*syri.filtered.out')
            syri_file_list = glob.glob(syri_pattern)

            if not syri_file_list:
                self.logger.error(f"未找到SyRI输出文件|No SyRI output files found: {syri_pattern}")
                return False

            # 按照基因组顺序排序SyRI文件|Sort SyRI files according to genome order
            # 比对顺序是 genome[0]_vs_genome[1], genome[1]_vs_genome[2], ...
            # Alignment order is genome[0]_vs_genome[1], genome[1]_vs_genome[2], ...
            syri_file_list_sorted = []
            for i in range(len(self.config.genomes) - 1):
                ref_name = self.config.names[i]
                query_name = self.config.names[i + 1]
                expected_file = os.path.join(
                    self.config.output_dir,
                    'syri',
                    f"{ref_name}_vs_{query_name}syri.filtered.out"
                )
                if expected_file in syri_file_list:
                    syri_file_list_sorted.append(expected_file)
                else:
                    self.logger.error(f"未找到期望的SyRI文件|Expected SyRI file not found: {expected_file}")

            if len(syri_file_list_sorted) != len(syri_file_list):
                self.logger.warning(f"SyRI文件数量不匹配|SyRI file count mismatch: expected {len(syri_file_list)}, found {len(syri_file_list_sorted)}")

            syri_file_list = syri_file_list_sorted

            self.logger.info(f"找到|Found {len(syri_file_list)} 个SyRI文件|SyRI files")

            # 构建PlotSR命令|Build PlotSR command
            cmd = [
                'plotsr',
                '--genomes', os.path.join(self.config.output_dir, 'plotsr', 'genomes.txt'),
                '-o', plot_output,
                '-s', str(self.config.min_sr_size),
                '-f', str(self.config.font_size),
                '-d', str(self.config.dpi),
                '-S', str(self.config.space_ratio)
            ]

            # 添加每个SyRI文件|Add each SyRI file
            for syri_file in syri_file_list:
                cmd.extend(['--sr', syri_file])

            # 添加染色体过滤参数|Add chromosome filter parameters
            if self.config.chromosomes:
                resolved_chrs = self._resolve_chromosome_names_for_extraction(self.config.chromosomes)
                if resolved_chrs:
                    self.logger.info(f"指定染色体|Specified chromosomes: {', '.join(resolved_chrs)}")
                    for chr_name in resolved_chrs:
                        cmd.extend(['--chr', chr_name])

            if self.config.vertical:
                cmd.append('-v')
            if self.config.itx:
                cmd.append('--itx')
            if self.config.nosyn:
                cmd.append('--nosyn')
            if self.config.noinv:
                cmd.append('--noinv')
            if self.config.notr:
                cmd.append('--notr')
            if self.config.nodup:
                cmd.append('--nodup')

            self.logger.info(f"运行PlotSR|Running PlotSR...")
            self.logger.debug(f"命令|Command: {' '.join(cmd)}")

            subprocess.run(
                cmd,
                check=True,
                capture_output=True,
                text=True
            )

            self.logger.info("PlotSR可视化完成|PlotSR visualization completed")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"PlotSR执行失败|PlotSR execution failed: {e.stderr}")
            return False
        except Exception as e:
            self.logger.error(f"PlotSR异常|PlotSR error: {e}")
            return False
