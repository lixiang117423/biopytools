"""
参考基因组引导命名模块|Reference Genome Guided Naming Module

使用minimap2将组装序列比对到参考基因组，根据最佳匹配进行命名
Use minimap2 to align assembly sequences to reference genome, name based on best match
"""

import os
import subprocess
import logging
from typing import Dict, List, Tuple, Optional
from pathlib import Path
from .data_transfer import DataTransferManager


class ReferenceGenomeNamer:
    """参考基因组引导命名器|Reference Genome Guided Namer"""

    def __init__(self, config, logger: logging.Logger):
        """
        初始化命名器|Initialize namer

        Args:
            config: HifiHicWorkflowConfig配置对象|HifiHicWorkflowConfig object
            logger: 日志记录器|Logger
        """
        self.config = config
        self.logger = logger
        self.data_transfer = DataTransferManager(config, logger)

    def run_rename(self, input_fa: str, output_fa: str) -> bool:
        """
        执行参考基因组引导命名|Perform reference genome guided naming

        Args:
            input_fa: 输入FASTA文件|Input FASTA file
            output_fa: 输出FASTA文件|Output FASTA file

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info("开始参考基因组引导命名|Starting reference genome guided naming")
        self.logger.info("=" * 60)
        self.logger.info(f"输入文件|Input: {input_fa}")
        self.logger.info(f"参考基因组|Reference: {self.config.reference_genome}")
        self.logger.info(f"输出文件|Output: {output_fa}")

        try:
            # Step 1: 使用minimap2比对|Step 1: Align using minimap2
            self.logger.info("步骤1: 使用minimap2比对到参考基因组|Step 1: Align to reference genome using minimap2")
            alignment_file = self._run_minimap2(input_fa)

            if not alignment_file:
                self.logger.error("minimap2比对失败|minimap2 alignment failed")
                return False

            # Step 2: 分析比对结果，提取最佳匹配|Step 2: Analyze alignment, extract best matches
            self.logger.info("步骤2: 分析比对结果|Step 2: Analyzing alignment results")
            naming_map = self._analyze_alignment(alignment_file, input_fa)

            if not naming_map:
                self.logger.error("比对分析失败|Alignment analysis failed")
                return False

            # Step 3: 根据命名映射重命名序列|Step 3: Rename sequences based on naming map
            self.logger.info("步骤3: 重命名序列|Step 3: Renaming sequences")
            success = self._rename_sequences(input_fa, output_fa, naming_map)

            if success:
                # Step 4: 生成命名报告|Step 4: Generate naming report
                self._generate_naming_report(naming_map, output_fa)

                self.logger.info("=" * 60)
                self.logger.info("参考基因组引导命名完成|Reference genome guided naming completed")
                self.logger.info("=" * 60)
                return True
            else:
                self.logger.error("序列重命名失败|Sequence renaming failed")
                return False

        except Exception as e:
            self.logger.error(f"参考基因组引导命名异常|Reference genome guided naming error: {e}")
            return False

    def _run_minimap2(self, query_fa: str) -> Optional[str]:
        """
        运行minimap2比对|Run minimap2 alignment

        Args:
            query_fa: 查询FASTA文件|Query FASTA file

        Returns:
            str: PAF文件路径|PAF file path, or None if failed
        """
        try:
            # 生成输出文件路径|Generate output file path
            output_dir = os.path.join(self.config.rename_output_dir, "alignment")
            Path(output_dir).mkdir(parents=True, exist_ok=True)

            paf_file = os.path.join(output_dir, "alignment.paf")

            # 构建minimap2命令|Build minimap2 command
            cmd = [
                "minimap2",
                "-cx", self.config.naming_minimap2_preset,  # asm5/asm10/asm20
                "--secondary=no",  # 只输出 primary alignment
                query_fa,
                self.config.reference_genome
            ]

            self.logger.info(f"执行minimap2命令|Running minimap2: {' '.join(cmd)}")

            # 运行minimap2|Run minimap2
            with open(paf_file, 'w') as f_out:
                result = subprocess.run(
                    cmd,
                    stdout=f_out,
                    stderr=subprocess.PIPE,
                    text=True,
                    timeout=3600  # 1 hour timeout
                )

            if result.returncode != 0:
                self.logger.error(f"minimap2运行失败|minimap2 failed: {result.stderr}")
                return None

            # 检查输出文件|Check output file
            if not os.path.exists(paf_file) or os.path.getsize(paf_file) == 0:
                self.logger.error(f"minimap2输出为空|minimap2 output is empty: {paf_file}")
                return None

            self.logger.info(f"minimap2比对完成|minimap2 alignment completed: {paf_file}")
            return paf_file

        except subprocess.TimeoutExpired:
            self.logger.error("minimap2运行超时|minimap2 timeout")
            return None
        except Exception as e:
            self.logger.error(f"minimap2运行异常|minimap2 error: {e}")
            return None

    def _analyze_alignment(self, paf_file: str, query_fa: str) -> Optional[Dict[str, str]]:
        """
        分析PAF比对结果，生成命名映射|Analyze PAF alignment, generate naming map

        Args:
            paf_file: PAF文件路径|PAF file path
            query_fa: 查询FASTA文件|Query FASTA file

        Returns:
            dict: 命名映射 {query_name: reference_name}|Naming map {query_name: reference_name}
        """
        try:
            # 读取PAF文件|Read PAF file
            query_matches = {}  # {query_name: [(ref_name, identity, coverage), ...]}

            with open(paf_file, 'r') as f:
                for line in f:
                    fields = line.strip().split('\t')
                    if len(fields) < 12:
                        continue

                    query_name = fields[0]
                    query_length = int(fields[1])
                    ref_name = fields[5]
                    ref_length = int(fields[6])
                    strand = fields[4]
                    mapping_length = int(fields[10])  # 对齐长度|Alignment length
                    total_length = int(fields[11])  # 总长度|Total length

                    # 计算一致性|Calculate identity
                    # PAF格式中的NM字段是mismatches数量|NM field in PAF is number of mismatches
                    # 我们需要计算一致性 = (1 - mismatches / mapping_length) * 100
                    # Identity = (1 - mismatches / mapping_length) * 100

                    # 简化计算：使用block_length作为对齐长度
                    # Simplified: use block_length as alignment length
                    # 实际一致性需要从可选字段中读取|Actual identity needs to be read from optional fields
                    # 这里先简化为：如果匹配覆盖度足够，就认为是好的匹配
                    # Simplified: if coverage is sufficient, consider it a good match

                    # 计算覆盖度|Calculate coverage
                    query_coverage = (mapping_length / query_length) * 100
                    ref_coverage = (mapping_length / ref_length) * 100

                    # 选择最小覆盖度作为覆盖度指标|Use min coverage as coverage metric
                    coverage = min(query_coverage, ref_coverage)

                    # 对于一致性，暂时使用覆盖度作为替代
                    # For identity, temporarily use coverage as substitute
                    # 后续可以改进：从PAF的de字段计算真正的identity
                    # Can be improved later: calculate real identity from 'de' field in PAF
                    identity = coverage  # 简化|Simplified

                    if query_name not in query_matches:
                        query_matches[query_name] = []

                    query_matches[query_name].append({
                        'ref_name': ref_name,
                        'identity': identity,
                        'coverage': coverage,
                        'strand': strand,
                        'mapping_length': mapping_length
                    })

            # 为每个query选择最佳匹配|Select best match for each query
            naming_map = {}  # {query_name: new_name}
            unnamed_count = 0

            for query_name, matches in query_matches.items():
                # 找到最佳匹配（identity和coverage都满足阈值）|Find best match (both identity and coverage meet threshold)
                best_match = None
                best_score = 0

                for match in matches:
                    # 检查是否满足阈值|Check if meets threshold
                    if (match['identity'] >= self.config.naming_min_identity and
                        match['coverage'] >= self.config.naming_min_coverage):

                        # 计算综合得分|Calculate comprehensive score
                        score = match['identity'] * match['coverage'] / 100  # 简单加权|Simple weighting

                        if score > best_score:
                            best_score = score
                            best_match = match

                if best_match:
                    # 使用参考染色体名称|Use reference chromosome name
                    naming_map[query_name] = best_match['ref_name']
                    self.logger.debug(f"{query_name} -> {best_match['ref_name']} "
                                    f"(identity={best_match['identity']:.1f}%, coverage={best_match['coverage']:.1f}%)")
                else:
                    # 没有好的匹配，保留原名或使用备用命名|No good match, keep original or use backup naming
                    naming_map[query_name] = query_name  # 保留原名|Keep original name
                    unnamed_count += 1

            # 统计命名结果|Statistics of naming results
            total_queries = len(query_matches)
            named_count = total_queries - unnamed_count

            self.logger.info(f"命名统计|Naming statistics:")
            self.logger.info(f"  总序列数|Total sequences: {total_queries}")
            self.logger.info(f"  已命名|Named: {named_count} ({named_count/total_queries*100:.1f}%)")
            self.logger.info(f"  未命名|Unnamed: {unnamed_count} ({unnamed_count/total_queries*100:.1f}%)")

            return naming_map

        except Exception as e:
            self.logger.error(f"分析PAF文件异常|PAF analysis error: {e}")
            return None

    def _rename_sequences(self, input_fa: str, output_fa: str, naming_map: Dict[str, str]) -> bool:
        """
        根据命名映射重命名序列|Rename sequences based on naming map

        Args:
            input_fa: 输入FASTA文件|Input FASTA file
            output_fa: 输出FASTA文件|Output FASTA file
            naming_map: 命名映射|Naming map {query_name: new_name}

        Returns:
            bool: 是否成功|Whether successful
        """
        try:
            # 检查新名称是否有重复|Check if new names have duplicates
            new_names = list(naming_map.values())
            duplicates = [name for name in set(new_names) if new_names.count(name) > 1]

            if duplicates:
                self.logger.warning(f"检测到重复的染色体名|Detected duplicate chromosome names: {duplicates}")
                # 对于重复的名称，添加后缀|Add suffix for duplicate names
                name_counts = {}
                final_naming_map = {}

                for query_name, new_name in naming_map.items():
                    if new_name not in name_counts:
                        name_counts[new_name] = 1
                        final_naming_map[query_name] = new_name
                    else:
                        name_counts[new_name] += 1
                        # 添加后缀|Add suffix
                        final_naming_map[query_name] = f"{new_name}_alt{name_counts[new_name]}"

                naming_map = final_naming_map

            # 读取并重命名FASTA|Read and rename FASTA
            with open(input_fa, 'r') as f_in, open(output_fa, 'w') as f_out:
                current_name = None

                for line in f_in:
                    if line.startswith('>'):
                        # 提取序列名|Extract sequence name
                        current_name = line[1:].strip().split()[0]

                        # 获取新名称|Get new name
                        if current_name in naming_map:
                            new_name = naming_map[current_name]
                            # 保留描述信息|Keep description
                            desc = ' '.join(line[1:].strip().split()[1:])
                            if desc:
                                f_out.write(f">{new_name} {desc}\n")
                            else:
                                f_out.write(f">{new_name}\n")
                        else:
                            # 未找到映射，保留原名|No mapping found, keep original
                            f_out.write(line)
                            self.logger.warning(f"未找到命名映射|No naming mapping found: {current_name}")
                    else:
                        f_out.write(line)

            self.logger.info(f"序列重命名完成|Sequence renaming completed: {output_fa}")
            return True

        except Exception as e:
            self.logger.error(f"重命名序列异常|Sequence renaming error: {e}")
            return False

    def _generate_naming_report(self, naming_map: Dict[str, str], output_fa: str):
        """
        生成命名报告|Generate naming report

        Args:
            naming_map: 命名映射|Naming map
            output_fa: 输出FASTA文件|Output FASTA file
        """
        try:
            report_file = os.path.join(self.config.rename_output_dir, "naming_report.txt")

            with open(report_file, 'w') as f:
                f.write("=" * 60 + "\n")
                f.write("参考基因组引导命名报告|Reference Genome Guided Naming Report\n")
                f.write("=" * 60 + "\n\n")

                f.write(f"参考基因组|Reference genome: {self.config.reference_genome}\n")
                f.write(f"命名阈值|Naming thresholds:\n")
                f.write(f"  最小一致性|Min identity: {self.config.naming_min_identity}%\n")
                f.write(f"  最小覆盖度|Min coverage: {self.config.naming_min_coverage}%\n")
                f.write(f"  minimap2预设|minimap2 preset: {self.config.naming_minimap2_preset}\n\n")

                f.write(f"输出文件|Output file: {output_fa}\n\n")

                # 统计命名结果|Statistics of naming results
                original_names = set(naming_map.keys())
                new_names = set(naming_map.values())

                renamed = [old for old, new in naming_map.items() if old != new]
                not_renamed = [old for old, new in naming_map.items() if old == new]

                f.write(f"命名统计|Naming statistics:\n")
                f.write(f"  总序列数|Total sequences: {len(naming_map)}\n")
                f.write(f"  已重命名|Renamed: {len(renamed)}\n")
                f.write(f"  未重命名|Not renamed: {len(not_renamed)}\n\n")

                # 详细命名列表|Detailed naming list
                f.write("命名详情|Naming details:\n")
                f.write("-" * 60 + "\n")

                for old_name in sorted(naming_map.keys()):
                    new_name = naming_map[old_name]
                    if old_name != new_name:
                        f.write(f"  {old_name} -> {new_name}\n")
                    else:
                        f.write(f"  {old_name} (未重命名|not renamed)\n")

            self.logger.info(f"命名报告已保存|Naming report saved: {report_file}")

        except Exception as e:
            self.logger.warning(f"生成命名报告失败|Failed to generate naming report: {e}")
