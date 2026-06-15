"""
NCBI分类学注释核心处理模块|NCBI Taxonomy Annotation Core Processing Module
"""

import os
import sys
from collections import Counter
from typing import Dict, List, Set, Tuple
from .config import NCBITaxoConfig
from .utils import CommandRunner, detect_input_type, parse_lineage, format_number


class BlastTaxonomyProcessor:
    """BLAST分类学注释处理器|BLAST Taxonomy Annotation Processor"""

    def __init__(self, config: NCBITaxoConfig, logger, cmd_runner: CommandRunner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 输出文件路径|Output file paths
        self.accessions_file = f"{self.config.output_prefix}.accessions.txt"
        self.acc2taxid_file = f"{self.config.output_prefix}.acc2taxid.txt"
        self.taxonomy_file = f"{self.config.output_prefix}.taxonomy.txt"
        self.accession2title_file = f"{self.config.output_prefix}.accession2title.txt"

    def extract_accessions(self) -> List[str]:
        """步骤1: 提取accession列表|Step 1: Extract accession list"""
        self.logger.info("开始提取accession|Starting accession extraction")

        # 检测输入类型|Detect input type
        if self.config.input_type == 'auto':
            input_type = detect_input_type(self.config.input_file)
            self.logger.info(f"自动检测输入类型|Auto-detected input type: {input_type}")
        else:
            input_type = self.config.input_type
            self.logger.info(f"使用指定的输入类型|Using specified input type: {input_type}")

        accessions = set()
        total_lines = 0
        filtered_lines = 0

        if input_type == 'blast':
            # 从BLAST结果提取第2列|Extract column 2 from BLAST results
            col_idx = self.config.blast_column - 1  # 转换为0-based索引
            length_col_idx = 3  # 第4列是比对长度（0-based: 3）

            self.logger.info(f"过滤条件：比对长度 >= {self.config.min_alignment_length} bp|Filter: alignment length >= {self.config.min_alignment_length} bp")

            with open(self.config.input_file, 'r') as f:
                for line_num, line in enumerate(f, 1):
                    line = line.strip()
                    if not line:
                        continue

                    total_lines += 1
                    fields = line.split('\t')

                    # 检查列数是否足够|Check if enough columns
                    if len(fields) > max(col_idx, length_col_idx):
                        # 提取比对长度|Extract alignment length
                        try:
                            alignment_length = int(fields[length_col_idx].strip())
                        except (ValueError, IndexError):
                            # 如果无法解析长度，跳过该行|Skip if cannot parse length
                            continue

                        # 过滤：只保留比对长度 >= min_alignment_length的记录
                        if alignment_length >= self.config.min_alignment_length:
                            accession = fields[col_idx].strip()
                            if accession:
                                accessions.add(accession)
                                filtered_lines += 1

                    if line_num % 100000 == 0:
                        self.logger.debug(f"已处理|Processed: {line_num} 行|lines, {len(accessions)} 唯一accessions|unique accessions")

        elif input_type == 'accession':
            # 直接读取accession列表|Read accession list directly
            with open(self.config.input_file, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line:
                        accessions.add(line)

        # 转换为排序后的列表|Convert to sorted list
        accessions_list = sorted(list(accessions))

        if input_type == 'blast':
            self.logger.info(f"总行数|Total lines: {total_lines}")
            self.logger.info(f"满足长度条件的行|Lines passing length filter: {filtered_lines}")
            self.logger.info(f"提取到唯一accession数量|Extracted unique accessions: {len(accessions_list)}")
        else:
            self.logger.info(f"提取到唯一accession数量|Extracted unique accessions: {len(accessions_list)}")

        # 写入文件|Write to file
        with open(self.accessions_file, 'w') as f:
            for acc in accessions_list:
                f.write(f"{acc}\n")

        self.logger.info(f"accession列表已保存|Accession list saved: {self.accessions_file}")

        return accessions_list

    def lookup_taxids(self, accessions: List[str]) -> Dict[str, str]:
        """步骤2: 使用zgrep查找taxid|Step 2: Lookup taxids using zgrep"""
        self.logger.info("开始查找taxid|Starting taxid lookup")

        # 创建临时文件用于zgrep|Create temp file for zgrep
        temp_acc_file = self.accessions_file

        # 构建zgrep命令|Build zgrep command
        cmd = f"zgrep -F -f {temp_acc_file} {self.config.taxid_db} | cut -f2,3"

        success, output = self.cmd_runner.run(
            cmd,
            "查找taxid|Lookup taxids"
        )

        if not success:
            self.logger.error("taxid查找失败|Taxid lookup failed")
            return {}

        # 解析结果|Parse results
        acc2taxid = {}
        found_count = 0
        not_found = []

        for line in output.split('\n'):
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            if len(fields) >= 2:
                accession = fields[0].strip()
                taxid = fields[1].strip()
                acc2taxid[accession] = taxid
                found_count += 1

        # 检查哪些accession没有找到taxid|Check which accessions have no taxid
        for acc in accessions:
            if acc not in acc2taxid:
                not_found.append(acc)

        self.logger.info(f"找到taxid的accession|Accessions with taxid: {found_count}/{len(accessions)}")

        if not_found:
            self.logger.warning(f"未找到taxid的accession数量|Accessions without taxid: {len(not_found)}")
            if len(not_found) <= 10:
                self.logger.debug(f"未找到taxid的accession|Accessions without taxid: {', '.join(not_found)}")
            else:
                self.logger.debug(f"未找到taxid的accession（前10个）|Accessions without taxid (first 10): {', '.join(not_found[:10])}")

        # 写入acc2taxid文件|Write acc2taxid file
        with open(self.acc2taxid_file, 'w') as f:
            for acc in accessions:
                if acc in acc2taxid:
                    f.write(f"{acc}\t{acc2taxid[acc]}\n")
                else:
                    f.write(f"{acc}\t\n")  # 空taxid

        self.logger.info(f"acc2taxid映射已保存|acc2taxid mapping saved: {self.acc2taxid_file}")

        return acc2taxid

    def get_lineage(self, acc2taxid: Dict[str, str]) -> Dict[str, str]:
        """步骤3: 使用taxonkit获取lineage|Step 3: Get lineage using taxonkit"""
        self.logger.info("开始获取分类学谱系|Starting lineage retrieval")

        # 构建taxonkit命令|Build taxonkit command
        # 使用-R显示rank，-n显示name，-t显示taxid（tab分隔）
        # 输出格式：accession\ttaxid\ttaxid;taxid;...\tname;name;...\trank;rank;...
        # 我们使用taxonkit lineage -R来同时获取rank和name信息
        cmd = f"cat {self.acc2taxid_file} | taxonkit lineage -i 2 -t -R -n"

        success, output = self.cmd_runner.run(
            cmd,
            "获取分类学谱系|Get lineage"
        )

        if not success:
            self.logger.error("分类学谱系获取失败|Lineage retrieval failed")
            return {}

        # 解析结果|Parse results
        # 输出格式（使用-R和-n参数）：
        # accession\ttaxid\tname_lineage\ttaxid_lineage\tspecies_name\trank_lineage
        # 共6列
        accession_rank_name = {}  # 格式: "rank:name;rank:name;..."用于统计
        accession_full_lineage = {}  # 完整lineage（仅名称）用于保存
        success_count = 0

        for line in output.split('\n'):
            line = line.strip()
            if not line:
                continue

            fields = line.split('\t')
            # 使用-R -n参数时输出6列
            if len(fields) >= 6:
                accession = fields[0].strip()
                taxid = fields[1].strip()
                name_lineage = fields[2].strip()  # 第3列: name lineage
                rank_lineage = fields[5].strip()  # 第6列: rank lineage

                # 组合rank和name信息
                name_parts = name_lineage.split(';')
                rank_parts = rank_lineage.split(';')

                # 构建rank:name格式
                rank_name_pairs = []
                for i in range(min(len(name_parts), len(rank_parts))):
                    name = name_parts[i].strip()
                    rank = rank_parts[i].strip() if i < len(rank_parts) else ''
                    if name:
                        if rank and rank != 'clade':  # 跳过clade这种非正式rank
                            rank_name_pairs.append(f"{rank}:{name}")
                        else:
                            # 对于clade这种，也保留name
                            rank_name_pairs.append(name)

                rank_name_str = ";".join(rank_name_pairs)

                # 保存rank:name格式用于统计
                accession_rank_name[accession] = rank_name_str

                # 保存纯名称用于输出文件
                accession_full_lineage[accession] = name_lineage

                success_count += 1
            # 兼容：如果输出列数不够
            elif len(fields) >= 3:
                accession = fields[0].strip()
                full_lineage = fields[2].strip()
                accession_rank_name[accession] = full_lineage
                accession_full_lineage[accession] = full_lineage
                success_count += 1

        self.logger.info(f"成功获取lineage的accession|Accessions with lineage: {success_count}/{len(acc2taxid)}")

        # 写入最终文件|Write final taxonomy file
        with open(self.taxonomy_file, 'w') as f:
            for acc in sorted(acc2taxid.keys()):
                taxid = acc2taxid.get(acc, "")
                lineage = accession_full_lineage.get(acc, "")
                # 保存3列：accession, taxid, full_lineage（纯名称）
                f.write(f"{acc}\t{taxid}\t{lineage}\n")

        self.logger.info(f"分类学注释文件已保存|Taxonomy annotation saved: {self.taxonomy_file}")

        # 返回rank:name格式lineage用于统计
        return accession_rank_name

    def get_accession_titles(self, accessions: List[str]) -> Dict[str, str]:
        """步骤3.5: 使用Entrez获取accession的序列描述|Step 3.5: Get accession titles using Entrez

        Args:
            accessions: accession列表

        Returns:
            Dict[str, str]: {accession: title}
        """
        self.logger.info("开始获取accession序列描述|Starting accession title retrieval")

        # 检查是否需要获取描述
        if not self.config.fetch_titles:
            self.logger.info("未启用序列描述获取|Title fetching disabled")
            return {}

        try:
            from Bio import Entrez
        except ImportError:
            self.logger.warning("Biopython未安装，跳过序列描述获取|Biopython not installed, skipping title fetching")
            self.logger.warning("请安装: pip install biopython|Please install: pip install biopython")
            return {}

        # 设置邮箱（NCBI要求）
        Entrez.email = "anonymous@example.com"

        accession2title = {}
        batch_size = 500  # NCBI允许批量查询最多500个

        # 分批查询
        for i in range(0, len(accessions), batch_size):
            batch = accessions[i:i+batch_size]
            self.logger.info(f"正在查询|Querying batch {i//batch_size + 1}/{(len(accessions)-1)//batch_size + 1}: {len(batch)} accessions")

            try:
                # 使用epost + esummary批量获取
                search_handle = Entrez.epost(db="nuccore", id=",".join(batch))
                search_results = Entrez.read(search_handle)
                search_handle.close()

                webenv = search_results["WebEnv"]
                query_key = search_results["QueryKey"]

                # 获取摘要
                fetch_handle = Entrez.esummary(db="nuccore", webenv=webenv, query_key=query_key, retmax=batch_size)
                fetch_results = Entrez.read(fetch_handle)
                fetch_handle.close()

                # 解析结果
                for doc in fetch_results:
                    if "Caption" in doc and "Title" in doc:
                        acc = doc["Caption"]
                        title = doc["Title"]
                        accession2title[acc] = title

            except Exception as e:
                self.logger.warning(f"批次查询失败|Batch query failed: {str(e)}")
                continue

            # 添加延迟以遵守NCBI API规则（每秒最多3个请求）
            import time
            time.sleep(0.5)

        self.logger.info(f"获取到|Retrieved {len(accession2title)}/{len(accessions)} 个accession的序列描述|titles")

        # 写入文件
        with open(self.accession2title_file, 'w', encoding='utf-8') as f:
            for acc in sorted(accession2title.keys()):
                f.write(f"{acc}\t{accession2title[acc]}\n")

        self.logger.info(f"序列描述文件已保存|Accession titles saved: {self.accession2title_file}")

        # 生成分类统计报告
        self.summarize_sequence_types(accession2title)

        return accession2title

    def summarize_sequence_types(self, accession2title: Dict[str, str]):
        """根据序列描述统计序列类型|Summarize sequence types based on titles

        Args:
            accession2title: {accession: title}
        """
        if not accession2title:
            return

        self.logger.info("=" * 60)
        self.logger.info("序列类型统计分析|Sequence Type Analysis")
        self.logger.info("=" * 60)

        # 定义分类规则和关键词（按优先级从高到低）
        category_rules = {
            "Mitochondrial": {
                "keywords": ["mitochondrial", "mitochondrion", "mt genome", "mtDNA", "COX1", "CO1", "cytochrome oxidase", "cob", "cox", "nad", "atp"],
                "description": "线粒体基因组|Mitochondrial genome",
                "examples": []
            },
            "Chloroplast": {
                "keywords": ["chloroplast", "plastid", "cp genome", "cpDNA", "psb", "pet", "rbcl", "rbcL", "matK"],
                "description": "叶绿体基因组|Chloroplast genome",
                "examples": []
            },
            "Ribosomal RNA": {
                "keywords": ["18S", "28S", "16S", "5.8S", "5S", "ribosomal RNA", "rRNA gene", "small subunit", "large subunit", "SSU", "LSU"],
                "description": "核糖体RNA基因|Ribosomal RNA genes",
                "examples": []
            },
            "Internal Transcribed Spacer": {
                "keywords": ["ITS", "internal transcribed spacer", "ITS1", "ITS2"],
                "description": "内转录间隔区|Internal Transcribed Spacer",
                "examples": []
            },
            "Protein coding": {
                "keywords": ["mRNA", "protein coding", "CDS", "coding sequence", "gene for", "putative", "predicted"],
                "description": "蛋白编码基因|Protein coding genes",
                "examples": []
            },
            "Genomic": {
                "keywords": ["genomic", "genome", "chromosome", "chrom", "scaffold", "contig", "whole genome", "complete genome"],
                "description": "基因组序列|Genomic sequences",
                "examples": []
            }
        }

        # 未分类
        unclassified = []

        # 统计每个accession的类别
        acc2category = {}

        for acc, title in accession2title.items():
            title_lower = title.lower()
            categorized = False

            # 按优先级检查分类
            for category, rule in category_rules.items():
                for keyword in rule["keywords"]:
                    if keyword.lower() in title_lower:
                        acc2category[acc] = category
                        rule["examples"].append((acc, title))
                        categorized = True
                        break
                if categorized:
                    break

            if not categorized:
                acc2category[acc] = "Unclassified"
                unclassified.append((acc, title))

        # 统计各类别数量
        category_counts = {}
        for category in acc2category.values():
            category_counts[category] = category_counts.get(category, 0) + 1

        # 排序并输出统计结果
        sorted_categories = sorted(category_counts.items(), key=lambda x: x[1], reverse=True)

        self.logger.info("")
        self.logger.info("序列类型统计|Sequence Type Statistics:")
        self.logger.info("-" * 60)

        total = len(acc2category)
        for category, count in sorted_categories:
            percentage = (count / total) * 100 if total > 0 else 0
            if category in category_rules:
                desc = category_rules[category]["description"]
            else:
                desc = "未分类|Unclassified"
            self.logger.info(f"{desc:<40} {count:>6} ({percentage:>5.1f}%)")

        self.logger.info("-" * 60)
        self.logger.info(f"总计|Total: {total}")
        self.logger.info("")

        # 写入分类统计文件
        summary_file = f"{self.config.output_prefix}.sequence_types.txt"
        with open(summary_file, 'w', encoding='utf-8') as f:
            f.write("=" * 80 + "\n")
            f.write("序列类型统计分析|Sequence Type Analysis\n")
            f.write("=" * 80 + "\n")
            f.write(f"总Accession数量|Total accessions: {total}\n")
            f.write("\n")

            f.write("分类统计|Category Statistics:\n")
            f.write("-" * 80 + "\n")
            f.write(f"{'Category':<25}{'Description':<30}{'Count':>10}{'Percentage':>12}\n")
            f.write("-" * 80 + "\n")

            for category, count in sorted_categories:
                percentage = (count / total) * 100 if total > 0 else 0
                if category in category_rules:
                    desc = category_rules[category]["description"]
                else:
                    desc = "未分类|Unclassified"
                f.write(f"{category:<25}{desc:<30}{count:>10}{percentage:>11.1f}%\n")

            f.write("\n" + "=" * 80 + "\n")
            f.write("详细分类结果|Detailed Classification:\n")
            f.write("=" * 80 + "\n\n")

            # 按类别输出详细结果
            for category, rule in category_rules.items():
                if rule["examples"]:
                    f.write(f"\n{rule['description']} ({len(rule['examples'])} 条记录|records):\n")
                    f.write("-" * 80 + "\n")
                    for acc, title in rule["examples"][:10]:  # 只显示前10个
                        f.write(f"{acc:<20}\t{title}\n")
                    if len(rule["examples"]) > 10:
                        f.write(f"... 和其他|and {len(rule['examples']) - 10} 条|records\n")

            if unclassified:
                f.write(f"\n未分类|Unclassified ({len(unclassified)} 条记录|records):\n")
                f.write("-" * 80 + "\n")
                for acc, title in unclassified[:10]:
                    f.write(f"{acc:<20}\t{title}\n")
                if len(unclassified) > 10:
                    f.write(f"... 和其他|and {len(unclassified) - 10} 条|records\n")

        self.logger.info(f"序列类型统计文件已保存|Sequence type analysis saved: {summary_file}")
        self.logger.info("=" * 60)

    def read_blast_hits(self) -> Dict[str, int]:
        """读取BLAST结果，统计每个accession的hit次数|Read BLAST results and count hits per accession"""
        self.logger.info("读取BLAST结果统计hit次数|Reading BLAST results to count hits")

        hit_counts = Counter()
        col_idx = self.config.blast_column - 1

        with open(self.config.input_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                if not line:
                    continue

                fields = line.split('\t')
                if len(fields) > col_idx:
                    accession = fields[col_idx].strip()
                    if accession:
                        hit_counts[accession] += 1

                if line_num % 100000 == 0:
                    self.logger.debug(f"已统计|Counted: {line_num} 行|lines, {sum(hit_counts.values())} hits")

        self.logger.info(f"总BLAST hit次数|Total BLAST hits: {sum(hit_counts.values())}")
        self.logger.info(f"唯一accession数量|Unique accessions: {len(hit_counts)}")

        return dict(hit_counts)

    def run_pipeline(self) -> Tuple[Dict[str, str], Dict[str, int], List[str]]:
        """运行完整的处理流程|Run complete processing pipeline

        Returns:
            Tuple[Dict[str, str], Dict[str, int], List[str]]:
                (accession_lineage, blast_hits, accessions_list)
        """
        self.logger.info("=" * 60)
        self.logger.info("开始NCBI分类学注释流程|Starting NCBI Taxonomy Annotation Pipeline")
        self.logger.info("=" * 60)

        # 步骤1: 提取accession
        accessions_list = self.extract_accessions()

        # 步骤2: 查找taxid
        acc2taxid = self.lookup_taxids(accessions_list)

        # 步骤3: 获取lineage
        accession_lineage = self.get_lineage(acc2taxid)

        # 步骤3.5: 获取accession序列描述（可选）
        accession2title = self.get_accession_titles(accessions_list)

        # 步骤4: 如果是blast输入，读取hit次数
        blast_hits = {}
        input_type = self.config.input_type if self.config.input_type != 'auto' else detect_input_type(self.config.input_file)

        if input_type == 'blast':
            if self.config.stats_target in ['blast_hits', 'both']:
                blast_hits = self.read_blast_hits()

        self.logger.info("=" * 60)
        self.logger.info("NCBI分类学注释流程完成|NCBI Taxonomy Annotation Pipeline Completed")
        self.logger.info("=" * 60)

        return accession_lineage, blast_hits, accessions_list
