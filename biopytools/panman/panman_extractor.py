"""
PanMAN数据提取器模块|PanMAN Data Extractor Module
"""

import os
from .utils import PanMANParser


class PanMANExtractor:
    """PanMAN数据提取器|PanMAN Data Extractor"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.parser = PanMANParser(logger)

    def extract(self):
        """
        根据配置选项提取数据|Extract data based on configuration options

        Returns:
            dict: 提取结果字典，包含各类输出文件路径|Extraction results dict containing output file paths
        """
        self.logger.info(f"开始从PanMAN提取数据|Starting data extraction from PanMAN")
        self.logger.info(f"PanMAN文件|PanMAN file: {self.config.panman_file}")

        results = {}

        # 按优先级提取各类数据|Extract data types by priority
        if self.config.extract_summary:
            results['summary'] = self._extract_summary()

        if self.config.extract_newick:
            results['newick'] = self._extract_newick()

        if self.config.extract_extended_newick:
            results['extended_newick'] = self._extract_extended_newick()

        if self.config.extract_fasta:
            results['fasta'] = self._extract_fasta()

        if self.config.extract_msa:
            results['msa'] = self._extract_msa()

        if self.config.extract_vcf:
            results['vcf'] = self._extract_vcf()

        if self.config.extract_gfa:
            results['gfa'] = self._extract_gfa()

        if self.config.extract_maf:
            results['maf'] = self._extract_maf()

        if self.config.extract_aa:
            results['aa'] = self._extract_aa()

        # 新增功能|New features
        if self.config.extract_subnet:
            results['subnet'] = self._extract_subnet()

        if self.config.extract_annotate:
            results['annotate'] = self._extract_annotate()

        if self.config.extract_reroot:
            results['reroot'] = self._extract_reroot()

        if self.config.extract_create_network:
            results['create_network'] = self._extract_create_network()

        if self.config.extract_print_mutations:
            results['print_mutations'] = self._extract_print_mutations()

        # 范围查询功能|Range query feature
        if self.config.range_query_index:
            results['range_query'] = self._extract_range_query()

        # 对于conda后端，从输入文件目录移动输出文件|For conda backend, move output files from input dir
        if self.config.backend == "conda":
            results = self._move_extracted_files(results)

        # 记录提取结果|Log extraction results
        self.logger.info(f"数据提取完成|Data extraction completed")
        for data_type, file_path in results.items():
            if file_path:
                size = self.parser.get_file_size(file_path) if isinstance(file_path, str) and os.path.exists(file_path) else "N/A"
                self.logger.info(f"  {data_type.upper()}: {file_path} ({size})")

        return results

    def _move_extracted_files(self, results: dict) -> dict:
        """移动从输入文件目录生成的输出文件到目标目录|Move output files from input dir to target dir"""
        import shutil
        import os
        import glob

        input_dir = os.path.dirname(self.config.panman_file)
        moved_results = {}

        for data_type, expected_path in results.items():
            if not expected_path or not isinstance(expected_path, str):
                moved_results[data_type] = expected_path
                continue

            # 检查文件是否在预期位置|Check if file is at expected location
            if os.path.exists(expected_path):
                moved_results[data_type] = expected_path
                continue

            # 检查info/子目录（panmanUtils某些输出在info/目录）|Check info/ subdir
            expected_dir = os.path.dirname(expected_path)
            info_dir = os.path.join(expected_dir, "info")
            if os.path.exists(info_dir):
                # 查找匹配的文件|Find matching files
                basename = os.path.basename(expected_path)
                # 匹配模式：basename 或 basename_*|Match pattern: basename or basename_*
                pattern = os.path.join(info_dir, f"{basename}_*")
                matches = glob.glob(pattern)
                if matches:
                    # 使用第一个匹配文件|Use first matched file
                    source = matches[0]
                    os.makedirs(expected_dir, exist_ok=True)
                    shutil.move(source, expected_path)
                    self.logger.info(f"从info目录移动文件|Moved file from info dir: {expected_path}")
                    moved_results[data_type] = expected_path
                    continue

            # 检查输入文件目录|Check input file directory
            basename = os.path.basename(expected_path)
            input_path = os.path.join(input_dir, basename)

            if os.path.exists(input_path):
                # 移动到预期位置|Move to expected location
                os.makedirs(os.path.dirname(expected_path), exist_ok=True)
                shutil.move(input_path, expected_path)
                self.logger.info(f"从输入目录移动文件|Moved file from input dir: {expected_path}")
                moved_results[data_type] = expected_path
            else:
                # 都找不到，返回原路径|Not found anywhere, return original path
                moved_results[data_type] = expected_path
                self.logger.warning(f"未找到输出文件|Output file not found: {expected_path} or {input_path}")

        # 清理空的info目录|Cleanup empty info directory
        info_dir = os.path.join(self.config.output_dir, "info")
        if os.path.exists(info_dir) and not os.listdir(info_dir):
            os.rmdir(info_dir)

        return moved_results

    def _extract_summary(self) -> str:
        """提取摘要统计|Extract summary statistics"""
        self.logger.info("提取摘要统计|Extracting summary statistics")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_summary.txt"
        )

        args = [
            "-I", self.config.panman_file,
            "--summary"
        ]

        if output_file:
            args.extend(["--output-file", output_file])

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取摘要|Extract summary"
        )

        if success:
            self.logger.info(f"摘要提取成功|Summary extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("摘要提取失败|Summary extraction failed")
            return None

    def _extract_fasta(self) -> str:
        """提取FASTA序列|Extract FASTA sequences"""
        self.logger.info("提取FASTA序列|Extracting FASTA sequences")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}.fa"
        )

        args = [
            "-I", self.config.panman_file,
            "--fasta",
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取FASTA|Extract FASTA"
        )

        if success:
            self.logger.info(f"FASTA提取成功|FASTA extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("FASTA提取失败|FASTA extraction failed")
            return None

    def _extract_msa(self) -> str:
        """提取MSA比对|Extract MSA alignment"""
        self.logger.info("提取MSA比对|Extracting MSA alignment")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_msa.fa"
        )

        args = [
            "-I", self.config.panman_file,
            "--fasta-aligned",
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取MSA|Extract MSA"
        )

        if success:
            self.logger.info(f"MSA提取成功|MSA extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("MSA提取失败|MSA extraction failed")
            return None

    def _extract_vcf(self) -> str:
        """提取VCF变异|Extract VCF variants"""
        self.logger.info("提取VCF变异|Extracting VCF variants")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}.vcf"
        )

        args = [
            "-I", self.config.panman_file,
            "--vcf",
            "--reference", self.config.reference,
            "--output-file", output_file
        ]

        # 添加 treeID 参数（如果指定）
        if self.config.tree_id:
            args.extend(["--treeID", self.config.tree_id])

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取VCF|Extract VCF"
        )

        if success:
            self.logger.info(f"VCF提取成功|VCF extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("VCF提取失败|VCF extraction failed")
            return None

    def _extract_gfa(self) -> str:
        """提取GFA格式|Extract GFA format"""
        self.logger.info("提取GFA格式|Extracting GFA format")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}.gfa"
        )

        args = [
            "-I", self.config.panman_file,
            "--gfa",
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取GFA|Extract GFA"
        )

        if success:
            self.logger.info(f"GFA提取成功|GFA extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("GFA提取失败|GFA extraction failed")
            return None

    def _extract_newick(self) -> str:
        """提取Newick树|Extract Newick tree"""
        self.logger.info("提取Newick树|Extracting Newick tree")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}.nwk"
        )

        args = [
            "-I", self.config.panman_file,
            "--newick",
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取Newick|Extract Newick"
        )

        if success:
            self.logger.info(f"Newick提取成功|Newick extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("Newick提取失败|Newick extraction failed")
            return None

    def _extract_extended_newick(self) -> str:
        """提取扩展Newick格式|Extract extended Newick format"""
        self.logger.info("提取扩展Newick格式|Extracting extended Newick format")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_extended.nwk"
        )

        args = [
            "-I", self.config.panman_file,
            "--extended-newick",
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取扩展Newick|Extract extended Newick"
        )

        if success:
            self.logger.info(f"扩展Newick提取成功|Extended Newick extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("扩展Newick提取失败|Extended Newick extraction failed")
            return None

    def _extract_maf(self) -> str:
        """提取MAF格式|Extract MAF format"""
        self.logger.info("提取MAF格式|Extracting MAF format")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}.maf"
        )

        args = [
            "-I", self.config.panman_file,
            "--maf",
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取MAF|Extract MAF"
        )

        if success:
            self.logger.info(f"MAF提取成功|MAF extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("MAF提取失败|MAF extraction failed")
            return None

    def _extract_aa(self) -> str:
        """提取氨基酸翻译|Extract amino acid translations"""
        self.logger.info("提取氨基酸翻译|Extracting amino acid translations")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_aa.tsv"
        )

        args = [
            "-I", self.config.panman_file,
            "--aa-translations",
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取氨基酸翻译|Extract amino acid translations"
        )

        if success:
            self.logger.info(f"氨基酸翻译提取成功|Amino acid translations extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("氨基酸翻译提取失败|Amino acid translations extraction failed")
            return None

    def _extract_subnet(self) -> str:
        """提取子网络|Extract subnet"""
        self.logger.info("提取子网络|Extracting subnet")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_subnet.panman"
        )

        args = [
            "-I", self.config.panman_file,
            "--subnet",
            "--input-file", self.config.input_file,
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="提取子网络|Extract subnet"
        )

        if success:
            self.logger.info(f"子网络提取成功|Subnet extraction successful: {output_file}")
            return output_file
        else:
            self.logger.error("子网络提取失败|Subnet extraction failed")
            return None

    def _extract_annotate(self) -> str:
        """注释节点|Annotate nodes"""
        self.logger.info("注释节点|Annotating nodes")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_annotated.panman"
        )

        args = [
            "-I", self.config.panman_file,
            "--annotate",
            "--input-file", self.config.input_file,
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="注释节点|Annotate nodes"
        )

        if success:
            self.logger.info(f"节点注释成功|Node annotation successful: {output_file}")
            return output_file
        else:
            self.logger.error("节点注释失败|Node annotation failed")
            return None

    def _extract_reroot(self) -> str:
        """重新扎根|Reroot tree"""
        self.logger.info("重新扎根树|Rerooting tree")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_rerooted.panman"
        )

        args = [
            "-I", self.config.panman_file,
            "--reroot",
            "--reference", self.config.reference,
            "--output-file", output_file
        ]

        # 添加 ACR 方法
        if self.config.acr_method != "fitch":
            args.extend(["--acr", self.config.acr_method])

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="重新扎根|Reroot tree"
        )

        if success:
            self.logger.info(f"重新扎根成功|Reroot successful: {output_file}")
            return output_file
        else:
            self.logger.error("重新扎根失败|Reroot failed")
            return None

    def _extract_create_network(self) -> str:
        """创建网络|Create network"""
        self.logger.info("创建网络|Creating network")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_network.panman"
        )

        args = [
            "--create-network",
            "--input-file", self.config.input_file,
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="创建网络|Create network"
        )

        if success:
            self.logger.info(f"网络创建成功|Network creation successful: {output_file}")
            return output_file
        else:
            self.logger.error("网络创建失败|Network creation failed")
            return None

    def _extract_print_mutations(self) -> str:
        """打印突变|Print mutations"""
        self.logger.info("打印突变信息|Printing mutations")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_mutations.txt"
        )

        args = [
            "-I", self.config.panman_file,
            "--printMutations"
        ]

        # 输出到文件
        if output_file:
            args.extend(["--output-file", output_file])

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="打印突变|Print mutations"
        )

        if success:
            self.logger.info(f"突变信息打印成功|Mutations print successful: {output_file}")
            return output_file
        else:
            self.logger.error("突变信息打印失败|Mutations print failed")
            return None

    def _extract_range_query(self) -> str:
        """范围查询|Range query"""
        self.logger.info(f"执行范围查询|Executing range query: {self.config.range_start}-{self.config.range_end}")

        output_file = os.path.join(
            self.config.output_dir,
            f"{self.config.output_prefix}_range_{self.config.range_start}_{self.config.range_end}.fa"
        )

        args = [
            "-I", self.config.panman_file,
            "--index", self.config.range_query_index,
            "-x", str(self.config.range_start),
            "-y", str(self.config.range_end),
            "--reference", self.config.reference,
            "--output-file", output_file
        ]

        success, output = self.cmd_runner.run_panman_command(
            args,
            description="范围查询|Range query"
        )

        if success:
            self.logger.info(f"范围查询成功|Range query successful: {output_file}")
            return output_file
        else:
            self.logger.error("范围查询失败|Range query failed")
            return None
