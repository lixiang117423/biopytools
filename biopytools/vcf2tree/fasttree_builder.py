"""
FastTree建树模块|FastTree Tree Building Module
"""

import shutil
import subprocess


class FastTreeBuilder:
    """FastTree建树器|FastTree Tree Builder"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build(self) -> bool:
        """使用FastTree构建系统发育树|Build phylogenetic tree with FastTree

        FastTree默认氨基酸模式，对SNP数据必须使用-nt（核酸/GTR模式）
        |FastTree defaults to amino acid mode; -nt (nucleotide/GTR) is required for SNP data

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info(
            "步骤2: FastTree构建系统发育树|"
            "Step 2: FastTree phylogenetic tree construction"
        )
        self.logger.info("=" * 60)

        fasttree_bin = shutil.which(self.config.fasttree_path)
        if not fasttree_bin:
            fasttree_bin = self.config.fasttree_path

        # 构建命令：fasttree -nt alignment.fa > tree.nwk
        # FastTree 输出到stdout|FastTree outputs to stdout
        cmd = [fasttree_bin, '-nt', str(self.config.snps_fa)]

        # 添加额外参数|Add extra params
        if self.config.fasttree_params:
            cmd.extend(self.config.fasttree_params.split())

        # 记录命令|Log command
        self.logger.info(f"执行|Executing: FastTree建树|FastTree tree building")
        self.logger.info(
            f"命令|Command: {' '.join(cmd)} > {self.config.tree_nwk}"
        )
        self.logger.info(f"工作目录|Working directory: {self.config.step2_dir}")

        try:
            result = subprocess.run(
                cmd,
                shell=False,
                capture_output=True,
                text=True,
                check=True,
                cwd=self.config.step2_dir
            )

            # FastTree 将树输出到stdout|FastTree outputs tree to stdout
            tree_output = result.stdout.strip()
            if not tree_output:
                # 某些版本输出到stderr|Some versions output to stderr
                tree_output = result.stderr.strip()

            if not tree_output or not tree_output.endswith(';'):
                self.logger.error(
                    "FastTree输出不是有效的Newick格式|"
                    "FastTree output is not valid Newick format"
                )
                self.logger.debug(f"输出内容|Output content: {repr(tree_output)}")
                return False

            with open(self.config.tree_nwk, 'w') as f:
                f.write(tree_output + '\n')

            self.logger.info(
                f"系统发育树已保存|Phylogenetic tree saved: {self.config.tree_nwk}"
            )

            # 报告关键统计|Report key stats from stderr
            if result.stderr:
                for line in result.stderr.strip().split('\n'):
                    if any(kw in line for kw in ['Score', 'Rates', 'ML', 'Length']):
                        self.logger.info(f"FastTree: {line.strip()}")

            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"FastTree执行失败|FastTree execution failed")
            if e.stderr:
                self.logger.error(f"错误信息|Error message: {e.stderr}")
            return False
        except FileNotFoundError:
            self.logger.error(f"FastTree未找到|FastTree not found: {fasttree_bin}")
            return False
