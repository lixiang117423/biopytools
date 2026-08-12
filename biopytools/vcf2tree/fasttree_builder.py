"""
FastTree建树模块|FastTree Tree Building Module
"""

import subprocess

from .utils import build_conda_command


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

        # FastTree 不支持 ASC 校正: SNP 数据仅含可变位点(无恒定位点), 分支长度会
        # 系统性低估。若需 ASC 校正, 请改用 IQ-TREE 后端(默认即开启 +ASC)。
        # |FastTree has no ASC correction: SNP data holds only variable sites (no invariant
        # sites), so branch lengths are systematically underestimated. Use the IQ-TREE
        # backend (which enables +ASC by default) if ASC correction is required.
        self.logger.warning(
            "FastTree不支持ASC校正, SNP分支长度可能低估; 需ASC请用IQ-TREE后端|"
            "FastTree lacks ASC correction; SNP branch lengths may be underestimated. "
            "Use IQ-TREE backend if ASC is required"
        )

        # 构建参数: fasttree -nt alignment.fa > tree.nwk
        # FastTree 始终把树写到stdout(从不写stderr)
        # |FastTree always writes the tree to stdout (never to stderr)
        args = ['-nt', str(self.config.snps_fa)]

        # 添加额外参数|Add extra params
        if self.config.fasttree_params:
            args.extend(self.config.fasttree_params.split())

        # 经build_conda_command包装(与IQ-TREE一致): 在conda env则conda run, 否则直调
        # |Wrap via build_conda_command (same as IQ-TREE): conda run if in an env, else direct
        cmd = build_conda_command(self.config.fasttree_path, args)

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

            # FastTree 始终把Newick树写到stdout; stderr只有进度文本, 绝不是树。
            # 旧版在stdout为空时误把stderr当树, 会把进度文本误判为有效Newick, 此处移除该fallback。
            # |FastTree always writes the Newick tree to stdout; stderr holds only progress
            # text, never a tree. The old code mistook stderr for a tree when stdout was empty
            # (misreading progress text as valid Newick); that fallback is removed here.
            tree_output = result.stdout.strip()

            if not tree_output or not tree_output.endswith(';'):
                self.logger.error(
                    "FastTree未输出有效的Newick树(stdout为空)|"
                    "FastTree produced no valid Newick tree (empty stdout)"
                )
                self.logger.debug(f"stdout: {repr(result.stdout)}")
                return False

            with open(self.config.tree_nwk, 'w') as f:
                f.write(tree_output + '\n')

            self.logger.info(
                f"系统发育树已保存|Phylogenetic tree saved: {self.config.tree_nwk}"
            )

            # 报告关键统计(从stderr读取进度信息)|Report key stats from stderr
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
            self.logger.error(
                f"FastTree未找到|FastTree not found: {self.config.fasttree_path}"
            )
            return False
