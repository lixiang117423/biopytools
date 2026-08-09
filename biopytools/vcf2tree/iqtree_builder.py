"""
IQ-TREE建树模块|IQ-TREE Tree Building Module
"""

import shutil
from pathlib import Path
from .utils import build_conda_command


class IqtreeBuilder:
    """IQ-TREE建树器|IQ-TREE Tree Builder"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def build(self) -> bool:
        """使用IQ-TREE构建系统发育树|Build phylogenetic tree with IQ-TREE

        默认：ModelFinder自动模型选择 + UFBoot 1000
        |Default: ModelFinder auto-selection + UFBoot 1000

        Returns:
            bool: 是否成功|Whether successful
        """
        self.logger.info("=" * 60)
        self.logger.info(
            "步骤2: IQ-TREE构建系统发育树|"
            "Step 2: IQ-TREE phylogenetic tree construction"
        )
        self.logger.info("=" * 60)

        # IQ-TREE输出前缀|IQ-TREE output prefix
        prefix = str(self.config.step2_dir / self.config.base_name)

        # 构建参数|Build arguments
        args = [
            '-s', str(self.config.snps_fa),
            '--prefix', prefix,
            '-T', str(self.config.threads),
        ]

        # 保持静默|Keep quiet
        args.append('--quiet')

        # 模型选择|Model selection
        # SNP对齐仅含可变位点(无恒定位点), 需ASC校正(Lewis 2001)避免分支长度低估。
        # 默认iqtree_asc=True时, 若模型名未含ASC则追加+ASC; 用户显式指定模型时同样适用。
        # |SNP alignments have only variable sites (no invariant sites); ASC correction
        # (Lewis 2001) is needed to avoid underestimated branch lengths. When
        # iqtree_asc=True (default), append +ASC unless the model already specifies ASC.
        if self.config.iqtree_model:
            model = self.config.iqtree_model
        else:
            model = 'MFP'

        if self.config.iqtree_asc and 'ASC' not in model.upper():
            model = model + '+ASC'

        args.extend(['-m', model])
        if self.config.iqtree_model:
            self.logger.info(
                f"使用指定模型|Using specified model: {model}"
            )
        else:
            self.logger.info(
                f"使用自动模型选择(+ASC校正)|Automatic model selection with ASC: {model}"
            )

        # UFBoot|UFBoot
        if self.config.iqtree_bootstrap > 0:
            args.extend(['-B', str(self.config.iqtree_bootstrap)])
            self.logger.info(
                f"使用UFBoot|Using UFBoot: {self.config.iqtree_bootstrap} replicates"
            )

        # 外群|Outgroup
        if self.config.outgroup:
            args.extend(['-o', self.config.outgroup])
            self.logger.info(
                f"设置外群|Setting outgroup: {self.config.outgroup}"
            )

        # 构建conda包装命令|Build conda-wrapped command
        cmd = build_conda_command(self.config.iqtree_path, args)

        # 执行|Execute
        result = self.cmd_runner.run(cmd, "IQ-TREE建树|IQ-TREE tree building")

        if result:
            # IQ-TREE输出.treefile → 复制为统一命名
            # |IQ-TREE outputs .treefile → copy to unified name
            iqtree_tree = f"{prefix}.treefile"
            if Path(iqtree_tree).exists():
                shutil.copy(iqtree_tree, self.config.tree_nwk)
                self.logger.info(
                    f"系统发育树已保存|Phylogenetic tree saved: {self.config.tree_nwk}"
                )
            else:
                self.logger.warning(
                    f"IQ-TREE树文件未找到|IQ-TREE tree file not found: {iqtree_tree}"
                )
                # 尝试.contree|Try .contree (consensus tree)
                contree = f"{prefix}.contree"
                if Path(contree).exists():
                    shutil.copy(contree, self.config.tree_nwk)
                    self.logger.info(
                        f"使用一致树|Using consensus tree: {self.config.tree_nwk}"
                    )
                else:
                    return False

        return result
