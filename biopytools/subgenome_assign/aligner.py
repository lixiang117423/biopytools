"""
比对模块|Alignment Module

合并每个亲本的 hap FASTA，将目标基因组 minimap2 比对到每个合并后的亲本参考
|Concatenate each parent's hap FASTAs and align target genome vs each combined
parental reference using minimap2
"""

import os
import subprocess
from pathlib import Path

from ..common.conda_runner import build_conda_command


class ParentAligner:
    """亲本比对器|Parent Aligner"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

        # 亲本合并 FASTA 路径（每个亲本一份）|Combined FASTA per parent
        # {parent_name: combined.fa}
        self.combined_paths = {}

    def combine_parent_haps(self) -> bool:
        """
        合并每个亲本的 hap FASTA 到 01_alignment/<parent>.combined.fa
        |Concatenate each parent's hap FASTAs to 01_alignment/<parent>.combined.fa
        """
        self.logger.info("合并亲本 hap FASTA|Combining parental hap FASTAs")
        for name, haps in self.config.parents.items():
            out = self.config.alignment_dir / f"{name}.combined.fa"
            self.combined_paths[name] = str(out)

            # 简单 cat 即可，FASTA 格式宽容|Simple cat works for FASTA
            self.logger.info(f"合并|Combine {name}: {len(haps)} 个 hap 文件|hap files -> {out.name}")
            with open(out, 'w') as fh:
                for h in haps:
                    with open(h) as src:
                        for line in src:
                            fh.write(line)

            size = os.path.getsize(out)
            self.logger.info(f"  {name}: {size / 1e6:.1f} MB")

        return True

    def align_target_to_all_parents(self) -> dict:
        """
        将目标基因组比对到每个亲本合并 FASTA
        |Align target genome to each combined parental FASTA

        Returns:
            {parent_name: paf_path} 字典|dict of parent_name to PAF path
        """
        results = {}
        opts = self.config.get_minimap2_options()

        for name, combined in self.combined_paths.items():
            paf = self.config.alignment_dir / f"target_vs_{name}.paf"
            results[name] = str(paf)

            # 断点续传：PAF 已存在且非空则跳过
            # |Resume: skip if PAF exists and non-empty
            if paf.exists() and paf.stat().st_size > 0:
                self.logger.info(f"跳过已完成|Skipping existing: {paf.name}")
                continue

            # 单工具命令 + stdout 重定向(build_conda_command 自动conda包装)
            # |Single tool command with stdout redirect (auto conda wrap)
            cmd = build_conda_command(
                self.config.minimap2_path,
                opts + [combined, self.config.target],
            )
            self.logger.info(f"命令|Command: {' '.join(cmd)} > {paf}")

            try:
                with open(paf, 'w') as f:
                    result = subprocess.run(
                        cmd, stdout=f, stderr=subprocess.PIPE, text=True)
                if result.returncode != 0:
                    self.logger.error(f"比对失败|Alignment failed: target vs {name}")
                    self.logger.error(f"错误信息|Stderr: {result.stderr.strip()[:1000]}")
                    return {}
            except Exception as e:
                self.logger.error(f"比对异常|Alignment error: target vs {name}: {e}")
                return {}

        return results
