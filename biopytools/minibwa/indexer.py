"""
基因组索引管理模块|Genome Index Management Module

封装 `minibwa index` 命令，支持标准索引和BS-seq索引（--meth）。
|Wraps `minibwa index`, supporting standard and BS-seq (--meth) indices.
"""

import os
from pathlib import Path

from .utils import build_conda_command


class GenomeIndexer:
    """基因组索引器|Genome Indexer"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner

    def _index_files_exist(self) -> bool:
        """检查标准索引文件是否齐全|Check standard index files"""
        l2b = Path(f"{self.config.index_prefix}.l2b")
        mbw = Path(f"{self.config.index_prefix}.mbw")
        return l2b.exists() and mbw.exists()

    def _meth_index_exists(self) -> bool:
        """检查BS-seq索引是否存在|Check BS-seq index"""
        meth_mbw = Path(f"{self.config.index_prefix}.meth.mbw")
        return meth_mbw.exists()

    def check_and_build_index(self) -> bool:
        """
        检查并按需构建索引|Check and build index if missing

        BS-seq模式下，minibwa需要标准索引+额外的.meth.mbw
        |In BS-seq mode, minibwa requires both standard and .meth.mbw indices
        """
        self.logger.info("检查基因组索引|Checking genome index")

        need_standard = not self._index_files_exist()
        need_meth = self.config.is_meth_mode() and not self._meth_index_exists()

        if not need_standard and not need_meth:
            self.logger.info("基因组索引已存在|Genome index already exists")
            return True

        if need_standard:
            self.logger.info("缺少标准索引|Missing standard index")
            if not self.build_standard_index():
                return False

        if need_meth:
            self.logger.info("缺少BS-seq索引|Missing BS-seq index")
            if not self.build_meth_index():
                return False

        return True

    def build_standard_index(self) -> bool:
        """构建标准索引|Build standard index"""
        args = [
            '-t', str(self.config.threads),
            self.config.genome,
            str(self.config.index_prefix),
        ]
        cmd = build_conda_command(self.config.minibwa_path, ['index'] + args)

        success = self.cmd_runner.run(
            cmd,
            f"构建标准索引|Building standard index: {self.config.genome_name}"
        )
        if success:
            self.logger.info(
                f"标准索引构建成功|Standard index built: "
                f"{self.config.index_prefix}.l2b + .mbw"
            )
        else:
            self.logger.error("标准索引构建失败|Standard index build failed")
        return success

    def build_meth_index(self) -> bool:
        """构建BS-seq索引（--meth）|Build BS-seq index"""
        args = [
            '--meth',
            '-t', str(self.config.threads),
            self.config.genome,
            str(self.config.index_prefix),
        ]
        cmd = build_conda_command(self.config.minibwa_path, ['index'] + args)

        success = self.cmd_runner.run(
            cmd,
            f"构建BS-seq索引|Building BS-seq index: {self.config.genome_name}"
        )
        if success:
            self.logger.info(
                f"BS-seq索引构建成功|BS-seq index built: "
                f"{self.config.index_prefix}.meth.mbw"
            )
        else:
            self.logger.error("BS-seq索引构建失败|BS-seq index build failed")
        return success

    def get_index_prefix(self) -> str:
        """获取索引前缀（minibwa map用）|Get index prefix for minibwa map"""
        return str(self.config.index_prefix)
