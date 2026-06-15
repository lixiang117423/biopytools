"""
PanMAN构建器模块|PanMAN Builder Module
"""

import os
from .utils import PanMANValidator


class PanMANBuilder:
    """PanMAN构建器|PanMAN Builder"""

    def __init__(self, config, logger, cmd_runner):
        self.config = config
        self.logger = logger
        self.cmd_runner = cmd_runner
        self.validator = PanMANValidator(logger)

    def build(self):
        """
        根据输入类型构建PanMAN|Build PanMAN based on input type

        Returns:
            str: 输出的PanMAN文件路径|Output PanMAN file path
        """
        input_type = self.config.get_input_type()
        input_file = self.config.get_input_file()

        self.logger.info(f"开始构建PanMAN|Starting PanMAN build")
        self.logger.info(f"输入类型|Input type: {input_type}")
        self.logger.info(f"输入文件|Input file: {input_file}")
        self.logger.info(f"系统发育树文件|Phylogenetic tree file: {self.config.newick_file}")

        # 验证输入文件|Validate input files
        if not self._validate_inputs(input_type, input_file):
            return None

        # 构建PanMAN|Build PanMAN
        success = False
        if input_type == "pangraph":
            success = self._build_from_pangraph()
        elif input_type == "gfa":
            success = self._build_from_gfa()
        elif input_type == "msa":
            success = self._build_from_msa()
        else:
            self.logger.error(f"不支持的输入类型|Unsupported input type: {input_type}")
            return None

        if success:
            output_file = os.path.join(
                self.config.output_dir,
                f"{self.config.output_prefix}.panman"
            )
            self.logger.info(f"PanMAN构建成功|PanMAN build successful: {output_file}")
            return output_file
        else:
            self.logger.error("PanMAN构建失败|PanMAN build failed")
            return None

    def _validate_inputs(self, input_type: str, input_file: str) -> bool:
        """验证输入文件|Validate input files"""
        self.logger.info("验证输入文件|Validating input files")

        # 验证输入文件|Validate input file
        if input_type == "pangraph":
            if not self.validator.validate_json_file(input_file):
                return False
        elif input_type == "gfa":
            if not self.validator.validate_gfa_file(input_file):
                return False
        elif input_type == "msa":
            if not self.validator.validate_fasta_file(input_file):
                return False

        # 验证Newick文件|Validate Newick file
        if not self.validator.validate_newick_file(self.config.newick_file):
            return False

        self.logger.info("输入文件验证通过|Input files validation passed")
        return True

    def _build_from_pangraph(self) -> bool:
        """从PanGraph构建|Build from PanGraph"""
        self.logger.info("从PanGraph构建PanMAN|Building PanMAN from PanMAN")

        import tempfile
        import shutil

        # 对于Singularity后端，将输入文件复制到本地临时目录避免NFS读取卡住
        # Conda/Docker后端直接在宿主机运行，不需要复制|Conda/Docker run on host, no copy needed
        if self.config.backend == "singularity":
            temp_dir = tempfile.mkdtemp(prefix="panman_")
            self.logger.info(f"创建临时目录|Created temp dir: {temp_dir}")

            # 复制输入文件到本地临时目录
            local_pangraph = os.path.join(temp_dir, os.path.basename(self.config.pangraph_file))
            local_newick = os.path.join(temp_dir, os.path.basename(self.config.newick_file))

            self.logger.info(f"复制输入文件到本地|Copying input files to local storage")
            shutil.copy2(self.config.pangraph_file, local_pangraph)
            shutil.copy2(self.config.newick_file, local_newick)

            pangraph_arg = local_pangraph
            newick_arg = local_newick
            container_output = f"/tmp/{self.config.output_prefix}"
            host_output = os.path.join(self.config.output_dir, self.config.output_prefix)
            cleanup_temp_dir = True
        else:
            # Conda/Docker后端：直接使用原始路径|Conda/Docker: use original paths directly
            temp_dir = None
            pangraph_arg = self.config.pangraph_file
            newick_arg = self.config.newick_file
            container_output = os.path.join(self.config.output_dir, self.config.output_prefix)
            host_output = None
            cleanup_temp_dir = False

        try:
            args = [
                "-P", pangraph_arg,
                "-N", newick_arg,
                "-o", container_output
            ]

            success, output = self.cmd_runner.run_panman_command(
                args,
                description="从PanGraph构建|Build from PanGraph"
            )

            if not success:
                return False

            # 处理输出文件|Handle output files
            # panmanUtils会在输出目录下创建panman子目录|panmanUtils creates panman subdirectory
            if self.config.backend == "singularity":
                # Singularity: 从/tmp复制文件|Singularity: Copy from /tmp
                self.logger.info(f"移动文件到最终目录|Moving files to final directory: {host_output}")
                src_file = f"{container_output}.panman"
                dst_file = f"{host_output}.panman"

                if os.path.exists(src_file):
                    shutil.move(src_file, dst_file)
                    self.logger.info(f"文件已移动|File moved: {dst_file}")
                else:
                    # 检查panman子目录|Check panman subdirectory
                    src_file = f"{container_output}/panman/{os.path.basename(container_output)}.panman"
                    if os.path.exists(src_file):
                        os.makedirs(os.path.dirname(dst_file), exist_ok=True)
                        shutil.move(src_file, dst_file)
                        self.logger.info(f"从panman子目录移动文件|Moved from panman subdirectory: {dst_file}")
                    else:
                        self.logger.warning(f"预期的输出文件不存在|Expected output file not found: {src_file}")
                        return False
            else:
                # Conda/Docker: 检查输入文件目录中的panman子目录|Conda/Docker: Check panman subdirectory in input file dir
                # 由于使用了相对路径和工作目录，输出在输入文件目录|Due to relative paths and cwd, output is in input file dir
                input_dir = os.path.dirname(self.config.pangraph_file)
                output_in_cwd = os.path.join(input_dir, "panman", f"{self.config.output_prefix}.panman")
                expected_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}.panman")

                if os.path.exists(expected_file):
                    self.logger.info(f"输出文件已生成|Output file generated: {expected_file}")
                elif os.path.exists(output_in_cwd):
                    # 移动到正确位置|Move to correct location
                    self.logger.info(f"从输入文件目录的panman子目录移动文件|Moving from input dir panman subdirectory: {output_in_cwd}")
                    os.makedirs(self.config.output_dir, exist_ok=True)
                    shutil.move(output_in_cwd, expected_file)
                    self.logger.info(f"文件已移动|File moved: {expected_file}")
                    # 清理空的panman目录|Cleanup empty panman directory
                    panman_dir = os.path.join(input_dir, "panman")
                    if os.path.exists(panman_dir) and not os.listdir(panman_dir):
                        os.rmdir(panman_dir)
                else:
                    self.logger.error(f"未找到输出文件|Output file not found: {expected_file} or {output_in_cwd}")
                    return False

            return True

        finally:
            # 清理临时目录
            if cleanup_temp_dir and temp_dir and os.path.exists(temp_dir):
                self.logger.info(f"清理临时目录|Cleaning temp dir: {temp_dir}")
                shutil.rmtree(temp_dir)

    def _build_from_gfa(self) -> bool:
        """从GFA构建|Build from GFA"""
        self.logger.info("从GFA构建PanMAN|Building PanMAN from GFA")

        import tempfile
        import shutil

        # 对于Singularity后端，将输入文件复制到本地临时目录避免NFS读取卡住
        if self.config.backend == "singularity":
            temp_dir = tempfile.mkdtemp(prefix="panman_")
            self.logger.info(f"创建临时目录|Created temp dir: {temp_dir}")

            # 复制输入文件到本地临时目录
            local_gfa = os.path.join(temp_dir, os.path.basename(self.config.gfa_file))
            local_newick = os.path.join(temp_dir, os.path.basename(self.config.newick_file))

            self.logger.info(f"复制输入文件到本地|Copying input files to local storage")
            shutil.copy2(self.config.gfa_file, local_gfa)
            shutil.copy2(self.config.newick_file, local_newick)

            gfa_arg = local_gfa
            newick_arg = local_newick
            container_output = f"/tmp/{self.config.output_prefix}"
            host_output = os.path.join(self.config.output_dir, self.config.output_prefix)
            cleanup_temp_dir = True
        else:
            temp_dir = None
            gfa_arg = self.config.gfa_file
            newick_arg = self.config.newick_file
            container_output = os.path.join(self.config.output_dir, self.config.output_prefix)
            host_output = None
            cleanup_temp_dir = False

        try:
            args = [
                "-G", gfa_arg,
                "-N", newick_arg,
                "-o", container_output
            ]

            success, output = self.cmd_runner.run_panman_command(
                args,
                description="从GFA构建|Build from GFA"
            )

            if not success:
                return False

            # 处理输出文件（同pangraph）|Handle output files (same as pangraph)
            if self.config.backend == "singularity":
                self.logger.info(f"移动文件到最终目录|Moving files to final directory: {host_output}")
                src_file = f"{container_output}.panman"
                dst_file = f"{host_output}.panman"

                if os.path.exists(src_file):
                    shutil.move(src_file, dst_file)
                    self.logger.info(f"文件已移动|File moved: {dst_file}")
                else:
                    src_file = f"{container_output}/panman/{os.path.basename(container_output)}.panman"
                    if os.path.exists(src_file):
                        os.makedirs(os.path.dirname(dst_file), exist_ok=True)
                        shutil.move(src_file, dst_file)
                        self.logger.info(f"从panman子目录移动文件|Moved from panman subdirectory: {dst_file}")
                    else:
                        self.logger.warning(f"预期的输出文件不存在|Expected output file not found: {src_file}")
                        return False
            else:
                # Conda/Docker: 检查输入文件目录中的panman子目录|Conda/Docker: Check panman subdirectory in input file dir
                input_dir = os.path.dirname(self.config.gfa_file)
                output_in_cwd = os.path.join(input_dir, "panman", f"{self.config.output_prefix}.panman")
                expected_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}.panman")

                if os.path.exists(expected_file):
                    self.logger.info(f"输出文件已生成|Output file generated: {expected_file}")
                elif os.path.exists(output_in_cwd):
                    self.logger.info(f"从输入文件目录的panman子目录移动文件|Moving from input dir panman subdirectory: {output_in_cwd}")
                    os.makedirs(self.config.output_dir, exist_ok=True)
                    shutil.move(output_in_cwd, expected_file)
                    self.logger.info(f"文件已移动|File moved: {expected_file}")
                    panman_dir = os.path.join(input_dir, "panman")
                    if os.path.exists(panman_dir) and not os.listdir(panman_dir):
                        os.rmdir(panman_dir)
                else:
                    self.logger.error(f"未找到输出文件|Output file not found: {expected_file} or {output_in_cwd}")
                    return False

            return True

        finally:
            if cleanup_temp_dir and temp_dir and os.path.exists(temp_dir):
                self.logger.info(f"清理临时目录|Cleaning temp dir: {temp_dir}")
                shutil.rmtree(temp_dir)

    def _build_from_msa(self) -> bool:
        """从MSA构建|Build from MSA"""
        self.logger.info("从MSA构建PanMAN|Building PanMAN from MSA")

        import tempfile
        import shutil

        # 对于Singularity后端，将输入文件复制到本地临时目录避免NFS读取卡住
        if self.config.backend == "singularity":
            temp_dir = tempfile.mkdtemp(prefix="panman_")
            self.logger.info(f"创建临时目录|Created temp dir: {temp_dir}")

            # 复制输入文件到本地临时目录
            local_msa = os.path.join(temp_dir, os.path.basename(self.config.msa_file))
            local_newick = os.path.join(temp_dir, os.path.basename(self.config.newick_file))

            self.logger.info(f"复制输入文件到本地|Copying input files to local storage")
            shutil.copy2(self.config.msa_file, local_msa)
            shutil.copy2(self.config.newick_file, local_newick)

            msa_arg = local_msa
            newick_arg = local_newick
            container_output = f"/tmp/{self.config.output_prefix}"
            host_output = os.path.join(self.config.output_dir, self.config.output_prefix)
            cleanup_temp_dir = True
        else:
            temp_dir = None
            msa_arg = self.config.msa_file
            newick_arg = self.config.newick_file
            container_output = os.path.join(self.config.output_dir, self.config.output_prefix)
            host_output = None
            cleanup_temp_dir = False

        try:
            args = [
                "-M", msa_arg,
                "-N", newick_arg,
                "-o", container_output
            ]

            success, output = self.cmd_runner.run_panman_command(
                args,
                description="从MSA构建|Build from MSA"
            )

            if not success:
                return False

            # 处理输出文件（同pangraph）|Handle output files (same as pangraph)
            if self.config.backend == "singularity":
                self.logger.info(f"移动文件到最终目录|Moving files to final directory: {host_output}")
                src_file = f"{container_output}.panman"
                dst_file = f"{host_output}.panman"

                if os.path.exists(src_file):
                    shutil.move(src_file, dst_file)
                    self.logger.info(f"文件已移动|File moved: {dst_file}")
                else:
                    src_file = f"{container_output}/panman/{os.path.basename(container_output)}.panman"
                    if os.path.exists(src_file):
                        os.makedirs(os.path.dirname(dst_file), exist_ok=True)
                        shutil.move(src_file, dst_file)
                        self.logger.info(f"从panman子目录移动文件|Moved from panman subdirectory: {dst_file}")
                    else:
                        self.logger.warning(f"预期的输出文件不存在|Expected output file not found: {src_file}")
                        return False
            else:
                # Conda/Docker: 检查输入文件目录中的panman子目录|Conda/Docker: Check panman subdirectory in input file dir
                input_dir = os.path.dirname(self.config.msa_file)
                output_in_cwd = os.path.join(input_dir, "panman", f"{self.config.output_prefix}.panman")
                expected_file = os.path.join(self.config.output_dir, f"{self.config.output_prefix}.panman")

                if os.path.exists(expected_file):
                    self.logger.info(f"输出文件已生成|Output file generated: {expected_file}")
                elif os.path.exists(output_in_cwd):
                    self.logger.info(f"从输入文件目录的panman子目录移动文件|Moving from input dir panman subdirectory: {output_in_cwd}")
                    os.makedirs(self.config.output_dir, exist_ok=True)
                    shutil.move(output_in_cwd, expected_file)
                    self.logger.info(f"文件已移动|File moved: {expected_file}")
                    panman_dir = os.path.join(input_dir, "panman")
                    if os.path.exists(panman_dir) and not os.listdir(panman_dir):
                        os.rmdir(panman_dir)
                else:
                    self.logger.error(f"未找到输出文件|Output file not found: {expected_file} or {output_in_cwd}")
                    return False

            return True

        finally:
            if cleanup_temp_dir and temp_dir and os.path.exists(temp_dir):
                self.logger.info(f"清理临时目录|Cleaning temp dir: {temp_dir}")
                shutil.rmtree(temp_dir)
