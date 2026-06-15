"""
PanGraph构建器模块|PanGraph Builder Module
用于从FASTA文件生成PanGraph JSON格式
"""

import os
import subprocess
from typing import Tuple, Optional
from ..common.paths import expand_path


class PanGraphBuilder:
    """PanGraph构建器|PanGraph Builder"""

    def __init__(self, logger, pangraph_path: str = None,
                 pangraph_sif: str = None, singularity_path: str = None):
        """
        初始化PanGraph构建器|Initialize PanGraph Builder

        Args:
            logger: 日志记录器|Logger
            pangraph_path: PanGraph可执行文件路径|PanGraph executable path
            pangraph_sif: PanGraph Singularity SIF镜像路径|PanGraph Singularity SIF image path
            singularity_path: Singularity可执行文件路径|Singularity executable path
        """
        self.logger = logger

        # 保存原始参数以判断是否明确提供|Save original params to check if explicitly provided
        self.pangraph_sif_explicit = pangraph_sif is not None
        self.singularity_path_explicit = singularity_path is not None

        # 标准化路径|Normalize paths
        if pangraph_sif:
            self.pangraph_sif = os.path.normpath(os.path.abspath(pangraph_sif))
        else:
            # 设置默认 PanGraph SIF 镜像路径|Set default PanGraph SIF image path
            default_pangraph_sif = expand_path("~/software/singularity/pangraph_0.7.3.sif")
            if os.path.exists(default_pangraph_sif):
                self.pangraph_sif = default_pangraph_sif
            else:
                self.pangraph_sif = None

        if singularity_path:
            self.singularity_path = os.path.normpath(os.path.abspath(singularity_path))
        else:
            # 设置默认 singularity 路径|Set default singularity path
            default_singularity = expand_path("~/miniforge3/envs/singularity_v.3.8.7/bin/singularity")
            if os.path.exists(default_singularity):
                self.singularity_path = default_singularity
            else:
                self.singularity_path = None

        # 设置 pangraph_path，优先使用明确指定的路径或自动查找的路径
        # 优先级：明确指定的 pangraph_path > 明确指定的 SIF > 本地可执行文件 > 默认 SIF
        use_sif = False

        # 1. 如果明确指定了 pangraph_path，使用它
        if pangraph_path:
            self.pangraph_path = pangraph_path
        # 2. 如果明确指定了 SIF 镜像，使用 SIF
        elif self.pangraph_sif_explicit and self.pangraph_sif and self.singularity_path:
            if os.path.exists(self.pangraph_sif) and os.path.exists(self.singularity_path):
                use_sif = True
                self.pangraph_path = None
                self.logger.info(f"使用 SIF 镜像中的 PanGraph|Using PanGraph from SIF image: {self.pangraph_sif}")
            else:
                self.logger.warning("指定的 SIF 镜像或 singularity 可执行文件不存在|"
                                   "Specified SIF image or singularity executable does not exist")
                self.pangraph_path = self._find_pangraph()
        # 3. 查找本地可执行文件
        else:
            self.pangraph_path = self._find_pangraph()
            # 如果找到本地可执行文件，使用它；否则尝试默认 SIF
            if self.pangraph_path and os.path.exists(self.pangraph_path):
                self.logger.info(f"使用本地 PanGraph 可执行文件|Using local PanGraph executable: {self.pangraph_path}")
            elif self.pangraph_sif and self.singularity_path:
                if os.path.exists(self.pangraph_sif) and os.path.exists(self.singularity_path):
                    use_sif = True
                    self.pangraph_path = None
                    self.logger.info(f"使用默认 SIF 镜像中的 PanGraph|Using PanGraph from default SIF image: {self.pangraph_sif}")

    def _find_pangraph(self) -> str:
        """查找PanGraph可执行文件|Find PanGraph executable"""
        # 常见安装位置
        possible_paths = [
            expand_path("~/software/pangraph/pangraph-0.7.3/pangraph/bin/pangraph"),
            expand_path("~/software/pangraph/pangraph"),
            "/usr/local/bin/pangraph",
            "/usr/bin/pangraph",
        ]

        # 检查环境变量
        env_path = os.environ.get("PANGRAPH_PATH")
        if env_path and os.path.exists(env_path):
            return env_path

        # 检查常见位置
        for path in possible_paths:
            if os.path.exists(path):
                return path

        # 检查 PATH
        for path in os.environ.get("PATH", "").split(":"):
            full_path = os.path.join(path, "pangraph")
            if os.path.exists(full_path):
                return full_path

        return None

    def build_from_fasta(self, fasta_file: str, output_prefix: str,
                        output_dir: str, threads: int = 8) -> Tuple[bool, dict]:
        """
        从FASTA文件构建PanGraph|Build PanGraph from FASTA file

        Args:
            fasta_file: 输入FASTA文件路径|Input FASTA file path
            output_prefix: 输出文件前缀|Output file prefix
            output_dir: 输出目录|Output directory
            threads: 线程数|Number of threads

        Returns:
            (成功状态, 输出文件路径字典)|(Success status, Output file paths dict)
        """
        import tempfile
        import shutil

        self.logger.info("开始构建PanGraph|Starting PanGraph build")
        self.logger.info(f"输入FASTA|Input FASTA: {fasta_file}")
        self.logger.info(f"输出前缀|Output prefix: {output_prefix}")

        # 验证输入文件
        if not os.path.exists(fasta_file):
            self.logger.error(f"FASTA文件不存在|FASTA file does not exist: {fasta_file}")
            return False, {}

        # 检查使用哪种方式运行 PanGraph
        # pangraph_path 为 None 表示应该使用 SIF 镜像
        use_singularity = self.pangraph_path is None

        if use_singularity:
            if not self.pangraph_sif or not os.path.exists(self.pangraph_sif):
                self.logger.error(
                    f"找不到 PanGraph SIF 镜像|Cannot find PanGraph SIF image: {self.pangraph_sif}"
                )
                return False, {}
            if not self.singularity_path or not os.path.exists(self.singularity_path):
                self.logger.error(
                    f"找不到 Singularity 可执行文件|Cannot find Singularity executable: {self.singularity_path}"
                )
                return False, {}
        else:
            if not self.pangraph_path or not os.path.exists(self.pangraph_path):
                self.logger.error(
                    f"找不到 PanGraph 可执行文件|Cannot find PanGraph executable: {self.pangraph_path}\n"
                    f"请设置环境变量 PANGRAPH_PATH、使用 --pangraph-path 参数，或提供 --pangraph-sif 镜像|"
                    f"Please set PANGRAPH_PATH, use --pangraph-path, or provide --pangraph-sif image"
                )
                return False, {}

        # 对于singularity后端，复制FASTA到本地临时目录以避免NFS问题
        temp_dir = None
        actual_fasta_file = fasta_file
        actual_output_dir = output_dir

        if use_singularity:
            temp_dir = tempfile.mkdtemp(prefix="pangraph_")
            actual_fasta_file = os.path.join(temp_dir, os.path.basename(fasta_file))
            actual_output_dir = temp_dir
            shutil.copy2(fasta_file, actual_fasta_file)
            self.logger.info(f"复制FASTA到临时目录|Copied FASTA to temp dir: {actual_fasta_file}")

        # 设置输出文件路径
        json_output = os.path.join(actual_output_dir, f"{output_prefix}.json")
        newick_output = os.path.join(actual_output_dir, f"{output_prefix}.nwk")
        log_output = os.path.join(actual_output_dir, f"{output_prefix}_pangraph.log")

        # 构建命令
        if use_singularity:
            cmd = [
                self.singularity_path,
                "exec",
                self.pangraph_sif,
                "pangraph",
                "build",
                actual_fasta_file  # 使用实际的FASTA文件路径
            ]
            self.logger.info(f"Singularity路径|Singularity path: {self.singularity_path}")
            self.logger.info(f"PanGraph SIF镜像|PanGraph SIF image: {self.pangraph_sif}")
        else:
            cmd = [
                self.pangraph_path,
                "build",
                actual_fasta_file  # 使用实际的FASTA文件路径
            ]
            self.logger.info(f"PanGraph路径|PanGraph path: {self.pangraph_path}")

        # 设置环境变量
        env = os.environ.copy()
        env["JULIA_NUM_THREADS"] = str(threads)

        self.logger.info(f"执行命令|Executing command: {' '.join(cmd)}")

        try:
            # 运行 PanGraph
            with open(json_output, 'w') as json_out, open(log_output, 'w') as log_out:
                process = subprocess.run(
                    cmd,
                    stdout=json_out,  # JSON 输出到文件
                    stderr=log_out,  # 日志输出到文件
                    env=env,
                    check=False,
                    text=False
                )

            # 检查是否成功
            if process.returncode != 0:
                self.logger.error(f"PanGraph构建失败|PanGraph build failed with return code {process.returncode}")
                # 清理临时目录
                if temp_dir:
                    shutil.rmtree(temp_dir)
                return False, {}

            # 从日志中提取 Newick 树
            if not self._extract_newick(log_output, newick_output):
                self.logger.warning("未能从日志中提取Newick树|Failed to extract Newick tree from log")

            # 验证输出文件
            if not os.path.exists(json_output) or os.path.getsize(json_output) == 0:
                self.logger.error(f"JSON输出文件为空|JSON output file is empty: {json_output}")
                # 清理临时目录
                if temp_dir:
                    shutil.rmtree(temp_dir)
                return False, {}

            # 如果使用了临时目录，移动输出文件到目标目录
            final_json = json_output
            final_newick = newick_output if os.path.exists(newick_output) else None
            final_log = log_output

            if temp_dir:
                os.makedirs(output_dir, exist_ok=True)
                final_json = os.path.join(output_dir, os.path.basename(json_output))
                shutil.move(json_output, final_json)
                self.logger.info(f"移动JSON输出|Moved JSON output: {final_json}")

                if os.path.exists(newick_output):
                    final_newick = os.path.join(output_dir, os.path.basename(newick_output))
                    shutil.move(newick_output, final_newick)
                    self.logger.info(f"移动Newick输出|Moved Newick output: {final_newick}")

                final_log = os.path.join(output_dir, os.path.basename(log_output))
                shutil.move(log_output, final_log)
                self.logger.info(f"移动日志输出|Moved log output: {final_log}")

                # 清理临时目录
                shutil.rmtree(temp_dir)

            json_size = os.path.getsize(final_json)
            self.logger.info(f"PanGraph构建成功|PanGraph build successful")
            self.logger.info(f"  JSON文件|JSON file: {final_json} ({json_size} bytes)")

            if final_newick and os.path.exists(final_newick):
                newick_size = os.path.getsize(final_newick)
                self.logger.info(f"  Newick文件|Newick file: {final_newick} ({newick_size} bytes)")

            return True, {
                'json': final_json,
                'newick': final_newick,
                'log': final_log
            }

        except Exception as e:
            self.logger.error(f"PanGraph构建出错|Error during PanGraph build: {str(e)}")
            # 清理临时目录
            if temp_dir:
                shutil.rmtree(temp_dir)
            return False, {}

    def _extract_newick(self, log_file: str, newick_output: str) -> bool:
        """
        从PanGraph日志中提取Newick树|Extract Newick tree from PanGraph log

        Args:
            log_file: 日志文件路径|Log file path
            newick_output: Newick输出文件路径|Newick output file path

        Returns:
            是否成功提取|Whether extraction was successful
        """
        try:
            with open(log_file, 'r') as f:
                content = f.read()

            # 查找包含 "tree:" 的行
            for line in content.split('\n'):
                if 'tree:' in line:
                    # 提取 Newick 字符串
                    parts = line.split('tree:')
                    if len(parts) > 1:
                        newick_str = parts[1].strip()
                        with open(newick_output, 'w') as f:
                            f.write(newick_str + '\n')
                        self.logger.info("成功提取Newick树|Successfully extracted Newick tree")
                        return True

            return False

        except Exception as e:
            self.logger.warning(f"提取Newick树时出错|Error extracting Newick tree: {str(e)}")
            return False

    def validate_fasta(self, fasta_file: str) -> bool:
        """
        验证FASTA文件格式|Validate FASTA file format

        Args:
            fasta_file: FASTA文件路径|FASTA file path

        Returns:
            是否有效|Whether valid
        """
        try:
            with open(fasta_file, 'r') as f:
                first_line = f.readline().strip()
                if not first_line.startswith('>'):
                    self.logger.error(f"无效的FASTA文件|Invalid FASTA file: {fasta_file}")
                    return False

            # 统计序列数
            seq_count = 0
            with open(fasta_file, 'r') as f:
                for line in f:
                    if line.startswith('>'):
                        seq_count += 1

            self.logger.info(f"FASTA文件包含 {seq_count} 条序列|FASTA file contains {seq_count} sequences")
            return True

        except Exception as e:
            self.logger.error(f"读取FASTA文件失败|Failed to read FASTA file: {str(e)}")
            return False
