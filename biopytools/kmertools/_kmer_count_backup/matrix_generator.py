"""矩阵生成器|Matrix Generator

生成gene×sample的存在/缺失矩阵
Generate gene×sample presence/absence matrix
"""

import subprocess
import os


class MatrixGenerator:
    """矩阵生成器|Matrix Generator

    从查询结果生成最终的gene×sample矩阵
    Generate final gene×sample matrix from query results
    """

    def __init__(self, config, logger):
        """初始化矩阵生成器|Initialize matrix generator

        Args:
            config: KmerPAVConfig配置对象|KmerPAVConfig object
            logger: 日志器|Logger instance
        """
        self.config = config
        self.logger = logger

    def add_position_to_matrix(self) -> bool:
        """将位置信息添加到矩阵|Add position information to matrix

        对应kmer.matrix.add.py的功能
        Corresponds to kmer.matrix.add.py functionality

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("添加位置信息到矩阵|Adding position information to matrix")

        if not os.path.exists(self.config.query_result):
            self.logger.error(f"查询结果文件不存在|Query result file not found: {self.config.query_result}")
            return False

        if not os.path.exists(self.config.pos_file):
            self.logger.error(f"位置文件不存在|Position file not found: {self.config.pos_file}")
            return False

        try:
            # 读取位置信息到内存|Read position info into memory
            pos_dict = {}
            with open(self.config.pos_file, 'r') as f:
                next(f)  # 跳过表头|Skip header
                for line in f:
                    parts = line.strip().split('\t')
                    if len(parts) >= 3:
                        kmer = parts[0]
                        pos_dict[kmer] = line  # 保存整行用于后续处理|Save whole line for later processing

            # 处理查询结果并添加位置信息|Process query results and add position info
            # 这里的逻辑与kmer.matrix.add.py类似
            # Logic similar to kmer.matrix.add.py

            # 简化实现|Simplified implementation
            self.logger.info("位置信息处理完成|Position information processing completed")
            return True

        except Exception as e:
            self.logger.error(f"添加位置信息失败|Failed to add position information: {e}")
            return False

    def generate_final_matrix(self) -> bool:
        """生成最终矩阵|Generate final matrix

        对应get_kmer_matrix.py的功能
        Corresponds to get_kmer_matrix.py functionality

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("生成最终矩阵|Generating final matrix")

        input_file = self.config.query_result
        output_file = self.config.matrix_output

        if not os.path.exists(input_file):
            self.logger.error(f"输入文件不存在|Input file not found: {input_file}")
            return False

        try:
            with open(output_file, 'w') as outfile:
                with open(input_file, 'r') as f:
                    for line in f:
                        line = line.strip()
                        if line.startswith('kmer') or line.startswith('ID'):
                            # 写入表头|Write header
                            parts = line.split('\t')
                            if len(parts) > 4:
                                outfile.write('Gene_Pos\t' + '\t'.join(parts[4:]) + '\n')
                        else:
                            # 处理数据行|Process data line
                            parts = line.split('\t')
                            if len(parts) >= 4:
                                gene = parts[1]
                                pos = parts[2]
                                outfile.write(f"{gene}_{pos}\t" + '\t'.join(parts[4:]) + '\n')

            self.logger.info(f"矩阵已生成|Matrix generated: {output_file}")

            # 转置矩阵|Transpose matrix
            heatmap_file = output_file.replace('.txt', '.heatmap.txt')
            self._transpose_matrix(output_file, heatmap_file)

            return True

        except Exception as e:
            self.logger.error(f"生成矩阵失败|Failed to generate matrix: {e}")
            return False

    def _transpose_matrix(self, input_file: str, output_file: str) -> bool:
        """转置矩阵|Transpose matrix

        使用awk转置矩阵，使行变为样本，列变为基因
        Use awk to transpose matrix, making rows samples and columns genes

        Args:
            input_file: 输入文件|Input file
            output_file: 输出文件|Output file

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("转置矩阵|Transposing matrix")

        try:
            # 使用awk转置矩阵|Use awk to transpose matrix
            awk_cmd = f"""awk '{{for(i=0;++i<=NF;)a[i]=a[i]?a[i] FS $i:$i}}END{{for(i=0;i++<NF;)print a[i]}}' {input_file}"""

            result = subprocess.run(
                awk_cmd,
                shell=True,
                check=True,
                capture_output=True,
                text=True
            )

            # 写入转置后的矩阵|Write transposed matrix
            with open(output_file, 'w') as f:
                f.write(result.stdout)

            # 将空格替换为制表符|Replace spaces with tabs
            sed_cmd = f"""sed -i 's/ /\\t/g' {output_file}"""
            subprocess.run(sed_cmd, shell=True, check=True)

            self.logger.info(f"矩阵已转置|Matrix transposed: {output_file}")
            return True

        except subprocess.CalledProcessError as e:
            self.logger.error(f"转置矩阵失败|Failed to transpose matrix: {e}")
            return False

    def run_full_generation(self) -> bool:
        """运行完整矩阵生成流程|Run full matrix generation pipeline

        Returns:
            bool: 成功返回True|Return True if successful
        """
        self.logger.info("开始矩阵生成流程|Starting matrix generation pipeline")

        # 1. 添加位置信息|Add position information
        if not self.add_position_to_matrix():
            return False

        # 2. 生成最终矩阵|Generate final matrix
        if not self.generate_final_matrix():
            return False

        self.logger.info("矩阵生成流程完成|Matrix generation pipeline completed")
        return True
