"""
TASSEL GWAS分析核心模块 | TASSEL GWAS Analysis Core Module
"""

import os
import subprocess
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from .utils import TASSELLogger, check_dependencies, validate_vcf_file, validate_phenotype_file, prepare_phenotype_file, create_statistics_report


class TASSELGWASConfig:
    """TASSEL GWAS分析配置类 | TASSEL GWAS Analysis Configuration Class"""

    def __init__(self, **kwargs):
        # 基本参数 | Basic parameters
        self.vcf_file = Path(kwargs.get('vcf_file', '')).expanduser()
        self.pheno_file = Path(kwargs.get('pheno_file', '')).expanduser()
        self.output_prefix = kwargs.get('output_prefix', 'GWAS_Result')
        self.output_dir = Path(kwargs.get('output_dir', '.')).expanduser()

        # TASSEL路径 | TASSEL path
        self.tassel_path = kwargs.get('tassel_path', None)
        self.memory_min = kwargs.get('memory_min', '10g')
        self.memory_max = kwargs.get('memory_max', '100g')

        # 分析参数 | Analysis parameters
        self.model = kwargs.get('model', 'MLM').upper()
        self.threads = kwargs.get('threads', 4)
        self.maf_filter = kwargs.get('maf_filter', None)
        self.miss_filter = kwargs.get('miss_filter', None)
        self.skip_sort = kwargs.get('skip_sort', False)
        self.pca_components = kwargs.get('pca_components', 5)  # PCA主成分数量，默认为5

        # 矩阵文件 | Matrix files
        self.q_matrix = Path(kwargs.get('q_matrix', '')).expanduser() if kwargs.get('q_matrix') else None
        self.kinship = Path(kwargs.get('kinship', '')).expanduser() if kwargs.get('kinship') else None

        # 其他选项 | Other options
        self.keep_temp = kwargs.get('keep_temp', False)
        self.log_level = kwargs.get('log_level', 'INFO').upper()

    def validate(self):
        """验证配置参数 | Validate configuration parameters"""
        if not self.vcf_file:
            raise ValueError("VCF文件不能为空 | VCF file cannot be empty")

        if not self.pheno_file:
            raise ValueError("表型文件不能为空 | Phenotype file cannot be empty")

        if self.model not in ['GLM', 'MLM', 'BOTH']:
            raise ValueError("模型必须是 GLM、MLM 或 BOTH | Model must be GLM, MLM, or BOTH")

        # 确保输出目录存在 | Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)


class TASSELGWASAnalyzer:
    """TASSEL GWAS分析器 | TASSEL GWAS Analyzer"""

    def __init__(self, **kwargs):
        """初始化TASSEL GWAS分析器 | Initialize TASSEL GWAS analyzer"""
        self.config = TASSELGWASConfig(**kwargs)
        self.config.validate()

        # 设置日志 | Setup logging
        self.logger_manager = TASSELLogger(self.config.output_dir, f"{self.config.output_prefix}.pipeline.log")
        self.logger = self.logger_manager.get_logger()

        # 初始化TASSEL路径 | Initialize TASSEL path
        self.tassel_cmd = self._get_tassel_command()

        # 统计信息 | Statistics
        self.stats = {}
        self.temp_files = []

    def _get_tassel_command(self):
        """获取TASSEL命令路径 | Get TASSEL command path"""
        if self.config.tassel_path:
            tassel_path = Path(self.config.tassel_path).expanduser()
            if tassel_path.exists():
                return str(tassel_path)

        # 使用依赖检查找到的路径 | Use path found by dependency check
        success, tassel_path = check_dependencies(self.logger)
        if success and tassel_path:
            return str(tassel_path)

        raise RuntimeError("未找到有效的TASSEL安装 | No valid TASSEL installation found")

    def check_dependencies(self):
        """检查依赖软件 | Check dependencies"""
        return check_dependencies(self.logger)[0]

    def validate_input_files(self):
        """验证输入文件 | Validate input files"""
        self.logger.info("🔍 验证输入文件 | Validating input files")

        # 验证VCF文件 | Validate VCF file
        if not validate_vcf_file(self.config.vcf_file, self.logger):
            return False

        # 验证表型文件 | Validate phenotype file
        if not validate_phenotype_file(self.config.pheno_file, self.logger):
            return False

        # 验证Q矩阵文件 | Validate Q matrix file
        if self.config.q_matrix and self.config.q_matrix.exists():
            q_dims = self._count_matrix_columns(self.config.q_matrix)
            self.logger.info(f"✅ Q矩阵已提供 | Q matrix provided: {self.config.q_matrix} (K={q_dims})")

        # 验证Kinship矩阵文件 | Validate Kinship matrix file
        if self.config.kinship and self.config.kinship.exists():
            self.logger.info(f"✅ Kinship矩阵已提供 | Kinship matrix provided: {self.config.kinship}")

        return True

    def _count_matrix_columns(self, file_path: Path) -> int:
        """计算矩阵文件的列数 | Count matrix file columns"""
        try:
            with open(file_path, 'r') as f:
                first_line = f.readline().strip()
                return len(first_line.split('\t'))
        except Exception:
            return 0

    def _get_vcf_stats(self, vcf_file: Path) -> Tuple[int, int]:
        """获取VCF文件统计信息 | Get VCF file statistics (using bash logic)"""
        try:
            # 使用bash脚本的逻辑 | Use bash script logic
            # VCF_VARIANTS=$(zcat -f "$VCF_FILE" 2>/dev/null | grep -v "^#" | wc -l)
            if vcf_file.suffix.lower() == '.gz':
                variant_cmd = 'zcat -f "{}" 2>/dev/null | grep -v "^#" | wc -l'.format(str(vcf_file))
            else:
                variant_cmd = 'cat "{}" 2>/dev/null | grep -v "^#" | wc -l'.format(str(vcf_file))

            variant_result = subprocess.run(variant_cmd, shell=True, capture_output=True, text=True, timeout=120)

            if variant_result.returncode == 0:
                variant_count = int(variant_result.stdout.strip())
            else:
                variant_count = 0

            # VCF_SAMPLES=$(zcat -f "$VCF_FILE" 2>/dev/null | grep "^#CHROM" | awk '{print NF-9}')
            if vcf_file.suffix.lower() == '.gz':
                sample_cmd = 'zcat -f "{}" 2>/dev/null | grep "^#CHROM" | awk \'{{print NF-9}}\''.format(str(vcf_file))
            else:
                sample_cmd = 'cat "{}" 2>/dev/null | grep "^#CHROM" | awk \'{{print NF-9}}\''.format(str(vcf_file))

            sample_result = subprocess.run(sample_cmd, shell=True, capture_output=True, text=True, timeout=120)

            if sample_result.returncode == 0 and sample_result.stdout.strip():
                sample_count = int(sample_result.stdout.strip())
            else:
                sample_count = 0

            return variant_count, sample_count

        except Exception as e:
            self.logger.warning(f"⚠️ 无法获取VCF统计信息 | Cannot get VCF statistics: {e}")
            return 0, 0

    def _get_pheno_stats(self, pheno_file: Path) -> Tuple[int, int]:
        """获取表型文件统计信息 | Get phenotype file statistics"""
        try:
            with open(pheno_file, 'r') as f:
                lines = f.readlines()

            header_cols = len(lines[0].strip().split('\t'))
            data_rows = len([line for line in lines[1:] if line.strip() and not line.startswith('#')])

            return data_rows, header_cols - 1  # 样本数, 表型数

        except Exception as e:
            self.logger.warning(f"⚠️ 无法获取表型统计信息 | Cannot get phenotype statistics: {e}")
            return 0, 0

    def prepare_phenotype_file(self) -> Path:
        """准备表型文件 | Prepare phenotype file"""
        pheno_ready = self.config.output_dir / f"{self.config.output_prefix}.pheno_ready.txt"
        # 注意：不将pheno_ready添加到temp_files，因为TASSEL需要访问它
        # Note: don't add pheno_ready to temp_files as TASSEL needs to access it

        if prepare_phenotype_file(self.config.pheno_file, pheno_ready, self.logger):
            return pheno_ready
        else:
            raise RuntimeError("表型文件准备失败 | Phenotype file preparation failed")

    def sort_and_filter_vcf(self) -> Path:
        """VCF排序和过滤 | VCF sorting and filtering"""
        if self.config.skip_sort:
            self.logger.info("📄 跳过VCF排序 | Skipping VCF sorting (user requested)")
            return self.config.vcf_file

        vcf_sorted = self.config.output_dir / f"{self.config.output_prefix}.sorted.vcf"
        self.temp_files.append(vcf_sorted)

        try:
            # 检查是否有bcftools | Check if bcftools is available
            subprocess.run(['which', 'bcftools'], check=True, capture_output=True)

            filter_cmd = ['bcftools', 'view']
            if self.config.maf_filter:
                filter_cmd.extend(['-q', str(self.config.maf_filter)])
            if self.config.miss_filter:
                filter_cmd.extend(['-i', f'F_MISSING<{self.config.miss_filter}'])

            filter_cmd.append(str(self.config.vcf_file))

            # 使用bcftools排序 | Sort with bcftools
            sort_process = subprocess.Popen(filter_cmd, stdout=subprocess.PIPE)
            sort_cmd = ['bcftools', 'sort', '-O', 'v', '-T', str(self.config.output_dir), '-o', str(vcf_sorted)]
            subprocess.run(sort_cmd, stdin=sort_process.stdout, check=True)

            # 统计过滤后的变异位点数 | Count variants after filtering
            result = subprocess.run(['bcftools', 'stats', str(vcf_sorted)],
                                  capture_output=True, text=True)
            # 简单统计，假设每个非注释行都是一个变异
            variant_count = len([line for line in subprocess.run(['grep', '-v', '^#', str(vcf_sorted)],
                                                            capture_output=True, text=True).stdout.split('\n') if line.strip()])

            self.logger.info(f"✅ VCF排序完成 | VCF sorting completed: {vcf_sorted}")
            self.logger.info(f"   过滤后变异位点 | Variants after filter: {variant_count}")

            return vcf_sorted

        except subprocess.CalledProcessError as e:
            self.logger.warning(f"⚠️ bcftools不可用，使用原始VCF | bcftools not available, using original VCF")
            return self.config.vcf_file
        except FileNotFoundError:
            self.logger.warning(f"⚠️ bcftools未找到，使用原始VCF | bcftools not found, using original VCF")
            return self.config.vcf_file

    def calculate_pca(self, vcf_file: Path) -> Path:
        """计算PCA主成分作为协变量 | Calculate PCA eigenvectors as covariates"""
        # 在子目录中运行PCA，但VCF2PCACluster会在当前工作目录生成文件
        pca_prefix = f"{self.config.output_prefix}_pca"
        pca_eigenvec = self.config.output_dir / f"{pca_prefix}.eigenvec"

        self.logger.info("🧮 计算PCA主成分（MLM协变量）| Calculating PCA eigenvectors (for MLM covariates)...")
        self.logger.info(f"   使用工具 | Using tool: VCF2PCACluster")
        self.logger.info(f"   输入VCF | Input VCF: {vcf_file}")
        self.logger.info(f"   PCA输出 | PCA output: {pca_eigenvec}")

        try:
            # 检查VCF2PCACluster是否可用
            vcf2pca_path = "/share/org/YZWL/yzwl_lixg/software/VCF2PCACluster-1.42/bin/VCF2PCACluster"
            if not Path(vcf2pca_path).exists():
                self.logger.error(f"❌ VCF2PCACluster未找到 | VCF2PCACluster not found: {vcf2pca_path}")
                raise RuntimeError("VCF2PCACluster tool not available")

            # 构建VCF2PCACluster命令 - 使用绝对路径VCF，在子目录中执行
            cmd = [
                vcf2pca_path,
                '-InVCF', str(vcf_file.absolute()),  # 使用绝对路径
                '-OutPut', pca_prefix,  # 输出到子目录
                '-Threads', str(self.config.threads)
            ]

            shell_cmd = ' '.join(cmd)
            self.logger.info(f"   执行命令 | Executing command: {shell_cmd}")
            self.logger.info(f"   工作目录 | Working directory: {self.config.output_dir}")

            start_time = time.time()
            result = subprocess.run(shell_cmd, capture_output=True, text=True,
                                  shell=True, timeout=36000,  # 10小时超时
                                  cwd=str(self.config.output_dir))

            self.logger.info(f"   命令返回码 | Command return code: {result.returncode}")
            if result.stdout:
                self.logger.info(f"   标准输出 | Stdout: {result.stdout}")
            if result.stderr:
                self.logger.warning(f"   标准错误 | Stderr: {result.stderr}")

            end_time = time.time()
            pca_time = int(end_time - start_time)

            if pca_eigenvec.exists():
                file_size = pca_eigenvec.stat().st_size
                if file_size > 0:
                    self.logger.info(f"✅ PCA计算完成 | PCA calculation completed in {pca_time}s")
                    self.logger.info(f"   文件大小 | File size: {file_size} bytes")

                    # 统计样本和PC数量
                    with open(pca_eigenvec, 'r') as f:
                        lines = f.readlines()
                    sample_count = len(lines) - 1  # 减去表头
                    self.logger.info(f"   样本数 | Samples: {sample_count}")
                    self.logger.info(f"   PC数量 | PC components: 10 (default)")

                    # 格式化PCA文件作为TASSEL协变量
                    self.logger.info("🔄 格式化PCA文件作为TASSEL协变量 | Formatting PCA file as TASSEL covariate...")
                    pca_formatted = self._format_pca_for_tassel(pca_eigenvec)

                    return pca_formatted
                else:
                    self.logger.error(f"   PCA文件为空 | PCA file is empty: {pca_eigenvec}")
                    raise RuntimeError("PCA calculation failed, output file is empty")
            else:
                self.logger.error(f"   PCA文件不存在 | PCA file does not exist: {pca_eigenvec}")
                raise RuntimeError("PCA calculation failed, output file not generated")

        except subprocess.TimeoutExpired:
            self.logger.error(f"   PCA计算超时 | PCA calculation timed out")
            raise RuntimeError("PCA calculation timed out")
        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ PCA计算失败 | PCA calculation failed: {e.stderr}")
            raise RuntimeError(f"PCA calculation failed: {e.stderr}")
        except Exception as e:
            self.logger.error(f"❌ PCA计算过程中发生错误 | Error during PCA calculation: {e}")
            raise RuntimeError(f"PCA calculation error: {e}")

    def calculate_kinship(self, vcf_file: Path) -> Path:
        """计算Kinship矩阵 | Calculate Kinship matrix"""
        kinship_file = self.config.output_dir / f"{self.config.output_prefix}.kinship.txt"
        if not self.config.keep_temp:
            self.temp_files.append(kinship_file)

        # 检查是否已经有PCA生成的Kinship文件
        pca_kinship_file = self.config.output_dir / f"{self.config.output_prefix}_pca.Normalized_IBS.Kinship"
        if pca_kinship_file.exists() and pca_kinship_file.stat().st_size > 0:
            self.logger.info(f"   -> 发现PCA生成的Kinship矩阵 | Found PCA-generated Kinship matrix: {pca_kinship_file.name}")
            # 复制到期望的文件名
            import shutil
            shutil.copy2(pca_kinship_file, kinship_file)
            self.logger.info(f"   -> 复制到 | Copied to: {kinship_file.name}")
            return kinship_file
        elif self.config.kinship and self.config.kinship.exists():
            self.logger.info(f"   -> 使用提供的Kinship矩阵 | Using provided Kinship matrix: {self.config.kinship}")
            shutil.copy2(self.config.kinship, kinship_file)
            self.logger.info(f"   -> 复制到 | Copied to: {kinship_file.name}")
            return kinship_file

        self.logger.info("🧬 计算Kinship矩阵 | Calculating Kinship matrix (Centered_IBS)...")
        start_time = time.time()

        try:
            self.logger.debug(f"   运行目录 | Working directory: {self.config.output_dir}")
            self.logger.debug(f"   VCF文件 | VCF file: {vcf_file}")
            self.logger.debug(f"   Kinship输出 | Kinship output: {kinship_file}")

            # 使用用户指定的内存设置，不自动修改
            memory_max = self.config.memory_max
            memory_min = self.config.memory_min

            # 检查内存设置是否合理
            memory_gb = int(memory_max[:-1]) if memory_max.endswith('g') else 0
            if memory_gb < 20:
                self.logger.warning(f"⚠️ 内存设置可能不足 ({memory_max})，对于大规模数据建议至少20GB")

            cmd = [
                'perl', self.tassel_cmd,
                '-Xms' + memory_min,
                '-Xmx' + memory_max,
                '-importGuess', str(vcf_file.absolute()),
                '-KinshipPlugin', '-method', 'Centered_IBS', '-endPlugin',
                '-export', str(kinship_file.absolute()),
                '-exportType', 'SqrMatrix'
            ]

            self.logger.debug(f"   TASSEL命令 | TASSEL command: {' '.join(cmd)}")
            self.logger.info(f"   内存设置 | Memory settings: -Xms{memory_min} -Xmx{memory_max}")
            self.logger.info(f"   工作目录 | Working directory: {self.config.output_dir}")
            self.logger.info(f"   VCF文件 | VCF file: {vcf_file}")
            self.logger.info(f"   Kinship输出 | Kinship output: {kinship_file}")

            # 验证输入文件存在
            if not vcf_file.exists():
                raise RuntimeError(f"VCF文件不存在 | VCF file does not exist: {vcf_file}")

            # 确保输出目录存在
            self.config.output_dir.mkdir(parents=True, exist_ok=True)

            # 先检查基本命令是否可用
            try:
                test_cmd = ['perl', self.tassel_cmd, '-h']
                test_result = subprocess.run(test_cmd, capture_output=True, text=True, timeout=10)
                if test_result.returncode != 0:
                    self.logger.error(f"   TASSEL命令不可用 | TASSEL command not working: {test_result.stderr}")
                    raise RuntimeError(f"TASSEL command test failed: {test_result.stderr}")
            except Exception as e:
                self.logger.error(f"   TASSEL命令测试失败 | TASSEL command test failed: {e}")
                raise RuntimeError(f"TASSEL command test failed: {e}")

            self.logger.info(f"   开始执行TASSEL Kinship计算 | Starting TASSEL Kinship calculation...")

            try:
                # 使用shell=True和组合命令来更好地处理输出
                shell_cmd = ' '.join(cmd)
                self.logger.info(f"   执行命令 | Executing command: {shell_cmd}")
                self.logger.info(f"   工作目录 | Working directory: {self.config.output_dir}")

                result = subprocess.run(shell_cmd, capture_output=True, text=True,
                                        shell=True, timeout=1800,  # 30分钟超时
                                        cwd=str(self.config.output_dir))

                self.logger.info(f"   命令返回码 | Command return code: {result.returncode}")
                if result.stdout:
                    self.logger.info(f"   标准输出 | Stdout: {result.stdout[:500]}...")
                if result.stderr:
                    self.logger.warning(f"   标准错误 | Stderr: {result.stderr[:500]}...")
            except subprocess.TimeoutExpired:
                self.logger.error(f"   Kinship计算超时 | Kinship calculation timed out after 30 minutes")
                raise RuntimeError("Kinship计算超时 | Kinship calculation timed out")

            end_time = time.time()
            kinship_time = int(end_time - start_time)

            # 详细检查执行结果
            self.logger.info(f"   命令返回码 | Command return code: {result.returncode}")
            self.logger.info(f"   执行时间 | Execution time: {kinship_time}s")

            if result.returncode != 0:
                self.logger.error(f"   命令执行失败 | Command execution failed")
                self.logger.error(f"   标准错误 | Stderr: {result.stderr}")
                self.logger.error(f"   标准输出 | Stdout: {result.stdout[:500]}")  # 只显示前500字符
                raise RuntimeError(f"Kinship计算失败，返回码 | Kinship calculation failed, return code: {result.returncode}")

            # 检查输出文件
            self.logger.info(f"   检查输出文件 | Checking output file: {kinship_file}")
            self.logger.info(f"   文件存在 | File exists: {kinship_file.exists()}")

            if kinship_file.exists():
                file_size = kinship_file.stat().st_size
                self.logger.info(f"   文件大小 | File size: {file_size} bytes")

                if file_size > 0:
                    kinship_size = self._count_matrix_lines(kinship_file)
                    self.logger.info(f"✅ Kinship计算完成 | Kinship calculation completed in {kinship_time}s ({kinship_size}x{kinship_size})")
                    return kinship_file
                else:
                    self.logger.error(f"   Kinship文件为空 | Kinship file is empty")
                    raise RuntimeError("Kinship计算失败，输出文件为空 | Kinship calculation failed, output file is empty")
            else:
                self.logger.error(f"   Kinship文件不存在 | Kinship file does not exist")
                self.logger.error(f"   工作目录内容 | Working directory contents:")
                import os
                try:
                    for file in os.listdir(self.config.output_dir):
                        file_path = self.config.output_dir / file
                        if file_path.is_file():
                            stat = file_path.stat()
                            self.logger.error(f"     {file} - {stat.st_size} bytes")
                except Exception as e:
                    self.logger.error(f"   无法列出目录 | Cannot list directory: {e}")
                raise RuntimeError("Kinship计算失败，未生成输出文件 | Kinship calculation failed, no output file generated")

        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ Kinship计算失败 | Kinship calculation failed: {e}")
            self.logger.error(f"   返回码 | Return code: {e.returncode}")
            self.logger.error(f"   标准错误 | Stderr: {e.stderr}")
            self.logger.error(f"   标准输出 | Stdout: {e.stdout}")
            raise RuntimeError(f"Kinship计算失败 | Kinship calculation failed: {e.stderr}")

    def _count_matrix_lines(self, file_path: Path) -> int:
        """计算矩阵文件的行数 | Count matrix file lines"""
        try:
            with open(file_path, 'r') as f:
                return sum(1 for _ in f)
        except Exception:
            return 0

    def _format_pca_for_tassel(self, pca_file: Path) -> Path:
        """将VCF2PCACluster的PCA特征向量转换为TASSEL协变量格式 | Convert VCF2PCACluster PCA eigenvectors to TASSEL covariate format"""
        try:
            # 输出文件路径
            formatted_file = pca_file.parent / f"{pca_file.stem}.covariate.txt"

            with open(pca_file, 'r') as infile, open(formatted_file, 'w') as outfile:
                lines = infile.readlines()

                # 跳过表头
                if lines and lines[0].startswith('SampleName'):
                    lines = lines[1:]  # 跳过第一行表头

                # 处理每行数据，只保留样本名和PC列（跳过Group和Cluster列）
                pc_data = []
                for line in lines:
                    line = line.strip()
                    if line:
                        cols = line.split('\t')
                        if len(cols) >= 4:
                            # 格式：SampleName PC1 PC2 PC3...
                            sample_name = cols[0]
                            pc_values = cols[3:]  # 跳过Group和Cluster列
                            pc_data.append([sample_name] + pc_values)

                # 写入TASSEL表型文件格式的协变量
                if pc_data:
                    # 使用配置的PCA主成分数量
                    total_available_pc = len(pc_data[0]) - 1  # 可用的PC总数
                    requested_pc = min(self.config.pca_components, total_available_pc)

                    self.logger.info(f"   可用PC数量 | Available PC count: {total_available_pc}")
                    self.logger.info(f"   请求PC数量 | Requested PC count: {self.config.pca_components}")
                    self.logger.info(f"   实际使用PC数量 | Actually using PC count: {requested_pc}")

                    # 截取指定数量的PC
                    pc_data_limited = []
                    for row in pc_data:
                        pc_data_limited.append([row[0]] + row[1:1+requested_pc])

                    num_pc = requested_pc
                    header_cols = ['<Trait>'] + [f'<PC{i+1}>' for i in range(num_pc)]
                    outfile.write('\t'.join(header_cols) + '\n')

                    data_cols = ['<Data>'] + ['<Data>'] * num_pc
                    outfile.write('\t'.join(data_cols) + '\n')

                    # 写入数据行
                    for row in pc_data_limited:
                        outfile.write('\t'.join(row) + '\n')

            self.logger.info(f"✅ PCA协变量格式转换完成 | PCA covariate format conversion completed")
            self.logger.info(f"   输出文件 | Output file: {formatted_file}")
            self.logger.info(f"   最终PC数量 | Final PC count: {num_pc}")
            self.logger.info(f"   样本数量 | Sample count: {len(pc_data_limited)}")

            return formatted_file

        except Exception as e:
            self.logger.error(f"❌ PCA协变量格式转换失败 | PCA covariate format conversion failed: {e}")
            # 如果转换失败，返回原始文件
            return pca_file

    def run_glm(self, vcf_file: Path, pheno_file: Path) -> Optional[Dict]:
        """运行GLM模型 | Run GLM model"""
        if self.config.model not in ['GLM', 'BOTH']:
            return None

        self.logger.info("📊 运行GLM模型 | Running GLM (General Linear Model)...")
        start_time = time.time()

        glm_out = self.config.output_dir / f"{self.config.output_prefix}_GLM"
        glm_log = self.config.output_dir / f"{self.config.output_prefix}.glm.log"

        try:
            # 构建TASSEL命令 | Build TASSEL command
            # 所有文件都在工作目录中，直接使用文件名
            cmd = [
                'perl', self.tassel_cmd,
                '-Xms' + self.config.memory_min,
                '-Xmx' + self.config.memory_max,
                '-fork1', '-vcf', str(self.config.vcf_file.absolute()),
                '-fork2', '-t', str(pheno_file.absolute())
            ]

            if self.config.q_matrix and self.config.q_matrix.exists():
                self.logger.info("   -> 使用Q矩阵进行群体结构校正 | Using Q-matrix for population structure correction")
                # Q矩阵使用绝对路径（通常在外部）
                cmd.extend([
                    '-fork3', '-q', str(self.config.q_matrix.absolute()),
                    '-combine4', '-input1', '-input2', '-input3', '-intersect'
                ])
            else:
                self.logger.warning("   -> 未使用Q矩阵（可能增加假阳性）| Running without Q-matrix (may increase false positives)")
                cmd.extend(['-combine4', '-input1', '-input2', '-intersect'])

            # GLM输出使用文件名
            cmd.extend(['-FixedEffectLMPlugin', '-endPlugin', '-export', glm_out.name])

            glm_cmd = ' '.join(cmd)
            self.logger.info(f"   执行GLM命令 | Executing GLM command: {glm_cmd}")
            self.logger.info(f"   GLM工作目录 | GLM working directory: {self.config.output_dir}")

            result = subprocess.run(cmd, capture_output=True, text=True, check=True,
                                    cwd=str(self.config.output_dir))

            self.logger.info(f"   GLM命令返回码 | GLM command return code: {result.returncode}")
            if result.stdout:
                self.logger.info(f"   GLM标准输出 | GLM stdout: {result.stdout[:500]}...")
            if result.stderr:
                self.logger.warning(f"   GLM标准错误 | GLM stderr: {result.stderr[:500]}...")

            end_time = time.time()
            glm_time = int(end_time - start_time)

            # 提取结果 | Extract results
            # TASSEL GLM输出文件可能是目录形式或直接文件
            glm_result_file = None

            # 先检查目录形式: glm_out/1.txt
            dir_result = glm_out / "1.txt"
            if dir_result.exists():
                glm_result_file = dir_result
            else:
                # 检查直接文件形式: glm_out + "1.txt" 或 glm_out + "1.txt"
                direct_result1 = self.config.output_dir / f"{glm_out.name}1.txt"
                direct_result2 = self.config.output_dir / f"{glm_out.name}2.txt"
                if direct_result1.exists():
                    glm_result_file = direct_result1
                elif direct_result2.exists():
                    glm_result_file = direct_result2

            if glm_result_file and glm_result_file.exists():
                output_file = self.config.output_dir / f"{self.config.output_prefix}.glm.manht_input"
                self._extract_glm_results(glm_result_file, output_file)

                # 统计信息 | Statistics
                with open(output_file, 'r') as f:
                    lines = f.readlines()
                total_snps = len(lines) - 1 if lines else 0
                significant_snps = len([line for line in lines[1:] if line.strip() and
                                       float(line.strip().split('\t')[3]) < 1e-5])

                self.logger.info(f"✅ GLM完成 | GLM completed in {glm_time}s")
                self.logger.info(f"   • 总SNPs | Total SNPs: {total_snps}")
                self.logger.info(f"   • 显著SNPs (P<1e-5) | Significant (P<1e-5): {significant_snps}")
                self.logger.info(f"   • 输出文件 | Output file: {output_file}")

                return {
                    'runtime': glm_time,
                    'snps_tested': total_snps,
                    'significant': significant_snps,
                    'output_file': str(output_file)
                }
            else:
                self.logger.error(f"❌ GLM输出文件未找到 | GLM output file not found")
                self.logger.error(f"   期望文件 | Expected files: {glm_out}/1.txt or {glm_out.name}1.txt")
                # 列出实际存在的文件
                import os
                existing_files = [f for f in os.listdir(self.config.output_dir) if 'GLM' in f]
                self.logger.error(f"   实际GLM文件 | Actual GLM files: {existing_files}")
                return None

        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ GLM运行失败 | GLM execution failed: {e.stderr}")
            return None

    def _extract_glm_results(self, input_file: Path, output_file: Path):
        """提取GLM结果 | Extract GLM results"""
        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                outfile.write("Chr\tPos\tSNP\tP-value\n")  # 写入表头
                for line in infile:
                    if line.strip():
                        cols = line.strip().split()
                        if len(cols) >= 6 and cols[5].replace('.', '').replace('e-', '').replace('e+', '').isdigit():
                            # 提取Chr, Pos, SNP, P-value (假设格式为: Marker, Chr, Pos, Allele1, Allele2, P-value, ...)
                            outfile.write(f"{cols[1]}\t{cols[2]}\t{cols[0]}\t{cols[5]}\n")
        except Exception as e:
            self.logger.warning(f"⚠️ GLM结果提取失败 | GLM result extraction failed: {e}")

    def run_mlm(self, vcf_file: Path, pheno_file: Path, kinship_file: Path) -> Optional[Dict]:
        """运行MLM模型 | Run MLM model"""
        if self.config.model not in ['MLM', 'BOTH']:
            return None

        self.logger.info("📊 运行MLM模型 | Running MLM (Mixed Linear Model)...")
        start_time = time.time()

        mlm_out = self.config.output_dir / f"{self.config.output_prefix}_MLM"
        mlm_log = self.config.output_dir / f"{self.config.output_prefix}.mlm.log"

        try:
            # 处理PCA协变量 | Handle PCA covariate
            pca_covariate_file = None
            if self.config.q_matrix and self.config.q_matrix.exists():
                # 如果用户提供了q_matrix参数，实际上这是PCA协变量文件
                pca_covariate_file = self.config.q_matrix
                self.logger.info("   -> 使用提供的PCA协变量 + Kinship (PCA+K模型) | Using provided PCA covariate + Kinship (PCA+K model)")
            elif not self.config.skip_sort:
                self.logger.info("   -> 未提供PCA协变量，自动计算PCA | No PCA covariate provided, auto-computing PCA...")
                try:
                    # 使用原始VCF文件计算PCA，而不是排序后的文件
                    original_vcf = self.config.vcf_file
                    self.logger.info(f"   -> 使用原始VCF计算PCA | Using original VCF for PCA: {original_vcf}")
                    pca_covariate_file = self.calculate_pca(original_vcf)
                    self.logger.info("   -> 使用自动计算的PCA协变量 + Kinship (PCA+K模型) | Using auto-computed PCA covariate + Kinship (PCA+K model)")
                except Exception as e:
                    self.logger.warning(f"   -> PCA计算失败，使用仅Kinship (K模型) | PCA calculation failed, using Kinship only (K model): {e}")
                    pca_covariate_file = None
            else:
                self.logger.info("   -> 跳过PCA计算（批量处理模式）| Skipping PCA calculation (batch processing mode)")
                pca_covariate_file = None

            # 构建TASSEL命令 | Build TASSEL command
            # 由于在子目录中运行，需要使用绝对路径
            cmd = [
                'perl', self.tassel_cmd,
                '-Xms' + self.config.memory_min,
                '-Xmx' + self.config.memory_max,
                '-fork1', '-vcf', str(self.config.vcf_file.absolute()),
                '-fork2', '-t', str(pheno_file.absolute())
            ]

            # MLM输出使用文件名
            if pca_covariate_file and pca_covariate_file.exists():
                cmd.extend([
                    '-fork3', '-r', str(pca_covariate_file.absolute()),
                    '-fork4', '-k', str(kinship_file.absolute()),
                    '-combine5', '-input1', '-input2', '-input3', '-intersect',
                    '-combine6', '-input5', '-input4', '-mlm',
                    '-mlmVarCompEst', 'P3D',
                    '-mlmCompressionLevel', 'None',
                    '-export', mlm_out.name
                ])
            else:
                self.logger.info("   -> 仅使用Kinship (K模型) | Using Kinship only (K model)")
                cmd.extend([
                    '-fork4', '-k', str(kinship_file.absolute()),
                    '-combine5', '-input1', '-input2', '-intersect',
                    '-combine6', '-input5', '-input4', '-mlm',
                    '-mlmVarCompEst', 'P3D',
                    '-mlmCompressionLevel', 'None',
                    '-export', mlm_out.name
                ])

            mlm_cmd = ' '.join(cmd)
            self.logger.info(f"   执行MLM命令 | Executing MLM command: {mlm_cmd}")
            self.logger.info(f"   MLM工作目录 | MLM working directory: {self.config.output_dir}")

            result = subprocess.run(cmd, capture_output=True, text=True, check=True,
                                    cwd=str(self.config.output_dir))

            self.logger.info(f"   MLM命令返回码 | MLM command return code: {result.returncode}")
            if result.stdout:
                self.logger.info(f"   MLM标准输出 | MLM stdout: {result.stdout[:500]}...")
            if result.stderr:
                self.logger.warning(f"   MLM标准错误 | MLM stderr: {result.stderr[:500]}...")

            end_time = time.time()
            mlm_time = int(end_time - start_time)

            # TASSEL已完成，直接查找输出文件 | TASSEL completed, directly find output files
            self.logger.info(f"   🔍 搜索MLM输出文件... | Searching for MLM output files...")
            import glob

            # 根据你的实际输出，TASSEL生成的文件格式为：{output_prefix}12.txt 和 {output_prefix}13.txt
            patterns = [
                f"{mlm_out}12.txt",    # MLM statistics file (主要统计文件)
                f"{mlm_out}13.txt",    # MLM effects file (效应文件)
                f"{mlm_out}2.txt",
                f"{mlm_out}2.stats.txt",
                f"{mlm_out}*.txt",
                f"{mlm_out}*.stats*",
                f"{mlm_out}_2.txt",
                f"{mlm_out}_2.stats.txt",
                f"{mlm_out}_MLM.txt"
            ]

            export_patterns = [
                f"{mlm_out}*",
                f"*MLM*"
            ]

            # 列出当前目录的文件（调试用）
            try:
                current_files = list(self.config.output_dir.glob("*"))
                txt_files = [f for f in current_files if f.suffix == '.txt' and 'MLM' in f.name]
                self.logger.info(f"   📁 当前目录 .txt 文件 (包含MLM) | Current directory .txt files (containing MLM):")
                for f in txt_files[:10]:  # 只显示前10个
                    self.logger.info(f"     • {f.name} ({f.stat().st_size:,} bytes)")
                if len(txt_files) > 10:
                    self.logger.info(f"     ... 以及另外{len(txt_files)-10}个文件")
            except Exception as e:
                self.logger.warning(f"   ⚠️ 无法列出当前目录文件 | Cannot list current directory files: {e}")

            mlm_result_file = None

            # 直接查找输出文件，不需要循环 | Directly find output files, no loop needed
            self.logger.info(f"   🔍 尝试以下文件模式 | Trying file patterns:")
            for pattern in patterns:
                matches = glob.glob(str(pattern))
                if matches:
                    self.logger.info(f"   ✅ 模式匹配成功 | Pattern matched: {pattern}")
                    # 选择最新的文件
                    matches_with_time = [(Path(f), Path(f).stat().st_mtime) for f in matches]
                    matches_with_time.sort(key=lambda x: x[1], reverse=True)
                    candidate_file = matches_with_time[0][0]

                    if candidate_file.stat().st_size > 0:
                        mlm_result_file = candidate_file
                        self.logger.info(f"   🎯 找到MLM输出文件 | Found MLM output file: {mlm_result_file}")
                        self.logger.info(f"   📄 文件大小 | File size: {mlm_result_file.stat().st_size:,} bytes")
                        break

            # 如果主要模式没找到，尝试导出模式 | If main patterns not found, try export patterns
            if not mlm_result_file:
                self.logger.info(f"   🔍 尝试导出文件模式 | Trying export file patterns:")
                for pattern in export_patterns:
                    matches = glob.glob(str(pattern))
                    if matches:
                        self.logger.info(f"   ✅ 导出模式匹配成功 | Export pattern matched: {pattern}")
                        # 选择最大的文件
                        valid_matches = []
                        for f in matches:
                            try:
                                f_path = Path(f)
                                if f_path.stat().st_size > 0:
                                    valid_matches.append((f_path, f_path.stat().st_size))
                            except Exception as e:
                                self.logger.debug(f"   ⚠️ 无法访问文件 {f}: {e}")

                        if valid_matches:
                            valid_matches.sort(key=lambda x: x[1], reverse=True)
                            mlm_result_file = valid_matches[0][0]
                            self.logger.info(f"   🎯 找到MLM导出文件 | Found MLM export file: {mlm_result_file}")
                            self.logger.info(f"   📄 文件大小 | File size: {mlm_result_file.stat().st_size:,} bytes")
                            break

            if mlm_result_file and mlm_result_file.exists():
                output_file = self.config.output_dir / f"{self.config.output_prefix}.mlm.manht_input"
                self._extract_mlm_results(mlm_result_file, output_file)

                # 统计信息 | Statistics
                with open(output_file, 'r') as f:
                    lines = f.readlines()
                total_snps = len(lines) - 1 if lines else 0
                significant_snps = len([line for line in lines[1:] if line.strip() and
                                       float(line.strip().split('\t')[3]) < 1e-5])

                self.logger.info(f"✅ MLM完成 | MLM completed in {mlm_time}s")
                self.logger.info(f"   • 总SNPs | Total SNPs: {total_snps}")
                self.logger.info(f"   • 显著SNPs (P<1e-5) | Significant (P<1e-5): {significant_snps}")
                self.logger.info(f"   • 输出文件 | Output file: {output_file}")

                return {
                    'runtime': mlm_time,
                    'snps_tested': total_snps,
                    'significant': significant_snps,
                    'output_file': str(output_file)
                }
            else:
                self.logger.error(f"❌ MLM输出文件未找到 | MLM output file not found")
                self.logger.error(f"   搜索模式 | Search patterns: {patterns + export_patterns}")
                self.logger.error(f"   当前工作目录 | Current working directory: {self.config.output_dir}")
                self.logger.error(f"   当前目录文件 | Current directory files:")
                try:
                    cwd_files = list(self.config.output_dir.glob("*"))
                    for f in cwd_files[:30]:  # 显示前30个文件
                        file_info = f.name
                        if f.is_file():
                            file_info += f" ({f.stat().st_size} bytes)"
                        self.logger.error(f"     {file_info}")
                except Exception as e:
                    self.logger.error(f"     无法列出文件 | Cannot list files: {e}")
                return None

        except subprocess.CalledProcessError as e:
            self.logger.error(f"❌ MLM运行失败 | MLM execution failed: {e.stderr}")
            return None

    def _extract_mlm_results(self, input_file: Path, output_file: Path):
        """提取MLM结果 | Extract MLM results"""
        try:
            with open(input_file, 'r') as infile, open(output_file, 'w') as outfile:
                outfile.write("Chr\tPos\tSNP\tP-value\n")  # 写入表头
                for line in infile:
                    if line.strip():
                        cols = line.strip().split()
                        if len(cols) >= 7 and cols[6].replace('.', '').replace('e-', '').replace('e+', '').isdigit():
                            # 提取Chr, Pos, SNP, P-value (假设格式为: Marker, Chr, Pos, ..., P-value)
                            outfile.write(f"{cols[1]}\t{cols[2]}\t{cols[0]}\t{cols[6]}\n")
        except Exception as e:
            self.logger.warning(f"⚠️ MLM结果提取失败 | MLM result extraction failed: {e}")

    def run_analysis(self):
        """运行完整的TASSEL GWAS分析 | Run complete TASSEL GWAS analysis"""
        try:
            start_time = time.time()
            self.logger.info("🚀 开始TASSEL GWAS分析流程 | Starting TASSEL GWAS analysis pipeline")

            # 检查依赖 | Check dependencies
            if not self.check_dependencies():
                return False

            # 验证输入文件 | Validate input files
            if not self.validate_input_files():
                return False

            # 获取文件统计信息 | Get file statistics
            vcf_variants, vcf_samples = self._get_vcf_stats(self.config.vcf_file)
            pheno_samples, pheno_traits = self._get_pheno_stats(self.config.pheno_file)

            self.logger.info(f"📊 输入VCF | Input VCF: {self.config.vcf_file}")
            self.logger.info(f"   • 变异位点 | Variants: {vcf_variants}")
            self.logger.info(f"   • 样本数 | Samples: {vcf_samples}")
            self.logger.info(f"📊 表型文件 | Phenotype file: {self.config.pheno_file}")
            self.logger.info(f"   • 样本数 | Samples: {pheno_samples}")
            self.logger.info(f"   • 表型数量 | Traits: {pheno_traits}")

            # 准备表型文件 | Prepare phenotype file
            pheno_ready = self.prepare_phenotype_file()

            # VCF排序和过滤（如果未跳过）| VCF sorting and filtering (if not skipped)
            if self.config.skip_sort:
                # 使用原始VCF文件（假设已经预先处理过）| Use original VCF file (assume pre-processed)
                vcf_sorted = self.config.vcf_file
                self.logger.info("⏭️ 跳过VCF排序，使用原始文件 | Skipping VCF sorting, using original file")
            else:
                vcf_sorted = self.sort_and_filter_vcf()

            # 计算Kinship矩阵（如果需要且未预计算）| Calculate Kinship matrix (if needed and not pre-computed)
            kinship_file = None
            if self.config.kinship and self.config.kinship.exists():
                kinship_file = self.config.kinship
                self.logger.info(f"✅ 使用预计算的K矩阵 | Using pre-computed K matrix: {kinship_file}")
            elif self.config.model in ['MLM', 'BOTH']:
                if self.config.skip_sort:
                    # 跳过矩阵计算时，假设矩阵应该已经预计算 | When skipping matrix calculation, assume matrices are pre-computed
                    self.logger.warning("⚠️ 未提供K矩阵但跳过了排序，请确保已预计算K矩阵 | No K matrix provided but sorting skipped, please ensure K matrix is pre-computed")
                else:
                    kinship_file = self.calculate_kinship(vcf_sorted)

            # 运行分析模型 | Run analysis models
            results = {}
            output_files = []

            if self.config.model in ['GLM', 'BOTH']:
                glm_results = self.run_glm(vcf_sorted, pheno_ready)
                if glm_results:
                    results['glm_results'] = glm_results
                    output_files.append(glm_results['output_file'])

            if self.config.model in ['MLM', 'BOTH'] and kinship_file:
                mlm_results = self.run_mlm(vcf_sorted, pheno_ready, kinship_file)
                if mlm_results:
                    results['mlm_results'] = mlm_results
                    output_files.append(mlm_results['output_file'])

            # 生成统计报告 | Generate statistics report
            end_time = time.time()
            total_time = int(end_time - start_time)

            self.stats.update({
                'run_date': time.strftime('%Y-%m-%d %H:%M:%S'),
                'work_dir': str(Path.cwd()),
                'vcf_file': str(self.config.vcf_file),
                'vcf_variants': vcf_variants,
                'vcf_samples': vcf_samples,
                'pheno_file': str(self.config.pheno_file),
                'pheno_samples': pheno_samples,
                'pheno_traits': pheno_traits,
                'model': self.config.model,
                'memory': self.config.memory_max,
                'maf_filter': self.config.maf_filter or 'None',
                'miss_filter': self.config.miss_filter or 'None',
                'total_runtime': total_time,
                'log_file': str(self.logger_manager.log_file),
                'output_files': output_files
            })
            self.stats.update(results)

            stats_file = self.config.output_dir / f"{self.config.output_prefix}.stats.txt"
            if create_statistics_report(self.stats, stats_file):
                self.logger.info(f"✅ 统计报告已保存 | Statistics report saved: {stats_file}")

            # 清理临时文件 | Clean up temporary files
            if not self.config.keep_temp:
                self._cleanup_temp_files()

            self.logger.info(f"🎉 分析完成！总运行时间 | Analysis completed! Total runtime: {total_time}s")
            return True

        except KeyboardInterrupt:
            self.logger.info("\n🛑 用户中断操作 | User interrupted operation")
            return False
        except Exception as e:
            self.logger.error(f"❌ 分析过程中发生错误 | Error occurred during analysis: {str(e)}")
            return False

    def _cleanup_temp_files(self):
        """清理临时文件 | Clean up temporary files"""
        if self.temp_files:
            self.logger.info("🧹 清理临时文件 | Cleaning up temporary files...")
            for temp_file in self.temp_files:
                if temp_file.exists():
                    try:
                        temp_file.unlink()
                        self.logger.debug(f"已删除 | Removed: {temp_file}")
                    except Exception as e:
                        self.logger.warning(f"无法删除临时文件 | Cannot remove temp file: {temp_file} ({e})")