"""
TASSEL GWAS工具函数模块|TASSEL GWAS Utility Functions Module
"""

import logging
import subprocess
import sys
import os
from pathlib import Path


class TASSELLogger:
    """TASSEL GWAS日志管理器|TASSEL GWAS Logger Manager"""

    def __init__(self, output_dir, log_name: str = "tassel_gwas.log"):
        # 确保output_dir是Path对象|Ensure output_dir is Path object
        self.output_dir = Path(output_dir).expanduser() if not isinstance(output_dir, Path) else output_dir.expanduser()
        self.log_file = self.output_dir / log_name
        self.setup_logging()

    def setup_logging(self):
        """设置日志|Setup logging"""
        # 确保输出目录存在|Ensure output directory exists
        self.output_dir.mkdir(parents=True, exist_ok=True)

        if self.log_file.exists():
            self.log_file.unlink()

        logging.basicConfig(
            level=logging.INFO,
            format='%(asctime)s.%(msecs)03d - %(levelname)s - %(message)s',
            datefmt='%Y-%m-%d %H:%M:%S',
            handlers=[
                logging.FileHandler(self.log_file),
                logging.StreamHandler(sys.stdout)
            ]
        )
        self.logger = logging.getLogger(__name__)

    def get_logger(self):
        """获取日志器|Get logger"""
        return self.logger


def check_dependencies(logger=None):
    """检查TASSEL和相关依赖是否已安装|Check if TASSEL and dependencies are installed"""

    if logger is None:
        # 创建临时logger|Create temporary logger
        logging.basicConfig(level=logging.INFO)
        logger = logging.getLogger(__name__)

    # 检查 TASSEL|Check TASSEL
    tassel_paths = [
        "~/software/tassel-5-standalone/run_pipeline.pl",
        os.path.expanduser("~/software/tassel-5-standalone/run_pipeline.pl"),
        "run_pipeline.pl"
    ]

    tassel_found = False
    tassel_path = None

    for path in tassel_paths:
        if Path(path).exists():
            tassel_path = Path(path)
            tassel_found = True
            break

    if not tassel_found:
        # 在PATH中查找|Search in PATH
        try:
            which_cmd = ['which', 'run_pipeline.pl']
            logger.info(f"命令|Command: {' '.join(which_cmd)}")
            result = subprocess.run(which_cmd, capture_output=True, text=True, check=True)
            tassel_path = Path(result.stdout.strip())
            tassel_found = True
        except subprocess.CalledProcessError:
            pass

    if not tassel_found:
        logger.error("错误: 未找到TASSEL的run_pipeline.pl脚本|Error: TASSEL run_pipeline.pl script not found")
        logger.info("   请确保TASSEL已安装并且在系统PATH中|Please ensure TASSEL is installed and in your system PATH")
        logger.info("   或者设置正确的TASSEL_HOME路径|Or set the correct TASSEL_HOME path")
        logger.info("   安装方法: https://www.maizegenetics.net/tassel/|Installation: https://www.maizegenetics.net/tassel/")
        return False, None

    # 检查Perl|Check Perl
    try:
        perl_cmd = ['perl', '--version']
        logger.info(f"命令|Command: {' '.join(perl_cmd)}")
        subprocess.run(perl_cmd, check=True, capture_output=True)
        logger.info("Perl 已找到|Perl found")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("错误: 未找到Perl|Error: Perl not found")
        logger.info("   TASSEL需要Perl运行环境|TASSEL requires Perl runtime environment")
        return False, None

    # 检查Java|Check Java
    try:
        java_cmd = ['java', '-version']
        logger.info(f"命令|Command: {' '.join(java_cmd)}")
        result = subprocess.run(java_cmd, capture_output=True, text=True)
        java_output = result.stderr or result.stdout
        java_version = java_output.split('"')[1] if '"' in java_output else "Unknown"
        logger.info(f"Java 已找到|Java found: {java_version}")
    except (subprocess.CalledProcessError, FileNotFoundError):
        logger.error("错误: 未找到Java|Error: Java not found")
        logger.info("   TASSEL需要Java运行环境|TASSEL requires Java runtime environment")
        return False, None

    # 测试TASSEL是否正常工作|Test if TASSEL works properly
    try:
        test_cmd = [str(tassel_path), '-h']
        logger.info(f"命令|Command: {' '.join(test_cmd)}")
        result = subprocess.run(test_cmd, capture_output=True, text=True, timeout=10)
        if result.returncode == 0 or "Usage" in result.stdout or "help" in result.stdout:
            logger.info(f"TASSEL 可用: {tassel_path}|TASSEL available: {tassel_path}")
            return True, tassel_path
        else:
            logger.warning(f"TASSEL可能有问题|TASSEL may have issues: {result.stderr}")
            return True, tassel_path  # 仍然返回True，让用户尝试
    except subprocess.TimeoutExpired:
        logger.warning("TASSEL响应超时，但可能仍然可用|TASSEL response timeout, but may still be usable")
        return True, tassel_path
    except Exception as e:
        logger.error(f"TASSEL测试失败|TASSEL test failed: {e}")
        return False, None


def validate_vcf_file(vcf_file: Path, logger):
    """验证VCF文件是否有效|Validate VCF file"""
    if not vcf_file.exists():
        logger.error(f"VCF文件不存在|VCF file does not exist: {vcf_file}")
        return False

    if vcf_file.suffix.lower() not in ['.vcf', '.gz']:
        logger.error(f"文件格式错误，应为 .vcf 或 .vcf.gz|File format error, should be .vcf or .vcf.gz: {vcf_file}")
        return False

    # 检查VCF文件头|Check VCF file header
    try:
        # 检查是否有VCF格式头部
        if vcf_file.suffix.lower() == '.gz':
            import gzip
            header_cmd = 'zcat -f "{}" 2>/dev/null|head -10|grep "##fileformat=VCF"'.format(str(vcf_file))
        else:
            header_cmd = 'cat "{}" 2>/dev/null|head -10|grep "##fileformat=VCF"'.format(str(vcf_file))

        logger.info(f"命令|Command: {header_cmd}")
        header_result = subprocess.run(header_cmd, shell=True, capture_output=True, text=True, timeout=30)

        if header_result.returncode != 0 or not header_result.stdout.strip():
            logger.error(f"VCF文件格式错误，缺少VCF头部|VCF file format error, missing VCF header: {vcf_file}")
            return False

        # 使用bash逻辑获取统计信息|Use bash logic to get statistics
        # VCF_VARIANTS=$(zcat -f "$VCF_FILE" 2>/dev/null|grep -v "^#"|wc -l)
        if vcf_file.suffix.lower() == '.gz':
            variant_cmd = 'zcat -f "{}" 2>/dev/null|grep -v "^#"|wc -l'.format(str(vcf_file))
            sample_cmd = 'zcat -f "{}" 2>/dev/null|grep "^#CHROM"|awk \'{{print NF-9}}\''.format(str(vcf_file))
        else:
            variant_cmd = 'cat "{}" 2>/dev/null|grep -v "^#"|wc -l'.format(str(vcf_file))
            sample_cmd = 'cat "{}" 2>/dev/null|grep "^#CHROM"|awk \'{{print NF-9}}\''.format(str(vcf_file))

        # 增加超时时间以支持大规模VCF文件（600秒 = 10分钟）
        logger.info(f"命令|Command: {variant_cmd}")
        variant_result = subprocess.run(variant_cmd, shell=True, capture_output=True, text=True, timeout=600)
        logger.info(f"命令|Command: {sample_cmd}")
        sample_result = subprocess.run(sample_cmd, shell=True, capture_output=True, text=True, timeout=600)

        if variant_result.returncode == 0:
            variant_count = int(variant_result.stdout.strip())
        else:
            variant_count = 0

        if sample_result.returncode == 0 and sample_result.stdout.strip():
            sample_count = int(sample_result.stdout.strip())
        else:
            sample_count = 0

        logger.info(f"VCF文件验证通过|VCF file validation passed: {vcf_file}")
        logger.info(f"   • 变异位点|Variants: {variant_count}")
        logger.info(f"   • 样本数|Samples: {sample_count}")

        return True

    except subprocess.TimeoutExpired:
        logger.error(f"VCF文件读取超时|VCF file read timeout: {vcf_file}")
        return False
    except Exception as e:
        logger.error(f"VCF文件验证失败|VCF file validation failed: {e}")
        logger.error(f"   文件路径|File path: {vcf_file}")
        logger.error(f"   错误类型|Error type: {type(e).__name__}")
        return False


def validate_phenotype_file(pheno_file: Path, logger):
    """验证表型文件是否有效|Validate phenotype file"""
    if not pheno_file.exists():
        logger.error(f"表型文件不存在|Phenotype file does not exist: {pheno_file}")
        return False

    try:
        with open(pheno_file, 'r') as f:
            lines = f.readlines()

        if len(lines) < 2:
            logger.error(f"表型文件至少需要包含表头和一行数据|Phenotype file needs header and at least one data row: {pheno_file}")
            return False

        # 检查格式|Check format - 先尝试制表符分割，如果只有一列则尝试空格分割
        header_split = lines[0].strip().split('\t')
        if len(header_split) <= 1:
            # 尝试空格分割（支持多个连续空格）
            header_split = lines[0].strip().split()

        header_cols = len(header_split)
        if header_cols < 2:
            logger.error(f"表型文件至少需要2列（样本ID + 表型值）| Phenotype file needs at least 2 columns (Sample ID + Trait value): {pheno_file}")
            return False

        data_rows = 0
        for i, line in enumerate(lines[1:], 1):
            if line.strip() and not line.startswith('#'):
                # 使用相同的分隔符逻辑
                line_split = line.strip().split('\t')
                if len(line_split) <= 1:
                    line_split = line.strip().split()
                cols = len(line_split)
                if cols != header_cols:
                    logger.error(f"第{i+1}行列数({cols})与表头列数({header_cols})不匹配|Line {i+1} columns({cols}) don't match header columns({header_cols}): {pheno_file}")
                    return False
                data_rows += 1

        logger.info(f"表型文件验证通过|Phenotype file validation passed: {pheno_file}")
        logger.info(f"   • 数据行数|Data rows: {data_rows}")
        logger.info(f"   • 表型数量|Traits: {header_cols - 1}")

        return True

    except Exception as e:
        logger.error(f"表型文件验证失败|Phenotype file validation failed: {e}")
        return False


def prepare_phenotype_file(input_pheno: Path, output_pheno: Path, logger):
    """准备符合TASSEL格式的表型文件|Prepare TASSEL-compatible phenotype file"""
    try:
        logger.info(f"准备表型文件|Preparing phenotype file: {input_pheno}")

        with open(input_pheno, 'r') as infile, open(output_pheno, 'w') as outfile:
            # 处理表头|Process header
            header = infile.readline().strip()
            if not header:
                raise ValueError("表型文件为空|Phenotype file is empty")

            # 先尝试制表符分割，如果只有一列则尝试空格分割
            header_cols = header.split('\t')
            if len(header_cols) <= 1:
                # 尝试空格分割（支持多个连续空格）
                header_cols = header.split()

            # 第一列替换为<Trait>|Replace first column with <Trait>
            header_cols[0] = "<Trait>"
            outfile.write('\t'.join(header_cols) + '\n')

            # 添加数据类型行|Add data type line
            data_types = ["<Data>"] * len(header_cols)
            outfile.write('\t'.join(data_types) + '\n')

            # 处理数据行|Process data rows
            data_count = 0
            for line in infile:
                line = line.strip()
                if line and not line.startswith('#'):
                    outfile.write(line + '\n')
                    data_count += 1

        logger.info(f"表型文件准备完成|Phenotype file prepared: {output_pheno}")
        logger.info(f"   • 数据行数|Data rows: {data_count}")
        logger.info(f"   • 表型数量|Traits: {len(header_cols) - 1}")

        return True

    except Exception as e:
        logger.error(f"表型文件准备失败|Phenotype file preparation failed: {e}")
        return False


def create_statistics_report(stats: dict, output_file: Path):
    """创建统计报告|Create statistics report"""
    try:
        with open(output_file, 'w') as f:
            f.write("=" * 80 + "\n")
            f.write("                    TASSEL GWAS Pipeline Summary Report\n")
            f.write("=" * 80 + "\n")
            f.write(f"Pipeline Version:   v1.0\n")
            f.write(f"Run Date:           {stats.get('run_date', 'Unknown')}\n")
            f.write(f"Working Directory:  {stats.get('work_dir', 'Unknown')}\n\n")

            f.write("INPUT FILES\n")
            f.write("-" * 50 + "\n")
            f.write(f"VCF File:           {stats.get('vcf_file', 'Unknown')}\n")
            f.write(f"  • Total Variants: {stats.get('vcf_variants', 'N/A')}\n")
            f.write(f"  • Total Samples:  {stats.get('vcf_samples', 'N/A')}\n\n")

            f.write(f"Phenotype File:     {stats.get('pheno_file', 'Unknown')}\n")
            f.write(f"  • Total Samples:  {stats.get('pheno_samples', 'N/A')}\n")
            f.write(f"  • Total Traits:   {stats.get('pheno_traits', 'N/A')}\n\n")

            f.write("ANALYSIS SETTINGS\n")
            f.write("-" * 50 + "\n")
            f.write(f"Model:              {stats.get('model', 'Unknown')}\n")
            f.write(f"Memory:             {stats.get('memory', 'Unknown')}\n")
            f.write(f"MAF Filter:         {stats.get('maf_filter', 'None')}\n")
            f.write(f"Missing Filter:     {stats.get('miss_filter', 'None')}\n\n")

            f.write("RESULTS\n")
            f.write("-" * 50 + "\n")

            if 'glm_results' in stats:
                glm = stats['glm_results']
                f.write(f"GLM Analysis:\n")
                f.write(f"  • Runtime:        {glm.get('runtime', 'N/A')}s\n")
                f.write(f"  • SNPs Tested:    {glm.get('snps_tested', 'N/A')}\n")
                f.write(f"  • Significant:    {glm.get('significant', 'N/A')} (P < 1e-5)\n")
                f.write(f"  • Output File:    {glm.get('output_file', 'N/A')}\n\n")

            if 'mlm_results' in stats:
                mlm = stats['mlm_results']
                f.write(f"MLM Analysis:\n")
                f.write(f"  • Runtime:        {mlm.get('runtime', 'N/A')}s\n")
                f.write(f"  • SNPs Tested:    {mlm.get('snps_tested', 'N/A')}\n")
                f.write(f"  • Significant:    {mlm.get('significant', 'N/A')} (P < 1e-5)\n")
                f.write(f"  • Output File:    {mlm.get('output_file', 'N/A')}\n\n")

            f.write("PERFORMANCE\n")
            f.write("-" * 50 + "\n")
            total_time = stats.get('total_runtime', 0)
            hours = total_time // 3600
            minutes = (total_time % 3600) // 60
            seconds = total_time % 60
            f.write(f"Total Runtime:      {total_time}s ({hours:02d}:{minutes:02d}:{seconds:02d})\n\n")

            f.write("OUTPUT FILES\n")
            f.write("-" * 50 + "\n")
            f.write(f"Log File:           {stats.get('log_file', 'N/A')}\n")
            f.write(f"Statistics:         {output_file}\n")
            for output_file in stats.get('output_files', []):
                f.write(f"Results:           {output_file}\n")

            f.write("\n" + "=" * 80 + "\n")

        return True

    except Exception as e:
        return False


def extract_vcf_samples(vcf_file: Path, logger) -> set:
    """
    从VCF文件中提取样本ID（使用bcftools，快速）

    Args:
        vcf_file: VCF文件路径
        logger: 日志对象

    Returns:
        set: 样本ID集合
    """
    try:
        # 使用bcftools查询样本ID
        cmd = f"bcftools query -l {vcf_file}"
        logger.info(f"命令|Command: {cmd}")
        result = subprocess.run(
            cmd,
            shell=True,
            capture_output=True,
            text=True,
            timeout=120  # 2分钟超时
        )

        if result.returncode != 0:
            logger.error(f"Failed to extract samples from VCF: {result.stderr}")
            return set()

        samples = set()
        for line in result.stdout.strip().split('\n'):
            if line:
                samples.add(line)

        return samples

    except subprocess.TimeoutExpired:
        logger.error(f"VCF样本提取超时|VCF sample extraction timeout: {vcf_file}")
        return set()
    except Exception as e:
        if logger:
            logger.error(f"Failed to extract VCF samples: {str(e)}")
        return set()


def find_common_samples_and_filter(
    vcf_file: Path,
    pheno_file: Path,
    output_prefix: str,
    output_dir: Path,
    logger
) -> tuple:
    """
    找到VCF和表型的样本交集，并使用bcftools从VCF筛选样本

    Args:
        vcf_file: 原始VCF文件路径
        pheno_file: 表型文件路径
        output_prefix: 输出文件前缀
        output_dir: 输出目录
        logger: 日志对象

    Returns:
        tuple: (是否成功, 过滤后的VCF文件路径, 过滤后的表型文件路径, 样本数量)
    """
    try:
        logger.info(" 开始样本交集处理|Starting common sample filtering...")

        # 1. 从VCF文件读取样本ID（使用bcftools，快速）
        vcf_samples = extract_vcf_samples(vcf_file, logger)

        if len(vcf_samples) == 0:
            logger.error("VCF文件中未找到样本|No samples found in VCF file!")
            return False, None, None, 0

        logger.info(f"   VCF样本数|VCF samples: {len(vcf_samples)}")

        # 2. 从表型文件读取样本ID
        pheno_samples = set()
        with open(pheno_file, 'r') as f:
            # 跳过可能的表头行
            first_line = f.readline().strip()

            # 检查第一行是否是数据行还是表头
            # 如果第一行包含<Trait>或不是数字，很可能是表头
            is_header = False
            if '<Trait>' in first_line or first_line.startswith('#'):
                is_header = True
            else:
                # 检查第二列是否为数字
                cols = first_line.split('\t')
                if len(cols) >= 2:
                    try:
                        float(cols[1])
                    except ValueError:
                        is_header = True

            # 如果第一行是表头，已经读取了，需要处理
            if is_header:
                # 从数据行读取样本
                # 将第一行加回来，因为后续循环需要处理所有行
                lines = [first_line] + f.readlines()
            else:
                lines = [first_line] + f.readlines()

            # 从数据行提取样本ID（跳过可能的<Trait>/<Data>行）
            for line in lines:
                line = line.strip()
                if line and not line.startswith('#') and '<Trait>' not in line and '<Data>' not in line:
                    fields = line.split('\t')
                    if len(fields) >= 1:
                        pheno_samples.add(fields[0])

        logger.info(f"   表型样本数|Phenotype samples: {len(pheno_samples)}")

        # 3. 找到交集
        common_samples = vcf_samples.intersection(pheno_samples)
        logger.info(f"   共有样本数|Common samples: {len(common_samples)}")

        if len(common_samples) == 0:
            logger.error("VCF和表型文件之间没有共有样本|No common samples found between VCF and phenotype files!")
            return False, None, None, 0

        # 4. 使用bcftools从VCF筛选样本
        keep_file = output_dir / f"{output_prefix}_keep.txt"
        vcf_filtered = output_dir / f"{output_prefix}_filtered.vcf.gz"

        # 写入样本列表（每行一个样本ID）
        with open(keep_file, 'w') as f:
            for sample_id in sorted(common_samples):
                f.write(f"{sample_id}\n")

        logger.info(f"   从{len(vcf_samples)}个样本筛选到{len(common_samples)}个共有样本|Filtering from {len(vcf_samples)} to {len(common_samples)} common samples...")

        # 使用bcftools筛选样本
        cmd = f"bcftools view -S {keep_file} -Oz -o {vcf_filtered} {vcf_file}"

        logger.info(f"   命令|Command: {cmd}")
        result = subprocess.run(cmd, shell=True, capture_output=True, text=True)

        if result.returncode != 0:
            logger.error(f"bcftools筛选VCF失败|bcftools VCF filtering failed: {result.stderr}")
            return False, None, None, 0

        logger.info("   VCF筛选完成|VCF filtering completed")

        # 5. 过滤表型文件并确保格式符合TASSEL要求
        pheno_filtered = output_dir / f"{output_prefix}_pheno_filtered.txt"

        # 读取原始表型文件
        with open(pheno_file, 'r') as f_in:
            lines = f_in.readlines()

        # 处理表型文件
        has_data_line = False
        header_processed = False

        with open(pheno_filtered, 'w') as f_out:
            for line in lines:
                line_stripped = line.strip()

                if not line_stripped:
                    # 跳过空行
                    continue

                if line_stripped.startswith('#'):
                    # 保留注释行
                    f_out.write(line)
                    continue

                # 检查是否是表头行（第一行）
                if not header_processed:
                    header_processed = True

                    # 处理表头：确保第一列是<Trait>
                    fields = line_stripped.split('\t')
                    if len(fields) >= 1:
                        # 将第一列替换为<Trait>
                        fields[0] = '<Trait>'
                        f_out.write('\t'.join(fields) + '\n')
                    continue

                # 检查是否是<Data>行
                if '<Data>' in line_stripped:
                    has_data_line = True
                    f_out.write(line)
                    continue

                # 数据行：检查样本ID是否在交集中
                fields = line_stripped.split('\t')
                if len(fields) >= 1 and fields[0] in common_samples:
                    f_out.write(line)

            # 如果没有<Data>行，添加它
            if not has_data_line and header_processed:
                # 读取表头确定列数
                with open(pheno_filtered, 'r') as f_tmp:
                    header_line = f_tmp.readline().strip()
                    num_cols = len(header_line.split('\t'))

                # 在表头后插入<Data>行
                with open(pheno_filtered, 'r') as f_tmp:
                    content = f_tmp.read()

                with open(pheno_filtered, 'w') as f_out:
                    lines_content = content.split('\n', 1)  # 分割第一行和其余内容
                    if len(lines_content) == 2:
                        f_out.write(lines_content[0] + '\n')  # 写入表头
                        f_out.write('\t'.join(['<Data>'] * num_cols) + '\n')  # 写入<Data>行
                        f_out.write(lines_content[1])  # 写入其余内容

        # 6. 清理临时文件
        try:
            keep_file.unlink()
        except Exception:
            pass

        logger.info(f"   样本交集处理完成|Common sample filtering completed")
        logger.info(f"   过滤后的VCF|Filtered VCF: {vcf_filtered}")
        logger.info(f"   过滤后的表型|Filtered phenotype: {pheno_filtered}")
        logger.info(f"   最终样本数|Final sample count: {len(common_samples)}")

        return True, vcf_filtered, pheno_filtered, len(common_samples)

    except Exception as e:
        logger.error(f"样本交集处理失败|Common sample filtering failed: {str(e)}")
        return False, None, None, 0