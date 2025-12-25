"""
EGAPx批量配置类 | EGAPx Batch Configuration Class
"""

import os


class EGAPxBatchConfig:
    """EGAPx批量运行配置类 | EGAPx Batch Processing Configuration Class"""

    def __init__(
        self,
        genome: str,
        yaml_template: str = None,
        script_template: str = None,
        output_dir: str = None,
        egapx_path: str = "/share/org/YZWL/yzwl_lixg/software/EGAPX_v.0.4.1-alpha/egapx",
        chr_prefix: str = None,
        locus_tag_prefix: str = "",
        report_name: str = "EGAPx",
        short_reads: str = "",
        long_reads: str = ""
    ):
        """
        初始化配置 | Initialize configuration

        Args:
            genome: 基因组FASTA文件路径 | Genome FASTA file path
            yaml_template: YAML模板文件路径 | YAML template file path (默认: ~/software/scripts/62.EGAPx_temple.yaml)
            script_template: Shell脚本模板路径 | Shell script template path (默认: ~/software/scripts/63.EGAPx_temple.sh)
            output_dir: 输出目录路径 | Output directory path
            egapx_path: EGAPx安装路径 | EGAPx installation path
            chr_prefix: 染色体前缀过滤 | Chromosome prefix filter
            locus_tag_prefix: locus标签前缀 | Locus tag prefix (默认: 空)
            report_name: 报告名称 | Report name
            short_reads: 短读长测序数据文件路径 | Short reads file path
            long_reads: 长读长测序数据文件路径 | Long reads file path
        """
        # 设置默认模板路径 | Set default template paths
        default_yaml = os.path.expanduser("~/software/scripts/62.EGAPx_temple.yaml")
        default_script = os.path.expanduser("~/software/scripts/63.EGAPx_temple.sh")

        self.genome = genome
        self.yaml_template = yaml_template or default_yaml
        self.script_template = script_template or default_script
        self.output_dir = output_dir
        self.egapx_path = egapx_path
        self.chr_prefix = chr_prefix
        self.locus_tag_prefix = locus_tag_prefix
        self.report_name = report_name
        self.short_reads = short_reads
        self.long_reads = long_reads

        # 运行时变量 | Runtime variables
        self.chromosomes = []
        self.job_count = 0

    def validate(self):
        """
        验证配置参数 | Validate configuration parameters

        Raises:
            ValueError: 当配置参数无效时 | When configuration parameters are invalid
        """
        # 检查必需文件 | Check required files
        if not self.genome or not os.path.exists(self.genome):
            raise ValueError(f"基因组文件不存在 | Genome file not found: {self.genome}")

        # 检查模板文件 | Check template files
        if not os.path.exists(self.yaml_template):
            raise ValueError(f"YAML模板不存在 | YAML template not found: {self.yaml_template}")

        if not os.path.exists(self.script_template):
            raise ValueError(f"脚本模板不存在 | Script template not found: {self.script_template}")

        if not self.egapx_path or not os.path.exists(self.egapx_path):
            raise ValueError(f"EGAPx路径不存在 | EGAPx path not found: {self.egapx_path}")

        # 检查输出目录 | Check output directory
        if not self.output_dir:
            raise ValueError("输出目录不能为空 | Output directory cannot be empty")

        return True

    def __repr__(self):
        """配置的字符串表示 | String representation of configuration"""
        return (
            f"EGAPxBatchConfig(\n"
            f"  genome={self.genome!r},\n"
            f"  yaml_template={self.yaml_template!r},\n"
            f"  script_template={self.script_template!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  egapx_path={self.egapx_path!r},\n"
            f"  chr_prefix={self.chr_prefix!r},\n"
            f"  locus_tag_prefix={self.locus_tag_prefix!r},\n"
            f"  report_name={self.report_name!r}\n"
            f")"
        )
