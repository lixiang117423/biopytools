"""
TMHMM配置类|TMHMM Configuration Class
"""

import os

# 路径管理工具|Path management utilities
try:
    from common.paths import get_tool_path, expand_path
except ImportError:
    # 如果common模块不可用，使用简化版本|Fallback if common module unavailable
    def get_tool_path(tool_name, default_path, env_var=None):
        """获取工具路径|Get tool path"""
        if env_var and os.environ.get(env_var):
            return os.path.expandvars(os.path.expanduser(os.environ[env_var]))
        return os.path.expandvars(os.path.expanduser(default_path))

    def expand_path(path):
        """展开路径|Expand path"""
        return os.path.expandvars(os.path.expanduser(path))


class TmhmmConfig:
    """TMHMM配置类|TMHMM Configuration Class"""

    def __init__(
        self,
        input_file: str,
        output_dir: str,
        tmhmm_path: str = None,
        noplot: bool = True,
        output_prefix: str = None,
    ):
        """
        初始化配置|Initialize configuration

        Args:
            input_file: 输入蛋白质FASTA文件|Input protein FASTA file
            output_dir: 输出目录|Output directory
            tmhmm_path: tmhmm可执行文件路径|tmhmm binary path
            noplot: 不生成图形|No plot generation
            output_prefix: 输出文件前缀(默认使用输入文件名)|Output file prefix (default: input filename)
        """
        self.input_file = input_file
        self.output_dir = output_dir
        self.tmhmm_path = tmhmm_path
        self.noplot = noplot
        self.output_prefix = output_prefix

    def validate(self):
        """验证配置参数|Validate configuration parameters"""
        errors = []

        self.input_file = expand_path(self.input_file)
        self.output_dir = expand_path(self.output_dir)

        if not self.input_file:
            errors.append("输入文件不能为空|Input file cannot be empty")
        elif not os.path.exists(self.input_file):
            errors.append(f"输入文件不存在|Input file not found: {self.input_file}")

        if self.tmhmm_path:
            self.tmhmm_path = expand_path(self.tmhmm_path)
        else:
            self.tmhmm_path = get_tool_path(
                'tmhmm',
                '~/miniforge3/envs/tmmhmm_v.2.0c/bin/tmhmm',
                'TMHMM_PATH'
            )

        if not os.path.exists(self.tmhmm_path):
            errors.append(f"tmhmm不存在|tmhmm not found: {self.tmhmm_path}")

        if self.output_prefix is None:
            self.output_prefix = os.path.splitext(os.path.basename(self.input_file))[0]

        if errors:
            raise ValueError("\n".join(errors))

        os.makedirs(self.output_dir, exist_ok=True)
        return True

    def __repr__(self):
        return (
            f"TmhmmConfig(\n"
            f"  input_file={self.input_file!r},\n"
            f"  output_dir={self.output_dir!r},\n"
            f"  tmhmm_path={self.tmhmm_path!r}\n"
            f")"
        )
