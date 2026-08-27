"""splitstree6 配置|SplitsTree6 configuration"""
import os
from dataclasses import dataclass, field
from pathlib import Path
from typing import List, Optional

from ..common.paths import expand_path, get_tool_path

# 官方导出格式(RunWorkflow -e 合法值)|Official export formats
EXPORT_FORMATS = ('Newick', 'Nexus', 'GML', 'PlainText', 'NeXML',
                  'Phylip', 'FastA', 'Clustal')
# 默认输出格式|Default output formats
DEFAULT_EXPORT_FORMATS = ('Newick', 'Nexus', 'GML')
# 输入格式(RunWorkflow -f 合法值)|Legal input formats
INPUT_FORMATS = ('Unknown', 'CSV', 'FastA', 'GML', 'MSF', 'Newick',
                 'Nexml', 'Nexus', 'Phylip', 'Stockholm')

DEFAULT_WORKFLOW = '~/software/splitstree6/examples/large_DNA_phyml/' \
                   '10taxaExample.stree6'


@dataclass
class Splitstree6Config:
    """SplitsTree6 建网配置类|SplitsTree6 network/tree configuration class"""

    # 必需|Required
    input: str                        # 输入数据(VCF 时先转距离矩阵)|input data

    # 输出|Output
    output_dir: str = './splitstree6_output'
    export_formats: Optional[str] = None   # 逗号分隔;默认 Newick,Nexus,GML

    # 网络参数|Network parameters
    workflow: str = DEFAULT_WORKFLOW       # 官方 NeighborNet 工作流
    input_format: str = ''                 # 自动识别时留空|auto-detect if empty
    node_name: str = 'Splits'              # 导出节点(默认 Splits 网络)|export node

    # 运行|Runtime
    threads: int = 12
    log_level: str = 'INFO'

    # 工具路径|Tool paths
    splitstree_tools_dir: str = field(
        default_factory=lambda: get_tool_path(
            'splitstree6-tools', '~/software/splitstree6/splitstree6-tools-bin',
            'SPLITSTREE6_TOOLS_DIR'))
    xvfb_path: str = field(
        default_factory=lambda: get_tool_path(
            'Xvfb', '~/miniforge3/envs/xvfb_test/x86_64-conda-linux-gnu/'
                    'sysroot/usr/bin/Xvfb', 'XVFB_PATH'))

    def __post_init__(self):
        """展开路径、解析格式列表、建目录|Expand paths, parse formats, make dirs"""
        self.input = expand_path(self.input)
        self.output_dir = expand_path(self.output_dir)

        out = Path(self.output_dir)
        self.wf_dir = out / '01_workflow'
        self.info_dir = out / '00_pipeline_info'
        self.logs_dir = out / '99_logs'
        for d in (out, self.wf_dir, self.info_dir, self.logs_dir):
            d.mkdir(parents=True, exist_ok=True)

        self.report_file = out / '01_workflow/report.txt'
        self.log_file = self.logs_dir / 'splitstree6.log'

        # 导出格式解析|parse export formats
        if not self.export_formats:
            self.export_formats_list: List[str] = list(DEFAULT_EXPORT_FORMATS)
        else:
            parts = [p.strip() for p in self.export_formats.split(',') if p.strip()]
            bad = [p for p in parts if p not in EXPORT_FORMATS]
            if bad:
                raise ValueError(
                    f"非法导出格式|Illegal export formats: {bad} "
                    f"(可选|valid: {', '.join(EXPORT_FORMATS)})")
            self.export_formats_list = parts

        # 工作流与工具校验路径|resolve paths (expand later in validate)
        self.workflow = os.path.abspath(expand_path(self.workflow))
        self.splitstree_tools_dir = expand_path(self.splitstree_tools_dir)
        self.xvfb_path = expand_path(self.xvfb_path)

    def validate(self):
        """校验配置|Validate configuration"""
        errors = []
        wf = Path(self.workflow)
        tools_lib = Path(self.splitstree_tools_dir) / 'lib' / 'jars'
        # 输入存在性独立检查(VCF 等文件名含点,不能用"含点"判定)
        # |input existence checked independently (names like .vcf contain dots)
        if not Path(self.input).exists():
            errors.append(f"输入不存在|Input not found: {self.input}")
        if not wf.exists():
            errors.append(f"工作流文件不存在|Workflow file not found: {self.workflow}")
        if not tools_lib.is_dir():
            errors.append(f"SplitsTree6 jars 目录不存在|jars dir not found: "
                          f"{tools_lib}(需要 splitstree6-tools 安装包)")
        launcher = Path(self.splitstree_tools_dir) / 'tools' / 'workflow-run'
        if not launcher.exists():
            errors.append(f"workflow-run 未找到|workflow-run not found: {launcher}")
        java_bin = os.path.expanduser('~/.local/bin/java')
        if not (os.path.exists(java_bin) or os.path.exists('/usr/bin/java')):
            errors.append("java 不在 PATH 中|java not found on PATH")
        if errors:
            raise ValueError('\n'.join(errors))
        return True
