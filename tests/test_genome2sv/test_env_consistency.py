"""genome2sv 环境一致性单测|Env-consistency tests for genome2sv

比对管道(minimap2|samtools)必须单 conda 环境运行(§13.2.2 管道禁混 env);
六个工具全部固定 align 域环境全路径(svim-asm/SURVIVOR 2026-08-31 自
sv_calling 并入)。|Pipe tools must share one env (§13.2.2); all six tools
pin to align domain env full paths (svim-asm/SURVIVOR merged from
sv_calling on 2026-08-31).
"""
import os
import re
from types import SimpleNamespace

from biopytools.genome2sv.config import Genome2SVConfig
from biopytools.genome2sv.pipeline import Genome2SVPipeline


def _make_config(tmp_path) -> Genome2SVConfig:
    """默认工具路径的最小 config(fof 不存在则跳过样本解析)|Minimal config, default tool paths."""
    return Genome2SVConfig(input_fof="nonexistent.fof", ref_sample="ref",
                           output_dir=str(tmp_path / "out"))


def _env_of(tool_path: str):
    """从工具路径解析 conda 环境名|Parse conda env name from a tool path."""
    m = re.search(r"/envs/([^/]+)/", str(tool_path))
    return m.group(1) if m else None


def test_all_tools_pinned_to_align(tmp_path):
    """六工具默认全部固定 align 域环境|All six tools pin to align domain env."""
    config = _make_config(tmp_path)
    for path in (config.minimap2_path, config.samtools_path, config.svim_asm_path,
                 config.bcftools_path, config.bedtools_path, config.survivor_path):
        assert _env_of(path) == "align", path


def test_align_command_single_env_no_warning(tmp_path):
    """比对管道单 conda run 环境且无 mismatch WARNING|Single-env pipe without mismatch WARNING."""
    logs = []
    logger = SimpleNamespace(info=lambda m: logs.append(("info", m)),
                             warning=lambda m: logs.append(("warning", m)),
                             error=lambda m: logs.append(("error", m)))
    config = SimpleNamespace(
        alignment_dir=tmp_path / "01_alignment",
        reference_fasta=str(tmp_path / "ref.fasta"),
        threads=4, preset="asm5",
        minimap2_path=os.path.expanduser("~/miniforge3/envs/align/bin/minimap2"),
        samtools_path=os.path.expanduser("~/miniforge3/envs/align/bin/samtools"),
    )
    cmd = Genome2SVPipeline(config, logger)._align_command("s1", str(tmp_path / "q.fa"))
    # 同源conda绝对路径前缀: '<conda> run -p <env绝对前缀>'(解析不到回退 '-n align')|
    # Same-installation conda prefix: 'run -p <abs env prefix>' (or '-n align' fallback)
    assert cmd[1] == "run" and cmd[2] in ("-p", "-n")
    assert cmd[3].rstrip("/").endswith("align")
    assert "--no-capture-output" in cmd
    assert sum(t in ("-p", "-n") for t in cmd) == 1
    assert not any(k == "warning" for k, _ in logs)
