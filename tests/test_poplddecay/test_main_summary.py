"""PopLDdecayRunner汇总输出测试|Tests for PopLDdecayRunner summary output"""

from types import SimpleNamespace

import numpy as np
import pandas as pd

from biopytools.poplddecay.main import PopLDdecayRunner


class _RecorderLogger:
    """记录所有日志消息的假logger|Fake logger recording every message"""

    def __init__(self):
        self.messages = []

    def info(self, msg):
        self.messages.append(str(msg))

    def warning(self, msg):
        self.messages.append(str(msg))

    def error(self, msg):
        self.messages.append(str(msg))


class TestOutputThresholdSummary:

    def test_logs_full_decay_table(self, tmp_path):
        """汇总输出应把每个群体的完整衰减表写入日志|Summary should log each population's full decay table"""
        rec = {
            'population': 'pop1',
            'n': 60,
            'rho': 0.0105,
            'background_r2': 0.2351,
            'recommended': pd.DataFrame([{'threshold': 0.4, 'decay_kb': 27.0, 'bg_ratio': 1.7}]),
            'decay_table': pd.DataFrame([
                {'threshold': 0.5, 'decay_kb': np.nan, 'bg_ratio': 2.13,
                 'recommendation': 'not recommended: above fitted curve max'},
                {'threshold': 0.4, 'decay_kb': 27.0, 'bg_ratio': 1.7,
                 'recommendation': 'recommended'},
            ]),
            'gwas_window_kb': 50,
        }
        logger = _RecorderLogger()
        config = SimpleNamespace(
            get_output_files=lambda: {'threshold': str(tmp_path / 'threshold.tsv')})
        runner = object.__new__(PopLDdecayRunner)
        runner.config = config
        runner.logger = logger

        runner._output_threshold_summary([rec])

        joined = '\n'.join(logger.messages)
        assert 'pop1' in joined, "应标明群体名|Population name should appear"
        assert 'above fitted curve max' in joined, \
            "完整衰减表(含不可达标记)应写入日志|Full decay table incl. unattainable marks should be logged"
        assert (tmp_path / 'threshold.tsv').exists(), "汇总tsv仍应保存|Summary tsv should still be saved"
