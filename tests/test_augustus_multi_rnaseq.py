"""
AugustusMultiRnaseq模块测试 | AugustusMultiRnaseq Module Tests
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch

from biopytools.augustus_multi_rnaseq import AugustusMultiRnaseqProcessor, AugustusMultiRnaseqConfig
from biopytools.common.base import ConfigurationError

class TestAugustusMultiRnaseqConfig:
    """测试配置类 | Test configuration class"""
    
    @pytest.fixture
    def temp_dir(self):
        """创建临时目录 | Create temporary directory"""
        temp_dir = tempfile.mkdtemp()
        yield Path(temp_dir)
        shutil.rmtree(temp_dir)
    
    def test_config_creation(self, temp_dir):
        """测试配置创建 | Test config creation"""
        # 创建测试输入目录
        input_dir = temp_dir / "input"
        input_dir.mkdir()
        
        config = AugustusMultiRnaseqConfig(
            input_dir=str(input_dir),
            output_dir=str(temp_dir / "output")
        )
        
        assert str(config.input_dir) == str(input_dir)
        assert config.threads > 0

class TestAugustusMultiRnaseqProcessor:
    """测试处理器类 | Test processor class"""
    
    @pytest.fixture
    def temp_dir(self):
        """创建临时目录 | Create temporary directory"""
        temp_dir = tempfile.mkdtemp()
        yield Path(temp_dir)
        shutil.rmtree(temp_dir)
    
    def test_processor_creation(self, temp_dir):
        """测试处理器创建 | Test processor creation"""
        input_dir = temp_dir / "input"
        input_dir.mkdir()
        
        processor = AugustusMultiRnaseqProcessor(
            input_dir=str(input_dir),
            output_dir=str(temp_dir / "output")
        )
        
        assert processor.config.input_dir == str(input_dir)
        assert processor.logger is not None

if __name__ == "__main__":
    pytest.main([__file__])
