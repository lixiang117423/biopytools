"""
VcfLdHeatmap模块测试 | VcfLdHeatmap Module Tests
"""

import pytest
import tempfile
import shutil
from pathlib import Path
from unittest.mock import Mock, patch

from biopytools.vcf_ld_heatmap import VcfLdHeatmapProcessor, VcfLdHeatmapConfig
from biopytools.common.base import ConfigurationError

class TestVcfLdHeatmapConfig:
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
        
        config = VcfLdHeatmapConfig(
            input_dir=str(input_dir),
            output_dir=str(temp_dir / "output")
        )
        
        assert str(config.input_dir) == str(input_dir)
        assert config.threads > 0

class TestVcfLdHeatmapProcessor:
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
        
        processor = VcfLdHeatmapProcessor(
            input_dir=str(input_dir),
            output_dir=str(temp_dir / "output")
        )
        
        assert processor.config.input_dir == str(input_dir)
        assert processor.logger is not None

if __name__ == "__main__":
    pytest.main([__file__])
