"""
🔧 Ka/Ks Calculator配置管理模块
功能: 管理分析参数和默认设置 | Manage analysis parameters and default settings
"""

class KaKsConfig:
    """📋 Ka/Ks分析配置管理类 | Ka/Ks analysis configuration management"""
    
    # 🔧 默认参数 | Default parameters
    DEFAULT_METHOD = "GMYN"
    SUPPORTED_METHODS = [
        "GMYN", "MYN", "YN", "NG", "LWL", "LPB", 
        "MLWL", "MLPB", "GY", "MS", "MA", "GNG",
        "GLWL", "GLPB", "GMLWL", "GMLPB", "GYN"
    ]
    
    # 📁 文件扩展名 | File extensions
    FASTA_EXTENSIONS = ['.fasta', '.fa', '.fas', '.fna', '.cds']
    PAIR_EXTENSIONS = ['.txt', '.tsv', '.csv']
    AXT_EXTENSION = '.axt'
    
    # 🎯 输出文件名 | Output file names
    OUTPUT_FILES = {
        'summary': 'kaks_summary.xlsx',
        'detailed': 'kaks_detailed.tsv',
        'detailed_csv': 'kaks_detailed.csv', 
        'statistics': 'summary_stats.json',
        'report': 'analysis_report.html',
        'log': 'pipeline.log',
        'temp_axt': 'kaks_input.axt'
    }
    
    # ⚠️ 质量阈值 | Quality thresholds
    MIN_SEQUENCE_LENGTH = 50
    MAX_SEQUENCE_LENGTH = 50000
    MIN_SIMILARITY_THRESHOLD = 0.3
    MAX_MISSING_RATIO = 0.1
    
    # 🎨 选择压力分类 | Selection pressure classification
    SELECTION_THRESHOLDS = {
        'strong_negative': 0.5,
        'moderate_negative': 1.0,
        'neutral_lower': 0.95,
        'neutral_upper': 1.05,
        'weak_positive': 2.0
    }
    
    # 🧮 计算参数 | Calculation parameters
    CALCULATION_TIMEOUT = 300  # 5分钟超时 | 5 minutes timeout
    MAX_BATCH_SIZE = 1000     # 最大批处理大小 | Maximum batch size

    @classmethod
    def get_method_description(cls, method: str) -> str:
        """📖 获取计算方法描述 | Get calculation method description"""
        descriptions = {
            "GMYN": "🌟 推荐：Gamma分布修正的Yang-Nielsen方法 | Recommended: Gamma-corrected Yang-Nielsen",
            "MYN": "🔧 修正的Yang-Nielsen方法 | Modified Yang-Nielsen method",
            "YN": "📚 经典Yang-Nielsen方法 | Classic Yang-Nielsen method",
            "NG": "🏛️ Nei-Gojobori方法 | Nei-Gojobori method",
            "LWL": "⚡ Li-Wu-Luo方法 | Li-Wu-Luo method",
            "LPB": "🔬 Li-Pamilo-Bianchi方法 | Li-Pamilo-Bianchi method",
            "MLWL": "📈 修正的LWL方法 | Modified LWL method",
            "MLPB": "📊 修正的LPB方法 | Modified LPB method",
            "GY": "🎯 Goldman-Yang最大似然方法 | Goldman-Yang maximum likelihood",
            "MS": "🎯 基于AICc的模型选择 | Model Selection based on AICc",
            "MA": "📊 候选模型的模型平均 | Model Averaging on candidate models",
            "GNG": "🧬 Gamma-Nei-Gojobori方法 | Gamma-Nei-Gojobori method",
            "GLWL": "🔬 Gamma-LWL方法 | Gamma-LWL method",
            "GLPB": "📈 Gamma-LPB方法 | Gamma-LPB method",
            "GMLWL": "🧮 Gamma修正的MLWL方法 | Gamma-modified MLWL method",
            "GMLPB": "📊 Gamma修正的MLPB方法 | Gamma-modified MLPB method",
            "GYN": "🌟 Gamma-Yang-Nielsen方法 | Gamma-Yang-Nielsen method"
        }
        return descriptions.get(method, f"📊 {method}计算方法 | {method} calculation method")
