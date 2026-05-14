"""
HTML模板和样式定义|HTML Templates and Style Definitions
"""

def get_css_style(theme='modern'):
    """获取CSS样式|Get CSS style"""
    
    if theme == 'modern':
        return """
        * {
            margin: 0;
            padding: 0;
            box-sizing: border-box;
        }
        
        body {
            font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, 'Helvetica Neue', Arial, sans-serif;
            background: #f5f7fa;
            color: #2c3e50;
            line-height: 1.6;
            padding: 20px;
        }
        
        .container {
            max-width: 1400px;
            margin: 0 auto;
            background: white;
            border-radius: 8px;
            box-shadow: 0 2px 12px rgba(0,0,0,0.1);
            overflow: hidden;
        }
        
        header {
            background: linear-gradient(135deg, #667eea 0%, #764ba2 100%);
            color: white;
            padding: 30px;
        }
        
        header h1 {
            font-size: 28px;
            margin-bottom: 10px;
        }
        
        .analysis-info {
            display: flex;
            gap: 20px;
            font-size: 14px;
            opacity: 0.9;
            flex-wrap: wrap;
        }
        
        .analysis-info span {
            background: rgba(255,255,255,0.2);
            padding: 5px 12px;
            border-radius: 4px;
        }
        
        main {
            padding: 30px;
        }
        
        .search-filter {
            background: #f8f9fa;
            padding: 20px;
            border-radius: 6px;
            margin-bottom: 30px;
            display: flex;
            gap: 15px;
            flex-wrap: wrap;
            align-items: center;
        }
        
        .search-filter input, .search-filter select {
            padding: 10px 15px;
            border: 1px solid #ddd;
            border-radius: 4px;
            font-size: 14px;
            flex: 1;
            min-width: 200px;
        }
        
        .search-filter button {
            padding: 10px 20px;
            background: #667eea;
            color: white;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: background 0.3s;
        }
        
        .search-filter button:hover {
            background: #5568d3;
        }
        
        .sample-list, .alignment-list {
            display: grid;
            gap: 15px;
        }
        
        .sample-item, .alignment-item {
            background: white;
            border: 1px solid #e0e0e0;
            border-radius: 6px;
            padding: 20px;
            transition: transform 0.2s, box-shadow 0.2s;
        }
        
        .sample-item:hover, .alignment-item:hover {
            transform: translateY(-2px);
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
        }
        
        .sample-header, .alignment-header {
            display: flex;
            justify-content: space-between;
            align-items: center;
            margin-bottom: 15px;
        }
        
        .sample-header h3, .alignment-header h3 {
            color: #2c3e50;
            font-size: 18px;
        }
        
        .metrics {
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            margin-top: 10px;
        }
        
        .metric {
            font-size: 14px;
            color: #666;
        }
        
        .metric strong {
            color: #2c3e50;
        }
        
        .identity-badge {
            padding: 4px 12px;
            border-radius: 12px;
            font-size: 13px;
            font-weight: 600;
        }
        
        .identity-high { background: #d4edda; color: #155724; }
        .identity-medium { background: #fff3cd; color: #856404; }
        .identity-low { background: #f8d7da; color: #721c24; }
        
        .btn {
            padding: 8px 16px;
            border: none;
            border-radius: 4px;
            cursor: pointer;
            font-size: 14px;
            transition: all 0.3s;
        }
        
        .btn-primary {
            background: #667eea;
            color: white;
        }
        
        .btn-primary:hover {
            background: #5568d3;
        }
        
        .btn-secondary {
            background: #6c757d;
            color: white;
        }
        
        .btn-secondary:hover {
            background: #5a6268;
        }
        
        .alignment-content {
            margin-top: 15px;
            padding-top: 15px;
            border-top: 2px solid #f0f0f0;
        }
        
        .alignment-view {
            font-family: 'Courier New', Consolas, Monaco, monospace;
            font-size: 13px;
            line-height: 1.6;
            background: #f8f9fa;
            padding: 20px;
            border-radius: 4px;
            overflow-x: auto;
            white-space: pre;
        }
        
        .seq-line {
            margin-bottom: 1px;
        }
        
        .match { color: #28a745; }
        .mismatch { color: #fd7e14; }
        .gap { color: #dc3545; }
        
        .alignment-stats {
            margin-top: 15px;
            padding: 15px;
            background: #f8f9fa;
            border-radius: 4px;
            display: flex;
            gap: 20px;
            flex-wrap: wrap;
            font-size: 14px;
        }
        
        .statistics {
            margin-top: 30px;
            padding: 20px;
            background: #f8f9fa;
            border-radius: 6px;
        }
        
        .statistics h2 {
            color: #2c3e50;
            margin-bottom: 15px;
        }
        
        .back-link {
            display: inline-block;
            margin-bottom: 20px;
            color: #667eea;
            text-decoration: none;
            font-size: 14px;
        }
        
        .back-link:hover {
            text-decoration: underline;
        }
        
        @media (max-width: 768px) {
            body {
                padding: 10px;
            }
            
            .container {
                border-radius: 0;
            }
            
            header {
                padding: 20px;
            }
            
            main {
                padding: 20px;
            }
            
            .search-filter {
                flex-direction: column;
            }
            
            .search-filter input, .search-filter select {
                width: 100%;
            }
        }
        """
    
    return ""

def get_javascript():
    """获取JavaScript代码|Get JavaScript code"""
    return """
    function applyFilters() {
        const searchTerm = document.getElementById('search')?.value.toLowerCase() || '';
        const minIdentity = parseFloat(document.getElementById('identity-filter')?.value || 0);
        
        document.querySelectorAll('.sample-item, .alignment-item').forEach(item => {
            const name = (item.dataset.sampleName || item.dataset.alignmentName || '').toLowerCase();
            const identity = parseFloat(item.dataset.identity || 100);
            
            const matchesSearch = name.includes(searchTerm);
            const matchesIdentity = identity >= minIdentity;
            
            item.style.display = (matchesSearch && matchesIdentity) ? 'block' : 'none';
        });
    }
    
    function toggleAlignment(alignmentId) {
        const content = document.getElementById(alignmentId);
        if (content) {
            content.style.display = content.style.display === 'none' ? 'block' : 'none';
        }
    }
    
    function expandAll() {
        document.querySelectorAll('.alignment-content').forEach(el => {
            el.style.display = 'block';
        });
    }
    
    function collapseAll() {
        document.querySelectorAll('.alignment-content').forEach(el => {
            el.style.display = 'none';
        });
    }
    
    function copySequence(alignmentId) {
        const content = document.getElementById(alignmentId);
        if (content) {
            const text = content.querySelector('.alignment-view')?.innerText || '';
            navigator.clipboard.writeText(text).then(() => {
                showNotification(' 已复制到剪贴板');
            }).catch(err => {
                console.error('复制失败:', err);
                showNotification(' 复制失败');
            });
        }
    }
    
    function showNotification(message) {
        const notification = document.createElement('div');
        notification.textContent = message;
        notification.style.cssText = 'position: fixed; top: 20px; right: 20px; background: #333; color: white; padding: 15px 20px; border-radius: 4px; z-index: 1000; animation: slideIn 0.3s ease;';
        document.body.appendChild(notification);
        
        setTimeout(() => {
            notification.style.animation = 'slideOut 0.3s ease';
            setTimeout(() => notification.remove(), 300);
        }, 2000);
    }
    
    // 实时搜索
    document.getElementById('search')?.addEventListener('input', applyFilters);
    """

def get_index_template():
    """获取主页模板|Get index page template"""
    return """<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>BLAST Alignment Visualization</title>
    <style>
{css_content}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1> BLAST Alignment Visualization Report</h1>
            <div class="analysis-info">
                <span>分析日期|Date: {analysis_date}</span>
                <span> BLAST类型: {blast_type}</span>
                <span> 样品数: {sample_count}</span>
                <span> 总比对数: {total_alignments}</span>
            </div>
        </header>
        
        <main>
            <section class="search-filter">
                <input type="text" id="search" placeholder=" 搜索样品名称...">
                <select id="identity-filter">
                    <option value="0">所有相似度</option>
                    <option value="70">≥ 70%</option>
                    <option value="80">≥ 80%</option>
                    <option value="90">≥ 90%</option>
                    <option value="95">≥ 95%</option>
                </select>
                <button onclick="applyFilters()">应用筛选</button>
            </section>
            
            <section class="sample-list">
                <h2> 样品列表</h2>
                {sample_list_html}
            </section>
            
            <section class="statistics">
                <h2>快速统计|Quick Statistics</h2>
                {statistics_html}
            </section>
        </main>
    </div>
    
    <script>
{js_content}
    </script>
</body>
</html>
"""

def get_sample_template():
    """获取单样品页面模板|Get single sample page template"""
    return """<!DOCTYPE html>
<html lang="zh-CN">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{sample_name} - BLAST Alignments</title>
    <style>
{css_content}
    </style>
</head>
<body>
    <div class="container">
        <header>
            <h1> {sample_name} - Sequence Alignments</h1>
            <div class="analysis-info">
                <span> 输入文件: {file_name}</span>
                <span> 比对数量: {alignment_count}</span>
                <span> 平均相似度: {avg_identity:.2f}%</span>
            </div>
        </header>
        
        <main>
            <a href="index.html" class="back-link">← 返回主页</a>
            
            <section class="search-filter">
                <input type="text" id="search" placeholder=" 搜索目标序列ID...">
                <select id="identity-filter">
                    <option value="0">所有相似度</option>
                    <option value="70">≥ 70%</option>
                    <option value="80">≥ 80%</option>
                    <option value="90">≥ 90%</option>
                    <option value="95">≥ 95%</option>
                </select>
                <button onclick="applyFilters()">应用筛选</button>
                <button onclick="expandAll()" class="btn-secondary">展开全部</button>
                <button onclick="collapseAll()" class="btn-secondary">折叠全部</button>
            </section>
            
            <section class="alignment-list">
                {alignments_html}
            </section>
        </main>
    </div>
    
    <script>
{js_content}
    </script>
</body>
</html>
"""