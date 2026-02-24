# =====================================
# EZA (exa继任者) ZSH 配置
# 将此配置添加到你的 ~/.zshrc 文件中
# =====================================

# 添加 .local/bin 到 PATH（如果不存在）
if [[ -d "$HOME/.local/bin" ]] && [[ ":$PATH:" != *":$HOME/.local/bin:"* ]]; then
    export PATH="$HOME/.local/bin:$PATH"
fi

# 检查并设置 EZA 命令
if [[ -x "$HOME/.local/bin/eza" ]]; then
    EZA_CMD="$HOME/.local/bin/eza"
elif command -v eza &> /dev/null; then
    EZA_CMD="eza"
elif command -v exa &> /dev/null; then
    EZA_CMD="exa"
else
    EZA_CMD="ls"
fi

# =====================================
# 颜色和样式设置
# =====================================

# 设置 EZA 颜色主题
export EZA_COLORS="ur=0:uw=0:ux=0:ue=0:gr=0:gw=0:gx=0:tr=0:tw=0:tx=0:su=0:sf=0:xa=0:sn=38;5;244:sb=38;5;244:da=38;5;61:gm=38;5;203:gd=38;5;203:gv=38;5;203:gt=38;5;203"
export EZA_TIME_STYLE="long-iso"

# =====================================
# 基础别名 - 替换传统 ls 命令
# =====================================

alias ls="$EZA_CMD --icons --group-directories-first"
alias ll="$EZA_CMD -l --icons --group-directories-first --time-style=long-iso"
alias la="$EZA_CMD -la --icons --group-directories-first --time-style=long-iso"
alias l="$EZA_CMD -la --icons --group-directories-first --git --time-style=long-iso --total-size"
alias lf="$EZA_CMD -la --icons --group-directories-first --git --time-style=long-iso"

# =====================================
# 简短易用的别名
# =====================================
alias esa="$EZA_CMD --icons --group-directories-first"
alias esal="$EZA_CMD -la --icons --group-directories-first --git --time-style=long-iso"
alias esat="$EZA_CMD --tree --icons --git-ignore"

alias esa-help="eza_help"
alias esa-man="eza_man"
alias esa-options="eza_options"
alias esa-examples="eza_examples"
alias esa-h="eza_help"

alias esa-recent="lrecent"
alias esa-big="lbig"
alias esa-stat="lstat"
alias esa-git="lgit"
alias esa-type="ltype"
alias esa-proj="lproj"
alias esa-tree="ltree"

alias esa-dark="eza_dark"
alias esa-light="eza_light"

# =====================================
# 扩展别名 - 特定用途
# =====================================
alias lt2="$EZA_CMD --tree --level=2 --icons --git-ignore"
alias lt3="$EZA_CMD --tree --level=3 --icons --git-ignore"
alias lg="$EZA_CMD -la --git --icons --group-directories-first"
alias lgs="$EZA_CMD -la --git --git-repos --icons --group-directories-first"
alias lm="$EZA_CMD -la --sort=modified --icons --time-style=long-iso"
alias ls_size="$EZA_CMD -la --sort=size --icons"
alias lnn="$EZA_CMD -la --sort=name --icons"
alias le="$EZA_CMD -la --sort=extension --icons"
alias ld="$EZA_CMD -D --icons"
alias lh="$EZA_CMD -la --header --icons"

# =====================================
# 实用函数
# =====================================

# 智能树状显示函数
function ltree() {
    local level=${1:-2}
    local path=${2:-.}
    $EZA_CMD --tree --level=$level --icons --git-ignore "$path"
}

# 详细的项目浏览函数
function lproj() {
    local path=${1:-.}
    echo "📁 项目概览: $path"
    echo "===================="
    $EZA_CMD -la --git --icons --group-directories-first --time-style=long-iso "$path"
    echo ""
    echo "🌲 目录结构:"
    echo "===================="
    $EZA_CMD --tree --level=2 --icons --git-ignore "$path"
}

# 查找最近修改的文件
function lrecent() {
    local count=${1:-10}
    local path=${2:-.}
    echo "🕒 最近修改的 $count 个文件:"
    echo "========================="
    local current_line=0
    while IFS= read -r line && [[ $current_line -lt $((count + 1)) ]]; do
        echo "$line"
        ((current_line++))
    done < <($EZA_CMD -la --sort=modified --icons --time-style=long-iso "$path" 2>/dev/null)
}

# 查找最大的文件
function lbig() {
    local count=${1:-10}
    local path=${2:-.}
    echo "📊 最大的 $count 个文件:"
    echo "===================="
    local current_line=0
    while IFS= read -r line && [[ $current_line -lt $((count + 1)) ]]; do
        echo "$line"
        ((current_line++))
    done < <($EZA_CMD -la --sort=size --reverse --icons "$path" 2>/dev/null)
}

# 显示文件统计信息
function lstat() {
    local path=${1:-.}
    echo "📈 文件统计: $path"
    echo "=================="

    if [[ ! -d "$path" ]]; then
        echo "❌ 路径不存在: $path"
        return 1
    fi

    echo "正在统计文件..."

    local total_files=0
    local dir_count=0
    local hidden_count=0
    local visible_files=0

    while IFS= read -r line; do
        [[ -n "$line" ]] && ((total_files++))
    done < <($EZA_CMD -a "$path" 2>/dev/null)

    while IFS= read -r line; do
        [[ -n "$line" ]] && ((dir_count++))
    done < <($EZA_CMD -D "$path" 2>/dev/null)

    while IFS= read -r line; do
        if [[ -n "$line" && "$line" =~ ^\. ]]; then
            ((hidden_count++))
        fi
    done < <($EZA_CMD -a "$path" 2>/dev/null)

    visible_files=$((total_files - hidden_count))

    echo "总文件数: $total_files"
    echo "目录数量: $dir_count"
    echo "隐藏文件: $hidden_count"
    echo "可见文件: $visible_files"
    echo "普通文件: $((total_files - dir_count))"

    echo ""
    echo "📋 详细信息:"
    echo "============"
    $EZA_CMD -la --icons --group-directories-first "$path"
}

# Git 状态快速查看
function lgit() {
    local path=${1:-.}
    if [[ -d "$path/.git" ]]; then
        echo "🔄 Git 状态: $path"
        echo "================="
        $EZA_CMD -la --git --icons --group-directories-first "$path"
        echo ""
        echo "🌿 Git 分支信息:"
        echo "==============="
        git -C "$path" status --short
    else
        echo "❌ 不是 Git 仓库: $path"
        $EZA_CMD -la --icons --group-directories-first "$path"
    fi
}

# 智能搜索函数 (结合 fzf)
function lsearch() {
    if command -v fzf &> /dev/null; then
        $EZA_CMD -la --icons --group-directories-first | fzf --preview="echo {}" --preview-window=right:50%
    else
        echo "需要安装 fzf 才能使用搜索功能"
        $EZA_CMD -la --icons --group-directories-first
    fi
}

# 比较两个目录
function lcompare() {
    if [[ $# -ne 2 ]]; then
        echo "用法: lcompare <目录1> <目录2>"
        return 1
    fi

    echo "📂 目录比较:"
    echo "============"
    echo "目录1: $1"
    $EZA_CMD -la --icons --group-directories-first "$1"
    echo ""
    echo "目录2: $2"
    $EZA_CMD -la --icons --group-directories-first "$2"
}

# 按文件类型分类显示
function ltype() {
    local path=${1:-.}
    echo "📋 按类型分类: $path"
    echo "=================="

    echo "📁 目录:"
    $EZA_CMD -D --icons "$path" 2>/dev/null || echo "  (无目录)"
    echo ""

    echo "🖼️  图片文件:"
    $EZA_CMD -a --icons "$path" 2>/dev/null | grep -E '\.(jpg|jpeg|png|gif|bmp|svg|webp)' || echo "  (无图片文件)"
    echo ""

    echo "📄 文档文件:"
    $EZA_CMD -a --icons "$path" 2>/dev/null | grep -E '\.(txt|md|pdf|doc|docx|rtf)' || echo "  (无文档文件)"
    echo ""

    echo "💻 代码文件:"
    $EZA_CMD -a --icons "$path" 2>/dev/null | grep -E '\.(js|ts|py|java|cpp|c|h|css|html|php|rb|go|rs|sh|zsh)' || echo "  (无代码文件)"
    echo ""

    echo "📦 压缩文件:"
    $EZA_CMD -a --icons "$path" 2>/dev/null | grep -E '\.(zip|tar|gz|rar|7z|xz|bz2)' || echo "  (无压缩文件)"
}

# =====================================
# 动态主题切换
# =====================================

function eza_dark() {
    export EZA_COLORS="ur=0:uw=0:ux=0:ue=0:gr=0:gw=0:gx=0:tr=0:tw=0:tx=0:su=0:sf=0:xa=0:sn=38;5;244:sb=38;5;244:da=38;5;61:gm=38;5;203:gd=38;5;203:gv=38;5;203:gt=38;5;203"
    echo "🌙 切换到暗色主题"
}

function eza_light() {
    export EZA_COLORS="ur=0:uw=0:ux=0:ue=0:gr=0:gw=0:gx=0:tr=0:tw=0:tx=0:su=0:sf=0:xa=0:sn=38;5;100:sb=38;5;100:da=38;5;33:gm=38;5;196:gd=38;5;196:gv=38;5;196:gt=38;5;196"
    echo "☀️ 切换到亮色主题"
}

# =====================================
# 快捷键绑定
# =====================================

# Ctrl+L 清屏并显示当前目录内容
clear_and_list() {
    clear
    echo "📍 当前位置: $(pwd)"
    echo "=============="
    $EZA_CMD -la --icons --group-directories-first --git
}
zle -N clear_and_list
bindkey '^L' clear_and_list

# =====================================
# 帮助函数
# =====================================

function eza_man() {
    if [[ "$EZA_CMD" == *"eza"* ]] || [[ "$EZA_CMD" == *"exa"* ]]; then
        echo "📖 EZA/EXA 原始帮助文档:"
        echo "======================="
        $EZA_CMD --help
    else
        echo "📖 LS 帮助文档:"
        echo "==============="
        ls --help 2>/dev/null || man ls
    fi
}

function eza_options() {
    if [[ "$EZA_CMD" == *"eza"* ]] || [[ "$EZA_CMD" == *"exa"* ]]; then
        echo "🔧 EZA/EXA 常用选项详解:"
        echo "======================="
        echo ""
        echo "📋 基础选项:"
        echo "  -l, --long         使用长格式显示详细信息"
        echo "  -a, --all          显示隐藏文件和以.开头的文件"
        echo "  -A, --almost-all   显示隐藏文件，但不显示. 和 .."
        echo "  -1, --oneline      每行显示一个文件"
        echo "  -r, --reverse      反向排序"
        echo "  -s, --sort=WORD    排序方式: name|size|extension|modified|created|accessed|type|inode"
        echo ""
        echo "🎨 显示选项:"
        echo "  --icons            显示文件类型图标"
        echo "  --no-icons         不显示图标"
        echo "  --color=WHEN       颜色显示: auto|always|never"
        echo "  --color-scale      根据年龄/大小使用颜色渐变"
        echo "  --group-directories-first  目录优先显示"
        echo ""
        echo "📊 信息选项:"
        echo "  -b, --binary       以二进制前缀显示文件大小"
        echo "  -B, --bytes        以字节为单位显示文件大小"
        echo "  -h, --header       显示表头"
        echo "  -H, --links        显示硬链接数"
        echo "  -i, --inode        显示 inode 号"
        echo "  -m, --modified     显示修改时间"
        echo "  -S, --blocks       显示文件系统块数"
        echo "  -u, --accessed     显示访问时间"
        echo "  -U, --created      显示创建时间"
        echo ""
        echo "🌳 树状显示:"
        echo "  -T, --tree         以树状结构显示"
        echo "  -L, --level=NUM    限制树的层级深度"
        echo "  -I, --ignore-glob=GLOB  忽略匹配的文件"
        echo "  --git-ignore       忽略.gitignore中的文件"
        echo ""
        echo "🔄 Git 集成:"
        echo "  --git              显示每个文件的Git状态"
        echo "  --git-repos        显示Git仓库状态"
        echo ""
        echo "⏰ 时间格式:"
        echo "  --time-style=STYLE 时间显示格式:"
        echo "    default          默认格式"
        echo "    iso              ISO 8601格式"
        echo "    long-iso         长ISO格式 (YYYY-MM-DD HH:MM)"
        echo "    full-iso         完整ISO格式"
        echo "    relative         相对时间 (2 days ago)"
        echo ""
        echo "🔍 过滤选项:"
        echo "  -D, --only-dirs    只显示目录"
        echo "  -f, --only-files   只显示文件"
        echo "  --group=GROUP      按组过滤"
        echo "  --owner=USER       按用户过滤"
    else
        echo "当前使用的是系统 ls 命令，选项有限"
        ls --help 2>/dev/null || echo "请安装 eza 获得更多功能"
    fi
}

function eza_examples() {
    echo "🚀 EZA 使用示例:"
    echo "==============="
    echo ""
    echo "⚡ 推荐: 简短易输入的命令"
    echo "  esa                    # 基础列表 + 图标"
    echo "  esal                   # 详细列表 + Git状态 + 图标"
    echo "  esat                   # 树状结构 + 图标"
    echo "  esa-recent 5           # 最近5个文件"
    echo "  esa-big 10             # 最大10个文件"
    echo "  esa-git                # Git状态查看"
    echo "  esa-tree 3             # 3层树状结构"
    echo ""
    echo "📁 传统基础用法:"
    echo "  eza                    # 基本列表"
    echo "  eza -l                 # 长格式"
    echo "  eza -la                # 长格式 + 隐藏文件"
    echo "  eza -la --icons        # 添加图标"
    echo ""
    echo "🌳 树状显示:"
    echo "  eza --tree             # 树状结构"
    echo "  eza --tree -L 2        # 限制2层"
    echo "  eza --tree --git-ignore # 忽略Git文件"
    echo "  esat -L 3              # 简短版本: 3层树状"
    echo ""
    echo "📊 排序:"
    echo "  eza -l --sort=size     # 按大小排序"
    echo "  eza -l --sort=modified # 按修改时间"
    echo "  eza -l --sort=extension # 按扩展名"
    echo ""
    echo "🔄 Git 集成:"
    echo "  eza -la --git          # 显示Git状态"
    echo "  eza --git-repos        # 显示仓库状态"
    echo "  esal                   # 简短版本: 详细+Git"
    echo ""
    echo "⏰ 时间格式:"
    echo "  eza -l --time-style=long-iso    # 长ISO格式"
    echo "  eza -l --time-style=relative    # 相对时间"
    echo ""
    echo "🎨 颜色主题:"
    echo "  eza --color=always     # 强制颜色"
    echo "  eza --color-scale      # 颜色渐变"
    echo "  esa-dark               # 切换暗色主题"
    echo "  esa-light              # 切换亮色主题"
    echo ""
    echo "📋 实用组合 (推荐):"
    echo "  # 完整功能显示"
    echo "  esal"
    echo "  # 等价于:"
    echo "  eza -la --git --icons --group-directories-first --time-style=long-iso"
    echo ""
    echo "  # 项目结构浏览"
    echo "  esa-tree 3 src/"
    echo "  # 等价于:"
    echo "  eza --tree --level=3 --icons --git-ignore src/"
    echo ""
    echo "  # 查找大文件"
    echo "  esa-big 5"
    echo "  # 等价于:"
    echo "  eza -l --sort=size --reverse --icons | head -6"
    echo ""
    echo "💡 日常工作流建议:"
    echo "  esa                    # 快速浏览"
    echo "  esal                   # 详细信息"
    echo "  esa-tree 2             # 查看结构"
    echo "  esa-recent 10          # 最近修改"
    echo "  esa-git                # Git状态"
}

function eza_help() {
    echo "🚀 EZA 快捷命令帮助"
    echo "=================="
    echo ""
    echo "🔧 当前配置:"
    echo "  EZA命令路径: $EZA_CMD"
    echo ""
    echo "⚡ 简短别名 (推荐使用):"
    echo "  esa               - 基础列表显示"
    echo "  esal              - 详细列表 + Git状态"
    echo "  esat              - 树状显示"
    echo "  esa-help          - 显示此帮助 (你正在使用!)"
    echo "  esa-h             - 帮助简写"
    echo "  esa-man           - 原始帮助文档"
    echo "  esa-options       - 详细选项说明"
    echo "  esa-examples      - 使用示例"
    echo ""
    echo "📋 传统别名:"
    echo "  ls, ll, la, l     - 基本列表显示"
    echo "  lt, lt2, lt3      - 树状显示 (1-3层)"
    echo "  lg, lgs           - Git 状态显示"
    echo "  ld, lf, lh        - 目录/文件/带表头显示"
    echo ""
    echo "🔄 排序别名:"
    echo "  lm                - 按修改时间排序"
    echo "  ls_size           - 按文件大小排序"
    echo "  ln, le            - 按名称/扩展名排序"
    echo ""
    echo "⚡ 简短功能函数:"
    echo "  esa-recent [数量] [路径] - 最近修改的文件"
    echo "  esa-big [数量] [路径]    - 最大的文件"
    echo "  esa-stat [路径]         - 文件统计信息"
    echo "  esa-git [路径]          - Git 状态查看"
    echo "  esa-type [路径]         - 按类型分类显示"
    echo "  esa-proj [路径]         - 项目概览"
    echo "  esa-tree [层数] [路径]   - 智能树状显示"
    echo ""
    echo "⚡ 完整功能函数:"
    echo "  ltree [层数] [路径]   - 智能树状显示"
    echo "  lproj [路径]         - 项目概览"
    echo "  lrecent [数量] [路径] - 最近修改的文件"
    echo "  lbig [数量] [路径]    - 最大的文件"
    echo "  lstat [路径]         - 文件统计信息"
    echo "  lgit [路径]          - Git 状态查看"
    echo "  ltype [路径]         - 按类型分类显示"
    echo "  lcompare 目录1 目录2  - 比较两个目录"
    echo "  lsearch             - 交互式搜索(需要fzf)"
    echo ""
    echo "🎨 主题切换:"
    echo "  esa-dark            - 切换到暗色主题"
    echo "  esa-light           - 切换到亮色主题"
    echo "  eza_dark            - 暗色主题 (完整命名)"
    echo "  eza_light           - 亮色主题 (完整命名)"
    echo ""
    echo "📖 帮助命令 (推荐简短版本):"
    echo "  esa-help            - 显示此帮助信息 ⭐"
    echo "  esa-man             - 显示EZA原始帮助文档 ⭐"
    echo "  esa-options         - 显示详细选项说明 ⭐"
    echo "  esa-examples        - 显示使用示例 ⭐"
    echo ""
    echo "⌨️  快捷键:"
    echo "  Ctrl+L              - 清屏并显示当前目录"
    echo ""
    echo "💡 提示:"
    echo "  - 推荐使用 'esa-' 开头的简短命令，更容易输入"
    echo "  - 大部分函数支持路径参数，默认为当前目录"
    echo "  - 数量参数通常默认为10"
    echo "  - 使用 Tab 键可以自动补全路径"
    echo "  - 所有命令都支持颜色和图标显示"
    echo ""
    echo "🔗 更多信息:"
    echo "  官方文档: https://github.com/eza-community/eza"
    echo "  使用 'esa-man' 查看完整选项列表"
    echo "  使用 'esa-examples' 查看详细使用示例"
}

# 模块加载成功标记
export ZSH_MODULE_EZA_LOADED=1
