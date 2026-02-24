# -----------------------------------------------------------------------------
#  [9] 环境备份与管理函数 (Environment Backup & Management)
# -----------------------------------------------------------------------------

# === 9.1 Conda环境备份系统 ===
# 颜色定义（用于备份函数）
_BACKUP_RED='\033[0;31m'
_BACKUP_GREEN='\033[0;32m'
_BACKUP_YELLOW='\033[1;33m'
_BACKUP_BLUE='\033[0;34m'
_BACKUP_NC='\033[0m'

# 打印带颜色的消息（备份专用）
_backup_info() { echo -e "${_BACKUP_BLUE}ℹ️  $1${_BACKUP_NC}"; }
_backup_success() { echo -e "${_BACKUP_GREEN}✅ $1${_BACKUP_NC}"; }
_backup_error() { echo -e "${_BACKUP_RED}❌ $1${_BACKUP_NC}"; }
_backup_warning() { echo -e "${_BACKUP_YELLOW}⚠️  $1${_BACKUP_NC}"; }
_backup_progress() { echo -e "${_BACKUP_BLUE}▶️  $1${_BACKUP_NC}"; }

# 获取conda环境列表 (优化版)
_get_conda_envs() {
    # 检查conda命令是否存在，如果不存在则提前报错并退出
    if ! command -v conda &> /dev/null; then
        _backup_error "Conda command not found. Cannot list environments."
        return 1
    fi

    # 优先使用jq，因为它对JSON的解析最健壮
    if command -v jq &> /dev/null; then
        # 使用管道直接处理，更高效，避免了临时文件和磁盘I/O
        command conda env list --json | jq -r '.envs[]' 2>/dev/null
    else
        # 如果jq不存在，提供一个不依赖JSON的、更稳定的后备方案
        # 这个方案解析conda env list的文本输出，而不是JSON
        _backup_warning "jq command not found. Using a text-based fallback method."
        command conda env list | grep -v '^#' | awk 'NF>1 && $1 != "base" {print $NF}'
    fi
}

# 导出单个conda环境
_export_conda_env() {
    local env_path="$1"
    local output_dir="$2"
    
    if [ -z "$env_path" ]; then
        return
    fi
    
    local env_name=$(basename "$env_path")
    
    # 跳过base环境
    if [[ "$env_path" != *"/envs/"* ]]; then
        _backup_info "跳过 'base' 环境 ($env_name)..."
        return
    fi
    
    _backup_progress "正在导出环境: $env_name ..."
    
    local output_file="$output_dir/conda/${env_name}.yml"
    
    if conda env export -n "$env_name" --no-builds > "$output_file" 2>/dev/null; then
        _backup_success "成功导出到 $output_file"
    else
        _backup_error "导出环境 '$env_name' 失败！"
        [ -f "$output_file" ] && rm -f "$output_file"
    fi
}

# 备份.zshrc文件
_backup_zshrc() {
    local destination_dir="$1"
    local zshrc_path="$HOME/.zshrc"
    
    echo "-------------------------------------"
    _backup_progress "正在尝试备份 .zshrc 文件..."
    
    if [ -f "$zshrc_path" ]; then
        local destination_file="$destination_dir/zshrc"
        if cp "$zshrc_path" "$destination_file" 2>/dev/null; then
            _backup_success "成功备份 .zshrc 到 $destination_file"
        else
            _backup_error "备份 .zshrc 失败！"
        fi
    else
        _backup_info "未找到 ~/.zshrc 文件，跳过备份。"
    fi
}

# 备份zsh/functions目录（简化版）
_backup_zsh_functions() {
    local destination_dir="$1"
    local source_path="$HOME/zsh/functions"
    
    echo "-------------------------------------"
    _backup_progress "正在尝试备份 ~/zsh/functions 目录..."
    
    if [ ! -d "$source_path" ]; then
        _backup_info "未找到 ~/zsh/functions 目录，跳过备份。"
        return 0
    fi
    
    local dest_path="$destination_dir/zsh_functions"
    
    if ! mkdir -p "$dest_path" 2>/dev/null; then
        _backup_error "无法创建目标目录 $dest_path"
        return 1
    fi
    
    # 直接复制整个目录
    if cp -r "$source_path"/* "$dest_path/" 2>/dev/null; then
        # 保持可执行权限
        find "$dest_path" -type f -exec chmod --reference="$source_path"/{} {} \; 2>/dev/null || \
        find "$dest_path" -type f -executable -exec chmod +x {} \; 2>/dev/null || true
        
        local file_count=$(find "$dest_path" -type f | wc -l)
        _backup_success "成功备份 ~/zsh/functions 到 $dest_path"
        echo "  📊 备份了 $file_count 个函数文件"
        
        # 创建简单的文件列表
        echo "# ~/zsh/functions 备份文件列表 - $(date)" > "$dest_path/files_list.txt"
        find "$dest_path" -type f -printf "%P\n" 2>/dev/null | sort >> "$dest_path/files_list.txt" || \
        find "$dest_path" -type f | sed "s|$dest_path/||" | sort >> "$dest_path/files_list.txt"
        
    else
        _backup_error "备份 ~/zsh/functions 失败！"
        return 1
    fi
}

# =============================================================================
# 备份joblogs目录（跳过大于20M的文件）
# =============================================================================
_backup_joblogs() {
    local destination_dir="$1"
    local source_path="$HOME/joblogs"
    local max_size_mb=20
    local max_size_bytes=$((max_size_mb * 1024 * 1024))
    
    echo "-------------------------------------"
    _backup_progress "正在尝试备份 ~/joblogs 目录（跳过 >${max_size_mb}M 的文件）..."
    
    if [ ! -d "$source_path" ]; then
        _backup_info "未找到 ~/joblogs 目录，跳过备份。"
        return 0
    fi
    
    local dest_path="$destination_dir/joblogs"
    
    if ! mkdir -p "$dest_path" 2>/dev/null; then
        _backup_error "无法创建目标目录 $dest_path"
        return 1
    fi
    
    local backed_up_count=0
    local skipped_count=0
    local skipped_files=""
    
    # 遍历所有文件并根据大小选择性复制
    while IFS= read -r -d '' file; do
        local file_size=$(stat -f%z "$file" 2>/dev/null || stat -c%s "$file" 2>/dev/null)
        local relative_path="${file#$source_path/}"
        local target_file="$dest_path/$relative_path"
        
        if [ "$file_size" -gt "$max_size_bytes" ]; then
            # 文件超过10M，跳过
            ((skipped_count++))
            local size_mb=$(echo "scale=2; $file_size / 1024 / 1024" | bc 2>/dev/null || echo "$(($file_size / 1024 / 1024))")
            skipped_files="${skipped_files}\n  - $relative_path (${size_mb}M)"
        else
            # 文件小于等于10M，复制
            local target_dir=$(dirname "$target_file")
            mkdir -p "$target_dir" 2>/dev/null
            
            if cp "$file" "$target_file" 2>/dev/null; then
                chmod --reference="$file" "$target_file" 2>/dev/null || \
                chmod 644 "$target_file" 2>/dev/null || true
                ((backed_up_count++))
            else
                _backup_error "复制文件失败: $relative_path"
            fi
        fi
    done < <(find "$source_path" -type f -print0)
    
    if [ "$backed_up_count" -gt 0 ]; then
        local total_size=$(du -sh "$dest_path" 2>/dev/null | cut -f1 || echo "未知")
        
        _backup_success "成功备份 ~/joblogs 到 $dest_path"
        echo "  📊 备份了 $backed_up_count 个日志文件"
        echo "  ⏭️  跳过了 $skipped_count 个大文件 (>${max_size_mb}M)"
        echo "  💾 总大小: $total_size"
        
        # 创建备份信息文件
        echo "# ~/joblogs 备份文件列表 - $(date)" > "$dest_path/backup_info.txt"
        echo "# 备份文件数量: $backed_up_count" >> "$dest_path/backup_info.txt"
        echo "# 跳过文件数量: $skipped_count (大于 ${max_size_mb}M)" >> "$dest_path/backup_info.txt"
        echo "# 备份总大小: $total_size" >> "$dest_path/backup_info.txt"
        echo "# =================================" >> "$dest_path/backup_info.txt"
        echo "" >> "$dest_path/backup_info.txt"
        
        # 记录已备份的文件
        echo "## 已备份的文件:" >> "$dest_path/backup_info.txt"
        find "$dest_path" -type f -not -name "backup_info.txt" -printf "%P\n" 2>/dev/null | sort >> "$dest_path/backup_info.txt" || \
        find "$dest_path" -type f -not -name "backup_info.txt" | sed "s|$dest_path/||" | sort >> "$dest_path/backup_info.txt"
        
        # 记录跳过的文件
        if [ "$skipped_count" -gt 0 ]; then
            echo "" >> "$dest_path/backup_info.txt"
            echo "## 跳过的大文件 (>${max_size_mb}M):" >> "$dest_path/backup_info.txt"
            echo -e "$skipped_files" >> "$dest_path/backup_info.txt"
        fi
        
        return 0
    else
        _backup_error "备份 ~/joblogs 失败或没有文件被备份！"
        return 1
    fi
}

# 备份submitted_jobs.txt文件
_backup_submitted_jobs() {
    local destination_dir="$1"
    local submitted_jobs_path="$HOME/submitted_jobs.txt"

    echo "-------------------------------------"
    _backup_progress "正在尝试备份 submitted_jobs.txt 文件..."

    if [ -f "$submitted_jobs_path" ]; then
        local destination_file="$destination_dir/submitted_jobs.txt"
        if cp "$submitted_jobs_path" "$destination_file" 2>/dev/null; then
            _backup_success "成功备份 submitted_jobs.txt 到 $destination_file"
        else
            _backup_error "备份 submitted_jobs.txt 失败！"
        fi
    else
        _backup_info "未找到 ~/submitted_jobs.txt 文件，跳过备份。"
    fi
}

# 备份history_commands.txt文件
_backup_history_commands() {
    local destination_dir="$1"
    local history_commands_path="$HOME/history_commands.txt"

    echo "-------------------------------------"
    _backup_progress "正在尝试备份 history_commands.txt 文件..."

    if [ -f "$history_commands_path" ]; then
        local destination_file="$destination_dir/history_commands.txt"
        if cp "$history_commands_path" "$destination_file" 2>/dev/null; then
            _backup_success "成功备份 history_commands.txt 到 $destination_file"
        else
            _backup_error "备份 history_commands.txt 失败！"
        fi
    else
        _backup_info "未找到 ~/history_commands.txt 文件，跳过备份。"
    fi
}

# 备份scripts目录（简化版）
_backup_scripts() {
    local destination_dir="$1"
    local source_path="$HOME/software/scripts"
    
    echo "-------------------------------------"
    _backup_progress "正在尝试备份 ~/software/scripts 目录..."

    if [ ! -d "$source_path" ]; then
        _backup_info "未找到 ~/software/scripts 目录，跳过备份。"
        return 0
    fi
    
    local dest_path="$destination_dir/scripts"
    
    if ! mkdir -p "$dest_path" 2>/dev/null; then
        _backup_error "无法创建目标目录 $dest_path"
        return 1
    fi
    
    # 直接复制整个目录
    if cp -r "$source_path"/* "$dest_path/" 2>/dev/null; then
        # 保持原有权限
        find "$dest_path" -type f -exec chmod --reference="$source_path"/{} {} \; 2>/dev/null || \
        find "$dest_path" -type f -exec chmod 644 {} \; 2>/dev/null || true
        
        local file_count=$(find "$dest_path" -type f | wc -l)
        local total_size=$(du -sh "$dest_path" 2>/dev/null | cut -f1 || echo "未知")
        
        _backup_success "成功备份 ~/software/scripts 到 $dest_path"
        echo "  📊 备份了 $file_count 个脚本文件"
        echo "  💾 总大小: $total_size"
        
        # 创建简单的文件列表和统计信息
        echo "# ~/software/scripts 备份文件列表 - $(date)" > "$dest_path/backup_info.txt"
        echo "# 备份文件数量: $file_count" >> "$dest_path/backup_info.txt"
        echo "# 备份总大小: $total_size" >> "$dest_path/backup_info.txt"
        echo "# =================================" >> "$dest_path/backup_info.txt"
        echo "" >> "$dest_path/backup_info.txt"
        
        find "$dest_path" -type f -not -name "backup_info.txt" -printf "%P\n" 2>/dev/null | sort >> "$dest_path/backup_info.txt" || \
        find "$dest_path" -type f -not -name "backup_info.txt" | sed "s|$dest_path/||" | sort >> "$dest_path/backup_info.txt"
        
    else
        _backup_error "备份 ~/software/scripts 失败！"
        return 1
    fi
}

# 备份biopytools目录 - 使用统一路径配置
_backup_biopytools() {
    local destination_dir="$1"
    local source_path="${BIOPYTOOLS_DIR:-$HOME/software/biopytools}"

    echo "-------------------------------------"
    _backup_progress "正在尝试备份 biopytools 目录..."

    if [ ! -d "$source_path" ]; then
        _backup_info "未找到 $source_path 目录，跳过备份。"
        return 0
    fi

    local dest_path="$destination_dir/biopytools"

    # 直接复制整个目录
    if cp -r "$source_path" "$dest_path" 2>/dev/null; then
        local total_size=$(du -sh "$dest_path" 2>/dev/null | cut -f1 || echo "未知")
        _backup_success "成功备份 biopytools 到 $dest_path"
        echo "  💾 总大小: $total_size"
    else
        _backup_error "备份 biopytools 失败！"
        return 1
    fi
}

# 修复后的 _backup_local_bin 函数
_backup_local_bin() {
    local destination_dir="$1"
    local source_path="$HOME/.local/bin"
    
    echo "-------------------------------------"
    _backup_progress "正在尝试备份 ~/.local/bin 目录..."
    
    if [ ! -d "$source_path" ]; then
        _backup_info "未找到 $source_path 目录，跳过备份。"
        return 0
    fi
    
    local dest_path="$destination_dir/local/bin"
    
    if ! mkdir -p "$dest_path" 2>/dev/null; then
        _backup_error "无法创建目标目录 $dest_path"
        return 1
    fi
    
    # 创建详细信息文件
    local symlink_info_file="$dest_path/symlink_info.txt"
    local restore_script="$dest_path/restore_symlinks.sh"
    local ls_output_file="$dest_path/ls_output.txt"
    
    # 1. 生成 ls -la 风格的输出文件
    _backup_progress "生成详细链接信息..."
    echo "# ~/.local/bin 目录软链接详情 - $(date)" > "$ls_output_file"
    echo "# 生成时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$ls_output_file"
    # 修复：分别计算软链接和普通文件数量
    local total_symlinks=$(find "$source_path" -maxdepth 1 -type l 2>/dev/null | wc -l)
    local total_files=$(find "$source_path" -maxdepth 1 -type f 2>/dev/null | wc -l)
    echo "# 总文件数: $((total_symlinks + total_files))" >> "$ls_output_file"
    echo "#" >> "$ls_output_file"
    
    # 使用 ls -la 获取详细信息
    ls -la "$source_path" >> "$ls_output_file" 2>/dev/null
    
    # 2. 生成结构化的软链接信息文件
    echo "# ~/.local/bin 软链接信息" > "$symlink_info_file"
    echo "# 格式: 软链接名称 -> 目标路径" >> "$symlink_info_file"
    echo "# 生成时间: $(date '+%Y-%m-%d %H:%M:%S')" >> "$symlink_info_file"
    echo "#" >> "$symlink_info_file"
    
    # 3. 生成简化的恢复脚本
    cat > "$restore_script" << 'EOF'
#!/bin/bash
# ~/.local/bin 软链接恢复脚本

echo "🔗 开始恢复 ~/.local/bin 软链接..."
echo "=================================="

# 确保目标目录存在
mkdir -p ~/.local/bin

success_count=0
error_count=0

EOF
    
    # 4. 复制文件并处理软链接
    local symlink_count=0
    local regular_file_count=0
    
    # 处理软链接
    find "$source_path" -maxdepth 1 -type l 2>/dev/null | while IFS= read -r file; do
        if [ -n "$file" ]; then
            local filename=$(basename "$file")
            local target=$(readlink "$file")
            
            # 记录软链接信息
            echo "$filename -> $target" >> "$symlink_info_file"
            
            # 添加到恢复脚本
            cat >> "$restore_script" << SYMLINK_EOF
# 恢复软链接: $filename
if [ ! -e "\$HOME/.local/bin/$filename" ]; then
    if ln -s "$target" "\$HOME/.local/bin/$filename" 2>/dev/null; then
        echo "✅ 创建软链接: $filename"
        ((success_count++))
    else
        echo "❌ 创建软链接失败: $filename"
        ((error_count++))
    fi
else
    echo "⚠️  已存在: $filename"
fi

SYMLINK_EOF
            ((symlink_count++))
        fi
    done
    
    # 处理普通文件
    find "$source_path" -maxdepth 1 -type f 2>/dev/null | while IFS= read -r file; do
        if [ -n "$file" ]; then
            local filename=$(basename "$file")
            
            # 复制文件
            if cp "$file" "$dest_path/" 2>/dev/null; then
                # 保持可执行权限
                [ -x "$file" ] && chmod +x "$dest_path/$filename" 2>/dev/null || true
                
                # 记录文件信息
                local file_size=$(stat -c%s "$file" 2>/dev/null || echo "unknown")
                echo "# 普通文件: $filename ($file_size bytes)" >> "$symlink_info_file"
                
                # 添加到恢复脚本
                cat >> "$restore_script" << FILE_EOF
# 恢复普通文件: $filename  
if [ ! -f "\$HOME/.local/bin/$filename" ]; then
    if cp "\$(dirname "\$0")/$filename" "\$HOME/.local/bin/$filename" 2>/dev/null; then
        chmod +x "\$HOME/.local/bin/$filename" 2>/dev/null || true
        echo "✅ 复制文件: $filename"
        ((success_count++))
    else
        echo "❌ 复制文件失败: $filename"
        ((error_count++))
    fi
else
    echo "⚠️  文件已存在: $filename"
fi

FILE_EOF
                ((regular_file_count++))
            fi
        fi
    done
    
    # 完成恢复脚本
    cat >> "$restore_script" << 'EOF'
echo ""
echo "📊 恢复完成统计:"
echo "  成功: $success_count"
echo "  失败: $error_count"
echo "✅ 恢复脚本执行完毕"
EOF
    
    # 设置脚本可执行权限
    chmod +x "$restore_script" 2>/dev/null || true
    
    # 写入统计到信息文件
    cat >> "$symlink_info_file" << EOF

# 统计信息
# 软链接数量: $symlink_count
# 普通文件数量: $regular_file_count
# 总计: $((symlink_count + regular_file_count))
EOF
    
    # 显示结果
    _backup_success "成功备份 ~/.local/bin 到 $dest_path"
    echo "  📊 发现 $symlink_count 个软链接，$regular_file_count 个普通文件"
    echo "  📄 详细信息: $symlink_info_file"
    echo "  📋 ls输出: $ls_output_file"  
    echo "  🔧 恢复脚本: $restore_script"
    echo ""
    echo "  🔗 快速恢复命令:"
    echo "     bash $restore_script"
}

# 新增：备份 ~/.config 目录（简化版）
_backup_config() {
    local destination_dir="$1"
    local config_source="$HOME/.config"
    
    echo "-------------------------------------"
    _backup_progress "正在尝试备份 ~/.config 目录..."
    
    if [ ! -d "$config_source" ]; then
        _backup_info "未找到 ~/.config 目录，跳过备份。"
        return 0
    fi
    
    local dest_path="$destination_dir/config"
    
    if ! mkdir -p "$dest_path" 2>/dev/null; then
        _backup_error "无法创建目标目录 $dest_path"
        return 1
    fi
    
    # 简单复制
    if cp -r "$config_source"/* "$dest_path/" 2>/dev/null; then
        _backup_success "成功备份 ~/.config 到 $dest_path"
        echo "  📊 已备份配置目录"
        
        # 创建简单标记文件
        echo "Config backup completed - $(date)" > "$dest_path/.backup_info"
    else
        _backup_error "备份 ~/.config 失败！"
        return 1
    fi
}

# 修复后的 _backup_shell_history 函数
_backup_shell_history() {
    local destination_dir="$1"
    
    echo "-------------------------------------"
    _backup_progress "正在备份命令行历史记录..."
    
    local history_backed_up=false
    local history_dir="$destination_dir/shell_history"
    
    # 确保历史记录目录存在
    if ! mkdir -p "$history_dir" 2>/dev/null; then
        _backup_error "无法创建历史记录目录: $history_dir"
        return 1
    fi
    
    # 备份zsh历史记录
    if [ -f "$HOME/.zsh_history" ]; then
        local zsh_dest="$history_dir/zsh_history"
        
        if cp "$HOME/.zsh_history" "$zsh_dest" 2>/dev/null; then
            _backup_success "成功备份 zsh 历史记录到 $zsh_dest"
            history_backed_up=true
            
            # 创建可读性更好的格式
            local readable_dest="$history_dir/zsh_history_readable.txt"
            if command -v fc &> /dev/null && [ -n "$ZSH_VERSION" ]; then
                # 使用fc命令获取格式化的历史记录
                fc -l 1 > "$readable_dest" 2>/dev/null || \
                awk -F';' '{if(NF>1) print $2}' "$HOME/.zsh_history" > "$readable_dest" 2>/dev/null
                [ -f "$readable_dest" ] && _backup_success "创建可读历史记录: zsh_history_readable.txt"
            fi
        else
            _backup_error "备份 zsh 历史记录失败！"
        fi
    fi
    
    # 备份bash历史记录
    if [ -f "$HOME/.bash_history" ]; then
        local bash_dest="$history_dir/bash_history"
        
        if cp "$HOME/.bash_history" "$bash_dest" 2>/dev/null; then
            _backup_success "成功备份 bash 历史记录到 $bash_dest"
            history_backed_up=true
        else
            _backup_error "备份 bash 历史记录失败！"
        fi
    fi
    
    # 创建历史记录统计信息
    if [ "$history_backed_up" = true ]; then
        local stats_file="$history_dir/history_stats.txt"
        echo "命令行历史记录统计 - $(date)" > "$stats_file"
        echo "==============================" >> "$stats_file"
        
        # zsh历史统计
        if [ -f "$history_dir/zsh_history" ]; then
            local zsh_count=$(wc -l < "$history_dir/zsh_history" 2>/dev/null || echo "0")
            echo "ZSH 历史记录条数: $zsh_count" >> "$stats_file"
            
            # 最常用的命令（前10个）
            echo "" >> "$stats_file"
            echo "最常用的命令 (ZSH):" >> "$stats_file"
            awk -F';' '{if(NF>1) print $2}' "$history_dir/zsh_history" 2>/dev/null | \
                awk '{print $1}' | \
                sort | uniq -c | sort -nr | head -10 >> "$stats_file" 2>/dev/null
        fi
        
        # bash历史统计
        if [ -f "$history_dir/bash_history" ]; then
            local bash_count=$(wc -l < "$history_dir/bash_history" 2>/dev/null || echo "0")
            echo "" >> "$stats_file"
            echo "BASH 历史记录条数: $bash_count" >> "$stats_file"
            
            echo "" >> "$stats_file"
            echo "最常用的命令 (BASH):" >> "$stats_file"
            awk '{print $1}' "$history_dir/bash_history" 2>/dev/null | \
                sort | uniq -c | sort -nr | head -10 >> "$stats_file" 2>/dev/null
        fi
        
        _backup_success "创建历史记录统计: history_stats.txt"
    else
        _backup_info "未找到任何shell历史记录文件"
    fi
}

# 主备份函数（回到原版本加上config备份）
backup_envs() {
    local BASE_OUTPUT_DIR="conda_env_backups"
    local custom_dir=""
    local verbose=false
    local backup_config=true
    
    # 简化的参数解析
    while [[ $# -gt 0 ]]; do
        case $1 in
            -d|--dir)
                custom_dir="$2"
                shift 2
                ;;
            -v|--verbose)
                verbose=true
                shift
                ;;
            --no-config)
                backup_config=false
                shift
                ;;
            -h|--help)
                cat << 'EOF'
用法: backup_envs [选项]

Conda环境备份函数 - 导出所有conda环境并备份相关配置文件

选项:
  -h, --help     显示此帮助信息
  -d, --dir DIR  指定备份基础目录 (默认: conda_env_backups)
  -v, --verbose  详细输出模式
  --no-config    跳过 ~/.config 目录备份

功能:
  • 导出所有conda环境到YAML文件
  • 备份 ~/.zshrc 配置文件
  • 备份 ~/.local/bin 目录
  • 备份 ~/zsh/functions 目录
  • 备份 ~/software/scripts 目录
  • 备份 biopytools 工具目录
  • 备份 ~/.config 目录（可选）
  • 备份 shell 历史记录
  • 按日期组织备份文件

示例:
  backup_envs                    # 使用默认设置进行备份
  backup_envs -d my_backups      # 指定备份目录
  backup_envs -v                 # 详细输出模式
  backup_envs --no-config        # 跳过配置目录备份

备份文件将保存在: $BASE_OUTPUT_DIR/YYYY-MM-DD/ 目录中
EOF
                return 0
                ;;
            *)
                _backup_error "未知选项: $1"
                echo "使用 backup_envs -h 查看帮助信息"
                return 1
                ;;
        esac
    done
    
    # 使用自定义目录（如果指定）
    if [ -n "$custom_dir" ]; then
        BASE_OUTPUT_DIR="$custom_dir"
    fi
    
    # 启用详细输出
    if [ "$verbose" = true ]; then
        set -x
    fi
    
    echo "🔄 Conda环境备份开始..."
    echo "=================================="
    
    # 检查conda是否可用
    if ! command -v conda &> /dev/null; then
        _backup_error "未找到 conda 命令。请确保 Conda 已正确安装并配置。"
        return 1
    fi
    
    # 创建备份目录
    local today_str=$(date '+%Y-%m-%d')
    local dated_output_dir="$BASE_OUTPUT_DIR/$today_str"
    
    if [ ! -d "$dated_output_dir" ]; then
        if mkdir -p "$dated_output_dir"; then
            _backup_success "创建归档目录: $dated_output_dir"
        else
            _backup_error "无法创建备份目录: $dated_output_dir"
            return 1
        fi
    else
        _backup_info "使用现有归档目录: $dated_output_dir"
    fi

    # 创建conda子文件夹存放conda环境的yml文件
    local conda_files_dir="$dated_output_dir/conda"
    if ! mkdir -p "$conda_files_dir"; then
        _backup_error "无法创建Conda备份子目录: $conda_files_dir"
        return 1
    fi
    
    # 获取并导出conda环境
    _backup_progress "正在获取 Conda 环境列表..."
    local env_paths=()
    while IFS= read -r env_path; do
        if [ -n "$env_path" ]; then
            env_paths+=("$env_path")
        fi
    done < <(_get_conda_envs)
    
    if [ ${#env_paths[@]} -gt 0 ]; then
        echo ""
        _backup_info "发现 ${#env_paths[@]} 个 Conda 环境。开始导出..."
        _backup_info "💡 提示：如果某个环境导出很慢，可以按 Ctrl+C 跳过"
        echo "-------------------------------------"
        
        for env_path in "${env_paths[@]}"; do
            [ -n "$env_path" ] && _export_conda_env "$env_path" "$dated_output_dir"
            echo "-------------------------------------"
        done
    else
        _backup_warning "未找到任何 Conda 环境。"
    fi
    
    # 备份其他文件
    _backup_zshrc "$dated_output_dir"
    _backup_local_bin "$dated_output_dir"
    _backup_zsh_functions "$dated_output_dir"
    _backup_joblogs "$dated_output_dir"
    _backup_submitted_jobs "$dated_output_dir"
    _backup_history_commands "$dated_output_dir"
    _backup_scripts "$dated_output_dir"
    _backup_biopytools "$dated_output_dir"
    
    # 备份 ~/.config 目录（如果启用）
    if [ "$backup_config" = true ]; then
        _backup_config "$dated_output_dir"
    else
        _backup_info "跳过 ~/.config 目录备份（--no-config 选项）"
    fi
    
    _backup_shell_history "$dated_output_dir"
    
    echo ""
    echo "🎉 所有任务完成！"
    _backup_success "备份文件已保存在 '$dated_output_dir' 目录中。"
    
    # 显示备份摘要
    echo ""
    echo "📋 备份内容摘要:"
    echo "-------------------------------------"
    if [ -d "$dated_output_dir" ]; then
        local yml_count=$(find "$dated_output_dir" -name "*.yml" 2>/dev/null | wc -l)
        echo "  • Conda环境文件: $yml_count 个"
        
        [ -f "$dated_output_dir/zshrc" ] && echo "  • .zshrc 配置文件: ✓"
        [ -d "$dated_output_dir/local/bin" ] && echo "  • ~/.local/bin 目录: ✓"
        [ -d "$dated_output_dir/zsh_functions" ] && echo "  • ~/zsh/functions 目录: ✓"
        [ -d "$dated_output_dir/scripts" ] && echo "  • ~/software/scripts 目录: ✓"
        [ -d "$dated_output_dir/biopytools" ] && echo "  • biopytools 工具目录: ✓"
        
        # 配置目录备份信息
        if [ -d "$dated_output_dir/config" ]; then
            echo "  • ~/.config 目录: ✓"
        fi
        
        # 历史记录备份信息
        if [ -d "$dated_output_dir/shell_history" ]; then
            echo "  • Shell历史记录: ✓"
            [ -f "$dated_output_dir/shell_history/zsh_history" ] && \
                echo "    └─ ZSH历史: $(wc -l < "$dated_output_dir/shell_history/zsh_history" 2>/dev/null || echo "0") 条"
            [ -f "$dated_output_dir/shell_history/bash_history" ] && \
                echo "    └─ BASH历史: $(wc -l < "$dated_output_dir/shell_history/bash_history" 2>/dev/null || echo "0") 条"
            [ -f "$dated_output_dir/shell_history/history_stats.txt" ] && \
                echo "    └─ 统计文件: ✓"
        fi
        
        local size=$(du -sh "$dated_output_dir" 2>/dev/null | cut -f1)
        echo "  • 总文件大小: $size"
    fi
    
    # 关闭详细输出
    if [ "$verbose" = true ]; then
        set +x
    fi
}

# === 9.2 智能增量备份系统 ===
# 状态记录目录
BACKUP_STATE_DIR="$HOME/.conda_backup_state"

# 获取历史记录指纹（用于变化检测）
_get_history_fingerprint() {
    local fingerprint=""
    
    # zsh历史记录指纹
    if [ -f "$HOME/.zsh_history" ]; then
        local zsh_fp=$(tail -100 "$HOME/.zsh_history" 2>/dev/null | md5sum | cut -d' ' -f1)
        fingerprint="${fingerprint}zsh:${zsh_fp};"
    fi
    
    # bash历史记录指纹
    if [ -f "$HOME/.bash_history" ]; then
        local bash_fp=$(tail -100 "$HOME/.bash_history" 2>/dev/null | md5sum | cut -d' ' -f1)
        fingerprint="${fingerprint}bash:${bash_fp};"
    fi
    
    # fish历史记录指纹
    if [ -f "$HOME/.local/share/fish/fish_history" ]; then
        local fish_fp=$(tail -100 "$HOME/.local/share/fish/fish_history" 2>/dev/null | md5sum | cut -d' ' -f1)
        fingerprint="${fingerprint}fish:${fish_fp};"
    fi
    
    echo "$fingerprint"
}

# 检查历史记录是否有变化
_history_has_changed() {
    local state_file="$BACKUP_STATE_DIR/shell_history.state"
    
    # 获取当前历史记录指纹
    local current_fingerprint=$(_get_history_fingerprint)
    
    # 如果无法获取指纹，认为有变化
    if [ -z "$current_fingerprint" ]; then
        return 0  # 有变化
    fi
    
    # 如果状态文件不存在，认为是新的
    if [ ! -f "$state_file" ]; then
        return 0  # 有变化（新的）
    fi
    
    # 比较指纹
    local last_fingerprint=$(head -1 "$state_file" 2>/dev/null)
    if [ "$current_fingerprint" != "$last_fingerprint" ]; then
        return 0  # 有变化
    else
        return 1  # 无变化
    fi
}

# 更新历史记录状态
_update_history_state() {
    local state_file="$BACKUP_STATE_DIR/shell_history.state"
    
    # 确保状态目录存在
    mkdir -p "$BACKUP_STATE_DIR"
    
    # 保存当前指纹
    local current_fingerprint=$(_get_history_fingerprint)
    
    if [ -n "$current_fingerprint" ]; then
        echo "$current_fingerprint" > "$state_file"
        echo "$(date '+%Y-%m-%d %H:%M:%S')" >> "$state_file"
    fi
}

_get_env_fingerprint() {
    local env_name="$1"
    if [ -z "$env_name" ]; then
        return 1
    fi
    
    # 生成环境指纹：包含包名和版本的MD5哈希
    conda list -n "$env_name" --json 2>/dev/null | \
        jq -r '.[] | "\(.name)=\(.version)"' 2>/dev/null | \
        sort | \
        md5sum | \
        cut -d' ' -f1
}

# 获取简化的环境指纹（不依赖jq）
_get_env_fingerprint_simple() {
    local env_name="$1"
    if [ -z "$env_name" ]; then
        return 1
    fi
    
    # 不使用jq的简化方法
    conda list -n "$env_name" 2>/dev/null | \
        awk 'NR>3 && !/^#/ && NF>=2 {print $1"="$2}' | \
        sort | \
        md5sum | \
        cut -d' ' -f1
}

# 检查环境是否有变化
_env_has_changed() {
    local env_name="$1"
    local state_file="$BACKUP_STATE_DIR/${env_name}.state"
    
    # 获取当前指纹
    local current_fingerprint
    if command -v jq &> /dev/null; then
        current_fingerprint=$(_get_env_fingerprint "$env_name")
    else
        current_fingerprint=$(_get_env_fingerprint_simple "$env_name")
    fi
    
    # 如果无法获取指纹，认为有变化
    if [ -z "$current_fingerprint" ]; then
        return 0  # 有变化
    fi
    
    # 如果状态文件不存在，认为是新环境
    if [ ! -f "$state_file" ]; then
        return 0  # 有变化（新环境）
    fi
    
    # 比较指纹
    local last_fingerprint=$(cat "$state_file" 2>/dev/null)
    if [ "$current_fingerprint" != "$last_fingerprint" ]; then
        return 0  # 有变化
    else
        return 1  # 无变化
    fi
}

# 更新环境状态记录
_update_env_state() {
    local env_name="$1"
    local state_file="$BACKUP_STATE_DIR/${env_name}.state"
    
    # 确保状态目录存在
    mkdir -p "$BACKUP_STATE_DIR"
    
    # 保存当前指纹
    local current_fingerprint
    if command -v jq &> /dev/null; then
        current_fingerprint=$(_get_env_fingerprint "$env_name")
    else
        current_fingerprint=$(_get_env_fingerprint_simple "$env_name")
    fi
    
    if [ -n "$current_fingerprint" ]; then
        echo "$current_fingerprint" > "$state_file"
        echo "$(date '+%Y-%m-%d %H:%M:%S')" >> "$state_file"
    fi
}

# 检查是否有新环境或环境变化
check_env_changes() {
    local verbose=false
    local show_details=false
    
    # 参数解析
    while [[ $# -gt 0 ]]; do
        case $1 in
            -v|--verbose)
                verbose=true
                shift
                ;;
            -d|--details)
                show_details=true
                shift
                ;;
            -h|--help)
                cat << 'EOF'
用法: check_env_changes [选项]

检查conda环境和历史记录是否有变化

选项:
  -v, --verbose   详细输出
  -d, --details   显示变化详情
  -h, --help      显示帮助

检查内容:
  • Conda环境变化（新环境或包更新）
  • Shell历史记录变化

返回值:
  0 - 有变化
  1 - 无变化
EOF
                return 0
                ;;
            *)
                echo "未知选项: $1"
                return 1
                ;;
        esac
    done
    
    local has_changes=false
    local changed_envs=()
    local new_envs=()
    local history_changed=false
    
    # 获取当前所有环境（使用内置路径解析避免函数调用问题）
    local current_envs=()
    while IFS= read -r env_path; do
        if [ -n "$env_path" ]; then
            # 直接使用bash内置的路径解析，避免函数调用
            local env_name="${env_path##*/}"
            if [ -n "$env_name" ]; then
                current_envs+=("$env_name")
            fi
        fi
    done < <(_get_conda_envs)
    
    if [ "$verbose" = true ]; then
        echo "🔍 正在检查环境和历史记录变化..."
        echo "发现 ${#current_envs[@]} 个conda环境"
    fi
    
    # 检查每个环境
    for env_name in "${current_envs[@]}"; do
        # 跳过base环境
        if [ "$env_name" = "base" ] || [ -z "$env_name" ]; then
            continue
        fi
        
        local state_file="$BACKUP_STATE_DIR/${env_name}.state"
        
        if [ ! -f "$state_file" ]; then
            # 新环境
            new_envs+=("$env_name")
            has_changes=true
            [ "$verbose" = true ] && echo "🆕 新环境: $env_name"
        elif _env_has_changed "$env_name"; then
            # 环境有变化
            changed_envs+=("$env_name")
            has_changes=true
            [ "$verbose" = true ] && echo "🔄 环境有变化: $env_name"
        else
            [ "$verbose" = true ] && echo "✅ 环境无变化: $env_name"
        fi
    done
    
    # 检查历史记录变化
    if _history_has_changed; then
        history_changed=true
        has_changes=true
        [ "$verbose" = true ] && echo "📜 历史记录有变化"
    else
        [ "$verbose" = true ] && echo "✅ 历史记录无变化"
    fi
    
    # 显示摘要
    if [ "$show_details" = true ] || [ "$verbose" = true ]; then
        echo ""
        echo "📊 变化摘要:"
        echo "  新环境: ${#new_envs[@]} 个"
        echo "  变化环境: ${#changed_envs[@]} 个"
        echo "  历史记录变化: $([ "$history_changed" = true ] && echo "是" || echo "否")"
        
        if [ ${#new_envs[@]} -gt 0 ]; then
            echo "  新环境列表: ${new_envs[*]}"
        fi
        if [ ${#changed_envs[@]} -gt 0 ]; then
            echo "  变化环境列表: ${changed_envs[*]}"
        fi
    fi
    
    if [ "$has_changes" = true ]; then
        [ "$verbose" = true ] && echo "🎯 检测到变化，建议进行备份"
        return 0
    else
        [ "$verbose" = true ] && echo "😌 所有内容无变化"
        return 1
    fi
}

# 智能备份函数（仅在有变化时备份）
smart_backup_envs() {
    local force_backup=false
    local backup_args=()
    
    # 参数解析
    while [[ $# -gt 0 ]]; do
        case $1 in
            -f|--force)
                force_backup=true
                shift
                ;;
            -h|--help)
                cat << 'EOF'
用法: smart_backup_envs [选项]

智能备份系统 - 仅在检测到环境变化时进行备份

选项:
  -f, --force     强制备份（忽略变化检测）
  -d, --dir DIR   指定备份目录
  -v, --verbose   详细输出
  --no-config     跳过配置目录备份
  -h, --help      显示帮助

其他选项会传递给 backup_envs 函数

功能:
  • 自动检测环境变化
  • 仅在有变化时进行备份
  • 更新环境状态记录
  • 备份完成后自动压缩
  • 支持强制备份模式
  • 支持所有 backup_envs 的功能（包括biopytools备份）

输出:
  • YYYY-MM-DD/ - 原始备份目录
  • YYYY-MM-DD_backup.tar.gz - 压缩文件
EOF
                return 0
                ;;
            *)
                backup_args+=("$1")
                shift
                ;;
        esac
    done
    
    echo "🧠 智能备份检查开始..."
    echo "=========================="
    
    # 检查是否有变化
    if [ "$force_backup" = false ]; then
        if ! check_env_changes -v; then
            _backup_info "没有检测到变化，跳过备份"
            echo "💡 提示: 使用 'smart_backup_envs -f' 强制备份"
            return 0
        fi
        echo ""
    else
        _backup_warning "强制备份模式，跳过变化检测"
        echo ""
    fi
    
    # 执行备份
    _backup_progress "检测到变化，开始备份..."
    if backup_envs "${backup_args[@]}"; then
        # 备份成功后，更新所有环境和历史记录的状态记录
        _backup_progress "更新状态记录..."
        
        # 更新环境状态
        local env_paths=()
        while IFS= read -r env_path; do
            if [ -n "$env_path" ]; then
                env_paths+=("$env_path")
            fi
        done < <(_get_conda_envs)
        
        for env_path in "${env_paths[@]}"; do
            # 使用bash内置路径解析
            local env_name="${env_path##*/}"
            if [[ "$env_path" == *"/envs/"* ]] && [ -n "$env_name" ]; then
                _update_env_state "$env_name"
            fi
        done
        
        # 更新历史记录状态
        _update_history_state

        _backup_success "状态记录已更新"

        # 压缩备份目录
        _backup_progress "开始压缩备份目录..."
        local today_str=$(date '+%Y-%m-%d')
        local backup_dir="${BASE_OUTPUT_DIR:-$HOME/conda_env_backups}/$today_str"
        if _compress_backup_directory "$backup_dir"; then
            _backup_success "备份压缩完成"
        else
            _backup_warning "备份压缩失败，但备份已完成"
        fi

        return 0
    else
        _backup_error "备份失败"
        return 1
    fi
}

# 压缩备份目录
_compress_backup_directory() {
    local backup_dir="$1"

    # 检查是否提供了备份目录
    if [ -z "$backup_dir" ]; then
        _backup_error "未提供备份目录路径"
        return 1
    fi

    # 检查备份目录是否存在
    if [ ! -d "$backup_dir" ]; then
        _backup_error "备份目录不存在: $backup_dir"
        return 1
    fi

    # 确保使用绝对路径
    if [[ ! "$backup_dir" = /* ]]; then
        backup_dir="$(pwd)/$backup_dir"
    fi

    # 生成压缩文件名
    local base_dir=$(dirname "$backup_dir")
    local dir_name=$(basename "$backup_dir")
    local archive_name="${backup_dir}_backup.tar.gz"

    # 检查是否已经存在压缩文件
    if [ -f "$archive_name" ]; then
        _backup_warning "压缩文件已存在，将覆盖: $archive_name"
        rm -f "$archive_name"
    fi

    # 检查是否安装了tar
    if ! command -v tar &> /dev/null; then
        _backup_error "未找到tar命令，无法压缩"
        return 1
    fi

    # 获取原始大小
    local original_size=$(du -sh "$backup_dir" 2>/dev/null | cut -f1)

    _backup_progress "正在压缩备份目录: $backup_dir"

    # 执行压缩
    if tar -czf "$archive_name" -C "$base_dir" "$dir_name" 2>/dev/null; then
        # 获取压缩后大小
        local compressed_size=$(du -sh "$archive_name" 2>/dev/null | cut -f1)

        _backup_success "压缩成功: $archive_name"
        echo "  📊 原始大小: $original_size"
        echo "  📦 压缩后大小: $compressed_size"
        echo "  📁 压缩文件: $archive_name"

        return 0
    else
        _backup_error "压缩失败"
        return 1
    fi
}

# 便捷别名
alias backup-envs='backup_envs'
alias backup-conda='backup_envs'
alias smart-backup='smart_backup_envs'
alias check-changes='check_env_changes'