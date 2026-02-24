#!/bin/bash

# vbump - Version Bump Tool for pyproject.toml with Git integration and CHANGELOG support
# Add this to your .zshrc or .bashrc to use globally

# 配置默认值
VBUMP_EDITOR=${VBUMP_EDITOR:-${EDITOR:-vim}}
VBUMP_REMOTE=${VBUMP_REMOTE:-origin}
VBUMP_BRANCH=${VBUMP_BRANCH:-main}
VBUMP_COMMIT_TEMPLATE=${VBUMP_COMMIT_TEMPLATE:-"chore: bump version to v{version}"}
VBUMP_AUTO_OPEN_CHANGELOG=${VBUMP_AUTO_OPEN_CHANGELOG:-false}

# 主函数
vbump() {
    local version_type=""
    local custom_message=""
    local additional_files=()
    local no_push=false
    local no_changelog=false
    local changelog_edit=false
    local dry_run=false
    local help=false

    # 参数解析
    while [[ $# -gt 0 ]]; do
        case $1 in
            patch|minor|major)
                version_type="$1"
                shift
                ;;
            -m|--message)
                custom_message="$2"
                shift 2
                ;;
            -f|--files)
                shift  # 移除 -f 参数本身
                # 读取后续所有参数，直到遇到下一个选项（以-开头）或参数结束
                while [[ $# -gt 0 && ! "$1" =~ ^- ]]; do
                    additional_files+=("$1")
                    shift
                done
                ;;
            -e|--changelog-edit)
                changelog_edit=true
                shift
                ;;
            --no-push)
                no_push=true
                shift
                ;;
            --no-changelog)
                no_changelog=true
                shift
                ;;
            --dry-run)
                dry_run=true
                shift
                ;;
            -h|--help)
                help=true
                shift
                ;;
            *)
                echo "❌ 未知参数: $1"
                _vbump_help
                return 1
                ;;
        esac
    done

    # 显示帮助
    if [[ "$help" == true ]]; then
        _vbump_help
        return 0
    fi

    # 检查版本类型参数
    if [[ -z "$version_type" ]]; then
        echo "❌ 错误: 必须指定版本类型 (patch|minor|major)"
        _vbump_help
        return 1
    fi

    # 前置检查
    if ! _vbump_pre_checks; then
        return 1
    fi

    # 获取当前版本
    local current_version
    current_version=$(_vbump_get_current_version)
    if [[ $? -ne 0 ]]; then
        echo "❌ 无法获取当前版本"
        return 1
    fi

    # 计算新版本
    local new_version
    new_version=$(_vbump_increment_version "$current_version" "$version_type")
    if [[ $? -ne 0 ]]; then
        echo "❌ 无法计算新版本号"
        return 1
    fi

    echo "🔄 版本变更: $current_version → $new_version"

    # Dry run 模式
    if [[ "$dry_run" == true ]]; then
        echo "🔍 Dry run 模式，将要执行的操作："
        echo "  - 更新 pyproject.toml 版本号到 $new_version"
        [[ "$no_changelog" == false ]] && echo "  - 更新 CHANGELOG.md"
        echo "  - Git add 更改的文件:"
        echo "    • pyproject.toml"
        [[ "$no_changelog" == false ]] && echo "    • CHANGELOG.md"
        for file in "${additional_files[@]}"; do
            echo "    • $file"
        done
        local commit_msg
        if [[ -n "$custom_message" ]]; then
            if [[ "$custom_message" =~ version.*[0-9]+\.[0-9]+\.[0-9]+ ]]; then
                commit_msg="$custom_message"
            else
                commit_msg="version $new_version: $custom_message"
            fi
        else
            commit_msg=$(echo "$VBUMP_COMMIT_TEMPLATE" | sed "s/{version}/$new_version/g")
        fi
        echo "  - Git commit: $commit_msg"
        echo "  - Git tag: v$new_version"
        [[ "$no_push" == false ]] && echo "  - Git push 到 $VBUMP_REMOTE/$VBUMP_BRANCH"
        return 0
    fi

    # 更新版本号
    if ! _vbump_update_version "$new_version"; then
        echo "❌ 更新版本号失败"
        return 1
    fi

    # 处理 CHANGELOG
    if [[ "$no_changelog" == false ]]; then
        if ! _vbump_handle_changelog "$new_version" "$changelog_edit" "$custom_message" "${additional_files[@]}"; then
            echo "❌ 处理 CHANGELOG 失败"
            # 回滚版本更改
            _vbump_update_version "$current_version"
            return 1
        fi
    fi

    # Git 操作
    if ! _vbump_git_operations "$new_version" "$custom_message" "$no_push" "${additional_files[@]}"; then
        echo "❌ Git 操作失败"
        # 回滚更改
        _vbump_update_version "$current_version"
        [[ "$no_changelog" == false ]] && git checkout HEAD -- CHANGELOG.md 2>/dev/null
        return 1
    fi

    echo "✅ 版本更新完成: $current_version → $new_version"
    return 0
}

# 显示帮助信息
_vbump_help() {
    cat << 'EOF'
📦 vbump - Version Bump Tool

用法:
    vbump <version_type> [选项]

版本类型:
    patch       递增补丁版本 (x.y.Z)
    minor       递增次版本 (x.Y.0)  
    major       递增主版本 (X.0.0)

选项:
    -m, --message <msg>     自定义 git commit 信息 (会自动补充版本号前缀)
    -f, --files <files...>  指定额外要提交的文件/文件夹 (支持多个参数)
    -e, --changelog-edit    交互式编辑 CHANGELOG
    --no-push              不推送到远程仓库
    --no-changelog         跳过 CHANGELOG 更新
    --dry-run              预览模式，不执行实际操作
    -h, --help             显示此帮助信息

示例:
    vbump patch                                    # 递增patch版本
    vbump minor -m "add new feature"               # 自动补充版本前缀: "version 0.1.0: add new feature"
    vbump minor -m "version 0.1.0: add feature"   # 手动指定完整版本信息(不会重复添加)
    vbump major -e                                 # 递增major版本并编辑CHANGELOG
    vbump patch --dry-run                          # 预览patch版本更新
    vbump minor --no-push                          # 更新版本但不推送
    vbump patch -f dist/ docs/                     # 同时提交dist和docs文件夹
    vbump minor -f setup.py package.json           # 同时提交多个文件
    vbump patch -f dist/ setup.py docs/ -m "major release"  # 简化输入，自动补充版本号

环境变量配置:
    VBUMP_EDITOR              编辑器 (默认: $EDITOR 或 vim)
    VBUMP_REMOTE              远程仓库名 (默认: origin)
    VBUMP_BRANCH              分支名 (默认: main)
    VBUMP_COMMIT_TEMPLATE     提交信息模板 (默认: "chore: bump version to v{version}")
    VBUMP_AUTO_OPEN_CHANGELOG 是否自动打开CHANGELOG编辑 (默认: false)

智能功能:
    • 自动版本号补充: -m "fix bug" → "version 1.2.3: fix bug"
    • 多文件灵活添加: -f dist/ docs/ setup.py
    • 错误自动回滚: 操作失败时恢复原始状态
    • Git状态检查: 确保仓库状态安全

文件要求:
    - pyproject.toml  (必须存在且包含version字段)
    - CHANGELOG.md    (如不存在会自动创建)
EOF
}

# 前置检查
_vbump_pre_checks() {
    # 检查是否在git仓库中
    if ! git rev-parse --git-dir >/dev/null 2>&1; then
        echo "❌ 错误: 当前目录不是git仓库"
        return 1
    fi

    # 检查pyproject.toml是否存在
    if [[ ! -f "pyproject.toml" ]]; then
        echo "❌ 错误: 未找到 pyproject.toml 文件"
        return 1
    fi

    # 检查工作区是否干净
    if [[ -n "$(git status --porcelain)" ]]; then
        echo "⚠️  警告: 工作区有未提交的更改"
        echo "继续操作？(y/N)"
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            echo "操作已取消"
            return 1
        fi
    fi

    # 检查当前分支
    local current_branch
    current_branch=$(git branch --show-current)
    if [[ "$current_branch" != "$VBUMP_BRANCH" ]]; then
        echo "⚠️  当前分支: $current_branch (期望: $VBUMP_BRANCH)"
        echo "继续操作？(y/N)"
        read -r response
        if [[ ! "$response" =~ ^[Yy]$ ]]; then
            echo "操作已取消"
            return 1
        fi
    fi

    return 0
}

# 获取当前版本
_vbump_get_current_version() {
    local version
    local line
    
    # 获取版本行
    line=$(grep -E '^version\s*=' pyproject.toml)
    
    if [[ -z "$line" ]]; then
        echo "❌ 在 pyproject.toml 中未找到版本号" >&2
        return 1
    fi
    
    # 提取版本号
    version=$(echo "$line" | sed -E 's/version\s*=\s*["\x27]([^"\x27]+)["\x27]/\1/')
    
    if [[ -z "$version" ]]; then
        echo "❌ 无法从 pyproject.toml 中读取版本号: $line" >&2
        return 1
    fi
    
    # 验证版本号格式
    if [[ ! "$version" =~ ^[0-9]+\.[0-9]+\.[0-9]+([.-].*)?$ ]]; then
        echo "❌ 版本号格式无效: '$version' (期望格式: x.y.z)" >&2
        return 1
    fi
    
    # 提取主要版本号部分 (去除后缀如 -alpha, -beta 等)
    version=$(echo "$version" | sed -E 's/^([0-9]+\.[0-9]+\.[0-9]+).*/\1/')
    
    echo "$version"
}

# 递增版本号
_vbump_increment_version() {
    local current_version="$1"
    local increment_type="$2"
    
    # 使用cut命令解析版本号，更可靠
    local major minor patch
    major=$(echo "$current_version" | cut -d. -f1)
    minor=$(echo "$current_version" | cut -d. -f2)  
    patch=$(echo "$current_version" | cut -d. -f3)
    
    # 验证解析结果
    if [[ -z "$major" || -z "$minor" || -z "$patch" ]]; then
        echo "❌ 版本号解析失败: '$current_version'" >&2
        return 1
    fi
    
    # 验证是否都是数字
    if [[ ! "$major" =~ ^[0-9]+$ ]] || [[ ! "$minor" =~ ^[0-9]+$ ]] || [[ ! "$patch" =~ ^[0-9]+$ ]]; then
        echo "❌ 版本号包含非数字字符: '$current_version'" >&2
        return 1
    fi
    
    # 根据类型递增
    case "$increment_type" in
        major)
            major=$((major + 1))
            minor=0
            patch=0
            ;;
        minor)
            minor=$((minor + 1))
            patch=0
            ;;
        patch)
            patch=$((patch + 1))
            ;;
        *)
            echo "❌ 无效的递增类型: $increment_type" >&2
            return 1
            ;;
    esac
    
    echo "$major.$minor.$patch"
}

# 更新版本号
_vbump_update_version() {
    local new_version="$1"
    
    # 创建备份
    cp pyproject.toml pyproject.toml.backup
    
    # 使用perl进行替换，更可靠
    if command -v perl >/dev/null 2>&1; then
        # 使用perl替换，支持各种引号格式
        perl -i -pe "s/(version\\s*=\\s*[\"'])([^\"']+)([\"'])/\${1}$new_version\${3}/" pyproject.toml
    else
        # 如果没有perl，使用更简单的sed方式
        # 分别处理双引号和单引号，避免复杂的正则
        if grep -q 'version.*=.*"' pyproject.toml; then
            # 处理双引号情况
            sed -i.tmp 's/version.*=.*"[^"]*"/version = "'"$new_version"'"/' pyproject.toml
        elif grep -q "version.*=.*'" pyproject.toml; then
            # 处理单引号情况  
            sed -i.tmp "s/version.*=.*'[^']*'/version = '$new_version'/" pyproject.toml
        else
            echo "❌ 无法识别版本号格式" >&2
            mv pyproject.toml.backup pyproject.toml
            return 1
        fi
        # 清理临时文件
        command rm -f pyproject.toml.tmp
    fi
    
    # 验证更新是否成功
    local updated_version
    updated_version=$(_vbump_get_current_version 2>/dev/null)
    
    if [[ "$updated_version" == "$new_version" ]]; then
        command rm -f pyproject.toml.backup
        echo "✅ 版本号已更新到 $new_version"
        return 0
    else
        echo "❌ 版本号更新失败，恢复原文件" >&2
        mv pyproject.toml.backup pyproject.toml
        return 1
    fi
}

# 处理 CHANGELOG
_vbump_handle_changelog() {
    local new_version="$1"
    local edit_mode="$2" 
    local custom_message="$3"
    shift 3  # 移除前三个参数，剩下的都是额外文件
    local additional_files=("$@")
    local changelog_file="CHANGELOG.md"
    
    # 如果CHANGELOG不存在，创建一个
    if [[ ! -f "$changelog_file" ]]; then
        cat > "$changelog_file" << 'EOF'
# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

EOF
    fi
    
    local date_str=$(date +%Y-%m-%d)
    
    # 构建新的版本条目
    local new_entry=""
    new_entry+="## [$new_version] - $date_str"$'\n'
    new_entry+=""$'\n'
    new_entry+="### Changed"$'\n'
    
    # 添加commit信息（去除版本号前缀）
    if [[ -n "$custom_message" ]]; then
        local clean_message
        clean_message=$(echo "$custom_message" | sed -E 's/^version [0-9]+\.[0-9]+\.[0-9]+: //')
        new_entry+="- $clean_message"$'\n'
    fi
    
    # 添加文件变更信息
    if [[ ${#additional_files[@]} -gt 0 ]]; then
        new_entry+="- Updated files: $(IFS=', '; echo "${additional_files[*]}")"$'\n'
    fi
    
    # 如果没有任何信息，添加默认信息
    if [[ -z "$custom_message" && ${#additional_files[@]} -eq 0 ]]; then
        new_entry+="- Version bump"$'\n'
    fi
    
    new_entry+=""$'\n'
    
    if [[ "$edit_mode" == true ]] || [[ "$VBUMP_AUTO_OPEN_CHANGELOG" == true ]]; then
        # 交互式编辑模式
        local temp_file=$(mktemp)
        
        # 预填充编辑内容，包含新版本信息
        cat > "$temp_file" << EOF
## [$new_version] - $date_str

### Added
- 

### Changed
EOF
        
        # 添加commit信息（去除版本号前缀）
        if [[ -n "$custom_message" ]]; then
            local clean_message
            clean_message=$(echo "$custom_message" | sed -E 's/^version [0-9]+\.[0-9]+\.[0-9]+: //')
            echo "- $clean_message" >> "$temp_file"
        fi
        
        # 添加文件变更信息
        if [[ ${#additional_files[@]} -gt 0 ]]; then
            echo "- Updated files: $(IFS=', '; echo "${additional_files[*]}")" >> "$temp_file"
        fi
        
        cat >> "$temp_file" << EOF

### Fixed
- 

### Removed
- 

---
参考最近的提交记录:
EOF
        
        # 添加最近的commits作为参考
        git log --oneline -5 --pretty=format:"- %s" >> "$temp_file"
        
        echo -e "\n\n=== 请编辑上方内容，保存后关闭编辑器 ===" >> "$temp_file"
        
        # 打开编辑器
        "$VBUMP_EDITOR" "$temp_file"
        
        # 提取用户编辑的内容（去除参考信息）
        local user_content
        user_content=$(sed '/^---$/q' "$temp_file" | head -n -1)
        
        # 重新组装CHANGELOG，新版本在最前面
        _vbump_insert_changelog_entry "$changelog_file" "$user_content"

        command rm -f "$temp_file"
        echo "✅ CHANGELOG 已更新 (交互式编辑)"
    else
        # 自动生成模式 - 新版本插入到最前面
        _vbump_insert_changelog_entry "$changelog_file" "$new_entry"
        echo "✅ CHANGELOG 已自动更新"
    fi
    
    return 0
}

# 插入CHANGELOG条目到正确位置（最新版本在前）
_vbump_insert_changelog_entry() {
    local changelog_file="$1"
    local new_entry="$2"
    local temp_changelog=$(mktemp)
    
    # 找到第一个版本条目的行号
    local first_version_line
    first_version_line=$(grep -n "^## \[" "$changelog_file" | head -1 | cut -d: -f1)
    
    if [[ -n "$first_version_line" ]]; then
        # 有现有版本条目：插入到第一个版本条目之前
        head -n $((first_version_line - 1)) "$changelog_file" > "$temp_changelog"
        echo "$new_entry" >> "$temp_changelog"
        tail -n +$first_version_line "$changelog_file" >> "$temp_changelog"
    else
        # 没有现有版本条目：插入到文件末尾
        cat "$changelog_file" > "$temp_changelog"
        echo "$new_entry" >> "$temp_changelog"
    fi
    
    mv "$temp_changelog" "$changelog_file"
}

# Git 操作
_vbump_git_operations() {
    local new_version="$1"
    local custom_message="$2"
    local no_push="$3"
    shift 3  # 移除前三个参数，剩下的都是额外文件
    local additional_files=("$@")
    
    # 准备提交信息
    local commit_message
    if [[ -n "$custom_message" ]]; then
        # 检查用户消息是否已经包含版本号
        if [[ "$custom_message" =~ version.*[0-9]+\.[0-9]+\.[0-9]+ ]]; then
            # 用户已经包含版本号，直接使用
            commit_message="$custom_message"
        else
            # 自动添加版本号前缀
            commit_message="version $new_version: $custom_message"
        fi
    else
        commit_message=$(echo "$VBUMP_COMMIT_TEMPLATE" | sed "s/{version}/$new_version/g")
    fi
    
    echo "📝 提交信息: $commit_message"
    
    # Git add - 首先添加必要的文件
    git add pyproject.toml
    [[ -f "CHANGELOG.md" ]] && git add CHANGELOG.md
    
    # 添加用户指定的额外文件
    if [[ ${#additional_files[@]} -gt 0 ]]; then
        for file in "${additional_files[@]}"; do
            if [[ -e "$file" ]]; then
                git add "$file"
                echo "✅ 已添加文件: $file"
            else
                echo "⚠️  文件不存在，跳过: $file"
            fi
        done
    fi
    
    # Git commit
    if ! git commit -m "$commit_message"; then
        echo "❌ Git commit 失败" >&2
        return 1
    fi
    
    echo "✅ Git commit 成功: $commit_message"
    
    # 创建tag
    local tag_name="v$new_version"
    if ! git tag "$tag_name"; then
        echo "⚠️  创建标签失败，可能已存在: $tag_name"
    else
        echo "✅ 创建标签: $tag_name"
    fi
    
    # Git push
    if [[ "$no_push" == false ]]; then
        if git push "$VBUMP_REMOTE" "$VBUMP_BRANCH" && git push "$VBUMP_REMOTE" --tags; then
            echo "✅ 推送到远程仓库成功"
        else
            echo "❌ 推送到远程仓库失败" >&2
            return 1
        fi
    else
        echo "⏭️  跳过推送到远程仓库"
    fi
    
    return 0
}

# 创建别名以便更好的使用体验
alias vb='vbump'