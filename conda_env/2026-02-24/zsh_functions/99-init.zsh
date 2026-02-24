# =============================================================================
#  99-init.zsh - 应用程序初始化模块 (必须最后加载)
#  Application Initialization Module (Must Load Last)
# =============================================================================

# =============================================================================
#  模块加载验证函数 (Module Loading Verification)
# =============================================================================

# 检查所有模块加载状态
check_zsh_modules() {
  echo "🔍 ZSH 模块加载状态检查"
  echo "=========================="

  local modules=(
    "ZSH_PATH_CONFIG:00-path-config.zsh:路径配置"
    "ZSH_MODULE_ENV_PATH:00-env-path.zsh:环境路径"
    "ZSH_MODULE_CORE_LOADED:01-core.zsh:核心配置"
    "ZSH_MODULE_PLUGINS_LOADED:02-plugins.zsh:插件系统"
    "ZSH_MODULE_APPEARANCE_LOADED:03-appearance.zsh:外观主题"
    "ZSH_MODULE_ALIASES_LOADED:04-aliases.zsh:别名系统"
    "ZSH_MODULE_UTILS_LOADED:05-utils.zsh:工具函数"
    "ZSH_MODULE_CLUSTER_TOOLS_LOADED:06-cluster-tools.zsh:集群工具"
    "ZSH_MODULE_DATA_PROCESSING_LOADED:07-data-processing.zsh:数据处理"
    "ZSH_MODULE_BACKUP_SYSTEM_LOADED:08-backup-system.zsh:备份系统"
    "ZSH_MODULE_BIO_TOOLS_LOADED:09-bio-tools.zsh:生物工具"
    "ZSH_MODULE_EZA_LOADED:10.eza.zsh:EZA配置"
    "ZSH_MODULE_INIT_LOADED:99-init.zsh:初始化"
  )

  local loaded_count=0
  local total_count=${#modules[@]}
  local failed_modules=()

  for module_info in "${modules[@]}"; do
    local var_name="${module_info%%:*}"
    local rest="${module_info#*:}"
    local file_name="${rest%%:*}"
    local desc="${rest#*:}"

    if [[ -n "${(P)var_name}" ]]; then
      echo "✅ $desc ($file_name)"
      ((loaded_count++))
    else
      echo "❌ $desc ($file_name) - 未加载"
      failed_modules+=("$file_name")
    fi
  done

  echo ""
  echo "📊 统计: $loaded_count / $total_count 个模块已加载"

  if [[ $loaded_count -eq $total_count ]]; then
    echo "🎉 所有模块加载成功！"
    return 0
  else
    echo "⚠️  以下模块未加载:"
    for module in "${failed_modules[@]}"; do
      echo "   • $module"
    done
    echo ""
    echo "💡 建议: 检查 ~/.zshrc 中的 source 顺序"
    return 1
  fi
}

# 显示模块加载摘要（简化版）
zsh_modules_summary() {
  local loaded=0
  local modules=(
    "ZSH_PATH_CONFIG" "ZSH_MODULE_ENV_PATH" "ZSH_MODULE_CORE_LOADED"
    "ZSH_MODULE_PLUGINS_LOADED" "ZSH_MODULE_APPEARANCE_LOADED"
    "ZSH_MODULE_ALIASES_LOADED" "ZSH_MODULE_UTILS_LOADED"
    "ZSH_MODULE_CLUSTER_TOOLS_LOADED" "ZSH_MODULE_DATA_PROCESSING_LOADED"
    "ZSH_MODULE_BACKUP_SYSTEM_LOADED" "ZSH_MODULE_BIO_TOOLS_LOADED"
    "ZSH_MODULE_EZA_LOADED" "ZSH_MODULE_INIT_LOADED"
  )

  for mod in "${modules[@]}"; do
    [[ -n "${(P)mod}" ]] && ((loaded++))
  done

  echo "ZSH: $loaded/14 modules loaded"
}

# 快捷别名
alias zsh-check='check_zsh_modules'
alias zsh-modules='check_zsh_modules'
alias zsh-status='zsh_modules_summary'

# -----------------------------------------------------------------------------
#  Conda/Mamba 初始化 (Conda/Mamba Initialization)
# -----------------------------------------------------------------------------
# >>> conda initialize >>>
__conda_setup="$('${CONDA_EXE:-$HOME/miniforge3/bin/conda}' 'shell.zsh' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "${MINIFORGE3_DIR:-$HOME/miniforge3}/etc/profile.d/conda.sh" ]; then
        . "${MINIFORGE3_DIR:-$HOME/miniforge3}/etc/profile.d/conda.sh"
    else
        export PATH="${MINIFORGE3_DIR:-$HOME/miniforge3}/bin:$PATH"
    fi
fi
unset __conda_setup
# <<< conda initialize <<<

# >>> mamba initialize >>>
export MAMBA_EXE="${MAMBA_EXE:-$HOME/miniforge3/bin/mamba}"
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$HOME/miniforge3}"
__mamba_setup="$("$MAMBA_EXE" shell hook --shell zsh --root-prefix "$MAMBA_ROOT_PREFIX" 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__mamba_setup"
else
    alias mamba="$MAMBA_EXE"
fi
unset __mamba_setup
# <<< mamba initialize <<<

# -----------------------------------------------------------------------------
#  Terminal 集成 (Terminal Integration)
# -----------------------------------------------------------------------------
# Tabby Terminal 集成
precmd() {
  echo -n "\x1b]1337;CurrentDir=$(pwd)\x07"
}

# -----------------------------------------------------------------------------
#  提示符初始化 (Prompt Initialization) - 必须在最后
# -----------------------------------------------------------------------------
# Starship 提示符 - 必须是文件的最后一行
if command -v starship >/dev/null 2>&1; then
    eval "$(starship init zsh)"
fi

# 模块加载成功标记
export ZSH_MODULE_INIT_LOADED=1