# =============================================================================
#  01-core.zsh - Shell核心配置模块
#  Shell Core Configuration Module
# =============================================================================

# -----------------------------------------------------------------------------
#  历史记录配置 (History Configuration)
# -----------------------------------------------------------------------------
export HISTSIZE=100000
export SAVEHIST=100000
export HISTFILE=~/.zsh_history

setopt APPEND_HISTORY          # 追加历史，而不是覆盖
setopt SHARE_HISTORY           # 在所有shell间共享历史记录
setopt HIST_IGNORE_ALL_DUPS    # 忽略连续重复的命令

# -----------------------------------------------------------------------------
#  Shell选项 (Shell Options)
# -----------------------------------------------------------------------------
setopt AUTO_CD                 # 自动跳转目录，输入目录名即可cd

# -----------------------------------------------------------------------------
#  补全系统 (Completion System)
# -----------------------------------------------------------------------------
autoload -Uz compinit
compinit

# 补全系统高级配置
zstyle ':completion:*' menu select
zstyle ':completion:*' list-colors "${(s.:.)LS_COLORS}"
zstyle ':completion:*' matcher-list 'm:{a-zA-Z}={A-Za-z}' 'r:|[._-]=* r:|=*' 'l:|=* r:|=*'

# 模块加载成功标记
export ZSH_MODULE_CORE_LOADED=1