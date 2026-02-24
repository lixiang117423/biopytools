# =============================================================================
#  02-plugins.zsh - 插件与扩展功能模块
#  Plugins & Extensions Module
# =============================================================================

# -----------------------------------------------------------------------------
#  自动建议插件 (Auto-suggestions Plugin)
# -----------------------------------------------------------------------------
if [ -f ~/.oh-my-zsh/custom/plugins/zsh-autosuggestions/zsh-autosuggestions.zsh ]; then
  source ~/.oh-my-zsh/custom/plugins/zsh-autosuggestions/zsh-autosuggestions.zsh
  ZSH_AUTOSUGGEST_HIGHLIGHT_STYLE="fg=#666666"
  ZSH_AUTOSUGGEST_STRATEGY=(history completion)
fi

# -----------------------------------------------------------------------------
#  历史子串搜索 (History Substring Search)
# -----------------------------------------------------------------------------
if [ -f ~/.oh-my-zsh/custom/plugins/zsh-history-substring-search/zsh-history-substring-search.zsh ]; then
  source ~/.oh-my-zsh/custom/plugins/zsh-history-substring-search/zsh-history-substring-search.zsh
  bindkey '^[[A' history-substring-search-up
  bindkey '^[[B' history-substring-search-down
fi

# -----------------------------------------------------------------------------
#  语法高亮 (必须最后加载) (Syntax Highlighting - Must Load Last)
# -----------------------------------------------------------------------------
if [ -f ~/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting/zsh-syntax-highlighting.zsh ]; then
  source ~/.oh-my-zsh/custom/plugins/zsh-syntax-highlighting/zsh-syntax-highlighting.zsh
fi

# 模块加载成功标记
export ZSH_MODULE_PLUGINS_LOADED=1